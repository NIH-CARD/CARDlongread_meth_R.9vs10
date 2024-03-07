import sys
import pysam
import os
import numpy as np
import pandas as pd
import subprocess
from time import time
import pickle
from math import ceil

#Yields on pair at a time allowing for pairwise iteration
def yield_slice(mod_array, pairs, avg_array):
    for pair in pairs:
        yield (mod_array[pair, :], avg_array)

#Load the cpg_df created through modkit motif-bed for CpG
def load_bed_file(path, chrom, strand):
    cpg_df = pd.read_csv(path, sep='\t', header=None)
    cpg_df.columns = ['chrom', 'start', 'stop', 'score1', 'score2', 'strand']
    cpg_df = cpg_df[(cpg_df['chrom'] == chrom) & (cpg_df['strand'] == strand)]
    cpg_df.reset_index(drop=True, inplace = True)
    return cpg_df


#Write all reads to a bed file, then intersect bedfile with itself to get all pairwise intersections
def get_read_pairs(sam_path):
    tmp_fp = f"{os.path.basename(sam_path)}_tmp_read_bed.bed"
    with open(tmp_fp, 'w') as out_file:
        index = 0
        for read in pysam.AlignmentFile(sam_path, "rb").fetch():
            if (read.is_secondary or 
                read.is_supplementary or 
                read.is_unmapped or 
                read.is_duplicate or 
                read.is_qcfail):
                continue
            out_file.write(f"{read.reference_name}\t{read.reference_start}\t{read.reference_end}\t{index}\n")
            index += 1
    with open(f"self_{tmp_fp}", 'w') as self_out:
        subprocess.run(['bedtools','intersect','-wa','-wb','-a',tmp_fp,'-b',tmp_fp], 
                       encoding='utf-8',
                       stdout=self_out
                      )
    read_pairs = set()
    with open(f"self_{tmp_fp}", "r") as read_file:
        for line in read_file:
            split_line = line.strip().split('\t')
            if split_line[3] != split_line[7]:
                read_pairs.add((int(split_line[3]), int(split_line[7])))
    
    subprocess.run(['rm', f"self_{tmp_fp}"])
    subprocess.run(['rm', tmp_fp])
    
    return read_pairs

#Create an array of mod confidence for each read, using these arrays we will
#decide if any given pair of reads has sufficient overlap (20 CpGs) for 
#covariance calculation as well as computing average confidence
def create_mod_array(sam_path, bed_cpg_df):
    
    index_dict = {k:v for k, v in zip(bed_cpg_df["start"],bed_cpg_df.index)}
    samfile = pysam.AlignmentFile(sam_path, "rb")
    x = np.full([samfile.count(read_callback='all'),bed_cpg_df.shape[0]], -1, dtype='short')
    samfile = pysam.AlignmentFile(sam_path, "rb")
    count = 0
    t0 = time()
    for read in samfile.fetch():
        if (read.is_secondary or 
            read.is_supplementary or 
            read.is_unmapped or 
            read.is_duplicate or 
            read.is_qcfail):
            continue
        if count % 10000 == 0 and count != 0:
            t1 = time()
            print(f"{(t1 - t0)/60} minutes to process 10000 reads")
            t0 = t1

        ref_pos = read.get_aligned_pairs()
        ref_index = [x for x in ref_pos if x[0] is not None]

        if len(read.modified_bases_forward) > 0:
            for pos, qual in read.modified_bases_forward[('C', 0, 'm')]:
                if ref_index[pos][1] is not None and ref_index[pos][1] in index_dict:
                    x[count][index_dict[ref_index[pos][1]]] = qual
        count += 1
        
    return x

#Calculate the covariance between mod confidence of two reads
def covariance(pair, mean):
    idx = np.argwhere(np.all(pair != -1, axis=0))
    if len(idx) == 0:
        return (None, None)
    running_sum = 0
    for index in idx:
        index = int(index)
        running_sum += (pair[0][index] - mean[index]) * (pair[1][index] - mean[index])
    return (float(running_sum/len(idx)), len(idx))

#Compute the average methylation at a given position for covariance calculation
def average_array(mod_array):
    #Calculate the average modification rate for each array
    col_mean_array = np.full(mod_array.shape[1], -1, dtype = float) 
    t0 = time()
    processed = 0
    for col_index in range(mod_array.shape[1]):
        if col_index % 100 == 0 and col_index != 0:
            processed += 100
            t1 = time()
            print(f"{(t1 - t0)/60} minutes to process 100 positions")
            print(f"{100 * processed/col_mean_array.shape[0]}% complete")
            t0 = t1
        col_mean_array[col_index] = mod_array[:,col_index][mod_array[:,col_index] != -1].mean()
        
    return col_mean_array

#Subset read pairs for multiprocessing
def subset_read_pairs(read_pair_set, n_subsets):
    subset_size = ceil(len(read_pair_set) / n_subsets)
    read_pair_list = list(read_pair_set)
    for i in range(0, len(read_pair_list), subset_size):
        yield read_pair_list[i:i + subset_size]


def main(sam_path, 
         bed_cpg_path, 
         out_dir, 
         scratch_dir, 
         prefix, 
         chrom='chr21', 
         strand = '+'):
    
    print("Loading bed CpG")
    bed_cpg = load_bed_file(bed_cpg_path, chrom, strand)
    
    print("Creating Mod Array")
    mod_array = create_mod_array(sam_path, bed_cpg)
    
    print("Writing Mod Array to Disk")
    
    mod_array_path = os.path.join(out_dir, "mod_array.npy")
    
    with open(mod_array_path, 'wb') as out_pkl:
        np.save(out_pkl, mod_array)
    
    print("Creating Mod Average Array")
    avg_array = average_array(mod_array)
    
    print("Writing Mod Average Array to Disk")
    
    avg_array_path = os.path.join(out_dir, "avg_mod.npy")
    
    with open(avg_array_path, 'wb') as out_pkl:
        np.save(out_pkl, avg_array)
    
    print("Calculating Read Pairs")
    read_pairs = get_read_pairs(sam_path)
    for pair in read_pairs:
        with open(os.path.join(scratch_dir, 
                               f"{pair[0]}_{pair[1]}_{os.path.basename(out_dir)}_slice.pkl"), 
                  "wb") as out_pkl:
            np.save(out_pkl, mod_array[pair, :])
        
    subset_read_pairs_20 = subset_read_pairs(read_pairs, 20)
    
    print("Subsetting and Processing Read Pairs")
    for i in range(20):
        with open(os.path.join(out_dir, f"{i}_read_pairs.pkl"), 'wb') as out_pkl:
            pickle.dump(next(subset_read_pairs_20), out_pkl)
    
    with open(os.path.join(out_dir, f"{prefix}_covariance.tsv"), 'w'):
        for i in range(20):
            with open(os.path.join(out_dir, f"{i}_read_pairs.pkl"), 'rb') as in_pkl:
                subset_list = pickle.load(in_pkl)
            for pair in subset_list:
                cov_tup = covariance(mod_array[pair, :],
                                     avg_mod_array)
                out_tsv.write(f"{cov_tup[0]}\t{cov_tup[1]}\t{pair[0]}\t{pair[1]}\n")
            
    

if __name__ == "__main__":
    sam_path = sys.argv[1]
    bed_cpg_path = sys.argv[2]
    out_dir = sys.argv[3]
    scratch_dir = sys.argv[4]
    prefix = sys.argv[5]
    main(sam_path, bed_cpg_path, out_dir, scratch_dir, prefix)