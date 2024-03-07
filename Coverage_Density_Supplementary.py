import pandas as pd
import numpy as np
import seaborn as sns
import time
import seaborn as sns
import argparse
import re
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator
plt.rcParams['pdf.fonttype'] = 42
plt.switch_backend('agg')
sns.set(rc={'figure.figsize':(11.7,8.27)})


#Create a dictionary with sequence for each reference chromosome
def sequence_dict(ref_path):
    #open reference and create a dictionary containing the sequences for 
    #each chromosome
    with open(ref_path, 'r') as f:
        chrom = None
        current_seq = []
        sequence_dict = {}
        
        #Iterate file line by line
        for line in f:
            if line[0] == '>' and '_' in line:
                continue
            #If this is the first chromosome set chrom and continue
            if line[0] == '>' and chrom is None:
                chrom = line[1:].split()[0].strip()
            
            #If this is the start of a new chromosome, add old chrom
            elif line[0] == '>':
                sequence_dict[chrom] = "".join(current_seq)
                chrom = line[1:].split()[0].strip()
                current_seq = []
            
            #If line is continuation of current chromosome, append seq
            else:
                current_seq.append(line.strip().upper())
                
        if chrom not in sequence_dict:
            sequence_dict[chrom] = "".join(current_seq)
    return sequence_dict

#Create a list of cpg arrays from the input reference sequence. One pair per
#chrom, an array for positive and an array for negative
def cpg_array(sequence_dict, chrom_index):
    list_np_cpg = []
    for chrom in chrom_index:
        print(f"{chrom} processing cpg")
        list_np_cpg.append(np.full((2, len(sequence_dict[chrom])), False))
        for m in re.finditer('CG', sequence_dict[chrom]):
            list_np_cpg[chrom_index[chrom]][0][m.start(0)] = True
        for m in re.finditer('CG', sequence_dict[chrom]):
            list_np_cpg[chrom_index[chrom]][1][m.start(0) + 1] = True
    return list_np_cpg

#Calculate the cpg density in a given window
def cpg_density_calc(start, flank_dist, cpg_array):
    return np.sum(cpg_array[start-flank_dist:start+flank_dist+1])


#Create a merged df for R9 and R10
def read_merge(R9_path, R10_path):
    
    chroms = ['chr1',
              'chr2',
              'chr3',
              'chr4',
              'chr5',
              'chr6',
              'chr7',
              'chr8',
              'chr9',
              'chr10',
              'chr11',
              'chr12',
              'chr13',
              'chr14',
              'chr15',
              'chr16',
              'chr17',
              'chr18',
              'chr19',
              'chr20',
              'chr21',
              'chr22',
              'chrX',
              'chrY',
              'chrM']
    
    R9_df = pd.read_csv(R9_path,
                    sep='\t',
                    header = None,
                    usecols = [0, 1, 9, 10, 11, 12]
                   )
    R9_df.columns = ['chrom', 'position', 'coverage', 'ratio', 'modified', 'canonical']
    R9_df = R9_df[R9_df['chrom'].isin(chroms)]
        
    R10_df = pd.read_csv(R10_path,
                        sep='\t',
                        header = None,
                        usecols = [0, 1, 9, 10, 11, 12]
                       )
    R10_df.columns = ['chrom', 'position', 'coverage', 'ratio', 'modified', 'canonical']
    R10_df = R10_df[R10_df['chrom'].isin(chroms)]
    
    R10_df = R10_df[(R10_df['coverage'] >= 20) & (R10_df['coverage'] <= 200)]
    R10_df['ratio'] = R10_df['ratio'] / 100
    
    R9_df = R9_df[(R9_df['coverage'] >= 20) & (R9_df['coverage'] <= 200)]
    R9_df['ratio'] = R9_df['ratio'] / 100
    
    merged_df = pd.merge(R10_df, 
                         R9_df, 
                         how = 'inner', 
                         left_on = ['chrom', 'position'], 
                         right_on=['chrom', 'position'])
    
    
    merged_df.columns = ['chrom', 
                         'position', 
                         'R10_coverage',
                         'R10_ratio', 
                         'R10_modified', 
                         'R10_canonical', 
                         'R9_coverage',
                         'R9_ratio', 
                         'R9_modified', 
                         'R9_canonical']
    return merged_df


def main(r9_hg002_path,
         r10_hg002_path, 
         ref_path):
    
        df = read_merge(r9_hg002_path, r10_hg002_path)
        df["delta_coverage"] = df["R10_coverage"] - df["R9_coverage"]
        df["delta_proportion"] = df["R10_ratio"] - df["R9_ratio"]
        
        
        print("calculating density")
        
        chrom_index = {
            'chr1': 0,
            'chr2': 1,
            'chr3': 2,
            'chr4': 3,
            'chr5': 4,
            'chr6': 5,
            'chr7': 6,
            'chr8': 7,
            'chr9': 8,
            'chr10': 9,
            'chr11': 10,
            'chr12': 11,
            'chr13': 12,
            'chr14': 13,
            'chr15': 14,
            'chr16': 15,
            'chr17': 16,
            'chr18': 17,
            'chr19': 18,
            'chr20': 19,
            'chr21': 20,
            'chr22': 21,
            'chrX': 22,
            'chrY': 23,
            'chrM': 24
        }
        
        cpg_density_list = cpg_array(sequence_dict(ref_path), chrom_index)
        
        i = 0
        
        #Calculate the cpg density in a +/- 500 nucleotide window
        cpg_density = [0] * df.shape[0]
        for chrom, center in zip(df["chrom"].to_list(), df["position"].to_list()):
            cpg_density[i] = cpg_density_calc(center, 
                                              500, 
                                              cpg_density_list[chrom_index[chrom]][0])
            i += 1
            if i % 1000000 == 0:
                print(i)

        
        df["density"] = cpg_density
             
        sns.set(rc={'figure.figsize':(11.7,8.27)})
        fig = plt.figure(dpi = 300)
        ax = fig.add_subplot(1, 1, 1)
        ax.set_facecolor('black')
        density = ax.hist2d(x = df['density'].to_list(), 
                            y = df['delta_proportion'].to_list(), 
                            cmap='viridis', 
                            norm=colors.LogNorm(), 
                            bins=20)
        # now plot both limits against eachother
        fig.colorbar(density[3], label='Number of points per square')
        ax.grid(False)
        plt.xlabel('CpG Density in 1000 base window')
        plt.ylabel('R10 Methylation Proportion - R9 Methylation Proportion')
        ax.get_figure().savefig("CpG_Density.pdf", dpi = 300)
        plt.clf()
        
        sns.set(rc={'figure.figsize':(11.7,8.27)})
        fig = plt.figure(dpi = 300)
        ax = fig.add_subplot(1, 1, 1)
        ax.set_facecolor('black')
        density = ax.hist2d(x = df['delta_coverage'].to_list(), 
                            y = df['delta_proportion'].to_list(), 
                            cmap='viridis', 
                            norm=colors.LogNorm(), 
                            bins=20)
        # now plot both limits against eachother
        fig.colorbar(density[3], label='Number of points per square')
        ax.grid(False)
        plt.xlabel('Delta Coverage')
        plt.ylabel('R10 Methylation Proportion - R9 Methylation Proportion')
        ax.get_figure().savefig("Delta_Coverage.pdf", dpi = 300)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    parser.add_argument(
        "--r9_hg002",
        type=str,
        required = True,
        help = "R9 HG002 modkit"
    )
    
    parser.add_argument(
        "--r10_hg002",
        type=str,
        required = True,
        help = "R10 HG002 modkit"
    )
    
    parser.add_argument(
        "--reference",
        type=str,
        required = True,
        help = "Genome to which reads were aligned"
    )

    FLAGS, unparsed = parser.parse_known_args()
    main(FLAGS.r9_hg002, 
         FLAGS.r10_hg002,
         FLAGS.reference)