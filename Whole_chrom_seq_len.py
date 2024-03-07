import sys
import subprocess
import argparse

#Get the sequence of the reference chromosome
def get_ref_seq(path_to_ref, ref_section):
    ref_seq = subprocess.run(['samtools', 
                              'faidx', 
                              f"{path_to_ref}", 
                              f"{ref_section}"],
                             stdout=subprocess.PIPE).stdout.decode('utf-8')
    ref_seq = "".join(ref_seq.split('\n')[1:])
    return ref_seq.upper()

#For the primary considered chromosomes calculate the length in bases
#This is used for creating numpy arrays that are initialized with the correct
#lenght in down stream processing.
def main(path_to_ref, out_path):
    seq_dict = {}
    for chrom in ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 
                  'chr8','chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 
                  'chr15','chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 
                  'chr22','chrX', 'chrY', 'chrM']:
        seq_dict[chrom] = get_ref_seq(path_to_ref, chrom)
    for chrom in seq_dict:
        print(f"{chrom}\t{len(seq_dict[chrom])}")
    with open(out_path, 'w') as f:
        for chrom in seq_dict:
            f.write(f"{chrom}\t{len(seq_dict[chrom])}\n")

    
        
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument(
        "--ref",
        "-r",
        type=str,
        required = True,
        help = "Reference genome to which alignments were made"
    )
    
    parser.add_argument(
        "--out",
        "-o",
        type=str,
        required = True,
        help = "R9 modbam2bed"
    )
    
    FLAGS, unparsed = parser.parse_known_args()

    main(FLAGS.ref, FLAGS.out)