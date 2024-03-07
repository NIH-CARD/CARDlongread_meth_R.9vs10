import pandas as pd
import numpy as np
import sys
from scipy.stats import mannwhitneyu
import time
import argparse
import os


#Print out the number of positions from input dataframes where the extremes (0% and 100%) are shared and where they are unique
def compare_proportions(df, alabel, blabel):
    print(f"{alabel} == 0.0, {blabel} == 0.0: {df[(df[alabel] == 0.0) & (df[blabel] == 0.0)].shape[0]}")
    print(f"{alabel} == 0.0, {blabel} > 0.0: {df[(df[alabel] == 0.0) & (df[blabel] > 0.0)].shape[0]}")
    print(f"{alabel} > 0.0, {alabel} < 1.0, {blabel} > 0.0, {blabel} < 1.0: {df[(df[alabel] > 0.0) & (df[alabel] < 1.0) & (df[blabel] > 0.0) & (df[blabel] < 1.0)].shape[0]}")
    print(f"{alabel} > 0.0, {alabel} < 1.0, {blabel} == 0.0 | {blabel} == 1.0: {df[(df[alabel] > 0.0) & (df[alabel] < 1.0) & ((df[blabel] == 0.0) | (df[blabel] == 1.0))].shape[0]}")
    print(f"{alabel} == 1.0, {blabel} == 1.0: {df[(df[alabel] == 1.0) & (df[blabel] == 1.0)].shape[0]}")
    print(f"{alabel} == 1.0, {blabel} < 1.0: {df[(df[alabel] == 1.0) & (df[blabel] < 1.0)].shape[0]}")


#Read in R9 and R10 datasets for comparisson
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
    
    print(f"R9 Pre-filter: {R9_df.shape}")
    
    R10_df = pd.read_csv(R10_path,
                        sep='\t',
                        header = None,
                        usecols = [0, 1, 9, 10, 11, 12]
                       )
    R10_df.columns = ['chrom', 'position', 'coverage', 'ratio', 'modified', 'canonical']
    R10_df = R10_df[R10_df['chrom'].isin(chroms)]
    
    print(f"R10 Pre-filter: {R10_df.shape}")
    
    R10_df = R10_df[(R10_df['coverage'] >= 20) & (R10_df['coverage'] <= 200)]
    R10_df['ratio'] = R10_df['ratio'] / 100
    print(f"R10 has {R10_df.shape} sites")
    
    R9_df = R9_df[(R9_df['coverage'] >= 20) & (R9_df['coverage'] <= 200)]
    R9_df['ratio'] = R9_df['ratio'] / 100
    print(f"R9 has {R9_df.shape} sites")
    
    merged_df = pd.merge(R10_df, 
                         R9_df, 
                         how = 'inner', 
                         left_on = ['chrom', 'position'], 
                         right_on=['chrom', 'position'])
    
    print(f"R9 and R10 share {merged_df.shape} sites")
    
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

def merge_with_illumina(illumina, merged_df):
    
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
    
    illumina_df = pd.read_csv(illumina, sep='\t', header = None, usecols = [0,1,4,5])
    illumina_df.columns = ['chrom', 'position', 'mod', 'canon']
    illumina_df['proportion'] = illumina_df['mod'] / (illumina_df['mod'] + illumina_df['canon'])
    illumina_df = illumina_df[(illumina_df['mod'] + illumina_df['canon'] >= 20) & 
                              (illumina_df['mod'] + illumina_df['canon'] <= 200) &
                              (illumina_df['chrom'].isin(chroms))]
    
    new_merged_df = pd.merge(merged_df,
                             illumina_df,
                             how='inner',
                             left_on = ['chrom', 'position'],
                             right_on = ['chrom', 'position'])
    new_merged_df.columns = ['chrom', 
                             'position', 
                             'R10_coverage',
                             'R10_ratio', 
                             'R10_modified', 
                             'R10_canonical', 
                             'R9_coverage',
                             'R9_ratio', 
                             'R9_modified', 
                             'R9_canonical',
                             'illumina_modified',
                             'illumina_canonical',
                             'illumina_ratio']
    return new_merged_df

def main(R9_hg002_path, 
         R10_hg002_path, 
         R9_blood_path, 
         R10_blood_path, 
         R9_brain_path, 
         R10_brain_path, 
         illumina_path):
    
    #HG002
    merged_df = read_merge(R9_hg002_path, R10_hg002_path)

    print(f"\nHG002 Proportion Comparisons R9 to R10")
    compare_proportions(merged_df, 'R9_ratio', 'R10_ratio')
    print(f"\nHG002 Proportion Comparisons R10 to R9")
    compare_proportions(merged_df, 'R10_ratio', 'R9_ratio')
    
    print("\nMann-Whitney test for R9 >= 0.9, testing if R10 > R9")
    print(mannwhitneyu(merged_df[merged_df['R9_ratio'] > 0.9]['R10_ratio'], merged_df[merged_df['R9_ratio'] > 0.9]['R9_ratio'], alternative='greater'))
    print("\nMann-Whitney test for R9 <= 0.9, testing if R10 < R9")
    print(mannwhitneyu(merged_df[merged_df['R9_ratio'] < 0.1]['R10_ratio'], merged_df[merged_df['R9_ratio'] < 0.1]['R9_ratio'], alternative='less'))

    #Illumina binning
    new_merged_df = merge_with_illumina(illumina_path, merged_df)
    print("\nMann-Whitney test for illumina >= 0.9, testing if R10 > R9")
    print(mannwhitneyu(new_merged_df[new_merged_df['illumina_ratio'] > 0.9]['R10_ratio'], new_merged_df[new_merged_df['illumina_ratio'] > 0.9]['R9_ratio'], alternative='greater'))
    print("\nMann-Whitney test for illumina <= 0.1, testing if R10 < R9")
    print(mannwhitneyu(new_merged_df[new_merged_df['illumina_ratio'] < 0.1]['R10_ratio'], new_merged_df[new_merged_df['illumina_ratio'] < 0.1]['R9_ratio'], alternative='less'))    


    
    #Blood
    merged_df = read_merge(R9_blood_path, R10_blood_path)
    print(f"\nBlood Proportion Comparisons R9 to R10")
    compare_proportions(merged_df, 'R9_ratio', 'R10_ratio')
    print(f"\nBlood Proportion Comparisons R10 to R9")
    compare_proportions(merged_df, 'R10_ratio', 'R9_ratio')
    
    print("\nMann-Whitney test for R9 >= 0.9, testing if R10 > R9")
    print(mannwhitneyu(merged_df[merged_df['R9_ratio'] > 0.9]['R10_ratio'], merged_df[merged_df['R9_ratio'] > 0.9]['R9_ratio'], alternative='greater'))
    print("\nMann-Whitney test for R9 <= 0.1, testing if R9 > R10")
    print(mannwhitneyu(merged_df[merged_df['R9_ratio'] < 0.1]['R10_ratio'], merged_df[merged_df['R9_ratio'] < 0.1]['R9_ratio'], alternative='less'))
    
    #Brain
    merged_df = read_merge(R9_brain_path, R10_brain_path)
    print(f"\nBrain Proportion Comparisons R9 to R10")
    compare_proportions(merged_df, 'R9_ratio', 'R10_ratio')
    print(f"\nBrain Proportion Comparisons R10 to R9")
    compare_proportions(merged_df, 'R10_ratio', 'R9_ratio')
    
    print("\nMann-Whitney test for R9 >= 0.9, testing if R10 > R9")
    print(mannwhitneyu(merged_df[merged_df['R9_ratio'] > 0.9]['R10_ratio'], merged_df[merged_df['R9_ratio'] > 0.9]['R9_ratio'], alternative='greater'))
    print("\nMann-Whitney test for R9 >= 0.9, testing if R10 > R9")
    print(mannwhitneyu(merged_df[merged_df['R9_ratio'] < 0.1]['R10_ratio'], merged_df[merged_df['R9_ratio'] < 0.1]['R9_ratio'], alternative='less'))
    
    
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
        "--r9_blood",
        type=str,
        required = True,
        help = "R9 Blood modkit"
    )
    
    parser.add_argument(
        "--r10_blood",
        type=str,
        required = True,
        help = "R10 Blood modkit"
    )
    parser.add_argument(
        "--r9_brain",
        type=str,
        required = True,
        help = "R9 HG002 modkit"
    )
    
    parser.add_argument(
        "--r10_brain",
        type=str,
        required = True,
        help = "R10 Blood modkit"
    )
    
    parser.add_argument(
        "--illumina",
        type=str,
        required = True,
        help = "Illumina data to use as binning for hg002 Mann Whitney Test"
    )
    

    FLAGS, unparsed = parser.parse_known_args()
    main(FLAGS.r9_hg002, 
         FLAGS.r10_hg002,
         FLAGS.r9_blood,
         FLAGS.r10_blood,
         FLAGS.r9_brain,
         FLAGS.r10_brain, 
         FLAGS.illumina)  
