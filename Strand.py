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

#Merge R9 and R10 datasets
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
                    usecols = [0, 1, 5, 9, 10, 11, 12]
                   )
    R9_df.columns = ['chrom', 'position', 'strand', 'coverage', 'ratio', 'modified', 'canonical']
    R9_df = R9_df[R9_df['chrom'].isin(chroms)]
        
    R10_df = pd.read_csv(R10_path,
                        sep='\t',
                        header = None,
                        usecols = [0, 1, 5, 9, 10, 11, 12]
                       )
    R10_df.columns = ['chrom', 'position', 'strand', 'coverage', 'ratio', 'modified', 'canonical']
    R10_df = R10_df[R10_df['chrom'].isin(chroms)]
    
    R10_df = R10_df[(R10_df['coverage'] >= 20) & (R10_df['coverage'] <= 200)]
    R10_df['ratio'] = R10_df['ratio'] / 100
    
    R9_df = R9_df[(R9_df['coverage'] >= 20) & (R9_df['coverage'] <= 200)]
    R9_df['ratio'] = R9_df['ratio'] / 100
    
    merged_df = pd.merge(R10_df, 
                         R9_df, 
                         how = 'inner', 
                         left_on = ['chrom', 'position', 'strand'], 
                         right_on=['chrom', 'position', 'strand'])    
    
    merged_df.columns = ['chrom', 
                         'position',
                         'strand',
                         'R10_coverage',
                         'R10_ratio', 
                         'R10_modified', 
                         'R10_canonical', 
                         'R9_coverage',
                         'R9_ratio', 
                         'R9_modified', 
                         'R9_canonical']
    return merged_df


#Bin data by strand and create violin plot of methylation proportions
def main(r9_hg002_path,
         r10_hg002_path):
    
        df = read_merge(r9_hg002_path, r10_hg002_path)
        df["delta_proportion"] = df["R10_ratio"] - df["R9_ratio"]
        
        strand = df['strand'].to_list()
        r10_bin = [round(x, 1) for x in df["R10_ratio"]]
        delta = df["delta_proportion"].to_list()
        
        strand_df = pd.DataFrame.from_dict({'strand':strand,
                                            'bin':r10_bin,
                                            'delta':delta})
        
        
        
        sns.set(rc={'figure.figsize':(11.7,8.27)})
        fig = plt.figure(dpi = 300)
        ax = fig.add_subplot(1, 1, 1)
        b = sns.violinplot(data=strand_df,
                           x='bin',
                           y='delta',
                           hue='strand', 
                           bw = 0.25)
        
        b.set_xlabel("R10 Methylation - 10% bins")
        b.set_ylabel("R10 Methylation Proportion - R9 Methylation Proportion")
        
        b.get_figure().savefig("strand_analysis.pdf")
        


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


    FLAGS, unparsed = parser.parse_known_args()
    main(FLAGS.r9_hg002, 
         FLAGS.r10_hg002)