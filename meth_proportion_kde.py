import pandas as pd
import numpy as np
import seaborn as sns
import time
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
from scipy.stats import pearsonr, describe
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator
plt.rcParams['pdf.fonttype'] = 42
plt.switch_backend('agg')
sns.set(rc={'figure.figsize':(11.7,80.27)})

#Plot the KDE of methylation proportions for each of the technologies
def main(r9, r10, pacbio, illumina, min_reads, out_dir, out_prefix):
    
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
    
    r9_df = pd.read_csv(r9, sep='\t', header=None, usecols=[0,11,12])
    r9_df.columns = ['chrom', 'mod', 'canon']
    r9_df['proportion'] = r9_df['mod'] / (r9_df['mod'] + r9_df['canon'])
    r9_df = r9_df[(r9_df['mod'] + r9_df['canon'] >= 20) & 
                  (r9_df['mod'] + r9_df['canon'] <= 200) &
                  (r9_df['chrom'].isin(chrom_index))]
    
    r10_df = pd.read_csv(r10, sep='\t', header=None, usecols=[0,11,12])
    r10_df.columns = ['chrom', 'mod', 'canon']
    r10_df['proportion'] = r10_df['mod'] / (r10_df['mod'] + r10_df['canon'])
    r10_df = r10_df[(r10_df['mod'] + r10_df['canon'] >= 20) & 
                    (r10_df['mod'] + r10_df['canon'] <= 200)&
                    (r10_df['chrom'].isin(chrom_index))]
    
    pacbio_df = pd.read_csv(pacbio, sep='\t', header = None, usecols = [0,11,12])
    pacbio_df.columns = ['chrom', 'mod', 'canon']
    pacbio_df['proportion'] = pacbio_df['mod'] / (pacbio_df['mod'] + pacbio_df['canon'])
    pacbio_df = pacbio_df[(pacbio_df['mod'] + pacbio_df['canon'] >= 20) & 
                          (pacbio_df['mod'] + pacbio_df['canon'] <= 200)&
                          (pacbio_df['chrom'].isin(chrom_index))]
    
    illumina_df = pd.read_csv(illumina, sep='\t', header = None, usecols = [0,4,5])
    illumina_df.columns = ['chrom', 'mod', 'canon']
    illumina_df['proportion'] = illumina_df['mod'] / (illumina_df['mod'] + illumina_df['canon'])
    illumina_df = illumina_df[(illumina_df['mod'] + illumina_df['canon'] >= 20) & 
                              (illumina_df['mod'] + illumina_df['canon'] <= 200)&
                              (illumina_df['chrom'].isin(chrom_index))]
    
    fig, ax = plt.subplots(1, 1)
    fig.set_figwidth(80)
    fig.set_figheight(15)
    sns.kdeplot(data = r9_df, x = 'proportion', color = 'g', ax = ax, label = 'R9', clip = [0,1], bw_adjust = 2)
    sns.kdeplot(data = r10_df, x = 'proportion', color = 'b', ax = ax, label = 'R10', clip = [0,1], bw_adjust = 2)
    sns.kdeplot(data = pacbio_df, x = 'proportion', color = 'm', ax = ax, label = 'Pacbio', clip = [0,1], bw_adjust = 2)
    sns.kdeplot(data = illumina_df, x = 'proportion', color = 'k', ax = ax, label = 'Illumina', clip = [0,1], bw_adjust = 2)
    
    ax.grid(False)
    fig.legend()
    fig.get_figure().savefig('meth_proportion_kdeplots.pdf')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument(
        "--r9_modkit",
        "-9",
        type=str,
        required = True,
        help = "R9 modkit"
    )
    
    parser.add_argument(
        "--r10_modkit",
        "-10",
        type=str,
        required = True,
        help = "R10 modkit"
    )
    
    parser.add_argument(
        "--pacbio_modkit",
        "-pb",
        type=str,
        required = True,
        help = "pacbio modkit"
    )
    
    parser.add_argument(
        "--illumina_5mC",
        "-i",
        type=str,
        required = True,
        help = "Illumina in Bismark format"
    )
    
    parser.add_argument(
        "--min_reads",
        "-m",
        type = int,
        default = 10,
        help = "The minimum number of reads from both R9 and R10 to consider a loci"
    )
    
    parser.add_argument(
        "--output_directory",
        "-od",
        type = str,
        required = True,
        help = "Directory where plots will be saved"
    )
    
    parser.add_argument(
        "--output_prefix",
        "-op",
        type = str,
        help = "Prefix for output figures",
        default = "out"
    )
    

    FLAGS, unparsed = parser.parse_known_args()
    main(FLAGS.r9_modkit, 
         FLAGS.r10_modkit,
         FLAGS.pacbio_modkit,
         FLAGS.illumina_5mC,
         FLAGS.min_reads,
         FLAGS.output_directory, 
         FLAGS.output_prefix)  
