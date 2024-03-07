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
sns.set(rc={'figure.figsize':(11.7,8.27)})

#Get lengths of each chrom for data structure initialization
def chrom_lengths_from_txt(txt_path):
    chrom_lengths = {}
    with open(txt_path, 'r') as f:
        for line in f:
            chrom_lengths[line.split('\t')[0]] = int(line.split('\t')[1])
    return chrom_lengths

#Lookup dict to translate chrom to index for numpy array
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

def main(r9_modkit, 
         r10_modkit,
         pacbio_modkit,
         illumina_bismark,
         chrom_lengths_path,
         min_reads,
         output_dir,
         output_prefix
        ):

    #Output stats fh
    outfile_fh = open(os.path.join(output_dir, f"{output_prefix}_comparisson_stats.txt"), 'w')
    
    
    #Data structure will be a numpy array #chroms x 8 x length of each chromosome
    list_np_methylation = []
    chrom_lengths = chrom_lengths_from_txt(chrom_lengths_path)
    
    print("Setting up data structure...")
    #Fill each cell w/ 0
    for chrom in chrom_index:
        list_np_methylation.append(np.full((8, chrom_lengths[chrom]), 0.0))

    r9_cpgs = 0
    r10_cpgs = 0
    pacbio_cpgs = 0
    illumina_cpgs = 0
        
    #Read in methylation information
    print("Calculating R9 Methylation")
    with open(r9_modkit, 'r') as f_R9:
        for line in f_R9:
            if line.split('\t')[0] not in chrom_index:
                continue
            c = chrom_index[line.split('\t')[0]]
            index = int(line.split('\t')[1])
            canon = int(line.split('\t')[12])
            mod = int(line.split('\t')[11])
            
            if canon + mod >= min_reads and canon + mod <= 200:
                r9_cpgs += 1
            
            if canon + mod == 0:
                continue
            else:
                list_np_methylation[c][0][index] = mod / (mod + canon)
                list_np_methylation[c][1][index] = mod + canon
    outfile_fh.write(f"R9 total cpg w/ {min_reads}+ coverage: {r9_cpgs}\n")

    print("Calculating R10 Methylation")
    with open(r10_modkit, 'r') as f_R10:
        for line in f_R10:
            if line.split('\t')[0] not in chrom_index:
                continue
            c = chrom_index[line.split('\t')[0]]
            index = int(line.split('\t')[1])
            canon = int(line.split('\t')[12])
            mod = int(line.split('\t')[11])
            
            if canon + mod >= min_reads and canon + mod <= 200:
                r10_cpgs += 1
                
            if canon + mod == 0:
                continue
            else:
                list_np_methylation[c][2][index] = mod / (mod + canon)
                list_np_methylation[c][3][index] = mod + canon
    outfile_fh.write(f"R10 total cpg w/ {min_reads}+ coverage: {r10_cpgs}\n")
    
    print("Calculating PacBio Methylation")
    with open(pacbio_modkit, 'r') as f_pacbio:
        for line in f_pacbio:
            if line.split('\t')[0] not in chrom_index:
                continue
            c = chrom_index[line.split('\t')[0]]
            index = int(line.split('\t')[1])
            canon = int(line.split('\t')[12])
            mod = int(line.split('\t')[11])

            if canon + mod >= min_reads and canon + mod <= 200:
                pacbio_cpgs += 1
                
            if canon + mod == 0:
                continue
            else:
                list_np_methylation[c][4][index] = mod / (mod + canon)
                list_np_methylation[c][5][index] = mod + canon
    outfile_fh.write(f"Pacbio total cpg w/ {min_reads}+ coverage: {pacbio_cpgs}\n")
    
    print("Calculating Illumina Methylation")
    with open(illumina_bismark, 'r') as f_illumina:
        for line in f_illumina:
            if line.split('\t')[0] not in chrom_index:
                continue
            c = chrom_index[line.split('\t')[0]]
            index = int(line.split('\t')[1])
            canon = int(line.split('\t')[5])
            mod = int(line.split('\t')[4])

            if canon + mod >= min_reads and canon + mod <= 200:
                illumina_cpgs += 1

            if canon + mod == 0:
                continue
            else:
                list_np_methylation[c][6][index] = mod / (mod + canon)
                list_np_methylation[c][7][index] = mod + canon    
    outfile_fh.write(f"Illumina total cpg w/ {min_reads}+ coverage: {illumina_cpgs}\n")
       
                
    print("Calculating Proportions R9 / R10")
    #create lists of proportion information
    r9_ratio = []
    r10_ratio = []
    for chrom in chrom_index:
        c = chrom_index[chrom]
        for i in range(0, len(list_np_methylation[c][0])):
            if (list_np_methylation[c][1][i] >= min_reads and 
                list_np_methylation[c][3][i] >= min_reads):
                r9_ratio.append(list_np_methylation[c][0][i])
                r10_ratio.append(list_np_methylation[c][2][i])
                
    r = pearsonr(r9_ratio, r10_ratio)
    outfile_fh.write(f"R9 / R10 Overlap cpg: {len(r9_ratio)}\n")
    outfile_fh.write(f"R9 / R10 pearson r: {round(r[0], 6)}\n")
    
    #Make scatter plot of R9 ratio vs R10 ratio
    df = pd.DataFrame.from_dict({'r10':r10_ratio, 'r9': r9_ratio})
    fig = plt.figure(dpi = 300)
    print("Making R9 vs. R10 Figures")
    ax = fig.add_subplot(1, 1, 1)
    ax.set_facecolor('black')
    density = ax.hist2d(x = df['r10'].to_list(), 
                        y = df['r9'].to_list(), 
                        cmap='viridis', 
                        #norm=colors.LogNorm(), 
                        bins=20)
    fig.colorbar(density[3], label='Number of points per square')
    plt.plot([0, 1],[0, 1], 'r-', alpha=1.0, zorder=10, linewidth=3)
    ax.grid(False)
    plt.xlabel('R10 Methylation Proportion')
    plt.ylabel('R9 Methylation Proportion')
    plt.title('R9 vs R10 Methylation Proportions')
    ax.get_figure().savefig(os.path.join(f'{output_dir}',f'{output_prefix}_R9_R10_methylation.pdf'), dpi = 300)
    plt.clf()

    #Create a measure of delta ratio by subtracting r9 ratio from r10 at every position
    print("Creating R9 / R10 Violin Plot")
    delta_ratio = [0] * len(df['r10'].to_list())
    r10_bin = [0] * len(df['r10'].to_list())
    for i in range(len(delta_ratio)):
        delta_ratio[i] = r10_ratio[i] - r9_ratio[i]
        #create bins for plotting (R10 rounded to nearst tenth)
        r10_bin[i] = round(r10_ratio[i], 1)
    delta_df = pd.DataFrame.from_dict({'delta':delta_ratio, 'r10_bin':r10_bin})
    
    plt.clf()
    sns.set(rc={"figure.figsize":(20, 6)}) #width=3, #height=4
    fig, ax = plt.subplots()
    b = sns.violinplot(data=delta_df, 
                       x = 'r10_bin', 
                       y = 'delta', 
                       ax=ax)
    plt.xlabel('R10 Proportion Methylation')
    plt.ylabel('R10 Proportion - R9 Proportion')
    plt.title('Delta Proportion of Methylation (R10 - R9) binned by R10 Proportion (nearest 0.1)')
    plt.ylim([-1.0, 1.0])
    b.axhline(0.0, color = 'r')
    b.get_figure().savefig(os.path.join(f'{output_dir}', f'{output_prefix}_Violin_Delta.R10.R9.pdf'), dpi = 300)
    
    #Create ratio and dataframes for r9 and pacbio
    print("Calculating Proportions R9 / Pacbio")
    r9_ratio = []
    pacbio_ratio = []
    for chrom in chrom_index:
        c = chrom_index[chrom]
        for i in range(0, len(list_np_methylation[c][0])):
            if (list_np_methylation[c][1][i] >= min_reads and 
                list_np_methylation[c][5][i] >= min_reads):
                r9_ratio.append(list_np_methylation[c][0][i])
                pacbio_ratio.append(list_np_methylation[c][4][i])

    r = pearsonr(r9_ratio, pacbio_ratio)
    outfile_fh.write(f"R9 / Pacbio Overlap cpg: {len(r9_ratio)}\n")
    outfile_fh.write(f"R9 / Pacbio pearson r: {round(r[0], 6)}\n")
    
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    df = pd.DataFrame.from_dict({'r9' : r9_ratio, 'pacbio' : pacbio_ratio})
    #Make scatter plot of R9 ratio vs Pacbio ratio
    fig = plt.figure(dpi = 300)
    print("Making R9 vs. Pacbio Figures")
    ax = fig.add_subplot(1, 1, 1)
    ax.set_facecolor('black')
    density = ax.hist2d(x = df['r9'].to_list(), 
                        y = df['pacbio'].to_list(), 
                        cmap='viridis', 
                        norm=colors.LogNorm(), 
                        bins=20)
    # now plot both limits against eachother
    fig.colorbar(density[3], label='Number of points per square')
    plt.plot([0, 1],[0, 1], 'r-', alpha=1.0, zorder=10, linewidth=3)
    ax.grid(False)
    plt.ylabel('Pacbio Methylation Proportion')
    plt.xlabel('R9 Methylation Proportion')
    plt.title('Pacbio vs R9 Methylation Proportions')
    ax.get_figure().savefig(os.path.join(f'{output_dir}',f'{output_prefix}_R9_Pacbio_methylation.pdf'), dpi = 300)
    plt.clf()

    #Create a measure of delta ratio by subtracting pacbio ratio from r9 at every position
    print("Creating R9 / Pacbio Violin Plots")
    delta_ratio = [0] * len(df['r9'].to_list())
    r9_bin = [0] * len(df['r9'].to_list())
    for i in range(len(delta_ratio)):
        delta_ratio[i] = r9_ratio[i] - pacbio_ratio[i]
        #create bins for plotting (R10 rounded to nearst tenth)
        r9_bin[i] = round(r9_ratio[i], 1)
    delta_df = pd.DataFrame.from_dict({'delta':delta_ratio, 'r9_bin':r9_bin})
    
    plt.clf()
    sns.set(rc={"figure.figsize":(20, 6)}) #width=3, #height=4
    fig, ax = plt.subplots()
    b = sns.violinplot(data=delta_df, 
                       x = 'r9_bin', 
                       y = 'delta', 
                       ax=ax)
    plt.xlabel('R9 Proportion Methylation')
    plt.ylabel('R9 Proportion - Pacbio Proportion')
    plt.title('Delta Proportion of Methylation (R9 - Pacbio) binned by R9 Proportion (nearest 0.1)')
    plt.ylim([-1.0, 1.0])
    b.axhline(0.0, color = 'r')
    b.get_figure().savefig(os.path.join(f'{output_dir}', f'{output_prefix}_Violin_Delta.R9.pacbio.pdf'), dpi = 300)
    
    #Make scatter plot of R10 ratio vs Pacbio ratio
    print("Calculating Proportions R10 / Pacbio")
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    r10_ratio = []
    pacbio_ratio = []
    for chrom in chrom_index:
        c = chrom_index[chrom]
        for i in range(0, len(list_np_methylation[c][0])):
            if (list_np_methylation[c][3][i] >= min_reads and 
                list_np_methylation[c][5][i] >= min_reads):
                r10_ratio.append(list_np_methylation[c][2][i])
                pacbio_ratio.append(list_np_methylation[c][4][i])
    
    r = pearsonr(r10_ratio, pacbio_ratio)
    outfile_fh.write(f"R10 / Pacbio Overlap cpg: {len(r10_ratio)}\n")
    outfile_fh.write(f"R10 / Pacbio pearson r: {round(r[0], 6)}\n")
    
    df = pd.DataFrame.from_dict({'r10' : r10_ratio, 'pacbio' : pacbio_ratio})
           
    fig = plt.figure(dpi = 300)
    print("Making R10 vs. Pacbio Figures")
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    ax = fig.add_subplot(1, 1, 1)
    ax.set_facecolor('black')
    density = ax.hist2d(x = df['r10'].to_list(), 
                        y = df['pacbio'].to_list(), 
                        cmap='viridis', 
                        norm=colors.LogNorm(), 
                        bins=20)
    # now plot both limits against eachother
    fig.colorbar(density[3], label='Number of points per square')
    plt.plot([0, 1],[0, 1], 'r-', alpha=1.0, zorder=10, linewidth=3)
    ax.grid(False)
    plt.xlabel('R10 Methylation Proportion')
    plt.ylabel('Pacbio Methylation Proportion')
    plt.title('Pacbio vs R10 Methylation Proportions')
    ax.get_figure().savefig(os.path.join(f'{output_dir}',f'{output_prefix}_R10_Pacbio_methylation.pdf'), dpi = 300)     

    #Create a measure of delta ratio by subtracting pacbio ratio from r10 at every position
    print("Creating R10 / Pacbio Violin Plot")
    delta_ratio = [0] * len(df['r10'].to_list())
    r10_bin = [0] * len(df['r10'].to_list())
    for i in range(len(delta_ratio)):
        delta_ratio[i] = r10_ratio[i] - pacbio_ratio[i]
        #create bins for plotting (R10 rounded to nearst tenth)
        r10_bin[i] = round(r10_ratio[i], 1)
    delta_df = pd.DataFrame.from_dict({'delta':delta_ratio, 'r10_bin':r10_bin})
    
    plt.clf()
    sns.set(rc={"figure.figsize":(20, 6)}) #width=3, #height=4
    fig, ax = plt.subplots()
    b = sns.violinplot(data=delta_df, 
                       x = 'r10_bin', 
                       y = 'delta', 
                       ax=ax)
    plt.xlabel('R10 Proportion Methylation')
    plt.ylabel('R10 Proportion - Pacbio Proportion')
    plt.title('Delta Proportion of Methylation (R10 - Pacbio) binned by R10 Proportion (nearest 0.1)')
    plt.ylim([-1.0, 1.0])
    b.axhline(0.0, color = 'r')
    b.get_figure().savefig(os.path.join(f'{output_dir}', f'{output_prefix}_Violin_Delta.R10.pacbio.pdf'), dpi = 300)

    #Make scatter plot of R10 ratio vs Bisulfite ratio
    print("Calculating Proportions R10 / Bisulfite")
    r10_ratio = []
    bisulfite_ratio = []
    for chrom in chrom_index:
        c = chrom_index[chrom]
        for i in range(0, len(list_np_methylation[c][0])):
            if (list_np_methylation[c][3][i] >= min_reads and 
                list_np_methylation[c][7][i] >= min_reads):
                r10_ratio.append(list_np_methylation[c][2][i])
                bisulfite_ratio.append(list_np_methylation[c][6][i])

    r = pearsonr(r10_ratio, bisulfite_ratio)
    outfile_fh.write(f"R10 / Bisulfite Overlap cpg: {len(r10_ratio)}\n")
    outfile_fh.write(f"R10 / Bisulfite pearson r: {round(r[0], 6)}\n")
                
    df = pd.DataFrame.from_dict({'r10' : r10_ratio, 'bisulfite' : bisulfite_ratio})
    
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    fig = plt.figure(dpi = 300)
    print("Making R10 vs. Bisulfite Figures")
    ax = fig.add_subplot(1, 1, 1)
    ax.set_facecolor('black')
    density = ax.hist2d(x = df['r10'].to_list(), 
                        y = df['bisulfite'].to_list(), 
                        cmap='viridis', 
                        norm=colors.LogNorm(), 
                        bins=20)
    # now plot both limits against eachother
    fig.colorbar(density[3], label='Number of points per square')
    plt.plot([0, 1],[0, 1], 'r-', alpha=1.0, zorder=10, linewidth=3)
    ax.grid(False)
    plt.xlabel('R10 Methylation Proportion')
    plt.ylabel('Bisulfite Methylation Proportion')
    plt.title('Bisulfite vs R10 Methylation Proportions')
    ax.get_figure().savefig(os.path.join(f'{output_dir}',f'{output_prefix}_R10_Bisulfite_methylation.pdf'), dpi = 300)     

    #Create a measure of delta ratio by subtracting pacbio ratio from r10 at every position
    print("Creating R10 / Bisulfite Violin Plots")
    delta_ratio = [0] * len(df['r10'].to_list())
    r10_bin = [0] * len(df['r10'].to_list())
    for i in range(len(delta_ratio)):
        delta_ratio[i] = r10_ratio[i] - bisulfite_ratio[i]
        #create bins for plotting (R10 rounded to nearst tenth)
        r10_bin[i] = round(r10_ratio[i], 1)
    delta_df = pd.DataFrame.from_dict({'delta':delta_ratio, 'r10_bin':r10_bin})
    
    plt.clf()
    sns.set(rc={"figure.figsize":(20, 6)}) #width=3, #height=4
    fig, ax = plt.subplots()
    b = sns.violinplot(data=delta_df, 
                       x = 'r10_bin', 
                       y = 'delta', 
                       ax=ax)
    plt.xlabel('R10 Proportion Methylation')
    plt.ylabel('R10 Proportion - Bisulfite Proportion')
    plt.title('Delta Proportion of Methylation (R10 - Bisulfite) binned by R10 Proportion (nearest 0.1)')
    plt.ylim([-1.0, 1.0])
    b.axhline(0.0, color = 'r')
    b.get_figure().savefig(os.path.join(f'{output_dir}', f'{output_prefix}_Violin_Delta.R10.Bisulfite.pdf'), dpi = 300)
    
    #R9 vs. Bisulfite Plots
    print("Calculating Proportions R9 / Bisulfite")
    r9_ratio = []
    bisulfite_ratio = []
    for chrom in chrom_index:
        c = chrom_index[chrom]
        for i in range(0, len(list_np_methylation[c][0])):
            if (list_np_methylation[c][1][i] >= min_reads and 
                list_np_methylation[c][7][i] >= min_reads):
                r9_ratio.append(list_np_methylation[c][0][i])
                bisulfite_ratio.append(list_np_methylation[c][6][i])

    r = pearsonr(r9_ratio, bisulfite_ratio)
    outfile_fh.write(f"R9 / Bisulfite Overlap cpg: {len(r9_ratio)}\n")
    outfile_fh.write(f"R9 / Bisulfite pearson r: {round(r[0], 6)}\n")
    
    df = pd.DataFrame.from_dict({'r9' : r9_ratio, 'bisulfite' : bisulfite_ratio})
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    fig = plt.figure(dpi = 300)
    print("Making R9 vs. Bisulfite Figures")
    ax = fig.add_subplot(1, 1, 1)
    ax.set_facecolor('black')
    density = ax.hist2d(x = df['r9'].to_list(), 
                        y = df['bisulfite'].to_list(), 
                        cmap='viridis', 
                        norm=colors.LogNorm(), 
                        bins=20)
    # now plot both limits against eachother
    fig.colorbar(density[3], label='Number of points per square')
    plt.plot([0, 1],[0, 1], 'r-', alpha=1.0, zorder=10, linewidth=3)
    ax.grid(False)
    plt.xlabel('R9 Methylation Proportion')
    plt.ylabel('Bisulfite Methylation Proportion')
    plt.title('Bisulfite vs R9 Methylation Proportions')
    ax.get_figure().savefig(os.path.join(f'{output_dir}',f'{output_prefix}_R9_Bisulfite_methylation.pdf'), dpi = 300)     

    #Create a measure of delta ratio by subtracting pacbio ratio from r10 at every position
    print("Creating R9 / Bisulfite Violin Plots")
    delta_ratio = [0] * len(df['r9'].to_list())
    r9_bin = [0] * len(df['r9'].to_list())
    for i in range(len(delta_ratio)):
        delta_ratio[i] = r9_ratio[i] - bisulfite_ratio[i]
        #create bins for plotting (R9 rounded to nearst tenth)
        r9_bin[i] = round(r9_ratio[i], 1)
    delta_df = pd.DataFrame.from_dict({'delta':delta_ratio, 'r9_bin':r9_bin})
    
    plt.clf()
    sns.set(rc={"figure.figsize":(20, 6)}) #width=3, #height=4
    fig, ax = plt.subplots()
    b = sns.violinplot(data=delta_df, 
                       x = 'r9_bin', 
                       y = 'delta', 
                       ax=ax)
    plt.xlabel('R9 Proportion Methylation')
    plt.ylabel('R9 Proportion - Bisulfite Proportion')
    plt.title('Delta Proportion of Methylation (R9 - Bisulfite) binned by R9 Proportion (nearest 0.1)')
    plt.ylim([-1.0, 1.0])
    b.axhline(0.0, color = 'r')
    b.get_figure().savefig(os.path.join(f'{output_dir}', f'{output_prefix}_Violin_Delta.R9.Bisulfite.pdf'), dpi = 300)
    
    #Pacbio vs. Bisulfite Plots
    print("Calculating Proportions Pacbio / Bisulfite")
    pacbio_ratio = []
    bisulfite_ratio = []
    for chrom in chrom_index:
        c = chrom_index[chrom]
        for i in range(0, len(list_np_methylation[c][0])):
            if (list_np_methylation[c][5][i] >= min_reads and 
                list_np_methylation[c][7][i] >= min_reads):
                pacbio_ratio.append(list_np_methylation[c][4][i])
                bisulfite_ratio.append(list_np_methylation[c][6][i])
                
    r = pearsonr(pacbio_ratio, bisulfite_ratio)
    outfile_fh.write(f"Pacbio / Bisulfite Overlap cpg: {len(pacbio_ratio)}\n")
    outfile_fh.write(f"Pacbio / Bisulfite pearson r: {round(r[0], 6)}\n")
    
    df = pd.DataFrame.from_dict({'pacbio' : pacbio_ratio, 'bisulfite' : bisulfite_ratio})  
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    fig = plt.figure(dpi = 300)
    print("Making R9 vs. Bisulfite Figures")
    ax = fig.add_subplot(1, 1, 1)
    ax.set_facecolor('black')
    density = ax.hist2d(x = df['pacbio'].to_list(), 
                        y = df['bisulfite'].to_list(), 
                        cmap='viridis', 
                        norm=colors.LogNorm(), 
                        bins=20)
    # now plot both limits against eachother
    fig.colorbar(density[3], label='Number of points per square')
    plt.plot([0, 1],[0, 1], 'r-', alpha=1.0, zorder=10, linewidth=3)
    ax.grid(False)
    plt.xlabel('Pacbio Methylation Proportion')
    plt.ylabel('Bisulfite Methylation Proportion')
    plt.title('Bisulfite vs Pacbio Methylation Proportions')
    ax.get_figure().savefig(os.path.join(f'{output_dir}',f'{output_prefix}_Pacbio_Bisulfite_methylation.pdf'), dpi = 300)     

    #Create a measure of delta ratio by subtracting pacbio ratio from r10 at every position
    print("Creating Pacbio / Bisulfite Violin Plots")
    delta_ratio = [0] * len(df['pacbio'].to_list())
    pacbio_bin = [0] * len(df['pacbio'].to_list())
    for i in range(len(delta_ratio)):
        delta_ratio[i] = pacbio_ratio[i] - bisulfite_ratio[i]
        #create bins for plotting (R9 rounded to nearst tenth)
        pacbio_bin[i] = round(pacbio_ratio[i], 1)
    delta_df = pd.DataFrame.from_dict({'delta':delta_ratio, 'pacbio_bin':pacbio_bin})   
    plt.clf()
    sns.set(rc={"figure.figsize":(20, 6)}) #width=3, #height=4
    fig, ax = plt.subplots()
    b = sns.violinplot(data=delta_df, 
                       x = 'pacbio_bin', 
                       y = 'delta', 
                       ax=ax)
    plt.xlabel('Pacbio Proportion Methylation')
    plt.ylabel('Pacbio Proportion - Bisulfite Proportion')
    plt.title('Delta Proportion of Methylation (Pacbio - Bisulfite) binned by Pacbio Proportion (nearest 0.1)')
    plt.ylim([-1.0, 1.0])
    b.axhline(0.0, color = 'r')
    b.get_figure().savefig(os.path.join(f'{output_dir}', f'{output_prefix}_Violin_Delta.Pacbio.Bisulfite.pdf'), dpi = 300)
       
    outfile_fh.close()
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument(
        "--chrom_lengths",
        "-c",
        type=str,
        required = True,
        help = "Length of reference chromosomes to which mod bam reads were aligned"
    )
    
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
        help = "Prefix for output figures"
    )
    

    FLAGS, unparsed = parser.parse_known_args()
    main(FLAGS.r9_modkit, 
         FLAGS.r10_modkit,
         FLAGS.pacbio_modkit,
         FLAGS.illumina_5mC,
         FLAGS.chrom_lengths,
         FLAGS.min_reads,
         FLAGS.output_directory, 
         FLAGS.output_prefix)  
