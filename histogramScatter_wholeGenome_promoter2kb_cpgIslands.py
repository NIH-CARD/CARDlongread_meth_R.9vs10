#!/usr/bin/env python3


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
import matplotlib.gridspec as gridspec
import pyranges as pr
import datetime
from multiprocessing import Pool


###############################################
### whole genome strand combined:           ###
###############################################
#  
modkit_names = ["chr",'start','end','modbase','score_Nvalid', 'strand', 'start2','end2','color',
                'Nvaid_cov','fracMod','Nmod','Ncanon','Nother_mod','Ndelete','Nfail','Ndiff','Nnocall']

r9_wg = pd.read_csv('HG002_R9_CARD_bam_GRCh38.strandCombined.bed', sep="\t", header=None,
                   names = modkit_names)

r10_wg = pd.read_csv('HG002_R10_GRCh38.strandCombined.bed', sep="\t", header=None,
                   names = modkit_names)

r9_r10_wg = pd.merge(r9_wg, r10_wg, on=['chr','start','end'], suffixes=('_r9', '_r10'))

r9_r10_wg['diff']=r9_r10_wg['fracMod_r10'].astype(float)-r9_r10_wg['fracMod_r9'].astype(float)


###############################################
### Whole genome 2kb Bins:                  ###
############################################### 

chrom_df_list = []
for chromosome in ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
       'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
       'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22',
       'chrX', 'chrY', 'chrM']:
    chrom_df_list.append(r9_r10_wg.loc[r9_r10_wg['chr'] == chromosome ] )
    
# Split the genome into 2kb windows and smooth the methylation rates within
print("starting date and time:", datetime.datetime.now())
# Create a Pool of processes
with Pool() as pool:
    # Map the task function to the Pool of processes
    results = pool.map(do_bin_and_smoothing_r9_r10, chrom_df_list)
    
print("ending date and time:", datetime.datetime.now())
# Concatenate the list of DataFrames
smoothed_r9_r10_wg_2kb_70cpg = pd.concat(results, ignore_index=True)



# Add fields and remove some 
def subtract_elements(lst):
    return lst[1] - lst[0]

smoothed_r9_r10_wg_2kb_70cpg['genomic_size'] = smoothed_r9_r10_wg_2kb_70cpg['end'] - smoothed_r9_r10_wg_2kb_70cpg['start']
smoothed_r9_r10_wg_2kb_70cpg['Ncpgs'] = smoothed_r9_r10_wg_2kb_70cpg['ungrouped_idx'].apply(subtract_elements)
smoothed_r9_r10_wg_2kb_70cpg['cpgDensity'] = smoothed_r9_r10_wg_2kb_70cpg.loc[smoothed_r9_r10_wg_2kb_70cpg['Ncpgs']<1000]['Ncpgs']/smoothed_r9_r10_wg_2kb_70cpg.loc[smoothed_r9_r10_wg_2kb_70cpg['Ncpgs']<1000]['genomic_size']

# 2Kb whole genome binned correlation
pearson_corr = smoothed_r9_r10_wg_2kb_70cpg['methylFq_r9'].corr(smoothed_r9_r10_wg_2kb_70cpg['methylFq_r10'])

print(smoothed_r9_r10_wg_2kb_70cpg.shape[0], 'sites',"\nPearson correlation coefficient:", pearson_corr)

print('RMSE:',np.sqrt(mean_squared_error(smoothed_r9_r10_wg_2kb_70cpg['methylFq_r10'], smoothed_r9_r10_wg_2kb_70cpg['methylFq_r9'])))


smoothed_r9_r10_wg_2kb_70cpg['diff'] = abs(smoothed_r9_r10_wg_2kb_70cpg['methylFq_r9']-smoothed_r9_r10_wg_2kb_70cpg['methylFq_r10'])
print('avd diff', smoothed_r9_r10_wg_2kb_70cpg['diff'].mean())
"""
2Kb Whole Genome Filtered Correlation:
48300 sites 
Pearson correlation coefficient: 0.9980392726556955
RMSE: 3.7484800508354565


avd diff 3.178550724637681

"""

# Filter by cpg Density of the window
df_2kb_filt = smoothed_r9_r10_wg_2kb_70cpg.loc[smoothed_r9_r10_wg_2kb_70cpg['cpgDensity']>0.025]




###############################################
### 2Kb Promoter CpGs:                      ###
###############################################
# 
headr = ["chr",'start','end','promoter_name','number','strand','ratio','avgMeth']
r9_mt_2kb = pd.read_csv('HG002_R9_CARD_bam.comb.regionMethyl.bed', sep="\t", header=None,
                names=headr)
r10_mt_2kb = pd.read_csv('HG002_R10.comb.regionMethyl.bed', sep='\t', header=None,
                names=headr)

for col in r9_mt_2kb.columns[6:]:
    r9_mt_2kb[col] = pd.to_numeric(r9_mt_2kb[col], errors='coerce')
for col in r10_mt_2kb.columns[6:]:
    r10_mt_2kb[col] = pd.to_numeric(r10_mt_2kb[col], errors='coerce')
                          

r9_r10_mt_2kb = pd.merge(r9_mt_2kb, r10_mt_2kb, on=['chr','start','end'], suffixes=('_r9', '_r10'))

r9_r10_mt_2kb['diff']=r9_r10_mt_2kb['avgMeth_r10'].astype(float)-r9_r10_mt_2kb['avgMeth_r9'].astype(float)
r9_r10_mt_2kb['genomic_size']= r9_r10_mt_2kb['end'] - r9_r10_mt_2kb['start']



###############################################
### CpG Island CpGs:                        ###
###############################################
# 
islandFields = ['chrom','start','end','name','length','cpgNum','gcNum','perCpg','perGc','obsExp','ratio','avgMeth']
r9_cpg = pd.read_csv('HG002_R9_CARD_bam.comb.regionMethyl.bed', sep="\t", header=None,names=islandFields)

r10_cpg = pd.read_csv('HG002_R10.comb.regionMethyl.bed', sep='\t', header=None,names=islandFields)


for col in r9_cpg.columns[6:]:
    r9_cpg[col] = pd.to_numeric(r9_cpg[col], errors='coerce')
for col in r10_cpg.columns[6:]:
    r10_cpg[col] = pd.to_numeric(r10_cpg[col], errors='coerce')


r9_r10_cpg = pd.merge(r9_cpg, r10_cpg, on=['chrom','start','end'], suffixes=('_r9', '_r10'))




###############################################
### Scatter plots of filtering :            ###
###############################################
## A scatter plot for 2kb Promoters after filtering for sites with coverage ratio > 0.9
# Create a subplots figure
fig, axs = plt.subplots(1, 2, figsize=(12, 5))
# plt.figure(figsize=(8, 8))


ratio_min = 0.9
alpha_val = 0.5
sc = axs[0].scatter(r9_r10_mt_2kb['avgMeth_r9'], r9_r10_mt_2kb['avgMeth_r10'],s=3, alpha=alpha_val) 
#                c=r9_r10_mt_2kb['num_cpgs_r10'], cmap='viridis')
filt_data_mt = r9_r10_mt_2kb.loc[(r9_r10_mt_2kb['ratio_r9'] > ratio_min) 
                  & (r9_r10_mt_2kb['ratio_r10'] > ratio_min) 
                  & (r9_r10_mt_2kb['avgMeth_r9'] >= 0)
                  & (r9_r10_mt_2kb['avgMeth_r10'] >= 0)]

sc = axs[1].scatter(filt_data_mt['avgMeth_r9'], filt_data_mt['avgMeth_r10'],s=3, alpha=alpha_val)
#               c=filt_data_isl_mt['num_cpgs_r10'], cmap='viridis')

axs[0].set_title('All Data: '+str(r9_r10_mt_2kb.shape[0])+' sites')
axs[0].set_xlabel('R9 Promoter Methylation Frequency')
axs[0].set_ylabel('R10 Promoter Methylation Frequency')

axs[1].set_title('Filtered Data : '+str(filt_data_mt.shape[0])+' sites')
axs[1].set_xlabel('R9 Promoter Methylation Frequency')
axs[1].set_ylabel('R10 Promoter Methylation Frequency')

# plt.tight_layout()
# Add colorbar to the right of the subplots
# cbar = fig.colorbar(sc, ax=axs.ravel().tolist(), label='Number of CpGs in R10')

fig.suptitle('HG002 2kb Promoter Average Methylation Frequency Regional Mean')


## A scatterplot for CpG Islands before and after filtering for sites with coverage ratio > 0.9
# Create a subplots figure
fig, axs = plt.subplots(1, 2, figsize=(12, 5))
# plt.figure(figsize=(8, 8))


ratio_min = 0.9
alpha_val = 0.5
sc = axs[0].scatter(r9_r10_cpg['avgMeth_r9'], r9_r10_cpg['avgMeth_r10'],s=3, alpha=alpha_val) 
#                c=r9_r10_mt_2kb['num_cpgs_r10'], cmap='viridis')
filt_data_cpg = r9_r10_cpg.loc[(r9_r10_cpg['ratio_r9'] > ratio_min) 
                  & (r9_r10_cpg['ratio_r10'] > ratio_min) 
                  & (r9_r10_cpg['avgMeth_r9'] >= 0)
                  & (r9_r10_cpg['avgMeth_r10'] >= 0)]

sc = axs[1].scatter(filt_data_cpg['avgMeth_r9'], filt_data_cpg['avgMeth_r10'],s=3, alpha=alpha_val)
#               c=filt_data_isl_mt['num_cpgs_r10'], cmap='viridis')

axs[0].set_title('All Data: '+str(r9_r10_cpg.shape[0])+' sites')
axs[0].set_xlabel('R9 CpG Island Methylation Frequency')
axs[0].set_ylabel('R10 CpG Island Methylation Frequency')

axs[1].set_title('Filtered Data : '+str(filt_data_cpg.shape[0])+' sites')
axs[1].set_xlabel('R10 CpG Island Methylation Frequency')
axs[1].set_ylabel('R10 CpG Island Methylation Frequency')

# plt.tight_layout()
# Add colorbar to the right of the subplots
# cbar = fig.colorbar(sc, ax=axs.ravel().tolist(), label='Number of CpGs in R10')

fig.suptitle('HG002 CpG Island Average Methylation Frequency Regional Mean')

###############################################
### Promoter and CpG Islands Correlations:  ###
###############################################

r9_r10_mt_2kb_dna = r9_r10_mt_2kb.dropna()

pearson_corr = r9_r10_mt_2kb['avgMeth_r9'].corr(r9_r10_mt_2kb['avgMeth_r10'])

print(r9_r10_mt_2kb.shape[0], 'sites',"\nPearson correlation coefficient:", pearson_corr)

print('RMSE:',np.sqrt(mean_squared_error(r9_r10_mt_2kb['avgMeth_r10'].dropna(), r9_r10_mt_2kb['avgMeth_r9'].dropna())))

pearson_corr = r9_r10_mt_2kb_dna['avgMeth_r9'].corr(r9_r10_mt_2kb_dna['avgMeth_r10'])
print('\n',r9_r10_mt_2kb_dna.shape[0], 'sites',"\nDNA Pearson correlation coefficient:", pearson_corr)
print('RMSE:',np.sqrt(mean_squared_error(r9_r10_mt_2kb_dna['avgMeth_r10'].dropna(), r9_r10_mt_2kb_dna['avgMeth_r9'].dropna())))



r9_r10_mt_2kb['diff'] = abs(r9_r10_mt_2kb['avgMeth_r9']-r9_r10_mt_2kb['avgMeth_r10'])
print('\n','avd diff', r9_r10_mt_2kb['diff'].mean())

r9_r10_mt_2kb_dna['diff'] = abs(r9_r10_mt_2kb_dna['avgMeth_r9']-r9_r10_mt_2kb_dna['avgMeth_r10'])
print('DNA avd diff', r9_r10_mt_2kb_dna['diff'].mean())

filt_pro_data = r9_r10_mt_2kb.loc[(r9_r10_mt_2kb['ratio_r9'] > ratio_min) 
                  & (r9_r10_mt_2kb['ratio_r10'] > ratio_min) 
                  & (r9_r10_mt_2kb['avgMeth_r9'] >= 0)
                  & (r9_r10_mt_2kb['avgMeth_r10'] >= 0)]
print('\n',filt_pro_data.shape[0], 'sites',"\nPearson correlation coefficient:", pearson_corr)

print('RMSE:',np.sqrt(mean_squared_error(filt_pro_data['avgMeth_r10'].dropna(), filt_pro_data['avgMeth_r9'].dropna())))

filt_pro_data_dna = r9_r10_mt_2kb_dna.loc[(r9_r10_mt_2kb_dna['ratio_r9'] > ratio_min) 
                  & (r9_r10_mt_2kb_dna['ratio_r10'] > ratio_min) 
                  & (r9_r10_mt_2kb_dna['avgMeth_r9'] >= 0)
                  & (r9_r10_mt_2kb_dna['avgMeth_r10'] >= 0)]
print('\n',filt_pro_data_dna.shape[0], 'sites',"\nPearson correlation coefficient:", pearson_corr)

print('RMSE:',np.sqrt(mean_squared_error(filt_pro_data_dna['avgMeth_r10'].dropna(), filt_pro_data_dna['avgMeth_r9'].dropna())))

"""
Promoter Raw Correlation:
29772 sites 
Pearson correlation coefficient: 0.998026744455362
RMSE: 24.089320830389386

# Filtered by dropping NA valies (DNA)
 29757 sites 
DNA Pearson correlation coefficient: 0.998026744455362
RMSE: 4.606664159909091

 avd diff 4.142792146210075
DNA avd diff 4.142792146210075

# Filtered by coverage min # of CpGs
 29124 sites 
Pearson correlation coefficient: 0.998026744455362
RMSE: 4.586001347260341

 29124 sites 
Pearson correlation coefficient: 0.998026744455362
RMSE: 4.586001347260341
"""

pearson_corr = r9_r10_cpg['avgMeth_r9'].corr(r9_r10_cpg['avgMeth_r10'])

print(r9_r10_cpg.shape[0], 'sites',"\nPearson correlation coefficient:", pearson_corr)
# print('RMSE:',np.sqrt(mean_squared_error(r9_r10_cpg['avgMeth_r10'].dropna(), r9_r10_cpg['avgMeth_r9'].dropna())))

r9_r10_cpg_dna = r9_r10_cpg.dropna()
pearson_corr = r9_r10_cpg_dna['avgMeth_r9'].corr(r9_r10_cpg_dna['avgMeth_r10'])

print('\n',r9_r10_cpg_dna.shape[0], 'sites',"\nDNA Pearson correlation coefficient:", pearson_corr)

print('RMSE:',np.sqrt(mean_squared_error(r9_r10_cpg_dna['avgMeth_r10'], r9_r10_cpg_dna['avgMeth_r9'])))


r9_r10_cpg['diff'] = abs(r9_r10_cpg['avgMeth_r9']-r9_r10_cpg['avgMeth_r10'])
print('\n','avd diff', r9_r10_cpg['diff'].mean())

filt_cpg_data = r9_r10_cpg.loc[(r9_r10_cpg['ratio_r9'] > ratio_min) 
                  & (r9_r10_cpg['ratio_r10'] > ratio_min) 
                  & (r9_r10_cpg['avgMeth_r9'] >= 0)
                  & (r9_r10_cpg['avgMeth_r10'] >= 0)]
print(filt_cpg_data.shape[0], 'sites',"\nPearson correlation coefficient:", pearson_corr)
print('RMSE:',np.sqrt(mean_squared_error(filt_cpg_data['avgMeth_r10'], filt_cpg_data['avgMeth_r9'])))

"""
CpG Island Raw Correlation:
27949 sites 
Pearson correlation coefficient: 0.9982378322387357

# Filtered by dropping NA valies (DNA)
 27767 sites 
DNA Pearson correlation coefficient: 0.9982378322387357
RMSE: 5.117415501526385

# Filtered by coverage min # of CpGs
 avd diff 4.6007682730546735
27579 sites 
Pearson correlation coefficient: 0.9982378322387357
RMSE: 5.119609441279698
"""

### HISTOGRAM WITH SCATTER PLOT ###

tab_colors = ['tab:blue', 'tab:olive', 'tab:cyan', 'tab:pink']
colors = tab_colors[0:2]
# Convert 'tab:blue' and 'tab:olive' to RGB
color1 = to_rgb('#1f77b4') 
color2 = to_rgb('#ff7f0e') 

# Blend the colors (average the RGB components)
blended_color = tuple((0.7 * c1 + 0.3 * c2) / 2 for c1, c2 in zip(color1, color2))

# probably make this a funciton as well
covMin = 5
wg_m_data = r9_r10_wg.loc[(r9_r10_wg['Nvaid_cov_r9'] > covMin) & (r9_r10_wg['Nvaid_cov_r10'] > covMin)]

cpg_data = r9_r10_cpg.loc[(r9_r10_cpg['read_depth_r9'] > covMin) & 
                   (r9_r10_cpg['read_depth_r10'] > covMin) &
                   (r9_r10_cpg['agg_region_mean_r9'] >= 0)&
                   (r9_r10_cpg['agg_region_mean_r10'] >= 0)]


pro_data = [
    r9_r10_mt.loc[(r9_r10_mt['read_depth_r9'] > 5) & (r9_r10_mt['read_depth_r10'] > 5) & (r9_r10_mt['agg_region_mean_r9'] >= 0)& (r9_r10_mt['agg_region_mean_r10'] >= 0)][
        'agg_region_mean_r9'],
    r9_r10_mt.loc[(r9_r10_mt['read_depth_r9'] > 5) & (r9_r10_mt['read_depth_r10'] > 5) & (r9_r10_mt['agg_region_mean_r10'] >= 0)& (r9_r10_mt['agg_region_mean_r9'] >= 0)][
        'agg_region_mean_r10']]

filt_pro_data = filt_r9_r10_mt_2kb.loc[(filt_r9_r10_mt_2kb['read_depth_r9'] > 5) 
                                         & (filt_r9_r10_mt_2kb['read_depth_r10'] > 5) 
                                         & (filt_r9_r10_mt_2kb['agg_region_mean_r9'] >= 0)
                                         & (filt_r9_r10_mt_2kb['agg_region_mean_r10'] >= 0)]

fig = plt.figure(figsize=(22, 19))

# Define the gridspec with 3 rows and 2 columns
# Specify the width ratios: [2, 1] means the first column is 2.5x the width of the second column
gs = gridspec.GridSpec(3, 2, width_ratios=[2.5, 1], wspace=0.2, hspace=0.3)

logbool = True
alpha_val = 0.2
axis_label_fontsize = 16
axis_title_fontsize = 18
figure_text_fontsize = 18
# CpG sites #################################
ax1 = fig.add_subplot(gs[0, 0])
# ax1.hist(data1, bins=30, alpha=0.7, color='blue')
n, bins, patches = ax1.hist([wg_m_data['fracMod_r9'],wg_m_data['fracMod_r10']], bins=np.floor(np.linspace(0, 101, 21)), 
                             alpha=0.25, log=True)
# Calculate the bin centers
bin_centers = np.floor(np.linspace(0, 101, 21))
ax1.legend(['R9', 'R10'], prop={'size': axis_title_fontsize})
ax1.set_xlabel('Methylation Frequency', fontsize=axis_label_fontsize)

# Set xticks and labels to align with the bars
ax1.set_xticks(bin_centers, range(0, 101, 5),rotation=90, fontsize=16)
# ax1.set_xticklabels(range(0, 101, 5))

ax1.tick_params(axis='y', labelsize=16) 
ax1.set_ylabel('Counts (cov>5)', fontsize=axis_label_fontsize)
ax1.set_title('HG002 Methylation Frequency R9 & R10 ('+str(r9_r10_wg.shape[0])+' sites)', fontsize=axis_title_fontsize)
ax1.text(-0.1, 1.1, 'a.', transform=ax1.transAxes, fontsize=figure_text_fontsize, va='top', ha='right')

##### CpG sites 2kb scatter: 
ax2 = fig.add_subplot(gs[0, 1])
scatter2 = ax2.scatter(df_2kb_filt['methylFq_r9'], df_2kb_filt['methylFq_r10'],
                       s=3, alpha=alpha_val, color=blended_color)
#               c=df_2kb_filt['cpgDensity'], cmap='viridis')
ax2.plot([0, 100],[0, 100],color='#00CED1',linewidth=2.5)

ax2.set_title('2Kb Whole Genome', fontsize=axis_title_fontsize)
ax2.set_xlabel('R9 2Kb Methylation Frequency', fontsize=axis_label_fontsize)
ax2.set_ylabel('R10 2Kb Methylation Frequency', fontsize=axis_label_fontsize)
ax2.tick_params(axis='x', labelsize=16) 
ax2.tick_params(axis='y', labelsize=16) 
# cbar2 = fig.colorbar(scatter2, ax=ax2)
# cbar2.set_label('CpG Density')
ax2.text(-0.2, 1.1, 'b.', transform=ax2.transAxes, fontsize=figure_text_fontsize, va='top', ha='right')


                       
# Promoters #################################
ax3 = fig.add_subplot(gs[1, 0])
n, bins, patches = ax3.hist([filt_pro_data['avgMeth_r9'],filt_pro_data['avgMeth_r10']], bins=np.floor(np.linspace(0, 101, 21)), log=logbool,
                             alpha=0.25)

ax3.legend(['R9', 'R10'], prop={'size': 20})
ax3.set_xlabel('Regional Methylation Frequency', fontsize=axis_label_fontsize)

# Set xticks and labels to align with the bars
ax3.set_xticks(bin_centers, range(0, 101, 5),rotation=90, fontsize=16)
ax3.tick_params(axis='y', labelsize=16) 
ax3.set_ylabel('Counts (cov>5)', fontsize=axis_label_fontsize)
ax3.set_title('HG002 Aggregated 2Kb Promoter Methylation Frequency ('+str(filt_pro_data.shape[0])+' sites)', fontsize=axis_title_fontsize)
ax3.text(-0.1, 1.1, 'c.', transform=ax3.transAxes, fontsize=figure_text_fontsize, va='top', ha='right')

## Promoters 2kb scatter: 
ax4 = fig.add_subplot(gs[1, 1])
scatter4 = ax4.scatter(filt_pro_data['avgMeth_r9'], filt_pro_data['avgMeth_r10'],
                       s=3, alpha=alpha_val, color=blended_color)
ax4.plot([0, 100],[0, 100],color='#00CED1',linewidth=2.5)

ax4.set_title('2Kb Promoters', fontsize=axis_title_fontsize)
ax4.set_xlabel('R9 Methylation Frequency', fontsize=axis_label_fontsize)
ax4.set_ylabel('R10 Methylation Frequency', fontsize=axis_label_fontsize)
ax4.set_xticks(np.arange(0, 101, 20), np.arange(0, 101, 20), fontsize=16)
ax4.set_yticks(np.arange(0, 101, 20), np.arange(0, 101, 20), fontsize=16)
# cbar4 = fig.colorbar(scatter4, ax=ax4)
# cbar4.set_label('R9 CpG Density')
ax4.text(-0.2, 1.1, 'd.', transform=ax4.transAxes, fontsize=figure_text_fontsize, va='top', ha='right')


# Islands #################################
ax5 = fig.add_subplot(gs[2, 0])
n, bins, patches = ax5.hist([cpg_data['avgMeth_r9'],cpg_data['avgMeth_r10']], bins=np.floor(np.linspace(0, 101, 21)), log=logbool,
                             alpha=0.25)

ax5.legend(['R9', 'R10'], prop={'size': 20})
ax5.set_xlabel('Regional Methylation Frequency', fontsize=axis_label_fontsize)


# Set xticks and labels to align with the bars
ax5.set_xticks(bin_centers, range(0, 101, 5),rotation=90, fontsize=16)
ax5.tick_params(axis='y', labelsize=16) 
ax5.set_ylabel('Counts (cov>5)', fontsize=axis_label_fontsize)
ax5.set_title('HG002 Aggregated CpG Island Methylation Frequency ('+str(cpg_data.shape[0])+' sites)', fontsize=axis_title_fontsize)
ax5.text(-0.1, 1.1, 'e.', transform=ax5.transAxes, fontsize=figure_text_fontsize, va='top', ha='right')

## Islands 2kb scatter: 
ax6 = fig.add_subplot(gs[2, 1])
scatter6 = ax6.scatter(cpg_data['avgMeth_r9'], cpg_data['avgMeth_r10'],
                       s=3, alpha=alpha_val, color=blended_color)
ax6.plot([0, 100],[0, 100],color='#00CED1',linewidth=2.5)              

ax6.set_title('CpG Islands', fontsize=axis_title_fontsize)
ax6.set_xlabel('R9 Methylation Frequency', fontsize=axis_label_fontsize)
ax6.set_ylabel('R10 Methylation Frequency', fontsize=axis_label_fontsize)
ax6.set_xticks(np.arange(0, 101, 20), np.arange(0, 101, 20), fontsize=16)
ax6.set_yticks(np.arange(0, 101, 20), np.arange(0, 101, 20), fontsize=16)
# cbar6 = fig.colorbar(scatter6, ax=ax6)
# cbar6.set_label('R9 CpG Density')
ax6.text(-0.2, 1.1, 'f.', transform=ax6.transAxes, fontsize=figure_text_fontsize, va='top', ha='right')


# plt.tight_layout()
# plt.show()
plt.savefig("HG002_methylation_2kbPromoter_scatter_log_06102024.png", dpi=300)
plt.savefig("HG002_methylation_2kbPromoter_scatter_log_06102024.svg", dpi=300)


