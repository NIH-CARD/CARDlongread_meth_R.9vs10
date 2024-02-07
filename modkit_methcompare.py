import pandas as pd
import numpy as np
import seaborn as sns
import argparse
import sys

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

'''SAMPLE COMMAND:

python modkit_methcompare.py \
--cov_min 20 \
--cov_max 200 \
--r9_modkit /Users/gennerrm/Desktop/modkit_files/HG002_R9.hg38.modkit.comb.bed \
--r10_modkit /Users/gennerrm/Desktop/modkit_files/HG002_R10.hg38.modkit.comb.bed \
--bis_modkit /Users/gennerrm/Desktop/modkit_files/CpG.gz.bismark.zero.merged.cov \
--interval 10 \
--sample_name HG002_bis \
--binning Bisulfite  \
--out_dir /Users/gennerrm/Desktop/figs/py_figs/test'''

parser = argparse.ArgumentParser()

parser.add_argument('--cov_min', type=int, default=20,
                    help='The minimum coverage, an integer, default 20, require cov_min <= cov_max.')

parser.add_argument('--cov_max', type=int, default=200,
                    help='The maximum coverage, an integer, default 200, require cov_min <= cov_max.')

parser.add_argument('--interval', type=int, default=10,
                    help='The binning interval, an integer, default 10.')

parser.add_argument('--custom_interval', type=list,
                    help='Custom binning intervals as a list, ex. [0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100]')
                    
parser.add_argument('--bw', type=float, default=0.1,
                    help='Factor (0.0-1.0) that scales the bandwidth to use more or less smoothing, default is 0.1') 
                    
parser.add_argument('--scale', type=str, default='width',
                    help='Method that normalizes each density to determine the violin’s width. Options are: area (each violin will have the same area), count (the width will be proportional to the number of observations), or width (default; each violin will have the same width).')

parser.add_argument('--r9_modkit', type=str, required=True,
                    help='The bedfile for r9_modkit, required.')

parser.add_argument('--r10_modkit', type=str, required=True,
                    help='The bedfile for r10_modkit, required.')

parser.add_argument('--bis_modkit', type=str, required=True,
                    help='The modkit file for bisulfite sequencing.')

parser.add_argument('--binning', type=str, required=True,
                   help='Dataset to bin by, must be either R9, R10 or Bisulfite, required')

parser.add_argument('--sample_name', type=str, required=True,
                    help='sample_name is the title and filename of your image file.')

parser.add_argument('--out_dir', type=str, required=False,
                    help='out_dir (optional) is the directory to put your image file into.')

args = vars(parser.parse_args())

# Extract variables from command line arguments.


cov_min = args['cov_min']
cov_max = args['cov_max']
interval = args['interval']
custom_interval = args['custom_interval']
r9_modkit = args['r9_modkit']
r10_modkit = args['r10_modkit']
bis_modkit = args['bis_modkit']
binning = args['binning']
bw = args['bw']
scale = args['scale']
scale = args['scale']
sample_name = args['sample_name']
out_dir = args['out_dir']

for key, value in args.items():
    print(key, '=', value)


if cov_min >= cov_max:
    print(f"Oops, cov_min <= cov_max is required, but you have cov_min={cov_min} and cov_max={cov_max}.",
          file=sys.stderr)
    exit(1)

if binning not in ('R9', 'R10', 'Bisulfite'):
    print(f'Oops, graph_binning has only two allowed values: "R9", "R10", "Bisulfite".',
          file=sys.stderr)
    exit(1)


select = {"chr1", "chr2", "chr3", "chr4",'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'}

print("Reading R9 data...")
R9_df = pd.read_csv(r9_modkit, sep="\t", header=None, engine="c",
    usecols=["chrom", "start", "end", "coverage", "freq", "mod", "canon"],
    dtype={
        "chrom":str, "start":int, "end":int, "name":str, "score":int, "strand":str,
        "start2":int, "end2":int, "color":str,
        "coverage":int, "freq":float, "mod":int, "canon":int, "other_mod":int, 'delete':int, 'fail':int, 'diff':int, 'no_call':int},
    names=[
        "chrom", "start", "end", "name", "score", "strand",
        "start2", "end2", "color",
        "coverage", "freq", "mod", "canon", "other_mod", "delete", "fail", 'diff', 'no_call'])

#filter for  chromosomes 1-22,X,Y,M and coverage > 20 and < 200
R9_df = R9_df.loc[
    (R9_df['chrom'].isin(select)) &
    (R9_df['coverage'] > cov_min) &
    (R9_df['coverage'] < cov_max)]

print("Reading R10 data...")
R10_df = pd.read_csv(r10_modkit,
    sep="\t", header=None, engine="c",
    usecols=["chrom", "start", "end", "coverage", "freq", "mod", "canon"],
    dtype={
        "chrom":str, "start":int, "end":int, "name":str, "score":int, "strand":str,
        "start2":int, "end2":int, "color":str,
        "coverage":int, "freq":float, "mod":int, "canon":int, "other_mod":int, 'delete':int, 'fail':int, 'diff':int, 'no_call':int},
    names=[
        "chrom", "start", "end", "name", "score", "strand",
        "start2", "end2", "color",
        "coverage", "freq", "mod", "canon", "other_mod", "delete", "fail", 'diff', 'no_call'])

#filter for coverage > 20 and < 200
R10_df = R10_df.loc[
    (R10_df['chrom'].isin(select)) &
    (R10_df['coverage'] > cov_min) &
    (R10_df['coverage'] < cov_max)]


print("Reading bisulfite data...")
bis_df = pd.read_csv(
    bis_modkit,
    sep="\t", header=None, engine="c",
    dtype={'chrom':str, 'start':int, 'end':int, 'freq.bis':float, 'mod':int, 'canon':int},
    names=['chrom', 'start', 'end', 'freq.bis', 'mod', 'canon']
)
bis_df["coverage.bis"] = bis_df["mod"] + bis_df["canon"]

#filter for coverage > 20 and < 200
bis_df = bis_df.loc[
    (bis_df['chrom'].isin(select)) &
    (bis_df['coverage.bis'] > cov_min) &
    (bis_df['coverage.bis'] < cov_max)]

# inner merge to get dataframe with CpG sites that are shared between R9 and R10
R9_R10 = pd.merge(R9_df, R10_df, how="inner", on=["chrom", "start", "end"], suffixes=[".R9", ".R10"])
inner = pd.merge(R9_R10, bis_df, how="inner", on=["chrom", "start", "end"])

def bin(interval):
    return list(range(0, 101, interval))    

if custom_interval is not None:
    bins = custom_interval
else:
    bins = bin(interval)
    
if binning == 'Bisulfite':
    binning_name = 'freq.bis'
else:
    binning_name = 'freq.' + binning

#make subset dataframe with R9 and R10 frequency columns 
freqs = inner[['freq.R9', 'freq.R10', 'freq.bis']]

# this adds a column that rounds up each R9 methylation frequency to the tenths place - ex. 87.3.5 -> 90.0, and assigns the corresponding R10 value to that same interval regardless of its value
# This is binning the intervals according to R9. Could also bin by R10 by switching it to pd.cut(freqs['freq.R10'])

freqs['ints'] = pd.cut(freqs[binning_name], bins=bins)

#convert df from int to str
freqs = freqs.astype({'ints':'string'})

#delete the start value of the interval so it's just one number - so (90.0, 100.0] becomes 100 - and then convert back to float 
freqs['ints'] = freqs['ints'].str[1:-1]
freqs[['int_start', 'int_end']] = freqs.ints.str.split(',', expand=True)
freqs['int_end'] = freqs['int_end'].astype(float)

#grab only these columns 
freqs = freqs[['freq.R9', 'freq.R10', 'freq.bis', 'int_end']]

# add two columns with R9 and R10 for when I separate them out later 
freqs.insert(1, 'R9', 'R9')
freqs.insert(2, 'R10', 'R10')
freqs.insert(3, 'bis', 'bis')

freqs_dict = {}
R9_freqs = {}
R10_freqs = {}
bis_freqs = {}

R9_R10_combined={}

for i in bins: 
    freqs_dict[i] = freqs.loc[freqs['int_end'] == float(i)]
    for row in freqs_dict:
        R9_freqs[i] = freqs_dict[i][['R9','freq.R9', 'int_end']] 
        R9_freqs[i] = R9_freqs[i].rename(columns={'R9':'R', 'freq.R9':'freq'})
        R10_freqs[i] = freqs_dict[i][['R10','freq.R10', 'int_end']]
        R10_freqs[i] = R10_freqs[i].rename(columns={'R10':'R', 'freq.R10':'freq'})
        bis_freqs[i] = freqs_dict[i][['bis','freq.bis', 'int_end']]
        bis_freqs[i] = bis_freqs[i].rename(columns={'bis':'R', 'freq.bis':'freq'})

#R9_R10_combined = pd.concat([R9_freqs[0], R9_freqs[5], R9_freqs[10], R9_freqs[20], R9_freqs[30], R9_freqs[40], R9_freqs[50], R9_freqs[60], R9_freqs[70], R9_freqs[80], R9_freqs[90], R9_freqs[95], R9_freqs[100], R10_freqs[0], R10_freqs[5], R10_freqs[10], R10_freqs[20], R10_freqs[30], R10_freqs[40], R10_freqs[50], R10_freqs[60], R10_freqs[70], R10_freqs[80], R10_freqs[90], R10_freqs[95], R10_freqs[100]], axis=0)
#R9_R10_combined

R9_combined = pd.concat({**R9_freqs})
R10_combined = pd.concat({**R10_freqs})
R9_R10_combined = pd.concat([R9_combined, R10_combined])


x_axis = []

bin_len = list(range(0, len(bins)-1, 1))

for i in bin_len:
    z = str([bins[i],bins[i+1]])
    x_axis.append(z)


R9_meds = []
R10_meds = []

bis_meds = []
for i in bins:
    R9_meds.append(R9_freqs[i]['freq'].median())
    R10_meds.append(R10_freqs[i]['freq'].median())
    bis_meds.append(bis_freqs[i]['freq'].median())


meds = pd.DataFrame({'bins':bins, 'R9': R9_meds, 'R10': R10_meds, 'Bisulfite':bis_meds})
meds["bins"] = meds["bins"].astype("str")


fig, ax = plt.subplots(figsize=(20, 15))
sns.set_style("whitegrid", {'axes.grid' : True})

sns.violinplot(data=R9_R10_combined, x="int_end", y="freq", hue='R', split='TRUE', cut=0, inner='quartile', linewidth=1, bw=bw, scale=scale, ax=ax, palette="pastel")


ax.set_yticks(bins)
ax.set_yticklabels(bins, fontsize=20)
ax.set_xticklabels(x_axis, fontsize=20)
ax.set_xlabel(f"{binning} methylation intervals (%)", fontsize=20, labelpad=25)
ax.set_ylabel("% Methylation", fontsize=20)
plt.grid()
plt.legend(fontsize=25)
plt.title(f"{sample_name} R9 vs. R10 methylation proportions (binned by {binning} methylation intervals)", fontsize=25, pad=30)
sns.despine(left=True)

n = 2 
[l.set_visible(False) for (i,l) in enumerate(ax.yaxis.get_ticklabels()) if i % n != 0]

#change color, thickness and style of quartile lines
for l in ax.lines:
    l.set_linestyle('--')
    l.set_linewidth(0.95)
    l.set_color('grey')
    l.set_alpha(0.8)
for l in ax.lines[1::3]:
    l.set_linestyle('-')
    l.set_linewidth(1.2)
    l.set_color('black')
    l.set_alpha(0.8)
    
plt.savefig(f'{out_dir}/{sample_name}_VP.jpg', dpi=300)
plt.savefig(f'{out_dir}/{sample_name}_VP.png', dpi=300)
plt.savefig(f'{out_dir}/{sample_name}_VP.svg', dpi=300)


# Code for lineplot connecting median values

fig, ax = plt.subplots(figsize=(20, 15))
sns.set_style("whitegrid",{'axes.grid' : True})

sns.lineplot(x='bins', y='R9', data=meds, marker='o', markersize=15, linewidth=8, color='tab:blue')
sns.lineplot(x='bins', y='R10', data=meds, marker='o', markersize=15, linewidth=8, color='tab:orange')
sns.lineplot(x='bins', y='Bisulfite', data=meds, marker='o', markersize=10, linewidth=5, color='tab:green')

ax.set_yticks(bins)
ax.set_yticklabels(bins, fontsize=20)
ax.set_xticklabels(x_axis, fontsize=20)
ax.set_xlabel(f"{binning} methylation intervals (%)", fontsize=20, labelpad=25)
ax.set_ylabel("% Methylation", fontsize=20)
ax.grid()
plt.grid()
plt.legend(fontsize=25)
plt.title(f"{sample_name} R9 vs. R10 methylation proportions (binned by {binning} methylation intervals)", fontsize=25, pad=30)

sns.despine(left=True)
[l.set_visible(False) for (i,l) in enumerate(ax.yaxis.get_ticklabels()) if i % n != 0]

plt.savefig(f'{out_dir}/{sample_name}_{binning}bins_lines.jpg', dpi=300)
plt.savefig(f'{out_dir}/{sample_name}_{binning}bins_lines.png', dpi=300)
plt.savefig(f'{out_dir}/{sample_name}_{binning}bins_lines.svg', dpi=300)

# Code for side graph with methylation distribution
# the graph has to be rotated since it will be added to the side of the violin plot, so axis labels are descending instead of ascending


sns.set_palette("pastel")
sns.set_style("whitegrid")

fig, ax = plt.subplots(figsize=(20, 3))

sns.kdeplot(inner['freq.R9'], ax=ax, label='R9', color='tab:blue', cut=0, fill=True)
sns.kdeplot(inner['freq.R10'], ax=ax, label='R10', color='tab:orange', cut=0, fill=True)
sns.kdeplot(inner['freq.bis'], ax=ax, label='bis', color='tab:green', linewidth=3, cut=0)

ax.set_yticks((0.00, 0.10))
ax.set_yticklabels(('0', '0.1'), fontsize=30)

# graph has to be rotated
ax.tick_params(axis='y', labelrotation = 90)

#set x axis limit
ax.set_xlim(105, -7)

plt.xticks(visible=False)
ax.set_ylabel('Density', fontsize=30, labelpad=20)
ax.set_xlabel(None)

plt.savefig(f'{out_dir}/{sample_name}_VP_cov.jpg', dpi=300)
plt.savefig(f'{out_dir}/{sample_name}_VP_cov.png', dpi=300)
plt.savefig(f'{out_dir}/{sample_name}_VP_cov.svg', dpi=300)




