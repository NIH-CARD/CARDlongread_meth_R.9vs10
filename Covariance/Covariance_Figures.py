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


#Plot the kdes of covariance
def kde_plots(files,
              out_dir):
    color_dict = {'R10':'b', 'R9':'r'}
    ax_dict = {'HG002':0, 'PPMI':1, 'NABEC':2}
    fig, ax = plt.subplots(1, 3)
    fig.set_figwidth(24)
    fig.set_figheight(8)
    for f in files:
        df = pd.read_csv(f[0], sep='\t', header=None).sample(n=100000, axis = 0, random_state=22)
        df.columns = ['Covariance', 'Length', 'Read 1', 'Read 2']
        sns.kdeplot(x=df['Length'], 
                    y=df['Covariance'], 
                    label = f[2], 
                    color=color_dict[f[2]], 
                    ax=ax[ax_dict[f[1]]], 
                    clip=[[0,5000],[-30000, 30000]])
        ax[ax_dict[f[1]]].set_title(f'{f[1]} (n={df.shape[0]})')
        ax[ax_dict[f[1]]].set_xlim(0,5000)
        ax[ax_dict[f[1]]].set_ylim(-30000,30000)
    fig.legend()
    fig.get_figure().savefig(os.path.join(out_dir, f"kde_covariance.pdf"))

#Plot 2d hist of covariance for one condition (R9 or R10)
def rx_2dhist(files,
               out_dir,
               chemistry):
    color_dict = {'R10':'b', 'R9':'r'}
    ax_dict = {'HG002':0, 'PPMI':1, 'NABEC':2}
    fig, ax = plt.subplots(1, 3)
    fig.set_figwidth(24)
    fig.set_figheight(8)
    for f in files:
        df = pd.read_csv(f[0], sep='\t', header=None)
        df.columns = ['Covariance', 'Length', 'Read 1', 'Read 2']
        graph = ax[ax_dict[f[1]]].hist2d(x=df['Length'], 
                                 y=df['Covariance'], 
                                 norm=colors.LogNorm(), 
                                 bins=25,
                                 cmap='viridis',
                                 density = True,
                                 range=[[0,5000],[-40000,40000]]
                                        )
        ax[ax_dict[f[1]]].set_title(f'{f[1]} {f[2]} (n={df.shape[0]})')
        ax[ax_dict[f[1]]].set_ylim(-40000,40000)
        ax[ax_dict[f[1]]].set_xlim(0,5000)
    fig.colorbar(graph[3], label='Distribution Density per Square')
    fig.get_figure().savefig(os.path.join(out_dir, f"{chemistry}_2dhist.pdf"))

#Plot scatter of covariance with linear regression, extrodinarily slow
def scatter(files,
            point_size,
            out_dir,
            chemistry
           ):
    color_dict = {'R10':'b', 'R9':'r'}
    ax_dict = {'HG002':0, 'PPMI':1, 'NABEC':2}
    fig, ax = plt.subplots(1, 3)
    fig.set_figwidth(24)
    fig.set_figheight(8)
    for f in files:
        df = pd.read_csv(f[0], sep='\t', header=None).sample(n=100000, axis = 0, random_state=22)
        df.columns = ['Covariance', 'Length', 'Read 1', 'Read 2']
        a = sns.scatterplot(x=df['Length'], 
                            y=df['Covariance'], 
                            label = f[2], 
                            s=point_size, 
                            alpha = 0.05, 
                            color=color_dict[f[2]], 
                            ax=ax[ax_dict[f[1]]])
        sns.regplot(x=df['Length'], 
                    y=df['Covariance'], 
                    color=color_dict[f[2]], 
                    ax=ax[ax_dict[f[1]]], 
                    scatter=False, 
                    lowess=True)
        ax[ax_dict[f[1]]].set_title(f'{f[1]} (n={df.shape[0]})')
        ax[ax_dict[f[1]]].set_xlim(0,5000)
        ax[ax_dict[f[1]]].set_ylim(-40000,40000)
    fig.get_figure().savefig(os.path.join(out_dir,
                                        f"covariance_scatter_{chemistry}_{point_size}.pdf"))
    
def main(ppmi_r9,
         ppmi_r10,
         nabec_r9,
         nabec_r10,
         hg002_r9,
         hg002_r10,
         out_dir
        ):
    r9_component_list = [(hg002_r9[0], hg002_r9[1], hg002_r9[2]),
                         (ppmi_r9[0], ppmi_r9[1], ppmi_r9[2]),
                         (nabec_r9[0], nabec_r9[1], nabec_r9[2])]

    r10_component_list = [(hg002_r10[0], hg002_r10[1], hg002_r10[2]),
                         (ppmi_r10[0], ppmi_r10[1], ppmi_r10[2]),
                         (nabec_r10[0], nabec_r10[1], nabec_r10[2])]
    
    #2dhist for each chemistry
    print("Starting 2dhist")
    rx_2dhist(r9_component_list,
              out_dir,
              'r9')
    print("R9 2dhist complete")
    rx_2dhist(r10_component_list,
              out_dir,
              'r10')
    print("R10 2dhist complete")
    
    #kdeplots for both chemistries
    print("Starting KDE plots")
    kde_plots(r9_component_list + r10_component_list,
              out_dir)
    print("KDE plots complete")
    
    print("Starting Scatter")
    for s in [1, 5, 10, 20]:
        print(f"Scatter {s} starting")
        scatter(r9_component_list,
                s,
                out_dir,
                "r9")
        scatter(r10_component_list,
                s,
                out_dir,
                "r10")
    
    print("Completed Scatter")
if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--ppmi_r9", 
                        "-pr9",
                        nargs=3,
                        required=True,
                        help="covariance_path sample(hg002/nabec/ppmi) chemistry(r9/r10)"
                       )
    parser.add_argument("--ppmi_r10", 
                        "-pr10",
                        nargs=3,
                        required=True,
                        help="covariance_path sample(hg002/nabec/ppmi) chemistry(r9/r10)"
                       )
    parser.add_argument("--nabec_r9", 
                        "-nr9",
                        nargs=3,
                        required=True,
                        help="covariance_path sample(hg002/nabec/ppmi) chemistry(r9/r10)"
                       )
    parser.add_argument("--nabec_r10", 
                        "-nr10",
                        nargs=3,
                        required=True,
                        help="covariance_path sample(hg002/nabec/ppmi) chemistry(r9/r10)"
                       )
    parser.add_argument("--hg002_r9",
                        "-hgr9",
                        nargs=3,
                        required=True,
                        help="covariance_path sample(hg002/nabec/ppmi) chemistry(r9/r10)"
                       )
    parser.add_argument("--hg002_r10", 
                        "-hgr10",
                        nargs=3,
                        required=True,
                        help="covariance_path sample(hg002/nabec/ppmi) chemistry(r9/r10)"
                       )
    parser.add_argument("--out_dir",
                        required=True)
    FLAGS, unparsed = parser.parse_known_args()
    
    main(FLAGS.ppmi_r9,
         FLAGS.ppmi_r10,
         FLAGS.nabec_r9,
         FLAGS.nabec_r10,
         FLAGS.hg002_r9,
         FLAGS.hg002_r10,
         FLAGS.out_dir
        )