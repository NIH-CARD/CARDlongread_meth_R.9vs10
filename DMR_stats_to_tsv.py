import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
import argparse

#Take the DMR stats dict and write it to a .tsv
def stats_dicts_to_tsv(stat_dict_list, label_list, outpath):
    d = {}
    for stat_dict, label in zip(stat_dict_list, label_list):
        d[label] = stat_dict
    df = pd.DataFrame.from_dict(d)
    df.T.to_csv(outpath, sep='\t')

#Calculate basic summary statistics for DMRs
def dmr_stats(inpath):
    df = pd.read_csv(inpath, sep = '\t')
    df['bases'] = df['end']-df['start']
    dmr_stat_dict = {'DMR Count': df.shape[0],
                     'Bases': df['bases'].sum(),
                     'nCG mean' : df['nCG'].mean(),
                     'nCG std' : df['nCG'].std(),
                     'nCG count' : df['nCG'].sum(),
                     'Methylation Difference Mean' : df['diff.Methy'].abs().mean(),
                     'Methylation Difference Std' : df['diff.Methy'].abs().std()
                    }
    return dmr_stat_dict

def main(dmr_path_list, dmr_label_list, outpath):
    dmr_stats_dicts_list = []
    dmr_stats_labels_list = []
    for path in dmr_path_list:
        dmr_stats_dicts_list.append(dmr_stats(path))
        
    stats_dicts_to_tsv(dmr_stats_dicts_list, dmr_label_list, outpath)

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--dmr_list",
                        "-i",
                        required=True,
                        nargs="+",
                        help="List of dmr files from Nanomethphase for parsing"
                       )
    parser.add_argument("--dmr_labels",
                        "-l",
                        required=True,
                        nargs="+",
                        help="List of dmr file labels for output tsv"
                       )
    parser.add_argument("--output",
                        "-o",
                        required = True,
                        type = str,
                        help = "Path to save tsv of DMR stats"
                       )
    
    args = parser.parse_known_args()[0]
    
    if len(args.dmr_list) != len(args.dmr_labels):
        print("The length of DMR files is different than the length of labels.")
        print("Exiting Program")
        
    else:
        main(args.dmr_list, args.dmr_labels, args.output)
