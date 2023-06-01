import sys
import os
import re
import gzip
import click
import subprocess
import pandas as pd
import numpy as np
from natsort import natsorted
from dbclients.colossus import ColossusApi

colossus_api = ColossusApi()

@click.command()
@click.argument("files", nargs = -1)
@click.argument("output_file", nargs = 1)
def merge_libs(**kwargs):
    '''
    (1) Given sample_ids and cached directory, filter raw reads and metrics for mappability >=0.99;
                                                                                read quality >=0.75;
                                                                                total mapped reads >= 500000;
                                                                                live cells;
                                                                                non S-phase cells;
                                                                                contaminants of SA501 & SA928;

    (2) Reshape to extract state, copy, and quantile transformed copy data

    Args:
        sample ids ("," sep str)
        cached dir path (path)
        good_cell_dir (path)
        output_dir (path)

    Returns:
        total combined state df
    '''

    best_loci=set()
    total_state_samples=pd.DataFrame()

    for filename in kwargs["files"]:
        try:
            filtered_state_data = pd.read_csv(filename)
        except:
            print("no filtered reads in", filename)
            continue

        if len(best_loci) == 0:
            best_loci = set(filtered_state_data.index.values)
        else:
            best_loci = best_loci.intersection(set(filtered_state_data.index.values))

        mask = pd.Series(filtered_state_data.index.values).isin(best_loci)

        if total_state_samples.empty:
            total_state_samples = filtered_state_data.loc[mask.values, :]
        else:
            total_state_samples = total_state_samples.join(filtered_state_data.loc[mask.values, :])


    total_state_samples = total_state_samples.reindex(index=natsorted(total_state_samples.index))
#     total_state_samples.to_csv(kwargs["output_file"], sep = ",", index_label = False)
    # get_noisy_cells
#     if not total_state_samples.empty:
#         pb = 0.9
#         get_noisy_cells(total_state_samples, kwargs["output_file"], pb)

    
def get_noisy_cells(cnv_df, output_file, pb=0.9):
    
    print("CNV mtx before removing noise: {0} {1}".format(cnv_df.shape[0],cnv_df.shape[1]))

    avg_ploidy = cnv_df.mean(axis=0)
    is_avg_ploidy = avg_ploidy < avg_ploidy.quantile(pb)
    
    change_cnv = []
    for c in range(cnv_df.shape[1]):
        count_chg = 0
        currv = cnv_df.iloc[0,c]
        for r in range(1, cnv_df.shape[0]):  
            if cnv_df.iloc[r,c] != currv:
                currv = cnv_df.iloc[r,c]
                count_chg = count_chg + 1


        change_cnv.append(count_chg) 
        
#     print(change_cnv[1:4]) 
    qtch = pd.Series(change_cnv).quantile(pb)
    is_change = change_cnv < qtch
    is_noise = is_avg_ploidy | is_change
    print("Nb cells as noise: {0}".format(sum(is_noise==False)))
    cnv_stat = pd.DataFrame({'avg_ploidy': avg_ploidy,
                            'is_avg_ploidy': is_avg_ploidy,
                            'nb_changes': change_cnv,
                            'is_noise': is_noise}, index=cnv_df.keys())
#     print("Statistic mtx: {0} {1}".format(cnv_stat.shape[0],cnv_stat.shape[1]))
    cnv_df = cnv_df.loc[:,cnv_stat['is_noise']]
    print("CNV mtx after removing noise: {0} {1}".format(cnv_df.shape[0],cnv_df.shape[1]))
    if not cnv_df.empty:
        cnv_df.to_csv(output_file, sep = ",", index_label = False)
        print("Save output into : {0}".format(output_file))

if __name__=='__main__':
    merge_libs()
