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
    total_state_samples.to_csv(kwargs["output_file"], sep = ",", index_label = False)


if __name__=='__main__':
    merge_libs()
