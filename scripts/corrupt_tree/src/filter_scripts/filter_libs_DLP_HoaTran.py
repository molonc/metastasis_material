import sys
import os
import re
import gzip
import argparse
import subprocess
import pandas as pd
import numpy as np
from natsort import natsorted
from distutils.version import StrictVersion



## control cells have the label as: *GM* or *NCC*, *NTC*, *gDNA*
def is_excluded_conditions(text):
    return bool(re.search(r'(?:GM|NCC|NTC|gDNA)', text, re.I))


def filter_libs(**kwargs):
    '''
    (1) Given sample_ids and cached directory, filter raw reads and metrics for mappability >=0.99;
                                                                                read quality >=0.75;
                                                                                total mapped reads >= 500000;
                                                                                live cells;
                                                                                non S-phase cells;
                                                                                contaminants cells;

    (2) Reshape to extract state, copy, and quantile transformed copy data

    Args:
        sample ids ("," sep str)
        cached dir path (path)
        good_cell_dir (path)
        output_dir (path)

    Returns:
        total combined state df
    '''
    
    ## Parameters setting
    
    
    ## Yaniv: from your filter: "non contaminated (!is_contaminated)". The contamination classification do not work well in some cases, 
    ## and you may cut off 1/3 of total cells depending on some specific library. So I usually do not include this filter in my analysis. 
    ## Shah lab also notice this part, and working on improving this classification. 
    ## I believe with the list of filters in the following scripts, almost contaminated cells are excluded. 
    ## Yaniv: if you want to apply this filter, you can check the number of cells before/after this filter, to avoid loosing lots of good cells. 
    
    
    ## Yaniv: mean %ideal >= 80% I don't really know about this filter. 
    ## From Daniel: %ideal is a boolean, true or false, not a percentage
    ## ideal is a column in HMMCopy output, it's basically anything that's valid (GC > 0, and mappability > 0.99) and also NOT an outlier
    ## outlier defined by "quantile" in R, filtering for 99% I think, it's mostly to ensure a good well-behaved set of data points enters the HMM segmentation function
    ## Yaniv: in my filter conditions, I do not include this condition. 
    
    
    ## Yaniv: Filter cells with outlier CN:  cell mean copy more than 1.5 sd of the population mean copy (up or down). 
    ## I used the threshold of 97.5% of 95%. 1.5 means 75%? You may filter out too many cells from this threshold. 
    
    
    save_dir = os.path.dirname(kwargs['output_file'])
    if not os.path.exists(save_dir): 
        os.makedirs(save_dir)
        print('Creating an output folder')
        
    best_loci = set()
    total_state_samples = pd.DataFrame()
    library_id = kwargs['library_id']
    # jira_id = find_latest_analysis_ticket(library_id)
    if 'samples' in kwargs:
        sample_id = list(kwargs['samples'].split(','))
        

    filtered_metrics = pd.DataFrame()
    filtered_reads = pd.DataFrame()
    tmp_metrics = pd.DataFrame()

    
    
    # excluded_cond = ['NCC','NTC','gDNA']
    sample_type = 'cell' # can be nucleus, depend on the nature of each sequencing library
    
    

    ## Read metrics file that contains all quality output for each cell
    if bool(re.search(".gz$", kwargs['metrics'])):
        metrics_df = pd.read_csv(kwargs['metrics'], header=0, compression='gzip',sep=',')
    else:
        metrics_df = pd.read_csv(kwargs['metrics'])
    
    
    
    ## Yaniv: I used the threshold of 0.75 to filter low quality cells, in your setting, you used 0.8. You may loose some cells from this threshold
    quality_thrs = 0.75
    filtered_metrics = metrics_df.loc[metrics_df['quality'] >= quality_thrs] #quality_thrs = 0.75
    # print('quality', filtered_metrics.shape)
    
    
    ## Yaniv: in dlp+ sequencing, different samples may be pool together, so we need to select cells from our observed sample. 
    if 'sample_id' not in filtered_metrics.columns:
        # print(filtered_metrics['cell_id'])
        filtered_metrics = filtered_metrics.assign(sample_id = filtered_metrics.apply(lambda row: row['cell_id'].split("-")[0], axis = 1))
        # print(filtered_metrics['sample_id'])
    for sample in sample_id:
        # print(sample)
        tmp_metrics = tmp_metrics.append(filtered_metrics[(filtered_metrics['sample_id'].str.startswith(sample))])
        # print(tmp_metrics.shape)
    filtered_metrics = tmp_metrics
    # print("sample filter", filtered_metrics.shape)
    
    
    ## Cell or nuclei as input???
    
    # if sample_type == "cell" or sample_type == "cells":
    #     cell_call_filter = "C1"
    # if sample_type == "nuclei":
    #     cell_call_filter = "C2"
    # cell_call_filter = "C1" # C1 is cell, C2 is nuclei
    # print(sample_type, filter_cond)
    # filtered_metrics = filtered_metrics[(filtered_metrics['cell_call'] == cell_call_filter)]
    # print('cell call filter', filtered_metrics.shape)
    # filtered_metrics = filtered_metrics[(filtered_metrics["experimental_condition"].isin(filter_cond))]
    
    
    ## experimental_condition
    ## control cells have the label as: *GM* or *NCC*, *NTC*, *gDNA*
    conds = np.unique(filtered_metrics["experimental_condition"])
    print(conds)
    conds = [st for st in conds if not is_excluded_conditions(st)]
    print("Experiment conditions that are included in downstream analysis: ")
    print(conds)
    
    
    ## Yaniv: From your filtering: "Experiment marked cells (!is_control if column exist, experimental_condition == ‘A’ otherwise)"
    ## Yaniv: there are some other cell types, not only 'A', so if you use experimental_condition == ‘A’, you loose lot of good cells, and other sub populations in your data
    print("** DEBUG before excluding control cells: **")
    print(filtered_metrics.shape)        
    filtered_metrics = filtered_metrics[filtered_metrics["experimental_condition"].isin(conds)]
    print("** DEBUG after excluding control cells: **")
    print(filtered_metrics.shape)  
    
    ## Old script: 
    # filtered_metrics = filtered_metrics[(~filtered_metrics["experimental_condition"].isin(excluded_cond))]
    # print('experimental_condition', filtered_metrics.shape)
    
    filtered_metrics = filtered_metrics['cell_id']


    ## Read the reads.csv.gz file from hmmcopy output
    if bool(re.search(".gz$", kwargs['reads'])):
        reads_df = pd.read_csv(kwargs['reads'], header=0, compression='gzip',sep=',')
    else:
        reads_df = pd.read_csv(kwargs['reads'])
    

    ## Yaniv: mappability threshold, 99% of mapping 
    mapping_thrs = 0.99
    filtered_reads = reads_df.loc[reads_df['map'] >= mapping_thrs]  #mapping_thrs=0.99
    # Hoa: modified, valid: TRUE, some bins have NA as reads count but assigned a HMMcopy state
#             filtered_reads = reads_df.loc[(reads_df['map'] >= 0.99) & (reads_df['valid'])]
    # print('mapping quality', filtered_reads.shape)

    
    ## cell predict matrix is metrics file here. 
    cell_predict_df = pd.read_csv(kwargs['metrics'], header=0, compression='gzip')
    # print(cell_predict_df.columns)
    cell_list_filtered = []
    for sample in sample_id:
        cell_list_filtered += list(filter(lambda x: re.search(sample, x), cell_predict_df['cell_id']))

    cell_predict_df = cell_predict_df[cell_predict_df['cell_id'].isin(cell_list_filtered)]
    
    ## Yaniv: Removing s-phase cells
    cell_predict_df = cell_predict_df[(cell_predict_df['is_s_phase'] == False)]
    
    
    
    ## Yaniv: from your comment: "No need to filter by number of reads - filtering by quality already makes sure the minimum reads per cell is around 125K, most of them with many more"
    ## Yaniv: I see many cases with quality > 0.75 but total mapped reads are very low, so better to include this filter. 
    
    # Hoa: add more filter condition to remove low mapped reads cells, default: 10000 
    print('Filtering low mapped reads cells with threshold: %s' % kwargs['thrs_map_reads'])
    cell_predict_df = cell_predict_df[cell_predict_df['total_mapped_reads'] >= kwargs['thrs_map_reads']]
    
    ## Yaniv: filter condition to remove mouse cells
    ## Yaniv: from your filter: Human (species == ‘grch37’), this is not correct, human cells can align to mouse reference, and have reads counts as well. 
    ## in pdx model, mouse cells can be mixed with human cells, so we need to exclude mouse cells from analysis. A cell is considered as mouse cells if the reads aligned to mouse reference greater than to human reference. 
    cell_predict_df = cell_predict_df[cell_predict_df['fastqscreen_mm10'] < cell_predict_df['fastqscreen_grch37']]
    
    not_s_phase = cell_predict_df['cell_id']

    if not (filtered_reads.empty and filtered_metrics.empty and not_s_phase.empty):
        filtered_cnv_data = filtered_reads[filtered_reads['cell_id'].isin(filtered_metrics)]
        filtered_cnv_data = filtered_cnv_data[filtered_cnv_data['cell_id'].isin(not_s_phase)]

        if filtered_cnv_data['state'].dtype == np.float64:
            filtered_cnv_data['state'] = filtered_cnv_data['state'].astype(int)
        ## Yaniv: drop NA bins, resulting in ~4375 good genomic bins for inference
        filtered_cnv_data['loci'] = filtered_cnv_data[['chr', 'start', 'end']].apply(lambda x: '_'.join(map(str, x)), axis=1)
        filtered_state_data = filtered_cnv_data.pivot(index='loci', columns='cell_id', values='state')
        filtered_state_data = filtered_state_data.dropna(axis='columns', how='all')
        filtered_state_data = filtered_state_data.dropna(axis='index', how='any')
        # print('**** Validate: Nb NaN values: {0}'.format(filtered_state_data.isnull().sum().sum()))

    filtered_state_data = filtered_state_data.reindex(index=natsorted(filtered_state_data.index))
    print("dimensions of total filtered data {0}, {1}".format(filtered_state_data.shape[0],filtered_state_data.shape[1]))
    # filtered_state_data.to_csv(kwargs['output_file'], sep=',', index_label=False)
    if not filtered_state_data.empty:
        filtered_state_data.to_csv(kwargs['output_file'], sep=',', index_label=False)
        print("Save filtered data -lib id:{0} into dir: {1}".format(library_id, kwargs['output_file']))
        # filtered_state_data_ds = downsampling(filtered_state_data, 0.2)
        # print("Downsampling filtered data {0}, {1}".format(filtered_state_data_ds.shape[0],filtered_state_data_ds.shape[1]))
        # filtered_state_data_ds.to_csv(kwargs['output_file'], sep=',', index_label=False)




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--metrics', required=True)
    parser.add_argument('--reads', required=True)
    # parser.add_argument('--cell_state_predict', required=False, default=None)
    parser.add_argument('--library_id', required=True)
    parser.add_argument('--samples', required=True)
    parser.add_argument('--output_file', required=True)
    parser.add_argument('--thrs_map_reads', required=False, default='10000') # you can change to 50K, 100K,..
    args = vars(parser.parse_args())


    library_id = args['library_id']
    samples = args['samples'] # sample id
    metrics = args['metrics'] # metric file directory
    reads = args['reads'] # read file directory
    # cell_state_predict = args['cell_state_predict']
    output_file = args['output_file']
    thrs_map_reads = int(args['thrs_map_reads']) # usually 250000, some special case: SA535: 10000, or SA919: 500000

    # print(output_file)
    print('library_id {0}'.format(library_id))
    print('samples {0}'.format(samples))
    print('thrs_map_reads {0}'.format(thrs_map_reads))
    filter_libs(library_id=library_id, samples=samples, metrics=metrics, reads=reads, output_file=output_file, thrs_map_reads=thrs_map_reads) #cell_state_predict=cell_state_predict,
    
    
                
                
## How to run this function from command line mode: 
## /home/htran/anaconda3/envs/sisyphus/bin/python /home/htran/Projects/hakwoo_project/corrupt_tree/src/filter_scripts/filter_libs_DLP_HoaTran.py --library_id A98166B --samples AT11391 --metrics /home/htran/storage/raw_DLP/metastasis_DLP/SA919_mixing_exp/A98166B/annotation/A98166B_metrics.csv.gz --reads /home/htran/storage/raw_DLP/metastasis_DLP/SA919_mixing_exp/A98166B/hmmcopy/A98166B_reads.csv.gz --output_file /home/htran/storage/datasets/metastasis_results/dlp_SA919_mixing_exp/A98166B_demo_filtered_states.csv --thrs_map_reads 10000
