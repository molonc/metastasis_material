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
# from dbclients.colossus import ColossusApi

# colossus_api = ColossusApi()

# def downsampling(data_df, percentds = 0.3):
#     total_cells = data_df.columns
#     cells_ext = np.random.choice(total_cells,size=int(len(total_cells)*percentds), replace=False)
#     print('Extract {} cells from {} total cells: '.format(len(cells_ext),len(total_cells)))
#     return data_df.loc[:,cells_ext]


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
    # quality_thrs = 0.75
    # mapping_thrs = 0.99
    quality_thrs = float(kwargs['quality_thrs']) #0.75
    mapping_thrs = float(kwargs['mapping_thrs'])#0.99
    print('Quality thrs: {0}, mapping_thrs: {1}'.format(quality_thrs, mapping_thrs))
    
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

    
    
    sample_type = 'cell' # can be nucleus, depend on the nature of each sequencing library
    
    # print(cell_call_filter)
    if bool(re.search(".gz$", kwargs['metrics'])):
        metrics_df = pd.read_csv(kwargs['metrics'], header=0, compression='gzip',sep=',')
    else:
        metrics_df = pd.read_csv(kwargs['metrics'])
    
    
    filtered_metrics = metrics_df.loc[metrics_df['quality'] >= quality_thrs] #quality_thrs = 0.75
    # print('quality', filtered_metrics.shape)
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
    cell_call_filter = "C1" # C1 is cell, C2 is nuclei
    # print(sample_type, filter_cond)
    # filtered_metrics = filtered_metrics[(filtered_metrics['cell_call'] == cell_call_filter)]
    # print('cell call filter', filtered_metrics.shape)
    # filtered_metrics = filtered_metrics[(filtered_metrics["experimental_condition"].isin(filter_cond))]
    
    
    ## experimental_condition
    conds = np.unique(filtered_metrics["experimental_condition"])
    print(conds)
    conds = [st for st in conds if not is_excluded_conditions(st)]
    print(conds)

    print("** DEBUG before excluding control cells: **")
    print(filtered_metrics.shape)        
    filtered_metrics = filtered_metrics[filtered_metrics["experimental_condition"].isin(conds)]
    print("** DEBUG after excluding control cells: **")
    print(filtered_metrics.shape)  
    
    ## Old script: 
    # filtered_metrics = filtered_metrics[(~filtered_metrics["experimental_condition"].isin(excluded_cond))]
    # print('experimental_condition', filtered_metrics.shape)
    
    filtered_metrics = filtered_metrics['cell_id']

    if bool(re.search(".gz$", kwargs['reads'])):
        reads_df = pd.read_csv(kwargs['reads'], header=0, compression='gzip',sep=',')
    else:
        reads_df = pd.read_csv(kwargs['reads'])
    

    # old version but quite ok
    filtered_reads = reads_df.loc[reads_df['map'] >= mapping_thrs]  #mapping_thrs=0.99
    # Hoa: modified, valid: TRUE, some bins have NA as reads count but assigned a HMMcopy state
#             filtered_reads = reads_df.loc[(reads_df['map'] >= 0.99) & (reads_df['valid'])]
    # print('mapping quality', filtered_reads.shape)

    
    # if StrictVersion(version.strip('v')) < StrictVersion('0.3.1'):
    #     cell_predict_df = pd.read_csv(kwargs['cell_state_predict'], compression='gzip')
    #     # print('DEBUG__version of pipeline: {0}'.format(version.strip('v')))
    #     print(kwargs['cell_state_predict'])
    # else:
    #     cell_predict_df = pd.read_csv(kwargs['metrics'], compression='gzip')
    #     print(kwargs['metrics'])
    cell_predict_df = pd.read_csv(kwargs['metrics'], header=0, compression='gzip')
    # print(cell_predict_df.columns)
    cell_list_filtered = []
    for sample in sample_id:
        cell_list_filtered += list(filter(lambda x: re.search(sample, x), cell_predict_df['cell_id']))

    cell_predict_df = cell_predict_df[cell_predict_df['cell_id'].isin(cell_list_filtered)]
    
    ## Removing s-phase cells
    cell_predict_df = cell_predict_df[(cell_predict_df['is_s_phase'] == False)]
    
    # Hoa: add more filter condition to remove low mapped reads cells, default: 10000 or 500000
    print('Filtering low mapped reads cells with threshold: %s' % kwargs['thrs_map_reads'])
    cell_predict_df = cell_predict_df[cell_predict_df['total_mapped_reads'] >= kwargs['thrs_map_reads']]
    
    # Hoa: add more filter condition to remove mouse cells
    cell_predict_df = cell_predict_df[cell_predict_df['fastqscreen_mm10'] < cell_predict_df['fastqscreen_grch37']]
    
    not_s_phase = cell_predict_df['cell_id']

    if not (filtered_reads.empty and filtered_metrics.empty and not_s_phase.empty):
        filtered_cnv_data = filtered_reads[filtered_reads['cell_id'].isin(filtered_metrics)]
        filtered_cnv_data = filtered_cnv_data[filtered_cnv_data['cell_id'].isin(not_s_phase)]

        if filtered_cnv_data['state'].dtype == np.float64:
            filtered_cnv_data['state'] = filtered_cnv_data['state'].astype(int)

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
    # parser.add_argument('--metrics', required=True)
    # parser.add_argument('--reads', required=True)
    # # parser.add_argument('--cell_state_predict', required=False, default=None)
    # parser.add_argument('--library_id', required=True)
    # parser.add_argument('--samples', required=True)
    # parser.add_argument('--output_file', required=True)
    parser.add_argument('--meta_fn', required=True) # meta sample csv file, containing library_id, and sample_id columns
    parser.add_argument('--input_dir', required=True)  #ex: yourDir where you keep hmmcopy files for a library id:YourDir/ inside this dir contains library id folder: ex: A98166B/
    parser.add_argument('--output_dir', required=True) # where you want to save output into
    parser.add_argument('--quality_thrs', required=False, default='0.75')
    parser.add_argument('--mapping_thrs', required=False, default='0.99')
    parser.add_argument('--thrs_map_reads', required=False, default='250000')
    args = vars(parser.parse_args())
    # 
    # 
    # library_id = args['library_id']
    # samples = args['samples'] # sample id
    # metrics = args['metrics'] # metric file directory
    # reads = args['reads'] # read file directory
    # # cell_state_predict = args['cell_state_predict']
    # output_file = args['output_file']
    thrs_map_reads = int(args['thrs_map_reads']) # usually 250000, some special case: SA535: 10000, or SA919: 500000
    quality_thrs = float(args['quality_thrs']) #0.75
    mapping_thrs = float(args['mapping_thrs'])#0.99
    input_dir = args['input_dir']
    output_dir = args['output_dir']
    meta_fn = args['meta_fn']
    # thrs_map_reads = 500000 # For SA919
    # quality_thrs = 0.75
    # mapping_thrs = 0.99    
    # 
    # # print(output_file)
    # print('library_id {0}'.format(library_id))
    # print('samples {0}'.format(samples))
    # print('thrs_map_reads {0}'.format(thrs_map_reads))
    # filter_libs(library_id=library_id, samples=samples, metrics=metrics, reads=reads, output_file=output_file, thrs_map_reads=thrs_map_reads) #cell_state_predict=cell_state_predict, 
    ## For many library ids, load library id, and sample id from command line
    # meta_fn = '/home/htran/Projects/hakwoo_project/metastasis_material/materials/dlp_trees/SA919/SA919_passageX8.csv'
    
    ## meta sample file with library_id, and sample_id as column name
    # meta_fn = '/home/htran/Projects/hakwoo_project/metastasis_material/materials/dlp_trees/SA919_mixing_experiment/library_groupings_extra.csv'
    # input_dir = '/home/htran/storage/raw_DLP/metastasis_DLP/SA919/'
    # save_dir = '/home/htran/storage/datasets/metastasis_results/SA919_X8/'
    # output_dir = '/home/htran/storage/datasets/metastasis_results/dlp_SA919_mixing_exp/'
    
    if not os.path.exists(output_dir): 
        os.makedirs(output_dir)
        
    
    
    lib_df = pd.read_csv(meta_fn, header=0)
    print(lib_df)
    # lib_df['library_id'] = lib_df['grouping'].values
    for i in range(lib_df.shape[0]): #
          library_id  = lib_df['library_id'][i]
          samples  = lib_df['sample_id'][i]
          print('library id: {0}'.format(library_id))
          print('samples: {0}'.format(samples))
          flag = True
          metrics_fn1 = os.path.join(input_dir,library_id,'annotation','metrics.csv.gz')
          metrics_fn2 = os.path.join(input_dir,library_id,'annotation',library_id + '_metrics.csv.gz')
          if os.path.isfile(metrics_fn1):
              metrics = metrics_fn1
          else:
             if os.path.isfile(metrics_fn2):
                 metrics = metrics_fn2
             else:
                 flag = False
                 print("Do not exist metrics file, check input file")
                 
          
          print(metrics)                  
          reads_fn1 = os.path.join(input_dir,library_id,'hmmcopy','reads.csv.gz')
          reads_fn2 = os.path.join(input_dir,library_id,'hmmcopy',library_id +'_reads.csv.gz')
          if os.path.isfile(reads_fn1):
              reads = reads_fn1
          else:
             if os.path.isfile(reads_fn2):
                 reads = reads_fn2
             else:
                 flag = False
                 print("Do not exist reads file, check input file")
                 
          print(reads)                                   
          
          if flag:
            output_file = os.path.join(output_dir,library_id+'_filtered_states.csv.gz')
            filter_libs(library_id=library_id, samples=samples, metrics=metrics,
                        reads=reads, output_file=output_file,
                        thrs_map_reads=thrs_map_reads, quality_thrs=quality_thrs, mapping_thrs=mapping_thrs)
                
                
                
## How to run this function from command line mode: 
# /home/htran/anaconda3/envs/sisyphus/bin/python /home/htran/Projects/hakwoo_project/metastasis_material/scripts/corrupt_tree/src/filter_scripts/filter_libraries_from_list.py --input_dir /home/htran/storage/raw_DLP/metastasis_DLP/SA919/ --thrs_map_reads 500000 --meta_fn /home/htran/Projects/hakwoo_project/metastasis_material/materials/dlp_trees/SA919_mixing_experiment/library_groupings_extra.csv --output_dir /home/htran/storage/datasets/metastasis_results/dlp_SA919_mixing_exp/
