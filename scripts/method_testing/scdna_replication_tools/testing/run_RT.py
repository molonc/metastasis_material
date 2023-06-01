import numpy as np
import pandas as pd
import os
import re
from datetime import datetime
from argparse import ArgumentParser
import statsmodels.api as sm
import logging
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
matplotlib.use('Agg') # for terminal mode


import scdna_replication_tools
## Some issues with import funcs from scdna_replication_tools
## to fix it temporarilly, I import all funcs here
from scdna_replication_tools.cncluster import kmeans_cluster
from scdna_replication_tools.compute_consensus_clone_profiles import compute_consensus_clone_profiles
from scdna_replication_tools.assign_s_to_clones import assign_s_to_clones
from scdna_replication_tools.bulk_gc_correction import bulk_g1_gc_correction
from scdna_replication_tools.normalize_by_cell import normalize_by_cell
from scdna_replication_tools.normalize_by_clone import normalize_by_clone
from scdna_replication_tools.binarize_rt_profiles import binarize_profiles
from scdna_replication_tools.compute_pseudobulk_rt_profiles import compute_pseudobulk_rt_profiles
from scdna_replication_tools.calculate_twidth import compute_time_from_scheduled_column, calculate_twidth
from scdna_replication_tools.pert_model import pyro_infer_scRT
from scdna_replication_tools.infer_scRT import scRT #main func

## utility functions
import RT_stat_utils



## cell_clones_fn: a dataframe that contains columns: 'cell_id','clone_id', 'library_id', 'cell_type_status': s or g cell type
## input_gfn: : a dataframe that contains columns: 'cell_id', 'chr', 'start', 'end', argv.gc_col, argv.cn_col, argv.input_col
## input_sfn: a dataframe that contains columns: 'cell_id', 'chr', 'start', 'end', argv.gc_col, argv.cn_col, argv.input_col
## 
def load_data(input_gfn, input_sfn, cell_clones_fn, argv=None):
    input_dir = os.path.dirname(input_gfn)
#     input_gfn = os.path.join(input_dir, 'hmmcopy/filtered_reads_RT_g_cells.csv.gz') 
#     input_sfn = os.path.join(input_dir, 'hmmcopy/filtered_reads_RT_s_cells.csv.gz') 
    if bool(re.search(".gz$", input_gfn)):
        cn_g = pd.read_csv(input_gfn, compression='gzip')
    else:
        cn_g = pd.read_csv(input_gfn) #header=0, index_col=0

    if bool(re.search(".gz$", input_sfn)):
        cn_s = pd.read_csv(input_sfn, compression='gzip')
    else:
        cn_s = pd.read_csv(input_sfn) #header=0, index_col=0

    print('G cells: ')
    print(cn_g.shape)
    print('S cells: ')
    print(cn_s.shape)
#     cell_clones_fn = os.path.join(input_dir, 'RT_input/A130854B_filtered_cell_clones.csv') 
    cell_clones = pd.read_csv(cell_clones_fn)
    print('Cell clones file: ')
    print(cell_clones.columns)
    print(cell_clones.shape)
#     select_values = cn_g['cell_id'].values
#     cells_clone_cut = cell_clones.query('cell_id in @select_values')  
# #     cells_clone_cut.index = cells_clone_cut['cell_id'].values
#     # cells_clone_cut.loc['AT13696-A130854B-R48-C07',:]
#     cells_clone_cut = cells_clone_cut.loc[:,['cell_id','clone_id', 'library_id']]
    cell_clones = cell_clones.loc[:,['cell_id','clone_id', 'library_id']]
    # included in input files
#     cn_s['library_id'] = np.repeat(libary_id,cn_s.shape[0])
#     cn_g['library_id'] = np.repeat(libary_id,cn_g.shape[0])
    cn_g = pd.merge(cn_g, cell_clones, on='cell_id')
    cn_s = pd.merge(cn_s, cell_clones, on='cell_id')
#     print('G cells: ')
#     print(cn_g.shape)
#     print('S cells: ')
#     print(cn_s.shape)
    
    if argv==None:# columns names definition
        argv = pd.Series(['reads', 'state', 'gc','clone_id','hmmcopy',
                 'model_cn_state','model_rep_state','reads'], 
                 index=['input_col','cn_col','gc_col','clone_col','cn_prior_method',
                       'frac_rt_col','rep_col','rpm_col'])
    # temp_cn_s = cn_s[['cell_id', 'chr', 'start', 'end', 'gc', 'state', 'library_id', 'reads']]
#     temp_cn_g = cn_g1[['cell_id', 'chr', 'start', 'end', 'gc', 'state', 'library_id', 'clone_id','reads']]
    temp_cn_s = cn_s[['cell_id', 'chr', 'start', 'end', argv.gc_col, argv.cn_col, 'library_id', argv.input_col]]
    temp_cn_g = cn_g[['cell_id', 'chr', 'start', 'end', argv.gc_col, argv.cn_col, 'library_id', argv.clone_col, argv.input_col]]
    print('G cells input: {0} {1}'.format(temp_cn_g.shape[0],temp_cn_g.shape[1]))
    print('S cells input: {0} {1}'.format(temp_cn_s.shape[0],temp_cn_s.shape[1]))
    
    return input_dir, temp_cn_g, temp_cn_s


def run_RT(datatag, input_dir, temp_cn_g, temp_cn_s, output_dir=None, nb_iterations = 500, argv=None):
    if not os.path.exists(input_dir): 
        os.makedirs(input_dir) 
    datatag = datatag.replace(' ', '_')
    if output_dir==None:
        output_dir = os.path.join(input_dir,'RT_results/')
    if not os.path.exists(output_dir): 
        os.makedirs(output_dir) 
        print('Creating an output dir: {0}'.format(output_dir))
        
#     nb_iterations = 10 #1500
    if argv==None:  # columns names definition
        argv = pd.Series(['reads', 'state', 'gc','clone_id','hmmcopy',
                 'model_cn_state','model_rep_state','reads'], 
                 index=['input_col','cn_col','gc_col','clone_col','cn_prior_method',
                       'frac_rt_col','rep_col','rpm_col'])
    print('Creating scrt object')
    # create SPF object with input
    scrt = scRT(temp_cn_s, temp_cn_g, input_col=argv.input_col, rt_prior_col=None, assign_col=argv.cn_col,
                cn_state_col=argv.cn_col, gc_col=argv.gc_col, cn_prior_method=argv.cn_prior_method, max_iter=nb_iterations)
    # run inference
    print('Running inference')
    cn_s_with_scrt, supp_s_output, cn_g_with_scrt, supp_g_output = scrt.infer_pyro_model()
    
    # cn_s_with_scrt, supp_s_output, cn_g_with_scrt, supp_g_output
    cn_s_with_scrt.to_csv(os.path.join(output_dir,datatag+'_cn_s_with_scrt.csv.gz'), index=False)
    supp_s_output.to_csv(os.path.join(output_dir,datatag+'_supp_s_output.csv.gz'), index=False)
    cn_g_with_scrt.to_csv(os.path.join(output_dir,datatag+'_cn_g_with_scrt.csv.gz'), index=False)
    supp_g_output.to_csv(os.path.join(output_dir,datatag+'_supp_g_output.csv.gz'), index=False)
    print('Save results into output dir: {0}'.format(output_dir))
    clonal_df = get_assigned_clone_ids(cn_s_with_scrt, cn_g_with_scrt, output_dir, datatag)
    return output_dir, clonal_df, cn_g_with_scrt, cn_s_with_scrt


## Assigned clone ids for s-phase cells, and clone ids from tree inference for g cells
def get_assigned_clone_ids(cn_s_with_scrt, cn_g_with_scrt, output_dir, datatag):
    
    cn_s_with_scrt_extracted = cn_s_with_scrt.loc[:,['cell_id','clone_id','library_id']] #,'cell_type_status'
    cn_g_with_scrt_extracted = cn_g_with_scrt.loc[:,['cell_id','clone_id','library_id']] #,'cell_type_status'
    cn_s_with_scrt_extracted = cn_s_with_scrt_extracted.groupby(['cell_id','clone_id','library_id'], as_index=False).size()
    cn_g_with_scrt_extracted = cn_g_with_scrt_extracted.groupby(['cell_id','clone_id','library_id'], as_index=False).size()
    
    cn_s_with_scrt_extracted['cell_type_status'] = np.repeat('cn_s',cn_s_with_scrt_extracted.shape[0])
    cn_g_with_scrt_extracted['cell_type_status'] = np.repeat('cn_g',cn_g_with_scrt_extracted.shape[0])
    
    cn_s_with_scrt_extracted = cn_s_with_scrt_extracted.drop(columns=['size'])
    cn_g_with_scrt_extracted = cn_g_with_scrt_extracted.drop(columns=['size'])
    
    clonal_df = pd.concat([cn_s_with_scrt_extracted, cn_g_with_scrt_extracted], ignore_index=True)
    clonal_df.to_csv(os.path.join(output_dir,datatag+'_clonal_RT.csv.gz'), index=False)
    print('Save assigned clones for s-phase cells, and g cells into output dir: {0}'.format(output_dir))
    return(clonal_df)
    
    
def save_images(basename, output_dir): 
#     import matplotlib.pyplot as plt
    if not os.path.exists(output_dir): os.makedirs(output_dir) 
    outname = os.path.join(output_dir, basename + '.png') 
    plt.savefig(outname, dpi=150) 
    plt.close()    

    
def summary_results(output_dir, datatag, cn_g_with_scrt, cn_s_with_scrt, cell_clones, argv=None):
    
    argv=None
    if argv==None:  # columns names definition
            argv = pd.Series(['reads', 'state', 'gc','clone_id','hmmcopy',
                     'model_cn_state','model_rep_state','reads'], 
                     index=['input_col','cn_col','gc_col','clone_col','cn_prior_method',
                           'frac_rt_col','rep_col','rpm_col'])

    
    cn = pd.concat([cn_s_with_scrt, cn_g_with_scrt], ignore_index=True)

    # compute the fraction of replicated bins within each cell
    # cn = compute_cell_frac(cn, frac_rt_col=argv.frac_rt_col, rep_state_col=argv.rep_col)

    # compute autocorrelation features to see which cells are truly low quality
    cell_metrics = RT_stat_utils.compute_quality_features(cn,  rpm_col=argv.rpm_col) #rep_state_col=argv.rep_col, cn_state_col=argv.cn_col,
    cell_metrics_RT_frac = RT_stat_utils.compute_cell_frac_v2(cn, frac_rt_col='cell_frac_rep', rep_state_col='model_rep_state')
    # cell_clones = cell_clones.loc[:,['cell_id','clone_id','cell_type_status']] #,'cell_type_status'
    
    # cell_metrics.to_csv(os.path.join(output_dir, datatag+'_RT_cell_metrics.csv.gz'), index=False)
    # cell_metrics_RT_frac.to_csv(os.path.join(output_dir, datatag+'_RT_cell_metrics_RT_frac.csv.gz'), index=False)
    print('Save summary results into output dir: {0}'.format(output_dir))
    
    cell_metrics = pd.merge(cell_metrics, cell_clones, on='cell_id')
    cell_metrics_RT_frac = pd.merge(cell_metrics_RT_frac, cell_clones, on='cell_id')
    cell_metrics.to_csv(os.path.join(output_dir, datatag+'_RT_cell_metrics_clones.csv.gz'), index=False)
    cell_metrics_RT_frac.to_csv(os.path.join(output_dir, datatag+'_RT_cell_metrics_RT_frac_clones.csv.gz'), index=False)
    
    p = sns.boxplot(data=cell_metrics, x="rep_bk", y="cell_type_status", hue="clone_id") #"rep_auto"
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=5)
    outname = os.path.join(output_dir, datatag + '_rep_bk.png') 
    plt.savefig(outname, dpi=150)
    plt.close()    

    p1 = sns.boxplot(data=cell_metrics_RT_frac, x="cell_frac_rep", y="cell_type_status", hue="clone_id") #"rep_auto"
    p1.axvline(0.05)
    p1.axvline(0.95)
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=5)
    outname1 = os.path.join(output_dir, datatag + '_cell_frac_rep.png') 
    plt.savefig(outname1, dpi=150)
    plt.close()   



## cell_clones_fn: a dataframe that contains columns: 'cell_id','clone_id', 'library_id'
## input_gfn: : a dataframe that contains columns: 'cell_id', 'chr', 'start', 'end', argv.gc_col, argv.cn_col, argv.input_col
## input_sfn: a dataframe that contains columns: 'cell_id', 'chr', 'start', 'end', argv.gc_col, argv.cn_col, argv.input_col
## 

if __name__ == '__main__':
    parser = ArgumentParser()  #argparse.ArgumentParser()
    parser.add_argument('--input_gfn', required=True)
    parser.add_argument('--input_sfn', required=True)
    parser.add_argument('--cell_clones_fn', required=True)
    parser.add_argument('--output_dir', required=False, default=None) # If None, using input directory from input files 
    parser.add_argument('--output_fn', required=False, default=None)
    parser.add_argument('--datatag', required=False, default='SA')
    parser.add_argument('--nb_iterations', required=False, default='5') # 5 for testing, in general, using 1500 iterations
    args = vars(parser.parse_args())

    start = datetime.now()
    print('input_gfn {0}'.format(args['input_gfn']))
    print('input_sfn {0}'.format(args['input_sfn']))
    print('cell_clones_fn {0}'.format(args['cell_clones_fn']))
    print('output_dir {0}'.format(args['output_dir']))
    print('datatag {0}'.format(args['datatag']))
    print('nb_iterations {0}'.format(args['nb_iterations']))
    input_gfn = args['input_gfn']
    input_sfn = args['input_sfn']
    cell_clones_fn = args['cell_clones_fn']
    nb_iterations = args['nb_iterations']
    output_dir = args['output_dir']
    datatag = args['datatag']
    
    input_dir, temp_cn_g, temp_cn_s = load_data(input_gfn, input_sfn, cell_clones_fn)
    
    output_dir, clonal_df, cn_g_with_scrt, cn_s_with_scrt = run_RT(datatag, input_dir, temp_cn_g, temp_cn_s, output_dir, int(nb_iterations))
    
    # summary_results(output_dir, datatag, cn_g_with_scrt, cn_s_with_scrt, clonal_df)

    end = datetime.now()
    td = (end - start).total_seconds() / 60 
    print(f"The time of execution of above program is : {td:.03f} mins")
    # log_fn=os.path.join(output_dir, datatag+'_log.txt')
    # print(log_fn)
    logging.basicConfig(filename='log.txt', encoding='utf-8', level=logging.DEBUG, filemode='w')
    logging.info(input_gfn)
    logging.info(input_sfn)
    logging.info(cell_clones_fn)
    logging.info(nb_iterations)
    logging.info(datatag)
    logging.info(output_dir)
    logging.info(f"The time of execution of replication timing program is : {td:.03f} mins")
    
    ##Ex:

    ## /home/htran/anaconda3/envs/myPython37/bin/python /home/htran/Projects/hakwoo_project/testing_methods/scdna_replication_tools/testing/run_RT.py --input_gfn /home/htran/storage/raw_DLP/metastasis_DLP/SA919/A130854B/RT_input/filtered_reads_RT_g_cells.csv.gz --input_sfn /home/htran/storage/raw_DLP/metastasis_DLP/SA919/A130854B/RT_input/filtered_reads_RT_s_cells.csv.gz --cell_clones_fn /home/htran/storage/raw_DLP/metastasis_DLP/SA919/A130854B/RT_input/A130854B_filtered_cell_clones_v2.csv --datatag SA919 --nb_iterations 3
    ## /home/htran/anaconda3/envs/myPython37/bin/python /home/htran/Projects/hakwoo_project/testing_methods/scdna_replication_tools/testing/run_RT.py --input_gfn /home/htran/storage/datasets/metastasis_results/replication_timing/filtered_reads_RT_g_cells.csv.gz --input_sfn /home/htran/storage/datasets/metastasis_results/replication_timing/filtered_reads_RT_s_cells.csv.gz --cell_clones_fn /home/htran/storage/datasets/metastasis_results/replication_timing/SA535X4XB05649_total_cells_clones.csv.gz --datatag SA535X4XB05649 --nb_iterations 500
    ## --cell_clones_fn /home/htran/storage/raw_DLP/metastasis_DLP/SA919/A130854B/RT_input/A130854B_filtered_cell_clones_v2.csv
    
    
    
