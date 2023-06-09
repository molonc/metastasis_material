import sys
import os
# import re
# import gzip
import argparse
# import subprocess

import pandas as pd
import numpy as np
from datetime import datetime
import logging
from logging.handlers import TimedRotatingFileHandler

import utils

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--expr_fn', required=True)
    parser.add_argument('--cnv_fn', required=True)
    parser.add_argument('--cell_clones_fn', required=False, default=None)
    parser.add_argument('--tree_fn', required=False, default=None)
    parser.add_argument('--TreeAlign_mode', required=False, default='clonal_infer') #TreeAlign_mode can be 'clonal_infer' or 'tree_infer'
    parser.add_argument('--datatag', required=False, default='SA')
    parser.add_argument('--assign_prob', required=False, default='0.6') # from 0.5 to 1, 0.6 means 60% of gene follow CN to assign a clonal label
    parser.add_argument('--iterations', required=False, default='500')
    parser.add_argument('--save_dir', required=False, default=None) # if None, save data into input folder
    
    
    args = vars(parser.parse_args())

    # (expr_fn, cnv_fn, cell_clones_fn=None, tree_fn=None)
    # min_clone_assign_prob=0.8, max_iter=400

    expr_fn = args['expr_fn']
    cnv_fn = args['cnv_fn'] 
    cell_clones_fn = args['cell_clones_fn'] 
    tree_fn = args['tree_fn']
    TreeAlign_mode = args['TreeAlign_mode']
    datatag = args['datatag']
    min_clone_assign_prob = float(args['assign_prob']) 
    iterations = int(args['iterations'])
    save_dir = args['save_dir']
    # print(output_file)
    print('expr_fn {0}'.format(expr_fn))
    print('cnv_fn {0}'.format(cnv_fn))
    print('cell_clones_fn {0}'.format(cell_clones_fn))
    print('TreeAlign_mode {0}'.format(TreeAlign_mode))
    print('min_clone_assign_prob {0} iterations: {1}'.format(min_clone_assign_prob, iterations))
    
    if save_dir==None: 
        save_dir = os.path.join(os.path.dirname(expr_fn),'TreeAlign_output')
        
    if not os.path.exists(save_dir): 
        os.makedirs(save_dir)
        print('Creating an output folder')
    
    print(save_dir)
    
    start = datetime.now()
    
    input_dic = utils.load_input_data(expr_fn, cnv_fn, cell_clones_fn, tree_fn)
        
    if cell_clones_fn is not None and TreeAlign_mode=='clonal_infer': 
        print('Run TreeAlign using clonal inference mode')
        utils.run_clone_TreeAlign(input_dic["expr"], input_dic["cnv"], input_dic["clone"], save_dir, datatag,
                        min_clone_assign_prob=min_clone_assign_prob, max_iter=iterations, min_clone_cell_count=1, cnv_cutoff=8)
    else:
        print('Run TreeAlign using tree phylogenetic as input for inference')
        utils.run_phylo_TreeAlign(expr, cnv, tree, save_dir, datatag)
        
    end = datetime.now()
    td = (end - start).total_seconds() / 60 
#     print(f"The time of execution of above program is : {td:.03f} mins")  
    log_time = f"The time of execution of above program is : {td:.03f} mins"
    print(log_time)
    mode_desc = 'TreeAlign_mode: ' + TreeAlign_mode
    prob_desc = 'min_clone_assign_prob: ' + str(min_clone_assign_prob)
    ite_desc = 'Iterations: ' + str(iterations)
    log_fn = os.path.join(save_dir, datatag+'_log_params.txt')
    h = TimedRotatingFileHandler(log_fn)
    
    logging.basicConfig(filename=os.path.join(save_dir, datatag+'_log_params.txt'), level=logging.DEBUG, filemode='w')#encoding='utf-8', 
    logging.info(expr_fn)
    logging.info(cnv_fn)
    logging.info(datatag)
    logging.info(mode_desc)
    logging.info(prob_desc)
    logging.info(ite_desc)
    logging.info(save_dir)
    logging.info(log_time)


# /home/htran/anaconda3/envs/sisyphus/bin/python exe_TreeAlign.py --expr_fn /home/htran/Projects/farhia_project/rnaseq/method_testing/TreeAlign/data/example_expr.csv --cnv_fn /home/htran/Projects/farhia_project/rnaseq/method_testing/TreeAlign/data/example_gene_cnv.csv --cell_clones_fn /home/htran/Projects/farhia_project/rnaseq/method_testing/TreeAlign/data/example_cell_clone.csv --datatag test --save_dir /home/htran/Projects/farhia_project/rnaseq/method_testing/TreeAlign/data/test_pipeline/


# /home/htran/anaconda3/envs/sisyphus/bin/python /home/htran/Projects/farhia_project/rnaseq/method_testing/TreeAlign/run_TreeAlign/exe_TreeAlign.py --expr_fn /home/htran/storage/rnaseq_datasets/testing_space/clonealign2/SA535X4XB05462/SA535X4XB05462_expr_introns.csv.gz --cnv_fn /home/htran/storage/rnaseq_datasets/testing_space/clonealign2/SA535X4XB05462/SA535X4XB05462_clones_cnv.csv.gz --cell_clones_fn /home/htran/storage/rnaseq_datasets/testing_space/clonealign2/SA535X4XB05462/SA535X4XB05462_cell_clones.csv.gz --datatag SA535X4XB05462
