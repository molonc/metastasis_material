import leidenalg as la
import igraph as ig
import os
import numpy as np
import pandas as pd 
import re
import argparse

def generate_clone_label(nbclones):
    test_list = [] 
    alpha = 'A'
    for i in range(0, nbclones): 
        test_list.append(alpha) 
        alpha = chr(ord(alpha) + 1) 
    return test_list
    


def get_cluster(output_fn, g, r, input_file, ptype=None, niter=20):
  
    input_dir = os.path.dirname(input_file)
    # There are many configurations - algos but RB is the best here
    # CPMVertexPartition, SignificanceVertexPartition, SurpriseVertexPartition, ModularityVertexPartition
    if ptype==None:
        ptype = la.RBConfigurationVertexPartition
        
    print(ptype)
    part = la.find_partition(g, partition_type=ptype, resolution_parameter=r, n_iterations=niter)
#     if tag=='MOD' or tag=='SUR':
#         part = la.find_partition(g, partition_type=ptype, n_iterations=niter)
#     else:
#         part = la.find_partition(g, partition_type=ptype, resolution_parameter=r, n_iterations=niter)
        
    # store output into adata.obs
    groups = np.array(part.membership)
    
    df = pd.DataFrame({'cell_id': part.graph.vs['name'],'clone_label': groups})
    locus = [l for l in df['cell_id'].values if re.search('locus_',l)]
    print(len(locus))
    root = [l for l in df['cell_id'].values if re.search('root',l)]
    df = df.loc[~df['cell_id'].isin(locus),:]
    df = df.loc[~df['cell_id'].isin(root),:]
    unique, counts = np.unique(df['clone_label'].values, return_counts=True)
    
    print('resolution %s and output of graph cut is: ' % r)
    print(dict(zip(unique,counts)))
    nbclones = len(unique)
    cid = generate_clone_label(nbclones)
    meta = pd.DataFrame({'clone_label': unique,'clone_id': cid})
    print(meta.shape)
    
    result = pd.merge(df, meta, on='clone_label',how='inner')
    print(result.shape)
    result.to_csv(os.path.join(input_dir,'cell_clones_r'+str(r)+'_nbclones_'+str(nbclones)+'.csv'), index=False)
    result.to_csv(output_fn, index=False)
#     print(result.head(10))


if __name__ == '__main__':
  
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', required=True)
    parser.add_argument('--resolution', required=True)
    parser.add_argument('--output_file', required=True)
    args = vars(parser.parse_args())
    
    input_file = args['input_file']
    resolution = np.float16(args['resolution'])
    output_file = args['output_file']
    g = ig.Graph()
    g = g.Read_GraphML(input_file)
    get_cluster(output_file, g, resolution, input_file)
    
    # You can try different grain - resolution from 0 to 1, depend on dataset
    # res = [0.01, 0.02, 0.03]
    # for r in res:
    #     get_cluster(g, r, input_dir)
