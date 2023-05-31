"""
Hoa Tran
Based on bbknn ridge regression demo and ridge class of sklearn
Ref 1: https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.Ridge.html
Ref 2: BBKNN: https://github.com/Teichlab/bbknn
"""
import argparse
import os
import pandas as pd
import numpy as np
import scipy
import sys
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
from packaging import version
from sklearn.linear_model import Ridge
import scanpy as sc 
import re
try:
	from scanpy import logging as logg
except ImportError:
	pass
try:
	import anndata
except ImportError:
	pass


def ridge_regression(adata, batch_key, confounder_key=[], chunksize=1e8, copy=False, **kwargs):
	'''
	Perform ridge regression on scaled expression data, accepting both technical and 
	biological categorical variables. The effect of the technical variables is removed 
	while the effect of the biological variables is retained. This is a preprocessing 
	step that can aid BBKNN integration `(Park, 2020) <https://science.sciencemag.org/content/367/6480/eaay3224.abstract>`_.
	
	Alters the object's ``.X`` to be the regression residuals, and creates ``.layers['X_explained']`` 
	with the expression explained by the technical effect.
	
	Input
	-----
	adata : ``AnnData``
		Needs scaled data in ``.X``.
	batch_key : ``list``
		A list of categorical ``.obs`` columns to regress out as technical effects.
	confounder_key : ``list``, optional (default: ``[]``)
		A list of categorical ``.obs`` columns to retain as biological effects.
	chunksize : ``int``, optional (default: 1e8)
		How many elements of the expression matrix to process at a time. Potentially useful 
		to manage memory use for larger datasets.
	copy : ``bool``, optional (default: ``False``)
		If ``True``, return a copy instead of writing to the supplied adata.
	kwargs
		Any arguments to pass to `Ridge <https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.Ridge.html>`_.
	'''
	start = logg.info('computing ridge regression')
	adata = adata.copy() if copy else adata
	#just in case the arguments are not provided as lists, convert them to such
	#as they need to be lists for downstream application
	if not isinstance(batch_key, list):
		batch_key = [batch_key]
	if not isinstance(confounder_key, list):
		confounder_key = [confounder_key]
	
	#construct a helper representation of the batch and biological variables
	#as a data frame with one row per cell, with columns specifying the various batch/biological categories
	#with values of 1 where the cell is of the category and 0 otherwise (dummy)
	#and subsequently identify which of the data frame columns are batch rather than biology (batch_index)
	#and subset the data frame to just those columns, in np.array form (dm)
	dummy = pd.get_dummies(adata.obs[batch_key+confounder_key],drop_first=False)
	if len(batch_key)>1:
		batch_index = np.logical_or.reduce(np.vstack([dummy.columns.str.startswith(x) for x in batch_key]))
	else:
		batch_index = np.vstack([dummy.columns.str.startswith(x) for x in batch_key])[0]
	dm = np.array(dummy)[:,batch_index]
	
	#compute how many genes at a time will be processed - aiming for chunksize total elements per
	chunkcount = np.ceil(chunksize/adata.shape[0])
	
	#make a Ridge with all the **kwargs passed if need be, and fit_intercept set to False
	#(as the data is centered). create holders for results
	LR = Ridge(fit_intercept=False, **kwargs)
	X_explained = []
	X_remain = []
	#loop over the gene space in chunkcount-sized chunks
	for ind in np.arange(0,adata.shape[1],chunkcount):
		#extract the expression and turn to dense if need be
		X_exp = adata.X[:,np.int(ind):np.int(ind+chunkcount)] # scaled data
		if scipy.sparse.issparse(X_exp):
			X_exp = X_exp.todense()
		#fit the ridge regression model, compute the expression explained by the technical 
		#effect, and the remaining residual
		LR.fit(dummy,X_exp)	
		X_explained.append(dm.dot(LR.coef_[:,batch_index].T))
		X_remain.append(X_exp - X_explained[-1])
	
	#collapse the chunked outputs and store them in the object
	X_explained = np.hstack(X_explained)
	X_remain = np.hstack(X_remain)
	adata.X = X_remain
	adata.layers['X_explained'] = X_explained
	logg.info('	finished', time=start,
		deep=('`.X` now features regression residuals\n'
		'	`.layers[\'X_explained\']` stores the expression explained by the technical effect'))
	return adata if copy else None

def save_images(output_dir, basename):
    outname = os.path.join(output_dir, basename + '.png')
#     plt.savefig(outname, dpi=200)
    plt.savefig(outname, dpi=200, bbox_inches='tight')
    plt.close()
    
# input_fn: a csv or csv.gz file, a matrix of genes x cells 
# output_fn: a csv file name to keep corrected data, genex x cells
def batch_effect_removal(save_dir, datatag, 
                         input_fn, meta_fn, output_fn, 
                         confounder_factor='clone', batch_factor='batch'):
    
    if not os.path.exists(save_dir): os.makedirs(save_dir)
#     input_fn = os.path.join(input_dir,'SA535_norm_untreated.csv.gz')
    
    if bool(re.search(".gz$", input_fn)):
        norm_df = pd.read_csv(input_fn, header=0,compression='gzip',sep=',', index_col=0)
    else:
        norm_df = pd.read_csv(input_fn, header=0, index_col=0)
        
    print(norm_df.shape)
    
#     meta_fn = os.path.join(input_dir,'SA535_meta_info_untreated.csv')
    meta_df = pd.read_csv(meta_fn, header=0)  #'meta_df should have batch','clone' columns
    print('Meta data: ')
    print(meta_df.shape)
    meta_df.index = meta_df['cell_id'].values
    print('Index: ')
    print(norm_df.index[1])
    # cells_use = np.intersect1d(norm_df.columns.values, meta_df['cell_id'].values)
    # print(len(cells_use))
    # norm_df = norm_df.loc[:,cells_use]
    # meta_df = meta_df.loc[cells_use,:]
    print('DEBUG')
    print(norm_df.shape)
    print(meta_df.shape)
    cells_use = norm_df.columns.values
    print(cells_use[1])
    var = pd.DataFrame(index=norm_df.index.values)
    print(var.shape)
    norm_df = np.transpose(norm_df)
    adata = sc.AnnData(np.asarray(norm_df), obs=meta_df, var=var)
    adata.obs_names = cells_use
    adata.obs['batch'] = adata.obs['batch'].astype("category")
    adata.obs['clone'] = adata.obs['clone'].astype("category")
    print(adata)

        
    print('Correcting technical effects...')
    corrected_adata = ridge_regression(adata, batch_key=[batch_factor], 
                 confounder_key=[confounder_factor],
                 copy=True)
    
    print(corrected_adata)
    print('Before correction: max is: {0}, after correction: max is: {1}'.format(np.round(np.max(adata.X),2),np.round(np.max(corrected_adata.X),2)))
    print('Visualizing corrected data...')
    sc.pp.highly_variable_genes(corrected_adata, n_top_genes=3000)
    sc.pp.pca(corrected_adata)
    sc.pp.neighbors(corrected_adata)
    sc.tl.umap(corrected_adata)
    sc.pl.umap(corrected_adata, color=[batch_factor,confounder_factor], show=False)
    save_images(save_dir, datatag+'_corrected_normalized_umap')
    adata_fn = os.path.join(save_dir,datatag+'_batch_correction.hdf5')
    corrected_adata.write(adata_fn)
    print('Save corrected data as: ')
    print(adata_fn)
    corrected_mtx = pd.DataFrame(corrected_adata.X.T, 
                            index=corrected_adata.var_names.values,
                            columns=corrected_adata.obs_names.values)
    print(corrected_mtx.shape)
#         output_fn=os.path.join(save_dir,datatag+'_untreated_batch_corrected_mtx.csv') 
    print('Save corrected matrix as: ')  
    print(output_fn)
    corrected_mtx.to_csv(output_fn, index=True, compression='gzip')
    
    print('Visualizing raw data...')
    sc.pp.highly_variable_genes(adata, n_top_genes=3000)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=[batch_factor,confounder_factor], show=False)
    save_images(save_dir, datatag+'_normalized_umap')

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_fn', required=True)  #csv fn of normalized mtx
    parser.add_argument('--meta_fn', required=True)  # cells metadata df
    parser.add_argument('--output_fn', required=True)
    parser.add_argument('--datatag', default='SA')
    parser.add_argument('--batch', default='batch')  # column in metafile that contain batch info
    parser.add_argument('--confounder', default='cell_type')# column in metafile that contain biological confounder info
    args = vars(parser.parse_args())
    print(args['input_fn'])
    print(args['meta_fn'])
    print(args['output_fn'])
    print(args['confounder'])
    print(args['batch'])
    batch_effect_removal(os.path.dirname(args['output_fn']), args['datatag'], args['input_fn'], 
                         args['meta_fn'], args['output_fn'], 
                         confounder_factor=args['confounder'], batch_factor=args['batch'])


# python $configPath/ridge_regression_batch_effect_removal.py --input_fn $input_fn --meta_fn $meta_fn --output_fn $input_fn --datatag $datatag --batch batch --confounder clone
