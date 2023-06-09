import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import re
from treealign import CloneAlignClone
from treealign import CloneAlignTree
# from Bio import Phylo


def load_input_data(expr_fn, cnv_fn, cell_clones_fn=None, tree_fn=None):
    # scRNA read count matrix where each row represents a gene, 
    # each column represents a cell
#     expr = pd.read_csv(expr_fn, index_col=0)
    if bool(re.search(".gz$", expr_fn)):
        expr = pd.read_csv(expr_fn, header=0, compression='gzip', sep=',', index_col=0)
    else:
        expr = pd.read_csv(expr_fn, index_col=0)

    print('Loading expression data: ')    
    print(expr.shape)    
    print(expr.head(2))
#     print(expr.index[2])
    # scDNA copy number matrix where each row represents a gene,
    # each column represents a cell
    # the numbers of the matrix reprents the copy number at given cells and genes
#     cnv = pd.read_csv("../data/example_gene_cnv.csv", index_col=0)
    if bool(re.search(".gz$", cnv_fn)):
        cnv = pd.read_csv(cnv_fn, header=0, compression='gzip', sep=',', index_col=0)
    else:
        cnv = pd.read_csv(cnv_fn, index_col=0)

    print('Loading CNV data: ')    
    print(cnv.shape)     
    print(cnv.head(2))
#     print(cnv.index[2])
    # clone labels for each cell in scDNA
#     clone = pd.read_csv("../data/example_cell_clone.csv")
    if cell_clones_fn is not None and bool(re.search(".gz$", cell_clones_fn)):
        clone = pd.read_csv(cell_clones_fn, header=0, compression='gzip', sep=',')
    elif cell_clones_fn is not None and not bool(re.search(".gz$", cell_clones_fn)):
        clone = pd.read_csv(cell_clones_fn) #, index_col=0
    else:
        clone = None

    print('Loading cell clones data: ')    
    print(clone.shape)                
    print(clone.head(2))
    
#     print(clone.head(3))
    if tree_fn is not None:
        tree = Phylo.read(tree_fn, "newick")
        print('Loading phylogenetic tree')  
    else:
        tree = None
        
    input_dic = {"expr":expr,"cnv":cnv,"clone":clone,"tree":tree}

    return input_dic


# class CloneAlignClone(treealign.clonealign.CloneAlign)
#  |  CloneAlignClone(clone, expr=None, cnv=None, hscn=None, snv_allele=None, snv=None, normalize_cnv=True, cnv_cutoff=10, infer_s_score=True, infer_b_allele=True, repeat=10, min_clone_assign_prob=0.8, min_clone_assign_freq=0.7, min_consensus_gene_freq=0.6, min_consensus_snv_freq=0.6, max_temp=1.0, min_temp=0.5, anneal_rate=0.01, learning_rate=0.1, max_iter=400, rel_tol=5e-05, cell_dirichlet_alpha=1, record_input_output=False, min_clone_cell_count=10, initialize_seed=False)
#  |      :param expr: expr read count matrix. row is gene, column is cell. (pandas.DataFrame)
#  |      :param cnv: cnv matrix. row is gene, column is cell. (pandas.DataFrame)
#  |      :param clone: groupings of cnv cells. (pandas.DataFrame)
#  |      :param normalize_cnv: whether to normalized cnv matrix by min or not. (bool)
#  |      :param cnv_cutoff: set cnv higher than cnv_cutoff to cnv_cutoff. (int)
#  |      :param model_select: "gene" for the extended clonealign model or "default" for the original clonelign model (str)
#  |      :param repeat: num of times to run clonealign to generate consensus results. (int)
#  |      :param min_clone_assign_prob: assign cells to a clone if clone assignment prob reaches min_clone_assign_prob (float)
#  |      :param min_clone_assign_freq: assign cells to a clone if a min proportion of runs generate the same results (float)
#  |      :param max_temp: starting temperature in Gumbel-Softmax reparameterization. (float)
#  |      :param min_temp: min temperature in Gumbel-Softmax reparameterization. (float)
#  |      :param anneal_rate: annealing rate in Gumbel-Softmax reparameterization. (float)
#  |      :param learning_rate: learning rate of Adam optimizer. (float)
#  |      :param max_iter: max number of iterations of elbo optimization during inference. (int)
#  |      :param rel_tol: when the relative change in elbo drops to rel_tol, stop inference. (float)

def run_clone_TreeAlign(expr, cnv, clone, save_dir, datatag,
                        min_clone_assign_prob=0.8, max_iter=500, min_clone_cell_count=1, cnv_cutoff=8):
#     if save_dir==None: 
#         save_dir = os.path.join(os.path.dirname(expr_fn),'test_TreeAlign/')

    if not os.path.exists(save_dir): 
        os.makedirs(save_dir)
        print('Creating an output folder')
        print(save_dir)
        
    # construct CloneAlignTree object for data preprocessing

    # `repeat` is set to 1 here for demonstration purposes. it would be better to set `repeat` larger than 5. 
    # obj = CloneAlignClone(clone=clone, expr=expr, cnv=cnv, hscn=hscn, snv_allele=snv_allele, snv=snv_total, repeat=1)

    # it is possible to run TreeAlign with total copy number data only
    obj = CloneAlignClone(clone=clone, expr=expr, cnv=cnv, repeat=1,
                         min_clone_assign_prob=min_clone_assign_prob, max_iter=max_iter, min_clone_cell_count=min_clone_cell_count, cnv_cutoff=cnv_cutoff)


#     obj = CloneAlignClone(clone=clone, expr=expr, cnv=cnv, repeat=1)

    # it is also possible to run TreeAlign with allele specific data only
    # obj = CloneAlignClone(clone=clone, hscn=hscn, snv_allele=snv_allele, snv=snv_total, repeat=1)

    # running TreeAlign to assign cells to phylogenetic subclades
    print('Holala')
    obj.assign_cells_to_clones()
    clone_assign_df, gene_type_score_df, allele_assign_prob_df = obj.generate_output()
    print('Store cells assignment into output folder: ')
    clone_assign_df.to_csv(os.path.join(save_dir, datatag + "_clone_assign.csv.gz"),index=False)
    gene_type_score_df.to_csv(os.path.join(save_dir, datatag + "_gene_type_score.csv.gz"),index=False)
    
    print("Output of clonal assignment is: ")
    cls = [str(c) for c in clone_assign_df['clone_id'].values]
    unique, counts = np.unique(cls, return_counts=True)
    print(dict(zip(unique, counts)))
    
    
def run_phylo_TreeAlign(expr, cnv, tree, save_dir, datatag):
    
    
    # construct CloneAlignTree object for data preprocessing

    # `repeat` is set to 1 here for demonstration purposes. it would be better to set `repeat` larger than 5. 
    # obj = CloneAlignClone(clone=clone, expr=expr, cnv=cnv, hscn=hscn, snv_allele=snv_allele, snv=snv_total, repeat=1)

    # it is possible to run TreeAlign with total copy number data only
    clonealignclone_obj = CloneAlignTree(expr, cnv, tree, repeat=1, min_gene_diff=50)
    clonealignclone_obj.assign_cells_to_tree()
    clone_assign_df, gene_type_score_df = clonealignclone_obj.generate_output()

    
    # Save data to folder
    clone_assign_df.to_csv(os.path.join(save_dir, datatag + "_clone_assign_test.csv.gz"))
    gene_type_score_df.to_csv(os.path.join(save_dir, datatag + "_gene_type_score_test.csv.gz"))

    
    
    
def save_images(basename, output_dir): 
    if not os.path.exists(output_dir): os.makedirs(output_dir) 
    outname = os.path.join(output_dir, basename + '.png') 
    plt.savefig(out)
    
    
    
