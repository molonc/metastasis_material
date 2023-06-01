results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535/SA535_total_v2/'

# Objective: map unassigned cells - non cells back to the tree
# 1. Get subtree which contain assigned cells only
# 2. Run corrupt-grow with subtree and unassigned cell state matrix

newick <- paste0(results_dir, 'tree.newick')
tree <- read.tree(newick)

cell_clones <- read.csv(paste0(results_dir,'cell_clones_v3.csv'), stringsAsFactors=F, check.names=F)
View(head(cell_clones))
summary(as.factor(cell_clones$clone_id))
dim(cell_clones)
cn <- read.csv(paste0(results_dir,'total_merged_filtered_states.csv'), row.names = 1, stringsAsFactors=F, check.names=F)
View(cn[1:4,1:4])
dim(cn)
none_cells <- colnames(cn)[colnames(cn) %in% cell_clones$cell_id]
length(none_cells)
# none_cells: cells that can not be assigned in the tree
# trim.internal =F, collapse.singles = F: avoid remove internal nodes from the tree 
# corrupt-grow input: tree and input matrix should have same cells or same loci list
none_cells <- cell_clones[cell_clones$clone_id=='None','cell_id']
length(none_cells)
summary(as.factor(cell_clones$clone_id))
none_cells <- paste0('cell_',none_cells)

none_cells <- paste0('cell_',cell_clones$cell_id)

tree_cg <- ape::drop.tip(tree, none_cells, trim.internal =T, collapse.singles = T)
tree_cg
root(tree_cg)
t <- grep('cell_',tree_cg$tip.label,value = T)
length(t)
# From Alexandre Bouchard: The method corrupt-grow assumes that the root is called 'ROOT'. 
# Adding a ROOT node fixed the issue.
root_name <- ''
is_root <- sum(tree_cg$node.label==root_name)
tree_cg$node.label[1] = 'root'

# verification
print('ROOT' %in% tree_cg$node.label)  
ggtree(tree_cg)


subtree_fn <- paste0(results_dir,'corrupt_grow/subtree.newick')
# https://www.rdocumentation.org/packages/ape/versions/5.3/topics/write.tree
write.tree(phy = tree_cg, file = subtree_fn, tree.names = F)
write.tree(phy = tree_cg, file = paste0(results_dir, 'tree.newick'), tree.names = F)


