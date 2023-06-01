

results_dir <- '/home/htran/storage/datasets/metastasis_results/SA535_new_encoding/SA535_wholedata_v2/'

cells_clone <- read.csv(paste0(results_dir,'cell_clones_v1.csv'), check.names=F, stringsAsFactors=F)
summary(as.factor(cells_clone$clone_id))
dim(cells_clone)
rownames(cells_clone) <- cells_clone$cell_id

cells_clone[outliers_cell,'clone_id'] <- 'None'
View(cells_clone[1:30,])

cells_clone[cell_clone_I$cell_id,'clone_id'] <- cell_clone_I$clone_id

cell_clone_I$clone_id <- paste0('I_',cell_clone_I$clone_id)
rownames(cell_clone_I) <- cell_clone_I$cell_id
outliers_cell <- cell_clone_I[cell_clone_I$clone_id=='I_A',]$cell_id
length(outliers_cell)

write.csv(cells_clone, paste0(paste0(results_dir,'cell_clones.csv')), quote = F, row.names = F)

total_cn <- data.table::fread(paste0(results_dir,'corrupt_grow/total_merged_filtered_states_original.csv'))
total_cn <- as.data.frame(total_cn)
rownames(total_cn) <- total_cn$V1
total_cn <- total_cn[, colnames(total_cn) != 'V1']
View(total_cn[1:3,1:3])
dim(total_cn)

cells_use <- cells_clone[cells_clone$clone_id=='I','cell_id']
length(cells_use)  
total_obs_cn <- total_cn[,colnames(total_cn) %in% cells_use]
dim(total_obs_cn)

length(med)
View(med)


cn_mat <- as.matrix(total_obs_cn)
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

compute_dist_mat_outliers <- function(cn_mat, use_hamming = TRUE) {
  
  
  med <- apply(cn_mat, 1, median)
  mode_val <- apply(cn_mat, 1, getmode)
  class(med)
  stat <- data.frame(med_val=med, mode_val=mode_val)
  View(stat)
  View()
  names(med) <- rownames(cn_mat)
  names(mode_val) <- rownames(cn_mat)
  
  intra_dist <- c()
  for(c in seq(ncol(cn_mat))){
      if (use_hamming) {
        intra_dist <- c(intra_dist, sum(cn_mat[,c]!=mode_val))
      } else {
        intra_dist <- c(intra_dist, sum(abs(cn_mat[,c]-mode_val)))
      }
  }
  names(intra_dist) <- colnames(cn_mat)
  outliers <- is_outlier(intra_dist, prob=0.7, type = "t") 
  
  outlier_cells <- names(outliers[outliers==T])
  length(outlier_cells)
}  


is_outlier <- function(x, prob, type = "t") {
  tt <- outliers::scores(x = x, type = type, prob = prob)
  ff <- outliers::scores(x = x, type = type) < 0
  tt & ff
}


outlier_cells <- cells_clone[cells_clone$clone_id !='I','cell_id']
print(length(outlier_cells))
outlier_cells <- paste0('cell_',outlier_cells)

tree <- read.tree(paste0(results_dir, 'tree.newick'))
tree_cg <- ape::drop.tip(tree, outlier_cells, trim.internal =T, collapse.singles = T)
all_cells <- grep('cell', tree_cg$tip.label, value = T)
print(length(all_cells))
# ggtree(tree_cg)
write.tree(phy = tree_cg, file = paste0(results_dir, 'corrupt_grow/clone_I_v2/cloneI.newick'), tree.names = F)

ggraph::ggraph(tree_cg)

# Get data from filtered cn
# Load cell clones, get clone I 
# Apply outlier functions 
# Get newick tree, cut clone I only 
# plot label of clone I 