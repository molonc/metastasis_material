# source_dir <- '/home/htran/Projects/farhia_project/rscript/dlp/visualize_tree/'
# source(paste0(source_dir,'trajectory_utils.R'))

# Farhia, SA1035
# results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_filtered/'
# cellclones <- paste0(results_dir, 'cell_clones.csv')
# datatag = 'SA919_Tyler'

# results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_filtered/'
# filter_df <- data.table::fread(paste0(results_dir,'filtered.csv'), stringsAsFactors = F)
# View(head(filter_df))
# dim(filter_df)
# loci_use <- unique(filter_df$loci)
# snv_2161_df <- data.table::fread(paste0(results_dir,'corrupt_grow/corrupt_tree_sync-cnv_features_2161.csv'), stringsAsFactors = F)
# dim(snv_2161_df)
# snv_2161_df <- snv_2161_df[snv_2161_df$loci %in% loci_use, ]
# data.table::fwrite(snv_2161_df, paste0(results_dir,'corrupt_grow/filtered_2161.csv'), row.names = F, quote = F)
state_mtx <- data.table::fread(paste0(results_dir,'total_merged_filtered_states_original.csv'), stringsAsFactors = F)
colnames(state_mtx)[1:3]
results_dir2 <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_filtered_v3/'
filter_df <- data.table::fread(paste0(results_dir2,'corrupt_tree_straightened_features.csv'), stringsAsFactors = F)
dim(filter_df)

library_ids <- get_library_labels_v2(colnames(state_mtx))
library_ids <- unique(library_ids)
'A108836B' %in% library_ids

cells_2161 <- unique(snv_2161_df$cells)
filter_df <- filter_df[filter_df$cells %in% cells_2161,]
filter_df <- filter_df[filter_df$loci %in% loci_use,]
length(unique(filter_df$loci))
length(loci_use)


get_clones_subgraph <- function(results_dir, cellclones){
  # Load tree
  newick <- paste0(results_dir, 'tree.newick')
  tree <- read.tree(newick)
  all_cells <- grep('cell', tree$tip.label, value = T)
  length(all_cells)
  cell_clones <- read.delim(cellclones, stringsAsFactors = FALSE, sep = ",")
  clones_cellidx <- all_cells[all_cells %in% paste0('cell_',cell_clones$cell_id)]
  length(clones_cellidx)
  none_cells <- all_cells[!all_cells %in% clones_cellidx]
  
  print(length(none_cells))
  cg_dir <- paste0(results_dir,'corrupt_grow/')
  if (!file.exists(cg_dir)){
    dir.create(cg_dir)
  }
  
  tree_cg <- ape::drop.tip(tree, none_cells, trim.internal =F, collapse.singles = F)
  tree_cg
  root_name <- ''
  root_name %in% tree_cg$node.label
  for(i in rep(1:length(tree_cg$node.label),1)){
    if(tree_cg$node.label[i]==root_name){
      tree_cg$node.label[i] = 'ROOT'
    }
  }
  print('ROOT' %in% tree_cg$node.label)  # verification
  # ggtree(tree_cg)
  subtree_fn <- paste0(cg_dir,'subtree.newick')
  # https://www.rdocumentation.org/packages/ape/versions/5.3/topics/write.tree
  write.tree(phy = tree_cg, file = subtree_fn, tree.names = F)
  # 
  if(file.exists(subtree_fn)){
    return(subtree_fn)
  } else{
    return(NULL)
  }
  
}


get_none_cells_cnstate <- function(results_dir, cellclones, tag='sohrab_version'){
  # Get tree with cells in clones only
  subtree_fn <- get_clones_subgraph(results_dir, cellclones)
  
  newick_path <- subtree_fn
  cg_dir <- paste0(results_dir,'corrupt_grow/')
  if (!file.exists(cg_dir)){
    dir.create(cg_dir)
  }
  filtered_df <- read.csv(paste0(results_dir, 'filtered.csv'), check.names = F, stringsAsFactors = FALSE)
  cell_clones <- read.delim(cellclones, stringsAsFactors = FALSE, sep = ",")
  filtered_nonecell_df <- filtered_df[!filtered_df$cells %in% paste0('cell_',cell_clones$cell_id),]
  length(unique(filtered_nonecell_df$cells))
  corrupt_cnv <- paste0(cg_dir, 'filtered_none_cells.csv')
  data.table::fwrite(filtered_nonecell_df, corrupt_cnv, row.names = F, quote = F)
  return(list(newick_path=newick_path, corrupt_cnv=corrupt_cnv))
  
}

# results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535/SA535_Tyler_v3/'
get_filtered_outliers <- function(results_dir){
  
  filtered_df <- read.csv(paste0(results_dir, 'filtered.csv'), check.names = F, stringsAsFactors = FALSE)
  dim(filtered_df)
  loci_ls <- unique(filtered_df$loci)
  length(loci_ls)
  filtered_df <- read.csv(paste0(results_dir, 'filtered.csv'), check.names = F, stringsAsFactors = FALSE)
  filtered_outliers <- data.table::fread(paste0(results_dir,'corrupt_grow/corrupt_tree_sync-cnv_features_outliers.csv'), stringsAsFactors = F)
  length(unique(filtered_outliers$loci))
  sum(unique(filtered_outliers$loci) %in% loci_ls)
  filtered_outliers <- filtered_outliers[filtered_outliers$loci %in% loci_ls,]
  
  length(unique(filtered_outliers$cells))
  data.table::fwrite(filtered_outliers, paste0(results_dir,'corrupt_grow/filtered_outliers.csv'), row.names = F, quote = F)
  
}

run_corrupt_grow <- function(results_dir,cellclones,
                             project_dir=NULL,
                             input_dir=NULL, output_dir = NULL,
                             ){
  # newick_path <- paste0(cg_dir,'subtree.newick')
  # corrupt_cnv <- paste0(cg_dir, 'filtered_none_cells.csv')
  corrupt_grow_path <- '/home/htran/storage/install_software/nowellpack/build/install/nowellpack/bin/'
  
  out_ls <- get_none_cells_cnstate(results_dir, cellclones)
  newick_path = out_ls[[newick_path]]
  corrupt_cnv = out_ls[[corrupt_cnv]]
  
  cmd_str <- sprintf('%scorrupt-grow --matrix NoisyBinaryCLMatrix --matrix.binaryMatrix %s --matrix.fpRate 0.01 --matrix.fnRate 0.5 --phylo file %s'
                     , corrupt_grow_path, corrupt_cnv, newick_path)
  print(cmd_str)
  print("Completed, Voila!!!")
 
  
}  




get_clones_subgraph <- function(results_dir, cellclones){
  # Load tree
  newick <- paste0(results_dir, 'tree.newick')
  tree <- read.tree(newick)
  all_cells <- grep('cell', tree$tip.label, value = T)
  length(all_cells)
  cell_clones <- read.delim(cellclones, stringsAsFactors = FALSE, sep = ",")
  clones_cellidx <- all_cells[all_cells %in% paste0('cell_',cell_clones$cell_id)]
  length(clones_cellidx)
  none_cells <- all_cells[!all_cells %in% clones_cellidx]
  
  print(length(none_cells))
  cg_dir <- paste0(results_dir,'corrupt_grow/')
  if (!file.exists(cg_dir)){
    dir.create(cg_dir)
  }
  
  tree_cg <- ape::drop.tip(tree, none_cells, trim.internal =F, collapse.singles = F)
  tree_cg
  root_name <- ''
  root_name %in% tree_cg$node.label
  for(i in rep(1:length(tree_cg$node.label),1)){
    if(tree_cg$node.label[i]==root_name){
      tree_cg$node.label[i] = 'ROOT'
    }
  }
  print('ROOT' %in% tree_cg$node.label)  # verification
  # ggtree(tree_cg)
  subtree_fn <- paste0(cg_dir,'subtree.newick')
  # https://www.rdocumentation.org/packages/ape/versions/5.3/topics/write.tree
  write.tree(phy = tree_cg, file = subtree_fn, tree.names = F)
  # 
  if(file.exists(subtree_fn)){
    return(subtree_fn)
  } else{
    return(NULL)
  }
  
}


get_none_cells_cnstate <- function(results_dir, cellclones, tag='sohrab_version'){
  # Get tree with cells in clones only
  subtree_fn <- get_clones_subgraph(results_dir, cellclones)
  
  newick_path <- subtree_fn
  cg_dir <- paste0(results_dir,'corrupt_grow/')
  if (!file.exists(cg_dir)){
    dir.create(cg_dir)
  }
  filtered_df <- read.csv(paste0(results_dir, 'filtered.csv'), check.names = F, stringsAsFactors = FALSE)
  cell_clones <- read.delim(cellclones, stringsAsFactors = FALSE, sep = ",")
  filtered_nonecell_df <- filtered_df[!filtered_df$cells %in% paste0('cell_',cell_clones$cell_id),]
  length(unique(filtered_nonecell_df$cells))
  corrupt_cnv <- paste0(cg_dir, 'filtered_none_cells.csv')
  data.table::fwrite(filtered_nonecell_df, corrupt_cnv, row.names = F, quote = F)
  return(list(newick_path=newick_path, corrupt_cnv=corrupt_cnv))
  
}

# results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535/SA535_Tyler_v3/'
get_filtered_outliers <- function(results_dir){
  
  filtered_df <- read.csv(paste0(results_dir, 'filtered.csv'), check.names = F, stringsAsFactors = FALSE)
  dim(filtered_df)
  loci_ls <- unique(filtered_df$loci)
  length(loci_ls)
  filtered_df <- read.csv(paste0(results_dir, 'filtered.csv'), check.names = F, stringsAsFactors = FALSE)
  filtered_outliers <- data.table::fread(paste0(results_dir,'corrupt_grow/corrupt_tree_sync-cnv_features_outliers.csv'), stringsAsFactors = F)
  length(unique(filtered_outliers$loci))
  sum(unique(filtered_outliers$loci) %in% loci_ls)
  filtered_outliers <- filtered_outliers[filtered_outliers$loci %in% loci_ls,]
  
  length(unique(filtered_outliers$cells))
  data.table::fwrite(filtered_outliers, paste0(results_dir,'corrupt_grow/filtered_outliers.csv'), row.names = F, quote = F)
  
}

run_corrupt_grow <- function(results_dir,cellclones,
                             project_dir=NULL,
                             input_dir=NULL, output_dir = NULL,
){
  # newick_path <- paste0(cg_dir,'subtree.newick')
  # corrupt_cnv <- paste0(cg_dir, 'filtered_none_cells.csv')
  corrupt_grow_path <- '/home/htran/storage/install_software/nowellpack/build/install/nowellpack/bin/'
  
  out_ls <- get_none_cells_cnstate(results_dir, cellclones)
  newick_path = out_ls[[newick_path]]
  corrupt_cnv = out_ls[[corrupt_cnv]]
  
  cmd_str <- sprintf('%scorrupt-grow --matrix NoisyBinaryCLMatrix --matrix.binaryMatrix %s --matrix.fpRate 0.01 --matrix.fnRate 0.5 --phylo file %s'
                     , corrupt_grow_path, corrupt_cnv, newick_path)
  print(cmd_str)
  print("Completed, Voila!!!")
  
  
}  
