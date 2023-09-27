## Input data
## tree phylogeny tree 
## cell_clones file cell_id, clone_id 
## copy number matrix, chr_start_end in rowname, cell_id in colnames 
library(dplyr)
library(signals)
## Reference: please cite the signals package if you use these scripts
## https://github.com/shahcompbio/signals
## I add some scripts to run this function but main functions are from signal packgage, Marc William
## Noted: comment main() function in this file
source('/home/htran/Projects/hakwoo_project/corrupt_tree/src/clustering/make_cell_copynumber_tree_heatmap.R')


## pctcells: percentage of minimum cells to consider as one clone
hdbscran_viz <- function(download_dir, results_dir, library_id, grouping_file, pctcells = 0.05){
  # pctcells = 0.05
  copynumber <- data.table::fread(paste0(results_dir,library_id,'_filtered_states.csv.gz'))%>% as.data.frame()
  rownames(copynumber) <- copynumber$V1
  copynumber$V1 <- NULL
  print(dim(copynumber))
  
  reads_df <- data.table::fread(paste0(download_dir,library_id,'/hmmcopy/','reads.csv.gz')) %>% as.data.frame()
  # reads_df <- data.table::fread(paste0(download_dir,library_id,'/hmmcopy/',library_id,'_reads.csv.gz')) %>% as.data.frame()
  print(dim(reads_df))
  # reads_df$copy[1:10]
  
  reads_df <- reads_df %>%
    dplyr::filter(cell_id %in% colnames(copynumber))
  
  reads_df$chr_desc <- paste0(reads_df$chr,'_',reads_df$start,'_',reads_df$end)
  reads_df <- reads_df %>%
    dplyr::filter(chr_desc %in% rownames(copynumber))
  
  print(dim(reads_df))
  ncells <- length(unique(reads_df$cell_id))
  ## Using field = "copy", copy number values as the input for cell clustering
  clusters <- signals::umap_clustering(reads_df,
                                        minPts = max(round(pctcells * ncells), 20),
                                        field = "copy")
  
  # clusters <- signals::umap_clustering(reads_df, 
  #                                      minPts = min(round(pctcells * ncells), 10), 
  #                                      field = "copy")
  
  tree <- clusters$tree
  clones <- clusters$clustering
  
  ## Save output of clones here
  
  # ## For A98166A 
  # clones$clone_id <- ifelse(clones$clone_id=='0','C',clones$clone_id)
  # clones$clone_id <- ifelse(clones$clone_id=='C','A',
  #                           ifelse(clones$clone_id %in% c('NA'),'C','NA'))
  # class(clones)
  # data.table::fwrite(clones, paste0(save_dir, library_id,'_cell_clones.csv'))
  # saveRDS(clones, paste0(save_dir, library_id,'_cell_clones.rds'))
  # ## For A98166B
  # clones$clone_id <- ifelse(clones$clone_id %in% c('A','0'),'C',
  #                           ifelse(clones$clone_id=='B','A','D'))
  # summary(as.factor(clones$clone_id))
  res_clones <- unique(clones$clone_id)
  print(res_clones)
  if('0' %in% res_clones & !'None' %in% res_clones){
    clones <- ifelse(clones=='0','None',clones)
  }
  data.table::fwrite(clones, paste0(save_dir,library_id,'_cell_clones.csv'))
  
  save_dir <- paste0(results_dir,'hdbscan_results/')
  if(!dir.exists(save_dir)){
    dir.create(save_dir)
  }
  tree_fn <- paste0(save_dir,library_id,'_tree.newick')
  ape::write.tree(phy = tree, file = tree_fn, tree.names = F)
  
  
  output_fn <- paste0(save_dir, library_id,'_hdbscan_cell_cn_tree_heatmap.png')
  png(output_fn, height = 2*700, width=2*1200, res = 2*72)
  make_cell_copynumber_tree_heatmap(
    tree, copynumber, clones, NULL, grouping_file
  )
  dev.off()
  
  
}

# Mixing experiment with primary cells library A98166A
# experiments: 1:0.25, expected: A: 1, B: 0.25
# Output: A: 250 cells, C: 13 cells. 


# results_dir <- '/home/htran/storage/gm_instability_results/HEK293/data_thrs_10000/'
# download_dir <- '/home/htran/storage/datasets/damian_DLP/'
# grouping_file <- '/home/htran/storage/gm_instability_results/HEK293/library_groupings.csv'
# libs_id <- c('A110616A','A110617A','A118346B','A118825A','A118866A')
# for(library_id in libs_id){
#   hdbscran_viz(results_dir, library_id, grouping_file)
# }

# results_dir <- '/home/htran/storage/datasets/metastasis_results/dlp_SA919_mixing_exp/'
# download_dir <- '/home/htran/storage/raw_DLP/metastasis_DLP/SA919/'
# grouping_file <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/dlp_trees/SA919_mixing_experiment/library_groupings.csv'
# libs_id <- c('A98166A','A98166B')
# libs_id <- c('A130841A','A130841B','A130839B')
# library_id <- 'A130839B'
# for(library_id in libs_id){
#   hdbscran_viz(download_dir, results_dir, library_id, grouping_file)
# }

# results_dir <- '/home/htran/storage/datasets/metastasis_results/SA1142/first_results/'
# download_dir <- '/home/htran/storage/raw_DLP/metastasis_DLP/SA1142/'
# grouping_file <- '/home/htran/storage/datasets/metastasis_results/SA1142/first_results/library_groupings.csv.gz'
# lib_df <- data.table::fread(grouping_file)
# 
# libs_id <- lib_df$grouping
# library_id <- lib_df$grouping[1]
# 
# for(library_id in lib_df$grouping[2:dim(lib_df)[1]]){
#   hdbscran_viz(download_dir, results_dir, library_id, grouping_file)
# }

results_dir <- '/home/htran/storage/datasets/metastasis_results/dlp_SA919_mixing_exp/'
download_dir <- '/home/htran/storage/raw_DLP/metastasis_DLP/SA919/'
grouping_file <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/dlp_trees/SA919_mixing_experiment/library_groupings_extra.csv'
lib_df <- data.table::fread(grouping_file)
libs_id <- unique(lib_df$grouping)
libs_id[2:3]

for(library_id in libs_id[2:3]){
  hdbscran_viz(download_dir, results_dir, library_id, grouping_file)
}



# results_dir <- '/home/htran/storage/datasets/metastasis_results/dlp_SA535_X8/'
# download_dir <- '/home/htran/storage/raw_DLP/metastasis_DLP/SA535_X8/'
# grouping_file <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/dlp_trees/SA535_X8/library_groupings.csv'
# lib_df <- data.table::fread(grouping_file)
# 
# libs_id <- lib_df$grouping
# library_id <- lib_df$grouping[2]
# library_id
# for(library_id in lib_df$grouping[2:dim(lib_df)[1]]){
#   hdbscran_viz(download_dir, results_dir, library_id, grouping_file)
# }

