#!/usr/bin/env Rscript

  # Hoa run tree cut with different sets of parameters
  suppressPackageStartupMessages({
    library(tidyr)
    library(ggtree)
    require("ape")
    require("gtools")
    require("igraph")
  })
  project_dir <- '/home/htran/Projects/hakwoo_project/corrupt_tree'
  source(paste0(project_dir, "/src/testing_only/make_cell_copynumber_tree_heatmap.R"))
  source(paste0(project_dir, "/src/tree_cut/utils.R"))
  source(paste0(project_dir, "/src/tree_cut/treecut_parameters.R"))
  
  # results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local_tyler_left_padding/'
  # results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_filtered/'
  
  # results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_new_encode/SA1035_filtered_v2/'
  # results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_new_encode/SA1035_filtered_FP001/'
  # results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_new_encode/SA1035_Tyler/'
  
  
  
  # results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_new_encode/SA1035_whole_local/'
  # results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local/'
run_treecut <- function(results_dir){
    save_dir <- paste0(results_dir,"tree_cut_out/")
    if (!file.exists(save_dir)){
      dir.create(save_dir)
    }
    newick <- paste0(results_dir, 'tree.newick')
    # newick <- paste0(results_dir, 'corrupt_grow/grown_nonecells.newick')
    
    copynumber <- paste0(results_dir, 'total_merged_filtered_states.csv')
    grouping_file <- paste0(results_dir, 'library_groupings.csv')
    
    
    minfrac_ls <- c(0.04, 0.08)
    # minfrac_ls <- c(0.1)
    # maxfrac_ls <- c(0.85)
    maxfrac_ls <- c(0.55,0.65,0.75)
    # maxfrac_ls <- c(0.45, 0.55, 0.65, 0.75)
    edge_df_ls <- c(0.04)
    prob_ls <- c(0.9)
    
    df <- expand.grid(minfrac_ls, maxfrac_ls, prob_ls, edge_df_ls)
    print(dim(df))
    clone_size_threshold <- 0.13
    # edge_diff_threshold <- 0.04
    adjust_ploidy <- 0
    # prob <- 0.75
    min_merge_back_fraction <- 0.03
    tree <- read.tree(newick)
    copy_number <- read.csv(copynumber, header=T, row.names=1, check.names = F,stringsAsFactors = FALSE)
    # copy_number <- read.delim(copynumber, check.names = FALSE, stringsAsFactors = FALSE, sep = ",", row.names = 1)
    # View(head(copy_number))
    for(i in rep(1:nrow(df),1)){
      minimum_fraction <- df[i,1]
      maximum_fraction <- df[i,2]
      prob <- df[i,3]
      edge_diff_threshold <- df[i,4]
      launch_run(tree,
                 copy_number,
                 grouping_file,
                 save_dir,
                 minimum_fraction,
                 maximum_fraction,
                 clone_size_threshold,
                 edge_diff_threshold,
                 adjust_ploidy,
                 prob,
                 min_merge_back_fraction)
      
      
    }
    
}

base_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535/'
fns <- c("SA535_Tyler")
# fns <- c("Tyler_2361","Tyler_2362","Tyler_2363","Tyler_2364")
for(subdir in fns){
  run_treecut(paste0(base_dir,subdir,'/'))
}

# base_dir2 <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/'
# fns2 <- c("Tyler_2251","Tyler_2252","Tyler_2253","Tyler_2254","Tyler_216","Tyler_2164")
# 
# for(subdir2 in fns2){
#   run_treecut(paste0(base_dir2,subdir2,'/'))
# }
  
# }
