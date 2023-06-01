#!/usr/bin/env Rscript


## Installation
## Install the following package to visualize data. 
# BiocManager::install('ape')
# BiocManager::install('argparse')
# BiocManager::install('dplyr')
# BiocManager::install('ggplot2')
# BiocManager::install('ggtree')
# BiocManager::install('gtools')
# BiocManager::install('RColorBrewer')
# BiocManager::install('devtools')
# BiocManager::install('igraph')
# library(devtools)
# install_github("jokergoo/ComplexHeatmap")

## Note: here are the package version that I used to plot figure
## but I think the script can work with almost other versions. 
## For ComplexHeatmap, I recommend to use the latest version 
# [1] "ape: 5.6.2"
# [1] "argparse: 2.1.6"
# [1] "ComplexHeatmap: 2.8.0"
# [1] "dplyr: 1.0.9"
# [1] "ggplot2: 3.3.6"
# [1] "ggtree: 3.0.4"
# [1] "gtools: 3.9.3"
# [1] "RColorBrewer: 1.1.3"
# [1] "igraph: 1.3.4"

## Original script is by Sohrab Salehi
## If you use this script, please cite sitka paper from link: https://github.com/UBC-Stat-ML/sitkatree
# Modified scripts by Hoa Tran

## How to run program, see main() function at run_treecut_parameters.R

suppressPackageStartupMessages({
    require("dplyr")
    require("ggtree")
    require("ape")
    require("gtools")
    require("igraph")
    require("data.table")
  })


  
run_treecut <- function(input_dir, datatag, save_dir=NULL,
                        newick_fn=NULL, copynumber_fn=NULL, grouping_file=NULL){
    if(is.null(save_dir)){
      save_dir <- paste0(input_dir,"tree_cut_out/")  
    }
    if (!file.exists(save_dir)){
      dir.create(save_dir)
    }
    if(is.null(newick_fn)){
      newick_fn <- paste0(input_dir, 'tree.newick')
    }
    # newick <- paste0(input_dir, 'corrupt_grow/grown_nonecells.newick')
  
    if(is.null(copynumber_fn)){
      copynumber_fn <- paste0(input_dir, 'total_merged_filtered_states.csv.gz')
    }
  
    if(is.null(grouping_file)){
      grouping_file <- paste0(input_dir, 'library_groupings.csv')
    }
    
    ## minimum and maximum fraction of cells in total cells to consider as a clone
    ## ex: a minimum clone should have size of 5% of total cells. 
    minimum_fraction <- 0.05
    maximum_fraction <- 0.3
    
    ## After do clustering, dividing cells into cluster, 
    ## looking into each cluster, and see if there are any outlier cells,
    ## then put these cells as None label. 
    ## Can set ex: 0.75; 0.975, 1. Here I used 1 means I don't want to remove any outlier cells
    prob <- 0.975
    prob <- 1
    ## Difference in term of median copy number profile to define two clones or one clone
    ## Not remember exactly, need to go into scripts...
    edge_diff_threshold <- 0.03
    clone_size_threshold <- 0.1
    adjust_ploidy <- 0
    
    ## If two clones are quite similar, can merge them back into one clone
    min_merge_back_fraction <- 0.03
    
    ## Reading input files
    tree <- ape::read.tree(newick_fn)
    copy_number <- data.table::fread(copynumber_fn)
    rownames(copy_number) <- copy_number$V1
    copy_number$V1 <- NULL
    desc <- paste0(datatag,'_if_',minimum_fraction,'_af_',maximum_fraction,'_p',prob,'_e',edge_diff_threshold)
    print(desc)
    
    ## Lumberjack tree cutting
    launch_run(desc, tree,
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
    
    
    ## In case you want to change clone labels and replotting
    ## Loading clone assignment results
    ## And change clone labels
    cell_clones <- data.table::fread(paste0(save_dir, 'cell_clones_',desc,'.csv.gz'))
    total_cells <- colnames(copy_number)
    unassigned_cells <- total_cells[!total_cells %in% cell_clones$cell_id]
    unassigned_clones <- tibble::tibble(cell_id=unassigned_cells, clone_id='Un')
    cell_clones <- dplyr::bind_rows(cell_clones, unassigned_clones)
    print(summary(as.factor(cell_clones$clone_id)))
    data.table::fwrite(cell_clones, paste0(save_dir, 'cell_clones_',desc,'_corrected.csv.gz'))
    
    ## And replotting using new labels: 
    print("Plot heatmap")
    output_hm <- paste0(save_dir, 'heatmap_',desc,'_corrected.png')
    png(output_hm, height = 2*800, width=2*1400, res = 2*72)
    make_cell_copynumber_tree_heatmap(
      tree, copy_number, cell_clones, NULL, grouping_file
    )
    dev.off()
    print("Holala, have fun!!")
    
    
    
    
    ## Playing with parameters to obtain a best cutting
    ## Run tree cut using different parameter sets and select the best one
    # minfrac_ls <- c(0.04, 0.08)
    # # minfrac_ls <- c(0.1)
    # # maxfrac_ls <- c(0.85)
    # maxfrac_ls <- c(0.55,0.65,0.75)
    # # maxfrac_ls <- c(0.45, 0.55, 0.65, 0.75)
    # edge_df_ls <- c(0.04)
    # prob_ls <- c(0.9)
    # 
    # df <- expand.grid(minfrac_ls, maxfrac_ls, prob_ls, edge_df_ls)
    # print(dim(df))
    # clone_size_threshold <- 0.13
    # # edge_diff_threshold <- 0.04
    # adjust_ploidy <- 0
    # # prob <- 0.75
    # min_merge_back_fraction <- 0.03
    # tree <- ape::read.tree(newick_fn)
    # copy_number <- data.table::fread(copynumber_fn)
    # rownames(copy_number) <- copy_number$V1
    # copy_number$V1 <- NULL
    # for(i in rep(1:nrow(df),1)){
    #   minimum_fraction <- df[i,1]
    #   maximum_fraction <- df[i,2]
    #   prob <- df[i,3]
    #   edge_diff_threshold <- df[i,4]
    #   launch_run(desc, tree,
                  # copy_number,
                  # grouping_file,
                  # save_dir,
                  # minimum_fraction,
                  # maximum_fraction,
                  # clone_size_threshold,
                  # edge_diff_threshold,
                  # adjust_ploidy,
                  # prob,
                  # min_merge_back_fraction)

    #   
    #   
    # }
    
}

main(){
  
  ## Yaniv: change to your input script directory
  project_dir <- '/home/htran/Projects/hakwoo_project/corrupt_tree/src/tree_cut_scripts/'
  
  source(paste0(project_dir, "graph_cut_lumberjack/utils.R"))
  source(paste0(project_dir, "graph_cut_lumberjack/treecut_parameters.R"))
  source(paste0(project_dir, "graph_cut_lumberjack/make_cell_copynumber_tree_heatmap.R"))
  
  ## Yaniv: change to your input testing data sitka directory
  input_dir <- '/home/htran/Projects/hakwoo_project/corrupt_tree/src/tree_cut_scripts/testing_data_sitka/'
  save_dir <- paste0(input_dir,'output_lumberjack/')
  newick_fn <- paste0(input_dir, 'tree.newick')
  copynumber_fn <- paste0(input_dir, 'A98207B_filtered_states.csv')
  grouping_file <- paste0(input_dir, 'library_groupings.csv')
  datatag <- 'A98207B'
  
  ## Execution time ~2 mins
  run_treecut(input_dir, datatag, save_dir, newick_fn, copynumber_fn, grouping_file)
    
}

