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


suppressPackageStartupMessages({
  library(dplyr)
  library(signals)
})  
options(dplyr.summarise.inform = FALSE)
options(tidyverse.quiet = TRUE)

## Reference: please cite the signals package if you use these scripts
## https://github.com/shahcompbio/signals
## Noted: original version from Tyler Funnel 
## Modified scripts by Hoa Tran
## Hoa: adding some scripts to run this function but main functions are from signal packgage, Marc William
## Note: comment main() function in this file








## Note: copy number file name with row name is chr regions, and column name is cell id
## filtered_CNV_fn file has this format: 
#                         AT11391-A98166B-R45-C08
# 1_2000001_2500000                       1
# 1_3000001_3500000                       1
# 1_4000001_4500000                       1

## Note: usually we have 4375 filtered genomic bins, but sometimes, you may exclude some noisy bins from a reference list, 
## and have ~3950 bins for analysis 

get_filtered_data <- function(filtered_CNV_fn){
  cnv <- data.table::fread(filtered_CNV_fn) %>% as.data.frame()
  rownames(cnv) <- cnv$V1 # chr regions are row names
  cnv$V1 <- NULL
  print(dim(cnv))
  ##  Note: format of cell id is: sampleId-libraryId-rowId-colId, ex: "AT11391-A98166B-R45-C08"
  filtered_cells <- colnames(cnv)
  
  ##  Note: format of chr region is: chr_start_end, ex: '1_2000001_2500000'
  filtered_chr_regions <- rownames(cnv)
  return(list(copynumber=cnv,
              filtered_cells=filtered_cells, 
              filtered_chr_regions=filtered_chr_regions))
}



## pctcells: percentage of minimum cells to consider as one clone
hdbscran_viz <- function(filtered_CNV_fn, reads_fn,
                         save_dir, datatag, grouping_file, pctcells = 0.05){
  res <- get_filtered_data(filtered_CNV_fn)
  copynumber <- res$copynumber
  if(!dir.exists(save_dir)){
    dir.create(save_dir)
  }
  
  ## Input hmmcopy reads file
  reads_df <- data.table::fread(reads_fn) %>% as.data.frame()
  # reads_df <- data.table::fread(paste0(download_dir,library_id,'/hmmcopy/',library_id,'_reads.csv.gz')) %>% as.data.frame()
  print(dim(reads_df))
  
  reads_df <- reads_df %>%
    dplyr::filter(cell_id %in% res$filtered_cells)
  
  reads_df$chr_desc <- paste0(reads_df$chr,'_',reads_df$start,'_',reads_df$end)
  reads_df <- reads_df %>%
    dplyr::filter(chr_desc %in% res$filtered_chr_regions)
  
  print(dim(reads_df))
  ncells <- length(unique(reads_df$cell_id))
  
  
  ##Note:  here I use field = "copy", copy number values as the input for cell clustering,
  ## instead of copy number state at 'state' column
  ## you can use different columns, and comparing the results of clustering
  clusters <- signals::umap_clustering(reads_df,
                                        minPts = max(round(pctcells * ncells), 30),
                                        field = "copy")
  
  # clusters <- signals::umap_clustering(reads_df, 
  #                                      minPts = min(round(pctcells * ncells), 30), 
  #                                      field = "copy")
  
  tree <- clusters$tree
  clones <- clusters$clustering
  print(unique(clones$clone_id))
  
  
  ## Save output of clones assignment here
  tree_fn <- paste0(save_dir,datatag,'_tree.newick')
  ape::write.tree(phy = tree, file = tree_fn, tree.names = F)
  
  data.table::fwrite(clones, paste0(save_dir,datatag,'_cell_clones.csv'))
  
  output_fn <- paste0(save_dir, datatag,'_hdbscan_cell_cn_tree_heatmap.png')
  png(output_fn, height = 2*700, width=2*1200, res = 2*72)
  make_cell_copynumber_tree_heatmap(
    tree, copynumber, clones, NULL, grouping_file
  )
  dev.off()
  
  
}



main <- function(){
  ##  Note: please change the directory to your script dir here. 
  script_dir <- '/home/htran/Projects/hakwoo_project/corrupt_tree/src/hdbscan_clustering/'
  source(paste0(script_dir,'make_cell_copynumber_tree_heatmap.R'))
  
  
  ## Note: you can change the input directory to point to hdbscan_clustering/testing_data/ 
  
  results_dir <- '/home/htran/Projects/hakwoo_project/corrupt_tree/src/hdbscan_clustering/testing_data/'
  save_dir <- paste0(results_dir,'output/')
  grouping_file <- paste0(results_dir,'library_grouping.csv')
  reads_fn <- paste0(results_dir,'A98166B_reads.csv.gz')
  datatag <- 'A98166B'
  filtered_CNV_fn <- paste0(results_dir,'A98166B_filtered_states.csv.gz')
  
  hdbscran_viz(filtered_CNV_fn, reads_fn,
               save_dir, datatag, grouping_file, pctcells = 0.05)
}
