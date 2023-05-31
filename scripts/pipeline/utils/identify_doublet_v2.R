# install.packages("remotes")
# remotes::install_github("huyuem/rscrublet")


# install.packages("/home/htran/storage/install_software/neldermead_1.0-11.tar", repos=NULL, type="source")
# install.packages("/home/htran/storage/install_software/optimbase_1.0-9.tar", repos=NULL, type="source")
# install.packages("/home/htran/storage/install_software/optimsimplex_1.0-7.tar", repos=NULL, type="source")
# devtools::install_local('/home/htran/storage/install_software/rscrublet')

suppressPackageStartupMessages({
  library(rscrublet)
  library(SingleCellExperiment)
  library(argparse)
  library(dplyr)
  library(tidyverse)
})
print("Run doublet analysis using rscrublet package version: ")
print(packageVersion("rscrublet"))

parser <- ArgumentParser(description = "Remove doublets from SingleCellExperiment")
parser$add_argument('--sce_file', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS")
parser$add_argument('--sce_output', type = 'character', metavar = 'FILE',
                    help="Output path for preprocessed SCE.")
parser$add_argument('--doublet_score', type = 'double', default = NULL,
                    help="doublet score threshold")
parser$add_argument('--output_csv', metavar='FILE', type = 'character',
                    help="Output dataframe with summary of # doublet cells")

parser$add_argument('--summary_csv', metavar='FILE', type = 'character',
                    help="Output dataframe metadata with summary of cells quality for total data")


args <- parser$parse_args()

# 
# doublet_score <- as.numeric(args$doublet_score)
# if(doublet_score==0){
#   doublet_score <- NULL
# }
# summary_csv <- args$summary_csv
# sce <- readRDS(args$sce_file)
# print(dim(sce))
# 
# doublets_df <- rscrublet::scrubDoublets(as.matrix(counts(sce)), n_neighbors = 15, 
#                                         doublet_score_threshold = doublet_score,
#                                         sim_doublet_ratio = 2, expected_doublet_rate = 0.06,
#                                         stdev_doublet_rate = 0.02, synthetic_doublet_umi_subsampling = 1,
#                                         use_approx_neighbors = TRUE, distance_metric = "euclidean",
#                                         min_counts = 2, min_cells = 3, min_gene_variability_pctl = 85,
#                                         log_transform = FALSE, z_score = TRUE, n_prin_comps = 30,
#                                         verbose = TRUE)
# 
# print(summary(doublets_df$scrubDoublets))
# print(summary(doublets_df$doublet_scores_obs))
# print(summary(doublets_df$doublet_scores_sim))
# 
# sce$is_doublet <- doublets_df$scrubDoublets
# sce$doublet_scores_obs <- round(doublets_df$doublet_scores_obs,2)
# 
# doublet_cells <- sce$Barcode[sce$is_doublet]
# print(paste0('Number of doublet cells: ', length(doublet_cells)))
# 
# tbl <- table(sce$is_doublet)
# as.data.frame(tbl) %>% 
#   dplyr::rename(is_doublet = Var1, n_cells = Freq) %>%
#   write_csv(args$output_csv)
# 
# 
# if(file.exists(summary_csv)){
#   summary_results <- read.csv(summary_csv, check.names=F, stringsAsFactors=F)
#   # There are low quality cells that we removed here, we suppose that they are not mouse cells
#   summary_results$is_doublet <- NULL
#   summary_results <- summary_results %>%
#     dplyr::mutate(is_doublet = case_when(Barcode %in% doublet_cells ~ T, 
#                                        !Barcode %in% doublet_cells ~ F)) %>% 
#     write_csv(summary_csv)
#   
# }else{
#   summary_results <- as.data.frame(colData(sce)) %>% 
#     write_csv(summary_csv)
# }
# sce <- sce[, !sce$is_doublet]
# print('Output of filtered data after removing doublets: ')
# print(dim(sce))
# saveRDS(sce, file = args$sce_output)
# cat("Completed.\n")


remove_doublet <- function(sce_file, sce_output, 
                           summary_csv=NULL, output_csv=NULL, doublet_score=0){
  doublet_score <- as.numeric(doublet_score)
  if(doublet_score==0){
    doublet_score <- NULL
  }
  lid <- basename(sce_file)
  lid <- gsub('_filtered','',lid)
  lid <- gsub('.rds','',lid)
  if(is.null(summary_csv)){
    summary_csv <- paste0(dirname(sce_file),'/',lid,'_qc.csv')
  }
  if(is.null(output_csv)){
    output_csv <- paste0(dirname(sce_file),'/',lid,'_doublet_filtered.csv')
  }
  # summary_csv <- args$summary_csv
  sce <- readRDS(sce_file)
  print(dim(sce))
  
  doublets_df <- rscrublet::scrubDoublets(as.matrix(counts(sce)), n_neighbors = 20, 
                                          doublet_score_threshold = doublet_score,
                                          sim_doublet_ratio = 2, expected_doublet_rate = 0.06,
                                          stdev_doublet_rate = 0.02, synthetic_doublet_umi_subsampling = 1,
                                          use_approx_neighbors = TRUE, distance_metric = "euclidean",
                                          min_counts = 2, min_cells = 3, min_gene_variability_pctl = 85,
                                          log_transform = FALSE, z_score = TRUE, n_prin_comps = 30,
                                          verbose = TRUE)
  
  print(summary(doublets_df$scrubDoublets))
  print(summary(doublets_df$doublet_scores_obs))
  print(summary(doublets_df$doublet_scores_sim))
  
  sce$is_doublet <- doublets_df$scrubDoublets
  sce$doublet_scores_obs <- round(doublets_df$doublet_scores_obs,2)
  
  doublet_cells <- sce$Barcode[sce$is_doublet]
  print(paste0('Number of doublet cells: ', length(doublet_cells)))
  
  tbl <- table(sce$is_doublet)
  as.data.frame(tbl) %>% 
    dplyr::rename(is_doublet = Var1, n_cells = Freq) %>%
    write_csv(output_csv)
  
  
  if(file.exists(summary_csv)){
    summary_results <- read.csv(summary_csv, check.names=F, stringsAsFactors=F)
    # There are low quality cells that we removed here, we suppose that they are not mouse cells
    summary_results$is_doublet <- NULL
    summary_results <- summary_results %>%
      dplyr::mutate(is_doublet = case_when(Barcode %in% doublet_cells ~ T, 
                                           !Barcode %in% doublet_cells ~ F)) %>% 
      write_csv(summary_csv)
    
  }else{
    summary_results <- as.data.frame(colData(sce)) %>% 
      write_csv(summary_csv)
  }
  sce <- sce[, !sce$is_doublet]
  print('Output of filtered data after removing doublets: ')
  print(dim(sce))
  saveRDS(sce, file = sce_output)
  cat("Completed.\n")
  
}

remove_doublet(args$sce_file, args$sce_output, args$summary_csv, args$output_csv, doublet_score=0)


# /usr/local/bin/Rscript /home/htran/Projects/farhia_project/rnaseq/pipeline/utils/identify_doublet_v2.R --sce_file /home/htran/storage/images_dataset/merfish_rnaseq/SITTA12_XP6873/analysis/filtered/SITTA12_XP6873_filtered.rds --sce_output /home/htran/Projects/farhia_project/rnaseq/pipeline/utils/identify_doublet_v2.R --sce_file /home/htran/storage/images_dataset/merfish_rnaseq/SITTA12_XP6873/analysis/filtered/SITTA12_XP6873_f_dt.rds --doublet_score 0 --output_csv /home/htran/Projects/farhia_project/rnaseq/pipeline/utils/identify_doublet_v2.R --sce_file /home/htran/storage/images_dataset/merfish_rnaseq/SITTA12_XP6873/analysis/filtered/SITTA12_XP6873_filtered.rds --sce_output /home/htran/Projects/farhia_project/rnaseq/pipeline/utils/identify_doublet_v2.R --sce_file /home/htran/storage/images_dataset/merfish_rnaseq/SITTA12_XP6873/analysis/filtered/SITTA12_XP6873__doublet_filtered.csv --summary_csv /home/htran/storage/images_dataset/merfish_rnaseq/SITTA12_XP6873/analysis/filtered/SITTA12_XP6873_qc.csv
# lids <- c('SCRNA10X_SA_CHIP0239_001','SCRNA10X_SA_CHIP0239_005',
#           'SCRNA10X_SA_CHIP0239_006',#'SCRNA10X_SA_CHIP0239_003',
#           'SCRNA10X_SA_CHIP0239_004','SCRNA10X_SA_CHIP0239_007','SCRNA10X_SA_CHIP0239_008')
# 
# output_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/filtered/'
# 
# for(library_id in lids){
#   sce_file <- paste0(output_dir,library_id,'_filtered.rds')
#   sce_output <- paste0(output_dir,library_id,'_f_dt.rds')
#   print(library_id)
#   remove_doublet(sce_file, sce_output, summary_csv=NULL, output_csv=NULL, doublet_score=0)
#   
# }

# sce_path <- '/home/htran/storage/datasets/hakwoo_metastasis_RNAseq/testing/human/filtered/SCRNA10X_SA_CHIP0146_004_mouse_filtered.rds'
# summary_csv <- '/home/htran/storage/datasets/hakwoo_metastasis_RNAseq/testing/human/filtered/SCRNA10X_SA_CHIP0146_004_qc.csv'

## example of command: 
# /usr/local/bin/Rscript /home/htran/Projects/farhia_project/rnaseq/pipeline/utils/identify_doublet_v2.R 
# --sce_file /home/htran/storage/datasets/metastasis_results/rnaseq_SA535/filtered/SCRNA10X_SA_CHIP0220_001_filtered.rds 
# --sce_output /home/htran/storage/datasets/metastasis_results/rnaseq_SA535/filtered/SCRNA10X_SA_CHIP0220_001_f_dt.rds --doublet_score 0 --output_csv /home/htran/storage/datasets/metastasis_results/rnaseq_SA535/filtered/SCRNA10X_SA_CHIP0220_001_doublet_filtered.csv --summary_csv /home/htran/storage/datasets/metastasis_results/rnaseq_SA535/filtered/SCRNA10X_SA_CHIP0220_001_qc.csv


# output_csv <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/filtered/SCRNA10X_SA_CHIP0220_001_doublet_filtered.csv'
# sce_output <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/filtered/SCRNA10X_SA_CHIP0220_001_f_dt.rds'
# summary_csv <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/filtered/SCRNA10X_SA_CHIP0220_001_qc.csv'
# sce_file <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/filtered/SCRNA10X_SA_CHIP0220_001_filtered.rds'
# 
# remove_doublet(sce_file, sce_output, summary_csv=NULL, output_csv=NULL, doublet_score=0)
