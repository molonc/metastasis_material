# install.packages("remotes")
# remotes::install_github("huyuem/rscrublet")

# htran
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


doublet_score <- as.numeric(args$doublet_score)
if(doublet_score==0){
  doublet_score <- NULL
}
summary_csv <- args$summary_csv
sce <- readRDS(args$sce_file)
print(dim(sce))

doublets_df <- rscrublet::scrubDoublets(as.matrix(counts(sce)), n_neighbors = 15, 
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
  write_csv(args$output_csv)


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
saveRDS(sce, file = args$sce_output)
cat("Completed.\n")




# sce_path <- '/home/htran/storage/datasets/hakwoo_metastasis_RNAseq/testing/human/filtered/SCRNA10X_SA_CHIP0146_004_mouse_filtered.rds'
# summary_csv <- '/home/htran/storage/datasets/hakwoo_metastasis_RNAseq/testing/human/filtered/SCRNA10X_SA_CHIP0146_004_qc.csv'
