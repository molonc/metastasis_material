# BiocManager::install("scMerge")
# Rscript compute_scSEG_stable_genes.R 
# -i /home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/SA609-v6/total_sce.rds 
# -o /home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/normalization_evaluation/segIndx_total.csv

suppressPackageStartupMessages({
  require("optparse")
  require("scater")
  require("argparse")
  require("SingleCellExperiment")
  require("scMerge")
  require("dplyr")
  # require("tidyverse")
  
})


option_list <- list(make_option(c("-i", "--input_file"), type="character", default=NULL, help="input_dir", metavar="character"),
                    make_option(c("-o", "--output_file"), type="character", default=NULL, help="output_file", metavar="character"),
                    make_option(c("-p", "--probIdx"), type="double", default=0.8, help="prob_SEG_Idx"),
                    make_option(c("-f", "--probFactor"), type="double", default=0.1, help="prob_Factor"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
print(opt$input_file)
print(opt$output_file)
print(opt$probIdx)
print(opt$probFactor)

source(paste0("/home/htran/Projects/farhia_project/rscript/pipeline/utils/normalize_utils.R"))

getSEGIndex <- function(input_file, output_file, probIdx=0.8, probFactor=0.1){
  sce <- readRDS(input_file)
  print("Original matrix: ")
  print(dim(sce))
  print(sce$series[1])
  sce$treatmentSt <- get_treatment_status(sce$series)
  metadata <- as.data.frame(colData(sce))
  metadata$cell_id <- colnames(sce)
  ext_cells <- downsampling_data(metadata, col_use='treatmentSt', downsample_ratio=0.6, thres_small_clone=600)
  sce <- sce[,ext_cells]
  output_dir <- paste0(dirname(output_file),'/')
  if (!file.exists(output_dir)){
    dir.create(output_dir)
  }
  print(output_dir)
  meta_genes <- data.frame(gene_ens=rownames(rowData(sce)),
                           gene_symb=rowData(sce)$Symbol, 
                           row.names = rownames(rowData(sce)), stringsAsFactors = F)
  print(dim(meta_genes))
  # get logcounts of sce
  # sce2 <- SingleCellExperiment(assays = list(counts = as.matrix(counts(sce)), logcounts=as.matrix(log2(counts(sce)+1))))
  exprs_mat = SummarizedExperiment::assay(sce, 'counts')
  # exprs_mat = SummarizedExperiment::assay(sce2, 'logcounts')
  print(dim(exprs_mat))
  # set.seed(1)
  param = BiocParallel::MulticoreParam(workers = 3, progressbar = FALSE)
  segIndx_df = scMerge::scSEGIndex(exprs_mat = exprs_mat, BPPARAM = param)
  ## Closing the parallelisation
  BiocParallel::register(BPPARAM = BiocParallel::SerialParam())
  segIndx_df$gene_ens <- rownames(segIndx_df)
  segIndx_df$gene_symb <- meta_genes[rownames(segIndx_df),'gene_symb']
  print("Select SEG genes total: ")
  print(dim(segIndx_df))
  saveRDS(segIndx_df, file=paste0(output_dir,'segIndx_total.rds'))
  write.csv(segIndx_df, file=paste0(output_dir,'segIndx_total.csv'), quote = F, row.names = F)
  segIndx_df1 <- segIndx_df %>% dplyr::filter(segIdx > quantile(x=segIdx, probs = probIdx) &
                                              rho > quantile(x=rho, probs = probFactor) &
                                              mu > quantile(x=mu, probs = probFactor) & 
                                              sigma > quantile(x=sigma, probs = probFactor))
  print("Select SEG genes with strict condition: ")
  print(dim(segIndx_df1))
  
  segIndx_df2 <- segIndx_df %>% dplyr::filter(segIdx > quantile(x=segIdx, probs = probIdx))
  print("Select SEG genes with threshold percentile of segIdx: ")
  print(dim(segIndx_df2))
  if(dim(segIndx_df1)[1]>0){
    write.csv(segIndx_df1, file=paste0(output_dir,'segIndx_strict_condition.csv'), quote = F, row.names = F)
  }
  if(dim(segIndx_df2)[1]>0){
    write.csv(segIndx_df2, file=paste0(output_dir,'segIndx_filtered_',(probIdx*100),'.csv'), quote = F, row.names = F)
  }  
  
}  
getSEGIndex(opt$input_file, opt$output_file, opt$probIdx, opt$probFactor)

#' #' ## Loading example data
# data('example_sce', package = 'scMerge')
# class(logcounts(example_sce))
#' ## subsetting genes to illustrate usage.

# df <- data("segList", package = "scMerge")
# data("segList_ensemblGeneID", package = "scMerge")
# length(segList_ensemblGeneID)
# length(segList_ensemblGeneID$human$human_scSEG)
# ref_genes <- segList_ensemblGeneID$human$human_scSEG
# ref_genes[1:3]
# input_dir = '/home/htran/storage/datasets/metastasis_results/rnaseq_SA919/'
# obs_seg <- read.csv(paste0(input_dir,'batch_effect/segIndx_df_filtered_90.csv'), stringsAsFactors = F, check.names = F)
# dim(obs_seg)
# res <- intersect(obs_seg$gene_ens, ref_genes)
# length(res)
# 
# obs_seg$gene_ens[1:3]
# dim(obs_seg)
# normalized_sce <- readRDS(paste0(input_dir,'normalized/SA919_normalized_output.rds'))
# dim(normalized_sce)
# res1 <- intersect(rownames(normalized_sce), ref_genes)
# length(res1)
# res2 <- intersect(rownames(normalized_sce), obs_seg$gene_ens)
# length(res2)
# res <- intersect(obs_seg$gene_ens, ref_genes)







