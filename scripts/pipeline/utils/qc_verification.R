suppressPackageStartupMessages({
  require("optparse")
  require("scater")
  require("argparse")
  require("SingleCellExperiment")
  require("stringr")
  require("tidyverse")
  require("scran")
  require("Seurat")
  
})


option_list <- list(make_option(c("-l", "--library_ids"), type="character", default=NULL, help="library_ids", metavar="character"),
                    make_option(c("-i", "--input_dir"), type="character", default=NULL, help="input_dir", metavar="character"),
                    make_option(c("-o", "--output_file"), type="character", default=NULL, help="output_file", metavar="character"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
print(opt$library_ids)
print(opt$input_dir)
print(opt$output_file)


get_feature_QC <- function(sce){
  
  mito_genes <- str_detect(rowData(sce)$Symbol, "^MT\\-")
  sum(mito_genes==TRUE)
  
  
  # For scater version >= 1.14.6
  per.cell <- perCellQCMetrics(sce, subsets=list(Mito=mito_genes))
  colData(sce) <- cbind(colData(sce), per.cell)
  
  retrieved_cells <- list()
  mito_thrs_ls <- c(20, 25, 30)
  for(max_mito in mito_thrs_ls){
    keep_mt <- sce$subsets_Mito_percent <= max_mito
    mt_max <- paste0('mt_',max_mito)
    retrieved_cells[[mt_max]] <- sum(keep_mt)
  }
  retrieved_cells[['median_genes']] <- median(sce$total_features_by_counts)
  retrieved_cells[['mean_reads']] <- mean(sce$sum)
  retrieved_cells[['max_sum']] <- max(sce$sum)
  return(retrieved_cells)
  
}
# results_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA919/'
# sce <- readRDS(paste0(results_dir,'normalized/SA919_raw.rds'))
# unique(sce$library_label)
# summary(as.factor((sce$library_label)))
# # sce <- readRDS(paste0(results_dir,'SCRNA10X_SA_CHIP0078_001_filtered.rds'))
# download_dir <- '/home/htran/storage/datasets/hakwoo_metastasis_RNAseq/'
# raw_sce <- readRDS(paste0(download_dir,'SCRNA10X_SA_CHIP0078_001/','SCRNA10X_SA_CHIP0078_001.rdata'))
# retrieved_cells <- get_feature_QC(raw_sce)
# max(raw_sce$sum)
# min(sce$sum)
# retrieved_cells$max_sum

verification_QC <- function(library_ids, input_dir, output_file){
  
  output_dir <- paste0(dirname(output_file),'/')
  if (!file.exists(output_dir)){
    dir.create(output_dir)
  }
  mito_20 <- c()
  mito_25 <- c()
  mito_30 <- c()
  mean_reads <- c()
  mean_reads_filtered <- c()
  median_genes <- c()
  median_genes_filtered <- c()
  our_out <- c()
  dev_out <- c()
  max_reads <- c()
  max_reads_filtered <- c()
  total_cells <- c()
  
  for (f in library_ids){
    # fn <- str_sub(f, 1, str_length(f)-19)  #-16
    print(paste0("Processing file:  ",f))
    sce <- readRDS(paste0(input_dir, '/', f,'/',f,'.rdata'))
    sce_dev <- readRDS(paste0(input_dir, '/human/dev/', f,'/',f,'.rdata'))
    sce_ours <- readRDS(paste0(output_dir, f,'_filtered.rds'))
    dev_out <- c(dev_out,dim(sce_dev)[2])
    retrieved_cells <- get_feature_QC(sce)
    mito_20 <- c(mito_20, retrieved_cells[['mt_20']])
    mito_25 <- c(mito_25, retrieved_cells[['mt_25']])
    mito_30 <- c(mito_30, retrieved_cells[['mt_30']])
    our_out <- c(our_out, dim(sce_ours)[2])
    total_cells <- c(total_cells, dim(sce)[2])
    median_genes <- c(median_genes, retrieved_cells[['median_genes']])
    median_genes_filtered <- c(median_genes_filtered, median(sce_ours$total_features_by_counts))
    mean_reads <- c(mean_reads, retrieved_cells[['mean_reads']])
    mean_reads_filtered <- c(mean_reads_filtered, mean(sce_ours$sum))
    max_reads <- c(max_reads, retrieved_cells[['max_sum']])
    max_reads_filtered <- c(max_reads_filtered, max(sce_ours$sum))
  }
  summary_df <- data.frame(library_id=library_ids,mito_20=mito_20,mito_25=mito_25,
                           mito_30=mito_30, our_out=our_out, dev_team_out=dev_out, total_cells, 
                           median_genes_in_wholedata=median_genes,
                           median_genes_filtered_data=median_genes_filtered,
                           mean_reads_in_wholedata=round(mean_reads),
                           mean_reads_filtered_data=round(mean_reads_filtered),
                           max_reads_in_wholedata=max_reads)
  write.csv(summary_df, file = output_file, quote=F, row.names = F)
  
}

library_ids = strsplit(opt$library_ids, ",")[[1]]
print(library_ids)
verification_QC(library_ids, opt$input_dir, opt$output_file)
