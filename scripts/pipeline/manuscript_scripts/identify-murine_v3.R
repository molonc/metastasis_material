## Based on script from Kieran Campbell
# Just put into function
# Hoa Tran
suppressPackageStartupMessages({
  library(tidyverse)
  library(SingleCellExperiment)
  library(scater)
  library(argparse)
  library(stringr)
  library(dplyr)
})

parser <- ArgumentParser(description = "Detecting mouse cells that are mixed with human cells")
parser$add_argument('--sce_human', metavar='FILE', type='character', 
                    help="Path to human SCE")

parser$add_argument('--sce_mouse', metavar='FILE', type='character', 
                    help="Path to mouse SCE")

parser$add_argument('--sample', type='character',
                    help="Sample name")

parser$add_argument('--output_figure', metavar='FILE', type = 'character',
                    help="Output figure file")

parser$add_argument('--output_csv', metavar='FILE', type = 'character',
                    help="Output dataframe with summary of # mouse cells")

parser$add_argument('--summary_results_csv', metavar='FILE', type = 'character',
                    help="Output dataframe with summary of # mouse cells for total data")

parser$add_argument('--output_sce', type = 'character', metavar = 'FILE',
                    help="Output path for SCE with mouse cells removed")

args <- parser$parse_args()

## Get rid of stuff that breaks rowRanges
unbreak_rowranges <- function(sce) {
  rd_to_rename <- c("start", "end", "chr","strand")
  for(r in rd_to_rename) {
    if(r %in% names(rowData(sce))) {
      rowData(sce)[[r]] <- NULL
    }
  }
  sce
}

get_metric_df <- function(sce) {
  as.data.frame(colData(sce)) %>% 
    as_tibble() %>% 
    dplyr::select(Barcode, total_counts)
}
qc_metrics <- function(sce) {
  mito_genes <- str_detect(rowData(sce)$Symbol, "^MT\\-")
  per.cell <- perCellQCMetrics(sce, subsets=list(Mito=mito_genes))
  colData(sce) <- cbind(colData(sce), per.cell)
  sce
}  
#args$output_figure
detect_mouse_cells <- function(sce_human, sce_mouse){
  obs_qc_features <- c('Barcode', 'total_counts')
  if(sum(obs_qc_features %in% colnames(colData(sce_human)))!=2){
    # sce_human <- calculateQCMetrics(sce_human)  #old version of scater package
    sce_human <- qc_metrics(sce_human)
  }
  
  if(sum(obs_qc_features %in% colnames(colData(sce_mouse)))!=2){
    # sce_mouse <- calculateQCMetrics(sce_mouse)
    sce_mouse <- qc_metrics(sce_mouse)
  }
  
  df_hm <- inner_join(
    get_metric_df(sce_human),
    get_metric_df(sce_mouse),
    by = "Barcode",
    suffix = c("_human", "_mouse")
  )
  
  df_hm <- dplyr::mutate(df_hm, is_mouse = total_counts_mouse > total_counts_human)
  # 
  # plt <- ggplot(df_hm, aes(x = total_counts_human, y = total_counts_mouse, colour = is_mouse)) + 
  #   geom_point() +
  #   geom_abline()
  # 
  # ggsave(output_figure, plt, width = 6, height = 5)
  
  mouse_bcs <- dplyr::filter(df_hm, is_mouse) %>% .$Barcode
  print(paste0('There are ',length(mouse_bcs),' mouse cells in this library'))
  return(mouse_bcs)
  
}

filter_mouse_cells <- function(datatag, sce_human_fn, sce_mouse_fn, output_dir){
  
  if (!file.exists(output_dir)){
    dir.create(output_dir)
  }
  if(file.exists(sce_human_fn) & file.exists(sce_mouse_fn)){
    sce_human <- readRDS(sce_human_fn)
    print(paste0('Human sce file: ', dim(sce_human)[1],' ',dim(sce_human)[2]))
    sce_mouse <- readRDS(sce_mouse_fn)
    print(paste0('Mouse sce file: ', dim(sce_mouse)[1],' ',dim(sce_mouse)[2]))
    
    sce_human <- unbreak_rowranges(sce_human)
    sce_mouse <- unbreak_rowranges(sce_mouse)
    
    # If cells are well aligned to mouse than human genome --> mouse cells
    mouse_bcs <- detect_mouse_cells(sce_human, sce_mouse)
    colnames(colData(sce_human))
    mouse_cells_df <- as.data.frame(colData(sce_human)) %>%
      dplyr::select(Barcode, Sample) %>%
      dplyr::filter(Barcode %in% mouse_bcs)
    print(paste0('There are ',dim(mouse_cells_df)[1],' mouse cells in the sample: ', datatag))
    data.table::fwrite(mouse_cells_df, paste0(output_dir, datatag, '_mouse_cells.csv.gz'))
  }else{
    print('---------------------------------------')
    print(paste0(datatag, ': do not exist data for this sample, double check!!! ',sce_human_fn))
  }
  
  
}

input_dir <- '/home/htran/storage/rnaseq_datasets/hakwoo_metastasis_RNAseq/'
fns <- basename(list.dirs(paste0(input_dir, 'SA535_human_exons/')))
fns <- fns[grepl('SCRNA', fns)]
# datatag <- 'SCRNA10X_SA_CHIP0176_001'
for(datatag in fns){
  sce_human_fn <- paste0(input_dir, 'SA535_human_exons/', datatag, '/', datatag, '.rdata')
  sce_mouse_fn <- paste0(input_dir, 'SA535_mouse_exons/', datatag, '/', datatag, '.rdata')
  output_dir <- paste0(input_dir, 'mouse_cells/')
  filter_mouse_cells(datatag, sce_human_fn, sce_mouse_fn, output_dir)
}


# sce_human_fn <- '/home/htran/storage/datasets/hakwoo_metastasis_RNAseq/testing/human/filtered/SCRNA10X_SA_CHIP0146_004_filtered.rdata'
# sce_mouse_fn <- '/home/htran/storage/datasets/hakwoo_metastasis_RNAseq/testing/mouse/SCRNA10X_SA_CHIP0146_004/SCRNA10X_SA_CHIP0146_004.rdata'
# output_figure <- '/home/htran/storage/datasets/hakwoo_metastasis_RNAseq/testing/human/filtered/SCRNA10X_SA_CHIP0146_004_mouse_filtered.png'
# output_csv <- '/home/htran/storage/datasets/hakwoo_metastasis_RNAseq/testing/human/filtered/SCRNA10X_SA_CHIP0146_004_mouse_filtered.csv'
# output_sce <- '/home/htran/storage/datasets/hakwoo_metastasis_RNAseq/testing/human/filtered/SCRNA10X_SA_CHIP0146_004_mouse_filtered.rds'
# summary_results_csv <- '/home/htran/storage/datasets/hakwoo_metastasis_RNAseq/testing/human/filtered/SCRNA10X_SA_CHIP0146_004_qc.csv'
# filter_mouse_cells(sce_human_fn, sce_mouse_fn, output_figure, output_csv, output_sce, summary_results_csv)
  
# sce_human <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/filtered/SCRNA10X_SA_CHIP0220_001_f_dt.rds'
# sce_mouse <- '/home/htran/storage/rnaseq_datasets/hakwoo_metastasis_RNAseq/SA535_mouse/SCRNA10X_SA_CHIP0220_001/SCRNA10X_SA_CHIP0220_001.rdata'
# output_figure <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/filtered/mouse_filter_info/SCRNA10X_SA_CHIP0220_001_mouse_filtered.png'
# output_csv <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/filtered/mouse_filter_info/SCRNA10X_SA_CHIP0220_001_mouse_filtered.csv'
# output_sce <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/filtered/SCRNA10X_SA_CHIP0220_001_f_dt_m.rds'
# summary_results_csv <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/filtered/SCRNA10X_SA_CHIP0220_001_qc.csv'
# 
# filter_mouse_cells(sce_human, sce_mouse, 
#                    output_figure, output_csv, 
#                    output_sce, summary_results_csv)