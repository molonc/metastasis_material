## Based on script from Kieran Campbell
# Just put into function
# Hoa Tran
suppressPackageStartupMessages({
  library(tidyverse)
  library(SingleCellExperiment)
  library(scater)
  library(argparse)
  library(stringr)
})

parser <- ArgumentParser(description = "Integrate SCEs")
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
  rd_to_rename <- c("start", "end", "strand")
  for(r in rd_to_rename) {
    if(r %in% names(rowData(sce))) {
      rowData(sce)[[r]] <- NULL
    }
  }
  sce
}

filter_mouse_cells <- function(sce_human_fn, sce_mouse_fn, 
                               output_figure, output_csv, output_sce, 
                               summary_results_csv){
  output_dir <- paste0(dirname(output_sce),'/')
  if (!file.exists(output_dir)){
    dir.create(output_dir)
  }
  sce_human <- readRDS(sce_human_fn)
  print(paste0('Human sce file: ', dim(sce_human)[1],' ',dim(sce_human)[2]))
  sce_mouse <- readRDS(sce_mouse_fn)
  print(paste0('Mouse sce file: ', dim(sce_mouse)[1],' ',dim(sce_mouse)[2]))
  ## Check if sce_mouse is missing and write output files if so
  if(is.character(sce_mouse)) {
    if(sce_mouse == "empty") {
      # Just write SCE human
      sce_human <- unbreak_rowranges(sce_human)
      
      sce_human$is_mouse <- FALSE
      
      ggplot()
      ggsave(output_figure, width = 6, height = 5)
      
      tibble() %>% 
        write_csv(output_csv)
      saveRDS(sce_human, file = output_sce)
      
      quit()
    }
  }
  
  if(nrow(sce_human) > 40000) {
    stop("Human SCE has > 40k genes: is this aligned to mouse hybrid?")
  }
  if(nrow(sce_mouse) > 40000) {
    stop("Mouse SCE has > 40k genes: is this aligned to mouse hybrid?")
  }
  
  sce_human <- unbreak_rowranges(sce_human)
  sce_mouse <- unbreak_rowranges(sce_mouse)
  
  # If cells are well aligned to mouse than human genome --> mouse cells
  mouse_bcs <- detect_mouse_cells(sce_human, sce_mouse, output_figure)
  
  sce_human <- readRDS(sce_human_fn)
  
  sce_human$is_mouse <- sce_human$Barcode %in% mouse_bcs
  
  tbl <- table(sce_human$is_mouse)
  as.data.frame(tbl) %>% 
    dplyr::rename(is_mouse = Var1, n_cells = Freq) %>%
    # dplyr::mutate(sample = args$sample) %>% 
    # dplyr::select(sample, everything()) %>% 
    write_csv(output_csv)
  
  if(file.exists(summary_results_csv)){
    summary_results <- read.csv(summary_results_csv, check.names=F, stringsAsFactors=F)
    # There are low quality cells that we removed here, we suppose that they are not mouse cells
    summary_results$is_mouse <- NULL
    summary_results <- summary_results %>%
      dplyr::mutate(is_mouse = case_when(Barcode %in% mouse_bcs ~ T, 
                                         !Barcode %in% mouse_bcs ~ F))%>% 
      write_csv(summary_results_csv)
    
  }else{
    summary_results <- as.data.frame(colData(sce_human)) %>% 
      write_csv(summary_results_csv)
  }
  
  sce_human <- sce_human[, !sce_human$is_mouse]
  print(dim(sce_human))
  # Write outputs
  saveRDS(sce_human, file = output_sce)
  cat("Completed.\n")
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
detect_mouse_cells <- function(sce_human, sce_mouse, output_figure){
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
  
  plt <- ggplot(df_hm, aes(x = total_counts_human, y = total_counts_mouse, colour = is_mouse)) + 
    geom_point() +
    geom_abline()
  
  ggsave(output_figure, plt, width = 6, height = 5)
  
  mouse_bcs <- dplyr::filter(df_hm, is_mouse) %>% .$Barcode
  print(paste0('There are ',length(mouse_bcs),' mouse cells in this library: ',basename(output_figure)))
  return(mouse_bcs)
  
}

# Params:
# args$sce_human
# args$sce_mouse
# args$output_figure
# args$output_csv
# args$output_sce
filter_mouse_cells(args$sce_human, args$sce_mouse, 
                   args$output_figure, args$output_csv, 
                   args$output_sce, args$summary_results_csv)



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