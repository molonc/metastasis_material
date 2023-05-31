suppressPackageStartupMessages({
  require("optparse")
  require("scater")
  require("argparse")
  require("SingleCellExperiment")
  require("stringr")
  require("tidyverse")
  require("scran")
  require("Seurat")
  # require("reticulate")
  
})

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

source(paste0(script.basename, "/utils/normalize_utils.R"))

option_list <- list(make_option(c("-i", "--input_file"), type="character", default=NULL, help="input_dir", metavar="character"),
                    make_option(c("-o", "--output_file"), type="character", default=NULL, help="output_file", metavar="character"),
                    make_option(c("-r", "--resolution"), type="double", default=0.2, help="grain_graph_cut"),
                    make_option(c("-d", "--datatag"), type="character", default='SA', help="basename", metavar="character"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
print(opt$input_file)
print(opt$output_file)
print(opt$resolution)
print(opt$datatag)

get_clusters <- function(input_file, output_file, resolution, datatag='SA'){
  output_dir <- paste0(dirname(output_file),'/')
  if (!file.exists(output_dir)){
    dir.create(output_dir)
  }
  # base_name <- basename(output_file)
  # base_name <- gsub('_normalized_clusters.rds','',base_name)
  # print(paste0("base_name is: ", base_name))
  
  sce <- readRDS(input_file)
  print(dim(sce))
  print(output_dir)
  
  # sce$pdxid <- gsub('^X0847_','',sce$pdxid)
  sce_cluster <- get_seurat_clusters(base_name=datatag, sce=sce, output_file=output_file, 
                                     save_dir=output_dir, res=resolution, 
                                     dims = 1:15, saveSrt = TRUE, plotRD=TRUE)
  run_QC_Cluster(datatag, sce_cluster, output_dir)
}




get_clusters(opt$input_file, opt$output_file, opt$resolution, opt$datatag)


# results_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA919/'
# sce_cluster <- readRDS(paste0(results_dir,'clustering/SA919_normalized_clusters.rds'))
# run_QC_Cluster(sce_cluster, output_dir = paste0(results_dir,'clustering/'))
# var_genes <- readRDS(paste0(results_dir,'clustering/var_genes.rds'))
# length(var_genes)
# var_genes[1:100]
# sce1 <- readRDS(paste0(results_dir, "clustering/SA919_sans_mito_ribo_genes.rds"))
# sce1 <- sce1[var_genes,]
# dim(sce1)
# rowData(sce1)$Symbol[1:3]
# saveRDS(sce1, file = paste0(results_dir, "clustering/SA919_var_genes.rds"))
# saveRDS(sce1, file = paste0(results_dir, "clustering/SA919_var_genes_sans5.rds"))
# srt <- readRDS(paste0(results_dir, "clustering/SA919_clustering_srt.rds"))
# sce1 <- sce1[,sce1$cluster_label!="5"]
# cells_use <- colnames(sce1)
# length(cells_use)
# srt <- srt[,cells_use]  # subset cells
