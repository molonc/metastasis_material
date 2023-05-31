suppressPackageStartupMessages({
  # require("scater")
  require("SingleCellExperiment")
  require("stringr")
  require("tidyverse")
  require("scran")
  require("scuttle") # remotes::install_github("LTLA/scuttle")
})


script_dir <- '/home/htran/Projects/farhia_project/rnaseq/pipeline/utils/'
source(paste0(script_dir, "normalize_utils.R"))

# rlog: log=TRUE, FALSE in logNormCounts function, Logical scalar indicating 
# whether normalized values should be log2-transformed.if rlog=FALSE - save data to normcounts, rlog=TRUE - save to logcounts
sce_normalize_size_factors <- function(sce, min_size=100, rlog=FALSE, exprs="counts", name){
  # library(scater)
  print("Quick clustering")
  if(min_size < dim(sce)[2]){
    qclust <- scran::quickCluster(sce, min.size = min_size, assay.type=exprs)
    print("Compute sum factors")
    sce <- scran::computeSumFactors(sce, clusters = qclust, assay.type=exprs)
    sce$size_factor <- sizeFactors(sce)
    print("Normalize data")
    # scater version 1.14.6
    # count --> normcounts
    # normcounts --> logcounts 
    # String containing an assay name for storing the output normalized values. 
    # Defaults to "logcounts" when log=TRUE and "normcounts" otherwise.
    sce_normalized <- scuttle::logNormCounts(sce, log=rlog, exprs_values=exprs, assay.type=exprs,
                                             name=name, size_factors=sce$size_factor)
    
  } else{
    print("Smaller nb cells than min size threshold in this library")
    sce_normalized <- scuttle::logNormCounts(sce, log=rlog, exprs_values=exprs, name=name, size_factors=NULL)
  }
  
  print(assayNames(sce_normalized))
  return(sce_normalized)
}

# sce_cbind_func <- function(sce_list, exprs = c("counts", "logcounts","normcounts"), 
#                               colData_names = NULL) {
#   # Combining expression
#   assay_list <- list()
#   for (i in seq_len(length(exprs))) {
#     assay_list[[i]] <- do.call(cbind, lapply(sce_list, 
#                                              function(y) assay(y, exprs[i])))
#   }
#   names(assay_list) <- exprs
#   
#   # Combining metadata
#   colData_list <- do.call(DelayedArray::rbind, 
#                           lapply(sce_list, function(y) colData(y)[, colData_names, drop = FALSE]))
#   sce_combine <- SingleCellExperiment::SingleCellExperiment(assay = assay_list, 
#                                                             colData = colData_list)
#   
#   print(paste0("Dim sce combine: ",dim(sce_combine)[1],' ',dim(sce_combine)[2]))
#   # print(colnames(colData(sce_combine)))
#   print(paste0("sce combine assay name: ",assayNames(sce_combine)))
#   return(sce_combine)
# }


twice_scran_normalize <- function(library_ids, input_dir, output_dir, datatag='SA'){
  #library_ids = strsplit(library_ids_ls, ",")[[1]]
  print(library_ids)
  #output_dir <- paste0(dirname(output_file),'/')
  if (!file.exists(output_dir)){
    dir.create(output_dir)
  }
  sce_list <- list()
  saved_rowdata <- list()
  for (f in library_ids){
    print(paste0("Processing file:  ",f))
    
    sce <- readRDS(paste0(input_dir, '/', f,'.rds'))
    saved_rowdata <- c(saved_rowdata, rowData(sce))
    print(dim(sce))
    if(dim(sce)[2]>0){
      # First time normalize data using scran, input is counts exp, normalized output are kept in normcounts exp values 
      sce$sample <- f   # adding the name so I can save the individual files later
      # rlog=FALSE means the normalized counts are NOT log transformed
      sce_normalized <- sce_normalize_size_factors(sce, min_size=100, rlog=FALSE, exprs="counts", name="normcounts")  ##  rlog=FALSE: count --> normcounts
      if(!is.null(sce_normalized)){  
        sce_list <- c(sce_list, sce_normalized)
      } else{
        print(paste0("Please double check library id: ",f))
      }
    }
  }
  # Combining all sces to 1 total sce
  # Take out the tenx column
  sce_normalized$tenx <- NULL
  sce_combine <- sce_cbind_func(sce_list, exprs = c("counts", "normcounts"),
                                colData_names = colnames(colData(sce_normalized)))
  
  # Twice scater normalization
  # normalized output are kept in logcounts exp values
  print(("Normalizing total data..."))
  sce_normalized_total <- sce_normalize_size_factors(sce_combine, min_size=300, rlog=TRUE, exprs="normcounts", name="logcounts") ##  rlog=TRUE: normcounts --> logcounts
  saveRDS(sce_normalized_total, file = paste0(output_dir,"/allsce_normalized.rds"))
  # Saving individual files
  j <- 1
  for (sample in library_ids) {
    print(sample)
    sample_level <- sce_normalized_total[,sce_normalized_total$sample==sample]
    rowData(sample_level) <- saved_rowdata[[j]]
    print(sample_level)
    saveRDS(sample_level, file=paste0(output_dir,"/",sample,".rdata"))
    j <- j+1
  }
}  

twice_scran_normalize_v2 <- function(sce, input_dir, output_dir, datatag, return_data=F){
  #library_ids = strsplit(library_ids_ls, ",")[[1]]
  # print(library_ids)
  #output_dir <- paste0(dirname(output_file),'/')
  if (!file.exists(output_dir)){
    dir.create(output_dir)
  }
  # sce <- sce_combine
  print(dim(sce))
  # colnames(colData(sce))
  # sce$sample[1:3]
  sce_list <- list()
  for (s in unique(sce$library_id)){
    print(paste0("Processing sample:  ",s))
    
    # sce <- readRDS(paste0(input_dir, '/', f,'.rds'))
    sce_tmp <- sce[,sce$library_id==s] 
    # saved_rowdata <- c(saved_rowdata, rowData(sce))
    print(dim(sce_tmp))
    if(dim(sce_tmp)[2]>0){
      # First time normalize data using scran, input is counts exp, normalized output are kept in normcounts exp values 
      # rlog=FALSE means the normalized counts are NOT log transformed
      sce_normalized <- sce_normalize_size_factors(sce_tmp, min_size=100, rlog=FALSE, exprs="counts", name="normcounts")  ##  rlog=FALSE: count --> normcounts
      if(!is.null(sce_normalized)){  
        sce_list <- c(sce_list, sce_normalized)
      } else{
        print(paste0("Please double check sample id: ",s))
      }
    }
  }
  # Combining all sces to 1 total sce
  # Take out the tenx column
  sce_normalized$tenx <- NULL
  # sce_combine <- sce_cbind_func(sce_list, rowData(sce), cut_off_overall = 0, exprs = c("counts", "normcounts"),
  #                               colData_names = colnames(colData(sce_normalized)), meta_data=NULL)
  sce_combine <- sce_cbind_func_v2(sce_list, cut_off_overall = 0.025, exprs = c("counts", "normcounts"), 
                                colData_names = colnames(colData(sce_normalized)), save_raw=F, save_dir='',tag='SA') 
  # colData(sce_combine)[1:3,1:3]
  # Twice scater normalization
  # normalized output are kept in logcounts exp values
  print(("Normalizing total data..."))
  # counts(sce_combine) <- NULL
  assayNames(sce_combine)
  sce_normalized_total <- sce_normalize_size_factors(sce_combine, min_size=500, rlog=TRUE, exprs="normcounts", name="logcounts") ##  rlog=TRUE: normcounts --> logcounts
  
  
  saveRDS(sce_normalized_total, file = paste0(output_dir,datatag,"_twice_scran_normalized_v3_logcounts.rds"))
  # sce_normalized_total <- sce_normalize_size_factors(sce_combine, min_size=300, rlog=F, exprs="normcounts", name="normcounts") ##  rlog=TRUE: normcounts --> logcounts
  # saveRDS(sce_normalized_total, file = paste0(output_dir,datatag,"_twice_scran_normalized_v3.rds"))
  if(return_data){
    return(sce_normalized_total)
  }  
}  

sce <- readRDS(paste0(output_dir,datatag,"_twice_scran_normalized_v3_logcounts.rds"))
dim(sce)
assayNames(sce)
get_umap <- function(sce){
  library(SingleCellExperiment)
  library(Seurat)
  library(ggplot2)
  library(sctransform)
  meta_genes <- SingleCellExperiment::rowData(sce)
  mito_genes <- meta_genes$ID[grepl("^MT-",meta_genes$Symbol)]
  genes_use <- rownames(sce)[!rownames(sce) %in% mito_genes]
  sce <- sce[genes_use,]
  dim(sce)
  pbmc <- Seurat::CreateSeuratObject(counts = as.matrix(counts(sce)), meta.data=as.data.frame(colData(sce)))
  
  pbmc <- SetAssayData(object = pbmc, slot = "data", new.data = as.matrix(logcounts(sce)))
  pbmc <- SetAssayData(object = pbmc, slot = "scale.data", new.data = as.matrix(logcounts(sce)))
  pbmc <- Seurat::FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 3000)
  pbmc <- RunPCA(pbmc, verbose = FALSE)
  pbmc <- RunUMAP(pbmc, dims = 1:25, verbose = FALSE)
  
  pbmc <- FindNeighbors(pbmc, dims = 1:25, verbose = FALSE)
  res = 0.3
  pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = res)
  p1 <- DimPlot(pbmc, label = TRUE, group.by = 'seurat_clusters') #+ NoLegend()
  p2 <- DimPlot(pbmc, label = TRUE, group.by = 'library_id') #+ NoLegend()
  pumap <- p1 + p2 + patchwork::plot_layout(ncol=2)
  png(paste0(output_dir,datatag,"_clusters_umap.png"), height = 2*400, width=2*1000,res = 2*72)
  print(pumap)
  dev.off()
  
  meta_cells <- as.data.frame(pbmc@meta.data)
  dim(meta_cells)
  meta_cells
  metasamples <- data.table::fread(paste0(dirname(input_dir),'/snakemake_10x_v2/SA535_10x_metadata.csv'))
  unique(metasamples$Grouping)
  unique(metasamples$library_id)
  colnames(metasamples)
  dim(meta_cells)
  meta_cells$cell_id <- colnames(pbmc)
  library(dplyr)
  meta_cells <- meta_cells %>%
    as.data.frame() %>%
    dplyr::left_join(metasamples, by='library_id')
  # sum(unique(meta_cells$library_id) %in% metasamples$library_id)
  sid_df <- meta_cells %>%
    select(cell_id, Grouping) %>%
    rename(MainSite=Grouping)
  rownames(sid_df) <- sid_df$cell_id
  sid_df$cell_id <- NULL
  pbmc <- AddMetaData(object = pbmc, metadata = sid_df, col.name = 'MainSite')
  
  lid_df <- data.frame(lid=gsub('SCRNA10X_SA_CHIP0','C',meta_cells$library_id), row.names = meta_cells$cell_id)
  pbmc <- AddMetaData(object = pbmc, metadata = lid_df, col.name = 'library_id')
  p3 <- DimPlot(pbmc, reduction = "umap", group.by = 'library_id')
  p4 <- DimPlot(pbmc, reduction = "umap", group.by = 'MainSite')
  pumap <- p1 + p3 + p4 + patchwork::plot_layout(ncol=3)
  png(paste0(output_dir,datatag,"_clusters_umap_twice_scran.png"), height = 2*350, width=2*1200,res = 2*72)
  print(pumap)
  dev.off()
  
  umap_SA535 <- Embeddings(object = pbmc, reduction = "umap") #_200123
  data.table::fwrite(as.matrix(umap_SA535), paste0(dirname(input_dir),'/normalized/', 'umap_SA535_200123.csv'))
  
}
input_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/filtered/'
output_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/normalized/'
suffixes <- '_f_dt_m.rds' ## filtered, doublet removed, mouse cells removed
tag <- 'SA535_200123'
datatag <- tag

sce <- readRDS(paste0(input_dir,"combined_",tag,".rds"))
dim(sce)


metagenes <- rowData(sce)
spike_genes <- metagenes$ID[grepl('ERCC',metagenes$Symbol)]



# Mirela
# twice_scran_normalize(gsub(".rds","", list.files("SA609_v6/sce_annotated/")), input_dir="SA609_v6/sce_annotated/", output_dir="SA609_v6/sce_twice_scran/")
# twice_scran_normalize(gsub(".rds","", list.files("SA1035_v6/sce_annotated/")), input_dir="SA1035_v6/sce_annotated/", output_dir="SA1035_v6/sce_twice_scran/")
# twice_scran_normalize(gsub(".rds","", list.files("SA535_v6/sce_annotated/")), input_dir="SA535_v6/sce_annotated/", output_dir="SA535_v6/sce_twice_scran/")
t1 <- as.matrix(counts(sce_combine))
sum(t1[,10])
assayNames(sce_combine) 

# meta_cells$library_id=='SCRNA10X_SA_CHIP0250_001'
sce_normalized_total <- sce_normalized_total[spike_genes,]
dim(sce_normalized_total)
stat <- DelayedArray::colSums(assay(sce_normalized_total, "logcounts"))
stat2 <- DelayedArray::colSums(assay(sce_normalized_total, "counts"))
length(stat)
head(stat)
names(stat)[1:3]
meta_cells <- as.data.frame(colData(sce_normalized_total))
meta_cells$cell_id <- rownames(meta_cells)
s <- data.frame(total_counts_exp=stat, cell_id=names(stat),total_logcounts_exp=stat2)

meta_cells <- meta_cells %>%
  left_join(s, by='cell_id')

meta_cells$library_id <- gsub('SCRNA10X_SA_CHIP0','C',meta_cells$library_id)
pc <- ggplot(data=meta_cells, aes(x=library_id, y=total_counts_exp, color=library_id)) + 
  geom_jitter(width = 0.1)
pc


pl <- ggplot(data=meta_cells, aes(x=library_id, y=total_logcounts_exp, color=library_id)) + 
  geom_jitter(width = 0.1)
pl

sce <- readRDS(paste0(input_dir,"SCRNA10X_SA_CHIP0249_001_f_dt_m.rds"))
dim(sce)
sce$library_id <- 'C249_001'
colnames(sce)[1:10]
obs_cell <- '_AACGGGACAAGCGGAT-1'
t <- as.matrix(counts(sce[, obs_cell]))
dim(t)
sum(t[,1])
t <- as.matrix(counts(sce_tmp))
sum(t[,1])

sce_combine <- readRDS(paste0(input_dir,"combined_",tag,".rds"))
colnames(sce_combine)[1:2]
t1 <- as.matrix(counts(sce_combine[, paste0('SCRNA10X_SA_CHIP0249_001',obs_cell)]))
dim(t1)
sum(t1[,1])

sce_normalized_total_2
t1 <- as.matrix(logcounts(sce_normalized_total_2[, paste0('SCRNA10X_SA_CHIP0249_001',obs_cell)]))
dim(t1)
sum(t1[,1])
