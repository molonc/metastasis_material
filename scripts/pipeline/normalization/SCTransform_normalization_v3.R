
library(SingleCellExperiment)
library(Seurat)
# library(DelayedArray) # need to install this package
sce_cbind_func_v3 <- function(sce_list, min_cells = 50, 
                              exprs = c("counts", "logcounts","normcounts"), 
                              colData_names = NULL, save_raw=T, save_dir='',tag='SA') { #, meta_data=NULL
  n_batch <- length(sce_list)
  # method = "intersect"
  
  assay_list <- list()
  for (i in seq_len(length(exprs))) {
    assay_list[[i]] <- do.call(cbind, lapply(sce_list, 
                                             function(y) assay(y, exprs[i])))
  }
  names(assay_list) <- exprs
  if(is.null(colData_names)){ 
    # print('Combining all meta cells features in colData(sce)')
    colData_names <- colnames(colData(sce_list[[1]]))
    for(i in rep(2:length(sce_list),1)){
      colData_names <- intersect(colData_names, colnames(colData(sce_list[[i]])))
    }
    # colData_list <- do.call(DelayedArray::rbind, 
    #                         lapply(sce_list, function(y) colData(y)))
  }
  colData_list <- do.call(DelayedArray::rbind, 
                          lapply(sce_list, function(y) colData(y)[, colData_names, drop = FALSE]))
  # genes_df <- as.data.frame(rowData(sce_list[[1]]))
  # genes_df <- genes_df %>%
  #   dplyr::rename(geo_strand=strand,
  #                 geo_chr=chr,
  #                 geo_start=start,
  #                 geo_end=end)%>%
  #   as.matrix()
  sce_combine <- SingleCellExperiment::SingleCellExperiment(assay = assay_list, 
                                                            colData = colData_list,
                                                            rowData=rowData(sce_list[[1]]))
  
  
  
  # sce_combine$mouse_id <- meta_data$mouse_id
  # sce_combine$passage <- meta_data$passage
  # sce_combine$treatmentSt <- meta_data$treatmentSt
  
  print(paste0("Dim sce combine: ",dim(sce_combine)))
  # print(colnames(colData(sce_combine)))
  print(paste0("sce combine assay name: ",assayNames(sce_combine)))
  if(save_raw){
    saveRDS(sce_combine, file = paste0(save_dir,"filtered_combined_",tag,".rds"))
  }
  if(min_cells > 0){
    nonzero_cbind <- DelayedArray::rowSums(assay(sce_combine, exprs[1]) > 0)
    sce_combine <- sce_combine[names(nonzero_cbind[nonzero_cbind >= min_cells]), ]
  }
  
  return(sce_combine)
}

sce_cbind_func_v2 <- function(sce_list, cut_off_overall = 0.01, exprs = c("counts", "logcounts","normcounts"), 
                              colData_names = NULL, save_raw=T, save_dir='',tag='SA') { #, meta_data=NULL
  # n_batch <- length(sce_list)
  # method = "intersect"
  
  assay_list <- list()
  for (i in seq_len(length(exprs))) {
    assay_list[[i]] <- do.call(cbind, lapply(sce_list, 
                                             function(y) assay(y, exprs[i])))
  }
  names(assay_list) <- exprs
  if(is.null(colData_names)){ 
    # print('Combining all meta cells features in colData(sce)')
    colData_names <- colnames(colData(sce_list[[1]]))
    for(i in rep(2:length(sce_list),1)){
      colData_names <- intersect(colData_names, colnames(colData(sce_list[[i]])))
    }
    # colData_list <- do.call(DelayedArray::rbind, 
    #                         lapply(sce_list, function(y) colData(y)))
  }
  colData_list <- do.call(DelayedArray::rbind, 
                          lapply(sce_list, function(y) colData(y)[, colData_names, drop = FALSE]))
  # genes_df <- as.data.frame(rowData(sce_list[[1]]))
  # genes_df <- genes_df %>%
  #   dplyr::rename(geo_strand=strand,
  #                 geo_chr=chr,
  #                 geo_start=start,
  #                 geo_end=end)%>%
  #   as.matrix()
  sce_combine <- SingleCellExperiment::SingleCellExperiment(assays = assay_list, 
                                                            colData = colData_list,
                                                            rowData=rowData(sce_list[[1]]))
  
  
  # sce_combine$mouse_id <- meta_data$mouse_id
  # sce_combine$passage <- meta_data$passage
  # sce_combine$treatmentSt <- meta_data$treatmentSt
  
  print(paste0("Dim sce combine: ",dim(sce_combine)))
  # print(colnames(colData(sce_combine)))
  print(paste0("sce combine assay name: ",assayNames(sce_combine)))
  if(save_raw){
    saveRDS(sce_combine, file = paste0(save_dir,"raw_combined_",tag,".rds"))
  }
  if(cut_off_overall > 0){
    zero_cbind <- DelayedArray::rowMeans(assay(sce_combine, exprs[1]) == 0)
    sce_combine <- sce_combine[names(zero_cbind[zero_cbind <= (1 - cut_off_overall)]), ]
  }
  
  return(sce_combine)
}

fns <- c('SCRNA10X_SA_CHIP0249_001','SCRNA10X_SA_CHIP0249_002','SCRNA10X_SA_CHIP0250_001')
download_dir <- '/home/htran/storage/rnaseq_datasets/hakwoo_metastasis_RNAseq/SA535_human/'
for(f in fns){
  sce <- readRDS(paste0(download_dir,f,'/',f, '.rdata'))
  print(dim(sce)[2])
  
}

input_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/filtered/'
suffixes <- '_f_dt_m.rds' ## filtered, doublet removed, mouse cells removed
tag <- 'SA535_200123'
load_data <- function(input_dir, suffixes, tag=''){
  
  fns <- list.files(input_dir)
  fns <- fns[grepl(suffixes, fns)]
  print(length(fns))
  stat <- tibble::tibble()
  
  sce_list <- list()
  c <- 0
  for (f in unique(fns)){
    fn <- paste0(input_dir,f)
    
    if(file.exists(fn)){
      sce_tmp <- readRDS(fn)
      c <- c + 1
      if(c==1){
        cols_use <- colnames(colData(sce_tmp))
      }else{
        cols_use <- intersect(cols_use, colnames(colData(sce_tmp)))
      }
      print(paste0('DEBUG: ',f,'  ',dim(sce_tmp)[1],' ',dim(sce_tmp)[2]))
      stat <- dplyr::bind_rows(stat, tibble::tibble(library=gsub(suffixes,'',f),nb_cells=dim(sce_tmp)[2]))
      sce_list[[f]] <- sce_tmp
    } else{
      warning('Error loading files')
      print(f)
    }  
    
    
  }  
  print(length(cols_use))
  data.table::fwrite(stat, paste0(input_dir, 'summary_nb_cells_200123.csv'))
  ## cut_off_overall = 0: do not filter any gene
  ## cut_off_overall = 0.025: filter genes with zero counts values in 97.5% of total cells
  sce_combine <- sce_cbind_func_v2(sce_list, cut_off_overall = 0.025, exprs = c("counts"), 
                                   colData_names = cols_use, save_raw=F, save_dir=input_dir, tag) 
  print(dim(sce_combine))
  # sce_combine <- readRDS(paste0(output_dir,'combined_total_genes_filtered.rds'))
  # Twice scater normalization
  # output are saved in logcounts exp values
  print(("Normalizing total data..."))
  
  ## cell name = library_id + barcode
  if('Sample' %in% colnames(colData(sce_combine))){
    sce_combine$library_id <- gsub('.cache/','',sce_combine$Sample)
    sce_combine$library_id <- gsub('/filtered_feature_bc_matrix','',sce_combine$library_id)
    print(unique(sce_combine$library_id))
    colnames(sce_combine) <- paste0(sce_combine$library_id,'_',sce_combine$Barcode)
  }
  # rownames(sce_combine)
  # print(rowData(sce_combine_raw)$Symbol[1:3])
  
  saveRDS(sce_combine, file=paste0(input_dir,"combined_",tag,".rds"))
  dim(sce_combine) 
  meta_cells <- colData(sce_combine)
  meta_cells$cell_id <- rownames(meta_cells)
  data.table::fwrite(as.matrix(meta_cells), paste0(input_dir, 'summary_meta_cells_200123.csv'))
  
  return(sce_combine)
  
}


normalize_scTransform_v3 <- function(sce, max_dim=25){
  library(Seurat)
  library(ggplot2)
  library(sctransform)
  rm(mtx)
  # counts: either a matrix-like object with unnormalized data with cells as columns and features as rows or an Assay-derived object
  pbmc <- Seurat::CreateSeuratObject(counts = mtx)
  
  # store mitochondrial percentage in object meta data
  # pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt") ## do not work because we don't use gene symbol as rownames
  ## Replace by the function above
  col.name = "percent.mt"
  total_counts <- pbmc[[paste0("nCount_", DefaultAssay(object = pbmc))]]
    
  meta_genes <- rowData(sce_combine)
  mito_genes <- meta_genes$ID[grepl("^MT-",meta_genes$Symbol)]
  percent.featureset <- colSums(x = GetAssayData(object = pbmc, 
                                                 assay = DefaultAssay(object = pbmc), slot = "counts")[mito_genes, , drop = FALSE])/total_counts * 100
  pbmc <- AddMetaData(object = pbmc, metadata = percent.featureset, 
                        col.name = col.name)
  
  
  # run sctransform
  pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
  
  
  ## Faster using GamPoi package  
  # BiocManager::install("glmGamPoi")
  # pbmc <- SCTransform(pbmc, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
  
  pbmc <- RunPCA(pbmc, verbose = FALSE)
  pbmc <- RunUMAP(pbmc, dims = 1:max_dim, verbose = FALSE)
  
  pbmc <- FindNeighbors(pbmc, dims = 1:max_dim, verbose = FALSE)
  res = 0.3
  pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = res)
  p1 <- DimPlot(pbmc, label = TRUE, group.by = 'seurat_clusters') #+ NoLegend()
  rownames(meta_cells)[1:3]
  sum(colnames(pbmc)==meta_cells$cell_id)
  dim(pbmc)
  is_lowquality <- ifelse(meta_cells$library_id=='SCRNA10X_SA_CHIP0250_001',
                          'CHIP0250_001','Others')
  summary(as.factor(is_lowquality))
  df <- data.frame(lowQuality_lib=is_lowquality)
  rownames(df) <- meta_cells$cell_id
  pbmc <- AddMetaData(object = pbmc, metadata = df, 
                      col.name = 'lowQuality_lib')
  p2 <- DimPlot(pbmc, reduction = "umap", group.by = 'lowQuality_lib')
  
  p2  
  output_dir <- input_dir
  datatag <- tag
  pumap <- p1 + p2 + patchwork::plot_layout(ncol=2)
  png(paste0(output_dir,datatag,"_clusters_umap.png"), height = 2*400, width=2*1000,res = 2*72)
  print(pumap)
  dev.off()
  
  metasamples <- data.table::fread(paste0(dirname(input_dir),'/snakemake_10x_v2/SA535_10x_metadata.csv'))
  unique(metasamples$Grouping)
  unique(metasamples$library_id)
  colnames(metasamples)
  dim(meta_cells)
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
  png(paste0(output_dir,datatag,"_clusters_umap.png"), height = 2*350, width=2*1200,res = 2*72)
  print(pumap)
  dev.off()
  
  umap_SA535 <- Embeddings(object = pbmc, reduction = "umap") #_200123
  data.table::fwrite(as.matrix(umap_SA535), paste0(dirname(input_dir),'/normalized/', 'umap_SA535_200123.csv'))
  
  saveRDS(pbmc, paste0(dirname(input_dir),'/normalized/', 'SA535_sctransform_normalized_srt_200123.rds'))
}
