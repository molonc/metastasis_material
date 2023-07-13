suppressPackageStartupMessages({
  require("SingleCellExperiment")
  # require("stringr")
  require("tidyverse")
  require("Seurat")
  # require("sctransform")
  require("dplyr")
  # require("inlmisc")
})

main(){
  
  meta <- data.table::fread('/home/htran/Projects/hakwoo_project/metastasis_material/materials/10x/SA535_10x_metadata_passage_X4_full_infos.csv')
  
  # meta$nb_cells_introns
  # meta <- meta %>% 
  #   dplyr::filter(mouse_id=='SA535X4XB05462')
  meta <- meta %>% 
    dplyr::filter(pdxid!='M2364')
  colnames(meta)
  dim(meta)
  input_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/filtered_introns/'
  output_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/normalized_introns/'
  datatag <- 'SA535' 
  sce_combine <- readRDS(paste0(output_dir,"combined_",tag,".rds"))
  sce <- sce_combine
}

normalize_scTransform_v3 <- function(sce, datatag='SA', max_dim=25){
  library(Seurat)
  library(ggplot2)
  library(sctransform)
  filtered_cells <- data.table::fread(paste0(output_dir, datatag, "_filtered_cells.csv.gz"))
  dim(filtered_cells)
  sce <- sce[,filtered_cells$cell_id]
  dim(sce)
  meta_cells <- as.data.frame(colData(sce))
  rownames(meta_cells)[1]
  colnames(sce)[1]
  # counts: either a matrix-like object with unnormalized data with cells as columns and features as rows or an Assay-derived object
  srt <- Seurat::CreateSeuratObject(counts = counts(sce),
                                     meta.data = meta_cells)
  
  # store mitochondrial percentage in object meta data
  # srt <- PercentageFeatureSet(srt, pattern = "^MT-", col.name = "percent.mt") ## do not work because we don't use gene symbol as rownames
  ## Replace by the function above
  col.name = "percent.mt"
  total_counts <- srt[[paste0("nCount_", DefaultAssay(object = srt))]]
  
  meta_genes <- rowData(sce)
  mito_genes <- meta_genes$ID[grepl("^MT-",meta_genes$Symbol)]
  percent.featureset <- colSums(x = GetAssayData(object = srt, 
                                                 assay = DefaultAssay(object = srt), slot = "counts")[mito_genes, , drop = FALSE])/total_counts * 100
  srt <- AddMetaData(object = srt, metadata = percent.featureset, 
                      col.name = col.name)
  
  
  # run sctransform
  srt <- SCTransform(srt, vars.to.regress = "percent.mt", verbose = FALSE)
  
  dim(srt)
  ## Faster using GamPoi package  
  # BiocManager::install("glmGamPoi")
  # srt <- SCTransform(srt, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
  
  srt <- RunPCA(srt, verbose = FALSE, npcs = max_dim)
  srt <- RunUMAP(srt, dims = 1:max_dim, verbose = FALSE)
  
  srt <- FindNeighbors(srt, dims = 1:max_dim, verbose = FALSE)
  res = 0.5
  srt <- FindClusters(srt, verbose = FALSE, resolution = res)
  
  meta_info <- srt@meta.data
  print(summary(as.factor(meta_info$seurat_clusters)))
  # pca_df <- as.data.frame(Seurat::Embeddings(object = srt, reduction = "pca"))
  umap_df <- as.data.frame(Seurat::Embeddings(object = srt, reduction = "umap"))
  # library(dplyr)
  dim(umap_df)
  
  meta_info$cell_id <- rownames(meta_info)
  # meta_info <- meta_info %>%
  #   dplyr::select(cell_id, Barcode, library_id, sample, treatmentSt, timepoint)
  
  umap_df$cell_id <- rownames(umap_df)
  umap_df <- umap_df %>% inner_join(meta_info, by=c("cell_id"))
  # umap_df$cell_id <- paste0(umap_df$library_id,'_',umap_df$Barcode)
  
  # pca_df$cell_id <- rownames(pca_df)
  # pca_df <- pca_df %>% inner_join(meta_info, by=c("cell_id"))
  # pca_df$cell_id <- paste0(pca_df$library_id,'_',pca_df$Barcode)
  print(dim(umap_df))
  colnames(umap_df)
  unique(umap_df$Sample)
  umap_df <- umap_df %>%  
    dplyr::select(cell_id, seurat_clusters,
                  Barcode, Sample, library_id,
                  UMAP_1, UMAP_2)
  
  # data.table::fwrite(pca_df, paste0(output_dir, datatag, "_norm_pca_",nb_hvg,".csv"))
  data.table::fwrite(umap_df, paste0(output_dir, datatag, "_scTransform_umap_filtered_cells.csv.gz"))

  ## Filtering small clusters of low quality data
  # filtered_umap_df <- umap_df %>% 
  #   dplyr::filter(!seurat_clusters %in% c(10, 11, 12)) %>% 
  #   dplyr::select(cell_id, seurat_clusters)
  # data.table::fwrite(filtered_umap_df, paste0(output_dir, datatag, "_filtered_cells.csv.gz"))
  # dim(filtered_umap_df)
  # 
  saveRDS(srt, paste0(output_dir, datatag, "_scTransform_srt.rds"))
  # p1 <- DimPlot(srt, label = TRUE, group.by = 'seurat_clusters') #+ NoLegend()
  # 
  # p2 <- DimPlot(srt, label = F, group.by = 'Sample') #+ NoLegend()
  # 
  # 
  # rownames(meta_cells)[1:3]
  # sum(colnames(srt)==meta_cells$cell_id)
  # dim(srt)
  # is_lowquality <- ifelse(meta_cells$library_id=='SCRNA10X_SA_CHIP0250_001',
  #                         'CHIP0250_001','Others')
  # summary(as.factor(is_lowquality))
  # df <- data.frame(lowQuality_lib=is_lowquality)
  # rownames(df) <- meta_cells$cell_id
  # srt <- AddMetaData(object = srt, metadata = df, 
  #                     col.name = 'lowQuality_lib')
  # p2 <- DimPlot(srt, reduction = "umap", group.by = 'lowQuality_lib')
  # 
  # p2  
  # output_dir <- input_dir
  # datatag <- tag
  # pumap <- p1 + p2 + patchwork::plot_layout(ncol=2)
  # png(paste0(output_dir,datatag,"_clusters_umap.png"), height = 2*400, width=2*1000,res = 2*72)
  # print(pumap)
  # dev.off()
  # 
  # metasamples <- data.table::fread(paste0(dirname(input_dir),'/snakemake_10x_v2/SA535_10x_metadata.csv'))
  # unique(metasamples$Grouping)
  # unique(metasamples$library_id)
  # colnames(metasamples)
  # dim(meta_cells)
  # library(dplyr)
  # meta_cells <- meta_cells %>%
  #   as.data.frame() %>%
  #   dplyr::left_join(metasamples, by='library_id')
  # # sum(unique(meta_cells$library_id) %in% metasamples$library_id)
  # sid_df <- meta_cells %>%
  #   select(cell_id, Grouping) %>%
  #   rename(MainSite=Grouping)
  # rownames(sid_df) <- sid_df$cell_id
  # sid_df$cell_id <- NULL
  # srt <- AddMetaData(object = srt, metadata = sid_df, col.name = 'MainSite')
  # 
  # lid_df <- data.frame(lid=gsub('SCRNA10X_SA_CHIP0','C',meta_cells$library_id), row.names = meta_cells$cell_id)
  # srt <- AddMetaData(object = srt, metadata = lid_df, col.name = 'library_id')
  # p3 <- DimPlot(srt, reduction = "umap", group.by = 'library_id')
  # p4 <- DimPlot(srt, reduction = "umap", group.by = 'MainSite')
  # pumap <- p1 + p3 + p4 + patchwork::plot_layout(ncol=3)
  # png(paste0(output_dir,datatag,"_clusters_umap.png"), height = 2*350, width=2*1200,res = 2*72)
  # print(pumap)
  # dev.off()
  # 
  # umap_SA535 <- Embeddings(object = srt, reduction = "umap") #_200123
  # data.table::fwrite(as.matrix(umap_SA535), paste0(dirname(input_dir),'/normalized/', 'umap_SA535_200123.csv'))
  # 
  # saveRDS(srt, paste0(dirname(input_dir),'/normalized/', 'SA535_sctransform_normalized_srt_200123.rds'))
}


get_umap <- function(sce, output_dir, nb_hvg=3000, datatag='SA', max_dim=30){
  # if(file.exists(sce_fn)){
  #   sce <- readRDS(sce_fn)  
  #   reducedDims(sce) <- NULL
  # }else{
  #   stop('Check input sce file, exit!!!')
  # }
  print(nb_hvg)
  print(datatag)
  print(dim(sce))
  
  dims <-  1:max_dim
  # dims <-  1:30, nb_hvg=3000 ## for trajectory analysis, take 30 PCs
  
  # logcounts(sce) <- as.matrix(logcounts(sce))
  # counts(sce) <- as.matrix(counts(sce))
  # library(Seurat)
  # sce$library_id <- gsub('SCRNA10X_SA_CHIP','',sce$library_id)
  srt <- Seurat::as.Seurat(sce, counts = "counts", data="logcounts") #set to NULL if only normalized data are present
  print(dim(srt)) #, , assay = "RNA", project = "SingleCellExperiment"
  srt[["originalexp"]]@scale.data <- as.matrix(logcounts(sce))
  # srt[["SCT"]]@scale.data <- as.matrix(logcounts(sce))
  srt <- Seurat::FindVariableFeatures(srt, selection.method = "vst", nfeatures = nb_hvg)
  # srt <- ScaleData(srt, features = all.genes) # scale for all genes, can scale by hvg genes
  srt <- Seurat::RunPCA(object = srt, verbose = FALSE, npcs = max_dim)
  # srt <- JackStraw(srt, num.replicate = 100)
  # srt <- ScoreJackStraw(srt, dims = 1:20)
  # JackStrawPlot(srt, dims = 1:15)
  # ElbowPlot(srt)
  srt <- Seurat::FindNeighbors(srt, dims = dims)
  # srt <- FindClusters(srt, resolution = res)
  srt <- Seurat::RunUMAP(object = srt, dims = dims, verbose=FALSE, n.components=3) #umap.method = "umap-learn", metric='correlation'
  
  pca_df <- as.data.frame(Seurat::Embeddings(object = srt, reduction = "pca"))
  umap_df <- as.data.frame(Seurat::Embeddings(object = srt, reduction = "umap"))
  # library(dplyr)
  meta_info <- srt@meta.data
  meta_info$cell_id <- rownames(meta_info)
  # meta_info <- meta_info %>%
  #   dplyr::select(cell_id, Barcode, library_id, sample, treatmentSt, timepoint)
  
  umap_df$cell_id <- rownames(umap_df)
  umap_df <- umap_df %>% inner_join(meta_info, by=c("cell_id"))
  # umap_df$cell_id <- paste0(umap_df$library_id,'_',umap_df$Barcode)
  
  pca_df$cell_id <- rownames(pca_df)
  pca_df <- pca_df %>% inner_join(meta_info, by=c("cell_id"))
  # pca_df$cell_id <- paste0(pca_df$library_id,'_',pca_df$Barcode)
  print(dim(umap_df))
  print(dim(pca_df))
  
  # data.table::fwrite(pca_df, paste0(output_dir, datatag, "_norm_pca_",nb_hvg,".csv"))
  data.table::fwrite(umap_df, paste0(output_dir, datatag, "_norm_umap_",nb_hvg,"_3UMAPS.csv"))
}
