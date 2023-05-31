# BiocManager::install('extrafont')
suppressPackageStartupMessages({
  # require("optparse")
  require("scater")
  # require("argparse")
  require("SingleCellExperiment")
  require("stringr")
  require("tidyverse")
  # require("scran")
  require("Seurat")
  # require("sctransform")
  require("dplyr")
  
})

# htran



format_theme_plot <- function(p, title='Seurat Clustering'){
  p <- p + labs(title)
  p <- p + theme(
    # legend.title = element_text(size=10), 
    legend.text = element_text(size=10), 
    plot.title = element_text(color="black", size=17, hjust = 0.5)
  )
  return(p)
}



# This function calls sctransform::vst. The sctransform package is available at https://github.com/ChristophH/sctransform. 
# Use this function as an alternative to the NormalizeData, FindVariableFeatures, ScaleData workflow. Results are saved 
# in a new assay (named SCT by default) with counts being (corrected) counts - ATTENTION, data being log1p(counts), 
# scale.data being pearson residuals; sctransform::vst intermediate results are saved in misc slot of new assay. 
# Here using sctransform function in Seurat3 package
normalize_SCTransform <- function(sce, output_dir, datatag, return_data=F, output_fn=NULL){
  if(is.null(output_fn)){
    output_fn <- paste0(output_dir, datatag,'_sctransform_normalized.rds')
  }
    # if(file.exists(output_fn)){
  #   print("Output exist in folder, reloading...")
  #   sce2 <- readRDS(output_fn)
  #   # return(sce2)
  # }
  # else{
  #   
  # }
  print("Removing mito, ribo genes prior to normalization...")
  print(rowData(sce)$Symbol[1])
  mito_genes <- str_detect(rowData(sce)$Symbol, "^MT\\-")
  print(sum(mito_genes==TRUE))
  ribo_genes <- str_detect(rowData(sce)$Symbol, "^RP(L|S)")  # or ^RP[L|S]?
  print(sum(ribo_genes==TRUE))
  sce <- sce[(!mito_genes) & (!ribo_genes),]
  print(dim(sce))
  
  
  meta_genes <- rowData(sce)
  rownames(meta_genes) <- meta_genes$ID
  print(meta_genes$ID[1])
  sce_counts <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=counts(sce)), 
                                                           colData = colData(sce), 
                                                           rowData = rowData(sce))
  print(dim(sce_counts))
  srt <- Seurat::as.Seurat(sce_counts, counts = "counts", data=NULL, assay = "RNA", project = "SingleCellExperiment") #set to NULL if only normalized data are present
  
  ncells_use <- 10000 # nb cells use to find residual
  if(dim(sce)[2]<ncells_use){
    ncells_use <- dim(sce)[2]
  }
  # ?Seurat::SCTransform
  srt <- Seurat::SCTransform(object = srt, verbose = T, 
                             do.scale=F, 
                             return.only.var.genes = T, 
                             assay = "RNA", new.assay.name = "SCT", ncells = ncells_use)
  # print(summary(srt[["RNA"]]@data))
  print(dim(srt[["SCT"]]@data))
  print(dim(srt[["SCT"]]@counts))
  # norm_mtx <- srt[["SCT"]]@data
  # # sc <- srt[["RNA"]]@scale.data
  # max(norm_mtx)
  # min(norm_mtx)
  # print('Scale sctransform normalized data...')
  # srt <- ScaleData(srt, features = rownames(srt), assay = 'SCT')
  sce2 <- Seurat::as.SingleCellExperiment(srt, assay = 'SCT')
  rownames(sce2) <- rownames(sce)
  # rowData(sce2) <- rowData(sce)
  print(dim(sce2))
  rowData(sce2) <- dplyr::bind_cols(as.data.frame(rowData(sce)),
                                         as.data.frame(rowData(sce2)))
  
  # normcounts(sce2) <- srt[["SCT"]]@scale.data  # save scale data to normcounts assay
  print(assayNames(sce2))
  # print(dim(srt[["SCT"]]@scale.data))
  # print(max(normcounts(sce2)))
  # # min(normcounts(sce2))
  # print('Save normalized data...')
  saveRDS(sce2, output_fn)
  
  # res = 0.3
  # dims = 1:15
  # srt <- RunPCA(object = srt, verbose = FALSE)
  # srt <- FindNeighbors(srt, dims = dims)
  # srt <- FindClusters(srt, resolution = res)
  # srt <- RunUMAP(object = srt, dims = dims, verbose = FALSE) #umap.method = "umap-learn", metric='correlation' 
  # 
  # p21 <- DimPlot(srt, reduction = "umap")
  # p21 <- format_theme_plot(p=p21, title='Seurat Clustering')
  # 
  # p22 <- DimPlot(srt, reduction = "umap", group.by = 'library_id')
  # p22 <- format_theme_plot(p=p22, title='library_id')
  # 
  # # p23 <- DimPlot(srt, reduction = "umap", group.by = 'Grouping')
  # # p23 <- format_theme_plot(p=p23, title='Grouping')
  # # 
  # # p24 <- DimPlot(srt, reduction = "umap", group.by = 'Site_origin')
  # # p24 <- format_theme_plot(p=p24, title='Site_origin')
  # # 
  # # pumap <- p21 + p22 + patchwork::plot_layout(ncol=2)
  # png(paste0(output_dir,datatag,"_clusters_umap.png"), height = 2*350, width=2*450,res = 2*72)
  # print(p21)
  # dev.off()
  # 
  # pumap2 <- p23 + p24 + patchwork::plot_layout(ncol=2)
  # png(paste0(output_dir,datatag,"_clusters_umap_2.png"), height = 2*350, width=2*800,res = 2*72)
  # print(pumap2)
  # dev.off()
  if(return_data){
    return(list(sce=sce, sctransform=sce2))
  }  
}  
viz_umap <- function(sce, output_dir, datatag='SA', res = 0.3, dims = 1:25, nb_hvg=3000,
                     return_seurat_obj=F){
  # library(Seurat)
  sce$library_id <- gsub('SCRNA10X_SA_CHIP','',sce$library_id)
  # sce$Grouping <- ifelse(sce$Site_origin=="Tumor_Recur","Primary",sce$Grouping)
  srt <- Seurat::as.Seurat(sce, counts = "counts", 
                           data="logcounts", assay = "RNA", project = "SingleCellExperiment") #set to NULL if only normalized data are present
  print(dim(srt))
  # rownames(srt)[1]
  # class(srt[["RNA"]]@data)
  # srt[["RNA"]]@scale.data <- normcounts(sce)
  srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = nb_hvg)
  
  
  all.genes <- rownames(srt)
  srt <- ScaleData(srt, features = all.genes) # scale for all genes, can scale by hvg genes
  srt <- RunPCA(object = srt, verbose = FALSE)
  srt <- JackStraw(srt, num.replicate = 100)
  srt <- ScoreJackStraw(srt, dims = 1:20)
  JackStrawPlot(srt, dims = 1:15)
  ElbowPlot(srt)
  srt <- FindNeighbors(srt, dims = dims)
  srt <- FindClusters(srt, resolution = res)
  srt <- RunUMAP(object = srt, dims = dims, verbose = FALSE) #umap.method = "umap-learn", metric='correlation'

  p21 <- DimPlot(srt, reduction = "umap", group.by ='seurat_clusters')
  p21 <- format_theme_plot(p=p21, title='Seurat Clustering')
  meta_df <- srt@meta.data
  # srt@meta.data <- meta_df
  unique(meta_df$passage.x)
  # meta_df$passage <- NULL
  # meta_df$passage.x <- NULL
  # meta_df$passage.y <- NULL
  meta_df$batch_info <- NULL
  # p22 <- DimPlot(srt, reduction = "umap", group.by = 'library_id')
  # p22 <- format_theme_plot(p=p22, title='library_id')
  # srt$treatment <- srt$treat
  # srt$clone <- 'unassigned'
  # srt$passage <- NULL
  
  # cells_use <- intersect(rownames(clonealign_df),colnames(srt))
  # srt@meta.data[cells_use,'clone'] <- clonealign_df[cells_use,'clone']
  # srt@meta.data[cells_use,'passage'] <- clonealign_df[cells_use,'passage']
  
  # srt@meta.data <- srt@meta.data %>% left_join(meta_cells, by=c('sample'='mouse_id'))
  srt@meta.data <- meta_df
  rownames(srt@meta.data) <- colnames(srt)
  p22 <- DimPlot(srt, reduction = "umap", group.by = 'treatment')
  p22 <- format_theme_plot(p=p22, title='Treatment')
  
  p22 <- DimPlot(srt, reduction = "pca", group.by = 'treatment')
  p22 <- format_theme_plot(p=p22, title='Treatment')
  # p22 <- DimPlot(srt, reduction = "umap", group.by = 'clone_desc')
  # p22 <- format_theme_plot(p=p22, title='Subclone')

  # p23 <- DimPlot(srt, reduction = "umap", group.by = 'Grouping')
  # p23 <- format_theme_plot(p=p23, title='Grouping')
  
  p23 <- DimPlot(srt, reduction = "umap", group.by = 'clone')
  p23 <- format_theme_plot(p=p23, title='Clone')

  # p24 <- DimPlot(srt, reduction = "umap", group.by = 'passage')
  # p24 <- format_theme_plot(p=p23, title='passage')
  # p24 <- DimPlot(srt, reduction = "umap", group.by = 'Site_origin')
  # p24 <- format_theme_plot(p=p24, title='Site_origin')

  # pumap <- p21 + p22 + patchwork::plot_layout(ncol=2)
  # png(paste0(output_dir,datatag,"_clusters_umap_1.png"), height = 2*350, width=2*800,res = 2*72)
  # print(pumap)
  # dev.off()
  # 
  # pumap2 <- p23 + p24 + patchwork::plot_layout(ncol=2)
  # png(paste0(output_dir,datatag,"_clusters_umap_2.png"), height = 2*350, width=2*800,res = 2*72)
  # print(pumap2)
  # dev.off()

  # pumap_total <- p21 + p24 + p22 + p23  + patchwork::plot_layout(ncol=2)
  pumap_total <- p22 + p23  + patchwork::plot_layout(ncol=2)
  out_fig <- paste0(output_dir,datatag,"_clusters_umap.png")
  png(out_fig, height = 2*300, width=2*800,res = 2*72)
  print(pumap_total)
  dev.off()
  # print('Save output as: ')
  # print(out_fig)
  if(return_seurat_obj){
    return(srt) #srt@meta.data
  }
}

normalize_Seurat <- function(sce, input_dir, output_dir, datatag, return_data=F){
  sce_counts <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=counts(sce)), 
                                                            colData = colData(sce), 
                                                           rowData = rowData(sce))
  print(dim(sce_counts))
  print('Normalizing data...')
  srt <- Seurat::as.Seurat(sce_counts, counts = "counts", data=NULL) #set to NULL if only normalized data are present
  srt <- NormalizeData(srt, normalization.method = "LogNormalize", scale.factor = 100000)
  srt <- ScaleData(srt, features = rownames(srt))
  # print(summary(srt[["RNA"]]@scale.data))
  # norm_mtx <- srt[["RNA"]]@data
  # sc <- srt[["RNA"]]@scale.data
  
  sce2 <- as.SingleCellExperiment(srt)
  normcounts(sce2) <- srt[["RNA"]]@scale.data  # save scale data to normcounts assay
  print(assayNames(sce2))
  print(dim(sce2))
  print('Save normalized data...')
  saveRDS(sce2, paste0(output_dir,datatag,'_seurat_normalized_v3.rds'))
  if(return_data){
    return(sce2)
  }  
  
  
  # max(normcounts(sce2))
  # min(normcounts(sce2))
  # max(logcounts(sce2))
  # srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = 2000)
  # srt <- RunPCA(srt, features = VariableFeatures(object = srt), verbose = F)
  # srt <- FindNeighbors(srt, dims = 1:15)
  # srt <- FindClusters(srt, resolution = 0.5)
  # srt <- RunUMAP(srt, dims = 1:15)
  
  # p <- DimPlot(srt, reduction = "umap")
  # p1 <- DimPlot(srt, reduction = "umap", group.by = 'series')
  # pumap <- p + p1 + patchwork::plot_layout(ncol=2)
  # png(paste0(output_dir, datatag, "_clusters_umap_Seurat_normalization.png"), height = 2*400, width=2*720,res = 2*72)
  # print(pumap)
  # dev.off()
}

