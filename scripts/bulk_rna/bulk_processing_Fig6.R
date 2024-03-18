# Figure 6
suppressPackageStartupMessages({
  library("dplyr")
  library("ggplot2")
  library("DESeq2")
})

## Loading utility function
script_dir <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/bulkRNAseq/'
source(paste0(input_dir, 'scripts/bulk_rna/bulk_utils.R'))

# Meta data for main + mixing, then extracting clone B + C 
load_metadata_SA919 <- function(){
  input_dir <- '/Users/hoatran/Documents/BCCRC_projects/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/'
  input_dir <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/bulkRNAseq/'
  meta_samples <- data.table::fread(paste0(input_dir, 'metadata/library_groupings_bulk_SA919_cloneIds.csv'))
  meta_samples$nb_cells <- NULL
  data.table::fwrite(meta_samples, paste0(input_dir, 'metadata/library_groupings_bulk_SA919_cloneIds.csv'))
  
  # meta_samples <- meta_samples %>%
  #   dplyr::select(mainsite, clone_id, bulk_sid) %>%
  #   dplyr::filter(bulk_sid %in% colnames(norm_df) & clone_id %in% c('B','C'))%>%
  #   dplyr::mutate(bulk_sid1 =stringr::str_sub(bulk_sid, nchar(bulk_sid)-8, nchar(bulk_sid)),
  #                 label=paste0(mainsite, '_',clone_id, '_',bulk_sid1))
  # # colnames(norm_df) %in% meta_samples$bulk_sid
  # meta_samples$label <- gsub('XB0','_',meta_samples$label)
  # # meta_samples$label
  # meta_samples <- as.data.frame(meta_samples)
  # rownames(meta_samples) <- meta_samples$label
  
  meta_samples_main_mixing <- data.table::fread(paste0(input_dir, 'metadata/metadata_Hakwoo_bulkRNA_mixing_main_exp.csv'))
  dim(meta_samples_main_mixing)
  
  meta_samples_main_mixing <- meta_samples_main_mixing %>%
    dplyr::filter(main_clone %in% c('B','C'))
  return(meta_samples_main_mixing)
}
# Load raw data B, C using txImport 
load_data <- function(){
  base_dir <- '/home/htran/storage/rnaseq_datasets/bulk_metastasis/example_bulkRNAseq_analysis/'
  input_dir <- paste0(base_dir, 'input_Kallisto_files/')
  datatag <- 'SA919Fig6'
  save_dir <- paste0(base_dir, 'preprocessed/')
  if(!dir.exists(save_dir)){
    dir.create(save_dir)
  }
  sample_ids=NULL
  df_counts <- load_raw_counts_kallisto(input_dir, datatag, save_dir, sample_ids)
  
  ## Do clustering, with, and without mixing experiment, and see output. 
  ## Loading counts - output of above function, and normalize data
  df_counts_fn <- paste0(save_dir, datatag,'_total_raw_counts.csv.gz')
  
  # Normalize data, get size factor for each sample
  res <- normalize_by_size_factor(df_counts_fn, datatag, save_dir)
  
  df_normalized <- res$norm_mtx
  size_factors <- res$size_factor
  
  ## Then do hierarchical clustering using normalized data
  normalized_fn <- paste0(save_dir, datatag, '_sizefactor_normalized.csv.gz')
  get_clustering(normalized_fn, save_dir, datatag)
  
  cts <- data.table::fread(df_counts_fn) %>% as.data.frame()
  print(dim(cts))
  
  # coldata <- meta_samples
  dds <- create_DESeq2_obj(coldata, cts, use_existing_size_factor=T)
  get_vst_clustering_results(dds)
  
  ## How to select samples here? 
  ddsBB <- dds[,sidBB]
  resBB <- get_DE_genes_DESeq2(ddsBB, DE_comp=c("Metastasis","Primary"),
                     filter_genes=T, min_total_exp_samples=10)
  
}
get_DE_genes_DESeq2 <- function(dds, DE_comp=c("Metastasis","Primary"),
                                filter_genes=T, min_total_exp_samples=10){
  
  print(dim(dds))
  print(sizeFactors(dds))
  print(colData(dds))
  ## ----prefilter----------------------------------------------------------------
  if(filter_genes){
    # smallestGroupSize <- 3 # group untreated with 3 samples is the smallest group
    # keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
    keep <- rowSums(counts(dds) >= min_total_exp_samples)
    dds <- dds[keep,]
    print(dim(dds))  
  }
  
  ## ----factorlvl----------------------------------------------------------------
  dds$condition <- factor(dds$condition, levels = DE_comp)
  
  ## ----relevel------------------------------------------------------------------
  # dds$condition <- relevel(dds$condition, ref = "untreated")
  
  ## ----droplevels---------------------------------------------------------------
  # dds$condition <- droplevels(dds$condition)
  
  ## ----deseq--------------------------------------------------------------------
  ## size factors are noted in colData column
  dds <- DESeq(dds)
  contrast_DE_comp <- c("condition")
  contrast_DE_comp <- c(contrast_DE_comp, DE_comp)
  res <- results(dds, contrast=contrast_DE_comp) #c("condition","treated","untreated")
  res$ensembl_gene_id <- rownames(res)
  return(res)
}
## Size factor is kept in coldata$size_factor
create_DESeq2_obj <- function(coldata, cts, use_existing_size_factor=T){
  # library("DESeq2")
  cts <- cts[, rownames(coldata)]
  print(all(rownames(coldata) == colnames(cts)))
  
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~ condition)
  print(dds)
  
  ## Supplying existing size factor for each sample to deseq2 object
  if(use_existing_size_factor){
    # is.null(sizeFactors(object)
    sizeFactors(dds) <- coldata$size_factor  
  }
  return(dds)
  
}
get_vst_clustering_results <- function(dds){
  library("pheatmap")
  library("RColorBrewer")
  vsd <- vst(dds, blind=FALSE)
  # rld <- rlog(dds, blind=FALSE)
  print(head(assay(vsd), 3))
  # this gives log2(n + 1)
  # ntd <- normTransform(dds)
  # # BiocManager::install("vsn")
  # library("vsn")
  # meanSdPlot(assay(ntd))
  # meanSdPlot(assay(vsd))
  # meanSdPlot(assay(rld))
  
  ## ----heatmap------------------------------------------------------------------
  
  select <- order(rowMeans(counts(dds, normalized=TRUE)),
                  decreasing=TRUE)[1:20]
  df <- as.data.frame(colData(dds)[,c("condition","type")])
  pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
           cluster_cols=FALSE, annotation_col=df)
  pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
           cluster_cols=FALSE, annotation_col=df)
  pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
           cluster_cols=FALSE, annotation_col=df)
  
  ## ----sampleClust--------------------------------------------------------------
  sampleDists <- dist(t(assay(vsd))) ## transform to the format in ML, samples in rows, features in cols
  
  ## ----figHeatmapSamples, fig.height=4, fig.width=6-----------------------------
  
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
  colnames(sampleDistMatrix) <- NULL
  # colnames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  p <- pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
  return(p)
  # sampleDistMatrix: full mtx
  # sampleDists: half mtx
  
}

# Run DESeq2 for B met vs B pri, C met vs B pri, C met mix vs B pri
# Get common cis/ trans genes 
# Comparing results with scRNA-seq at cluster level. 