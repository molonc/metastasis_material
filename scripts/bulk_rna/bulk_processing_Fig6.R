# Figure 6
suppressPackageStartupMessages({
  library("dplyr")
  library("ggplot2")
  library("DESeq2")
})

## Loading utility function
script_dir <- '/home/htran/Projects/hakwoo_project/metastasis_material/'
source(paste0(script_dir, 'scripts/bulk_rna/bulk_utils.R'))

# Meta data for main + mixing, then extracting clone B + C 
load_metadata_SA919 <- function(){
  # input_dir <- '/Users/hoatran/Documents/BCCRC_projects/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/'
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
  
  meta_samples <- data.table::fread(paste0(input_dir, 'metadata/metadata_Hakwoo_bulkRNA_mixing_main_exp.csv'))
  dim(meta_samples)
  
  meta_samples <- meta_samples %>%
    dplyr::filter(main_clone %in% c('B','C'))
  dim(meta_samples)
  # View(meta_samples)
  return(meta_samples)
}
# Load raw data B, C using txImport 
load_data <- function(meta_samples){
  input_dir <- '/home/htran/storage/rnaseq_datasets/bulk_metastasis/bulk_SA919/'
  datatag <- 'SA919Fig6'
  save_dir <- paste0(input_dir, 'preprocessed/')
  if(!dir.exists(save_dir)){
    dir.create(save_dir)
  }
  
  meta_samples_main <- meta_samples %>%
    dplyr::filter(experiment=="main_exp")
  df_counts_main <- load_raw_counts_kallisto(input_dir, datatag, save_dir, meta_samples_main$bulk_sid)
  dim(df_counts_main)
  data.table::fwrite(df_counts_main, paste0(save_dir, 'df_counts_main_exp_',datatag,'.csv.gz'))
  
  meta_samples_mixing <- meta_samples %>%
    dplyr::filter(experiment=="mixing_exp")
  df_counts_mixing <- load_raw_counts_kallisto(input_dir, datatag, save_dir, meta_samples_mixing$bulk_sid)
  dim(df_counts_mixing)
  data.table::fwrite(df_counts_mixing, paste0(save_dir, 'df_counts_mixing_exp_',datatag,'.csv.gz'))
  
  common_genes <- intersect(df_counts_main$ens_gene_id, df_counts_mixing$ens_gene_id)
  # df_counts_mixing <- df_counts_mixing %>%
  #   dplyr::filter(ens_gene_id %>% )
  counts_total <- df_counts_main  %>%
    inner_join(df_counts_mixing, by='ens_gene_id')
  dim(counts_total)
  head(counts_total)
  # ref_main <- annotables::grch38 %>%
  #   dplyr::select(ensgene,symbol,chr) %>%
  #   dplyr::rename(ens_gene_id=ensgene) %>%
  #   dplyr::filter(ens_gene_id %in% df_counts_main$ens_gene_id)
  # ref_main <- ref_main[!duplicated(ref_main$ens_gene_id),]
  # 
  # ref_mixing <- annotables::grch38 %>%
  #   dplyr::select(ensgene,symbol,chr) %>%
  #   dplyr::rename(ens_gene_id=ensgene) %>%
  #   dplyr::filter(ens_gene_id %in% df_counts_mixing$ens_gene_id)
  # ref_mixing <- ref_mixing[!duplicated(ref_mixing$ens_gene_id),]
  # 
  # sum(df_counts_mixing$ens_gene_id %in% df_counts_main$ens_gene_id)
  # sum(ref_main$symbol %in% ref_mixing$symbol)
  ## Do clustering, with, and without mixing experiment, and see output. 
  ## Loading counts - output of above function, and normalize data
  df_counts_fn <- paste0(save_dir, datatag,'_total_raw_counts.csv.gz')
  data.table::fwrite(counts_total, df_counts_fn)
  # colSums(is.na(counts_total))
  # Normalize data, get size factor for each sample
  res <- normalize_by_size_factor_v2(df_counts_fn, datatag, save_dir)
  
  df_normalized <- res$norm_mtx
  size_factors <- res$size_factor
  
  ## Then do hierarchical clustering using normalized data
  normalized_fn <- paste0(save_dir, datatag, '_sizefactor_normalized.csv.gz')
  df_normalized <- data.table::fread(normalized_fn)
  # head(df_normalized)
  sf_df <- data.table::fread(paste0(save_dir, datatag, '_sizefactors_Scran.csv.gz'))
  
  # get_clustering(normalized_fn, save_dir, datatag)
  df_counts_fn <- paste0(save_dir, datatag,'_total_raw_counts.csv.gz')
  cts <- data.table::fread(df_counts_fn) %>% as.data.frame()
  print(dim(cts))
  head(cts)
  colnames(meta_samples)
  
  cts <- cts %>%
    tibble::column_to_rownames('ens_gene_id')
  coldata <- meta_samples %>%
    dplyr::filter(bulk_sid %in% colnames(cts)) %>%
    dplyr::select(mainsite, main_clone, experiment, bulk_sid) %>%
    as.data.frame() %>%
    tibble::column_to_rownames('bulk_sid')
  coldata <- coldata[colnames(cts),]
  sfs <- sf_df$size_factor
  names(sfs) <- sf_df$sample
  coldata$size_factor <- sfs[rownames(coldata)]
  
  
  sids <- colnames(cts)
  for(s in sids){
    cts[,s] <- round(cts[,s],0)
  }
  # head(cts)
  # coldata <- meta_samples
  coldata$condition <- coldata$mainsite
  dds <- create_DESeq2_obj(coldata, cts, use_existing_size_factor=T)
  # dds <- DESeq(dds)
  get_vst_clustering_results(dds)
  
  
}
extract_DE_genes_DESeq2 <- function(dds){
  ## How to select samples here? 
  # backup_dds <- dds
  # dds <- backup_dds
  metaBB <- colData(dds) %>%
    as.data.frame() %>%
    filter(main_clone=='B' & experiment=='main_exp')
  sidBB <- rownames(metaBB)
  ddsBB <- dds[,sidBB]
  print(sizeFactors(ddsBB)) ## using existing size factors calculated using Scran method
  resBB <- get_DE_genes_DESeq2(ddsBB, DE_comp=c("Metastasis","Primary"),
                               filter_genes=F, min_total_exp_samples=10)
  dim(resBB)
  resBB1 <- resBB %>%
    as.data.frame() %>%
    filter(abs(log2FoldChange)>=0.25) #pvalue<0.05 & 
  dim(resBB1)
  colnames(resBB1)
  resBB2 <- resBB1 %>%
    dplyr::select(log2FoldChange, ensembl_gene_id) %>%
    dplyr::rename(log2FoldChange_BB=log2FoldChange)
  
  metaC <- colData(dds) %>%
    as.data.frame() %>%
    filter(main_clone=='C' & experiment=='main_exp')
  metaBpri <- colData(dds) %>%
    as.data.frame() %>%
    filter(main_clone=='B' & mainsite=='Primary' & experiment=='main_exp')
  sidCB <- c(rownames(metaC), rownames(metaBpri))
  ddsCB <- dds[,sidCB]
  print(sizeFactors(ddsCB)) ## using existing size factors calculated using Scran method
  resCB <- get_DE_genes_DESeq2(ddsCB, DE_comp=c("Metastasis","Primary"),
                               filter_genes=F, min_total_exp_samples=10)
  dim(resCB)
  
  resCB1 <- resCB %>%
    as.data.frame() %>%
    filter(abs(log2FoldChange)>=0.05) # 
  dim(resCB1)
  sum(resCB1$ensembl_gene_id %in% resBB1$ensembl_gene_id)
  resCB2 <- resCB1 %>%
    filter(ensembl_gene_id %in% resBB1$ensembl_gene_id)
  sum(resCB2$log2FoldChange>0)
  dim(resCB2)
  obs_genes <- resCB2$ensembl_gene_id
  ref <- annotables::grch38 %>%
    dplyr::select(ensgene,symbol,chr) %>%
    dplyr::rename(ensembl_gene_id=ensgene) %>%
    dplyr::filter(ensembl_gene_id %in% resCB2$ensembl_gene_id)
  ref <- ref[!duplicated(ref$ensembl_gene_id),]
  dim(ref)
  ref1 <- annotables::grch38 %>%
    dplyr::select(ensgene,symbol,chr) %>%
    dplyr::rename(ensembl_gene_id=ensgene) %>%
    dplyr::filter(ensembl_gene_id %in% resCB1$ensembl_gene_id)
  ref1 <- ref1[!duplicated(ref1$ensembl_gene_id),]
  dim(ref1)
  resCB1 <- resCB1 %>%
    dplyr::left_join(ref1, by='ensembl_gene_id')
  resCB2 <- resCB2 %>%
    dplyr::left_join(ref, by='ensembl_gene_id')
  
  resCB2 <- resCB2 %>%
    dplyr::left_join(resBB2, by='ensembl_gene_id')
  sum(resCB2$log2FoldChange>0 & resCB2$log2FoldChange_BB>0)
  sum(resCB2$log2FoldChange<0 & resCB2$log2FoldChange_BB<0)
  dim(resCB2)
  obs_genes_symb <- resCB2$symbol
  pw_df <- get_gprofiler_pathways_obsgenes(obs_genes_symb, save_dir, datatag, 
                                           custom_id="gp__RYib_WEFZ_gEw", pathway_fn=NULL, save_data=T)
  dim(pw_df)
  # View(pw_df)  
}

detect_cis_genes <- function(cnv){
  resCB3 <- resCB2 %>%
    dplyr::filter(symbol %in% cnv$gene_symbol)
  dim(resCB3)
  dim(resCB2)
  summary(resCB3$log2FoldChange)
  resCB2$ensembl_gene_id[1:3]
  sum(resCB1$symbol %in% cnv$gene_symbol)
  resCB11_up <- resCB1 %>%
    dplyr::filter(symbol %in% cnv$gene_symbol & log2FoldChange>0)
  resCB11_down <- resCB1 %>%
    dplyr::filter(symbol %in% cnv$gene_symbol & log2FoldChange<0)
  sum(resCB11$log2FoldChange>0)
  summary(resCB11$log2FoldChange)
  dim(resCB11)
  dim(resCB11_up)
  dim(resCB11_down)
  
  pw_set <- pw_df$reference_set[2]
  ref_genes <- ref_set[[pw_set]]
  length(ref_genes)
  sum(resCB11$symbol %in% ref_genes)
  resCB11$symbol[resCB11$symbol %in% ref_genes]
}
get_cis_genes <- function(save_dir){
  script_dir <- "/home/htran/Projects/hakwoo_project/metastasis_material/materials/bulkRNAseq/"
  cnv <- data.table::fread(paste0(script_dir, 'SA919_DE_analysis/mapping_gene_cnv_SA919.csv.gz'))  
  dim(cnv)
  cnv <- cnv[!duplicated(cnv$ensembl_gene_id),]
  # head(cnv)
  
  rownames(cnv) <- cnv$ensembl_gene_id
  rv <- rowVars(as.matrix(cnv[,c(5, 6)]))  # B, C median copy number profile
  # length(rv)
  cnv <- as.data.frame(cnv)
  genes_used <- cnv$ensembl_gene_id[rv>0]
  # length(genes_used)
  cnv <- cnv %>%
    dplyr::filter(ensembl_gene_id %in% genes_used)
  dim(cnv)
  data.table::fwrite(cnv, paste0(save_dir, 'cnv_variance_CB.csv.gz'))
  # cnv <- data.table::fread(paste0(save_dir, 'cnv_variance_CB.csv.gz'))
  # sum(cnv$B < cnv$C)
  return(cnv)
}
get_DE_genes_DESeq2 <- function(dds, DE_comp=c("Metastasis","Primary"),
                                filter_genes=F, min_total_exp_samples=10){
  
  print(dim(dds))
  print(sizeFactors(dds))
  print(colData(dds))
  ## ----prefilter----------------------------------------------------------------
  ## filter genes function pose an issue for next step, need to debug it later
  ## set filter_genes=F first
  if(filter_genes){
    # smallestGroupSize <- 3 # group untreated with 3 samples is the smallest group
    # keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
    keep <- rowSums(counts(dds) >= min_total_exp_samples)
    dds <- dds[keep,]
    print(dim(dds))  
  }
  print(colnames(dds))
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

# Run DESeq2 for B met vs B pri, C met vs B pri, C met mix vs B pri
# Get common cis/ trans genes 
# Comparing results with scRNA-seq at cluster level. 