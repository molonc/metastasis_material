suppressPackageStartupMessages({
  library(tidyverse)
  library(annotables)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  # library(scales)
  # library(ggrepel)
  library(stringr)
  # library(scran)
  library(SingleCellExperiment)
  # library(tximport)
  
})

## Loading utility functions
script_dir <- '/home/htran/Projects/hakwoo_project/metastasis_material/'
script_dir <- '/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/metastasis_material/'
source(paste0(script_dir, 'scripts/bulk_rna/bulk_utils.R'))


## To Do: loading data from all exps  
## First loading metadata, check full list of samples from clone A, B, C
load_input_data <- function(){
  input_dir <- "/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/metastasis_material/materials/bulkRNAseq/preprocessed_09April2024/"
  save_dir <- "/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/results_bulkRNAseq/"
  input_data_dir <- "/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/bulk_SA919/"
  
  datatag <- 'SA919_full'
  save_dir <- paste0(save_dir, datatag,'/')
  dir.create(save_dir, recursive=T)
  
  meta_samples <- data.table::fread(paste0(input_dir, '../metadata/metadata_Hakwoo_bulkRNA_mixing_main_exp.csv'))
  dim(meta_samples)
  # View(meta_samples)
  meta_samples_main <- meta_samples %>%
    dplyr::filter(experiment=='main_exp')
  meta_samples_mixing <- meta_samples %>%
    dplyr::filter(experiment=='mixing_exp')
  
  df_counts_main <- load_raw_counts_kallisto(input_data_dir, datatag, save_dir, 
                                             meta_samples_main$sample_id)
  
  df_counts_mixing <- load_raw_counts_kallisto(input_data_dir, datatag, save_dir, 
                                               meta_samples_mixing$bulk_sid)
  
  head(df_counts_main)
  head(df_counts_mixing)
  
  df <- dplyr::bind_rows(df_counts_main, df_counts_mixing)
  dim(df)
  df <- df %>%
    tidyr::pivot_wider(values_from = 'raw_counts', names_from = 'sample_id')
  # View(head(df))
  df1 <- df %>%
    na.omit()
  dim(df1)
  # View(head(df1))
  
  datatag <- 'SA919_full'
  df_counts_fn <- paste0(save_dir, datatag,'_total_raw_counts.csv.gz')
  data.table::fwrite(df1, df_counts_fn)
  
  ## Normalizing total 22 samples from all main, and mixing exp, and clone A,B,C
  df_normalized <- normalize_by_size_factor(df_counts_fn, datatag, save_dir)
  dim(df_normalized)
  sce <- readRDS(paste0(save_dir, datatag,'_sizefactor_normalized.rds'))
  dim(sce)
  df_normalized <- data.table::fread(paste0(save_dir, datatag, '_sizefactor_normalized.csv.gz'))
  ## Creating sce object for future use 
  ## Creating edgeR object for different cases: B vs. A, C vs. A
  ## 
  ## 
  ## 
  ## 
  ## 
  sf_df <- data.table::fread(paste0(save_dir, datatag, '_sizefactors_Scran.csv.gz'))
  sf_df <- colData(sce)
  sf_df$sample <- rownames(sf_df)
  head(sf_df)
  cts <- data.table::fread(df_counts_fn) %>% as.data.frame()
  print(dim(cts))
  # head(cts)
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
  
  head(df_normalized)
  df_normalized <- df_normalized %>%
    tibble::column_to_rownames('ens_gene_id') %>%
    dplyr::select(all_of(colnames(cts)))
  sum(colnames(df_normalized)==colnames(cts))
  dim(df_normalized)
  df_normalized <- df_normalized[rownames(cts),]
  sids <- colnames(df_normalized)
  for(s in sids){
    df_normalized[,s] <- round(df_normalized[,s],1)
  }
  
  meta_genes <- tibble::tibble(ens_gene_id=rownames(cts))
  ref <- annotables::grch38 %>%
    dplyr::select(ensgene,symbol,chr) %>%
    dplyr::rename(ens_gene_id=ensgene) %>%
    dplyr::filter(ens_gene_id %in% meta_genes$ens_gene_id)
  ref <- ref[!duplicated(ref$ens_gene_id),]
  dim(ref)
  meta_genes <- meta_genes %>%
    dplyr::left_join(ref, by='ens_gene_id')
  dim(meta_genes) ## some genes do not have gene symbol 
  sum(meta_genes$ens_gene_id==rownames(cts))
  
  
  ## Save data as SingleCellExperiment object for future uses
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts=as.matrix(cts), 
                                                         normcounts=as.matrix(df_normalized)),
                                                         colData=coldata,
                                                         rowData=meta_genes)
  
  print(sce)
  saveRDS(sce, paste0(save_dir, datatag, '_sce.rds'))
}  
create_DESeq2_object <- function(sce){
  coldata <- colData(sce) %>% as.data.frame()
  dim(coldata)
  cts <- counts(sce)
  dim(cts)
  dds <- create_DESeq2_obj(coldata, cts, use_existing_size_factor=T)
  return(dds)
}
extract_DE_genes_DESeq2 <- function(dds){
  library(DESeq2)
  ## How to select samples here? 
  # backup_dds <- dds
  # dds <- backup_dds
  input_dir <- "/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/metastasis_material/materials/bulkRNAseq/preprocessed_09April2024/"
  save_dir <- "/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/results_bulkRNAseq/SA919_full/"
  input_data_dir <- "/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/bulk_SA919/"
  datatag <- 'SA919_full'
  
  
  subtag <- 'Bmet_Apri'
  save_fig_dir <- paste0(save_dir, subtag, '/')
  # dir.create(save_fig_dir)
  meta_df <- colData(dds) %>%
    as.data.frame() %>%
    dplyr::mutate(desc=paste0(main_clone, '_',mainsite)) %>%
    filter(desc %in% c('B_Metastasis','A_Primary'))
  dim(meta_df)
  sids <- rownames(meta_df)
  dds_tmp <- dds[,sids]
  print(sizeFactors(dds_tmp)) ## using existing size factors calculated using Scran method
  res1 <- get_DE_genes_DESeq2(dds_tmp, DE_comp=c("Metastasis","Primary"),
                               filter_genes=F, min_total_exp_samples=10)
  dim(res1)
  res1 <- res1 %>%
    as.data.frame() %>%
    filter(abs(log2FoldChange)>=1) #pvalue<0.05 & 
  dim(res1)
  head(res1)
  data.table::fwrite(res1, paste0(save_fig_dir, subtag, '_DE_genes.csv.gz'))
  
  
  subtag <- 'Cmet_Apri'
  save_fig_dir <- paste0(save_dir, subtag, '/')
  dir.create(save_fig_dir)
  meta_df <- colData(dds) %>%
    as.data.frame() %>%
    dplyr::mutate(desc=paste0(main_clone, '_',mainsite)) %>%
    filter(desc %in% c('C_Metastasis','A_Primary'))
  dim(meta_df)
  sids <- rownames(meta_df)
  dds_tmp <- dds[,sids]
  print(sizeFactors(dds_tmp)) ## using existing size factors calculated using Scran method
  res2 <- get_DE_genes_DESeq2(dds_tmp, DE_comp=c("Metastasis","Primary"),
                              filter_genes=F, min_total_exp_samples=10)
  dim(res2)
  res2 <- res2 %>%
    as.data.frame() %>%
    filter(abs(log2FoldChange)>=1) #pvalue<0.05 & 
  dim(res2)
  head(res2)
  data.table::fwrite(res2, paste0(save_fig_dir, subtag, '_DE_genes.csv.gz'))
  
}  
create_edgeR_object <- function(counts_df, meta_df){
  
  groups_use <- c('Metastasis','Primary')
  head(meta_df)
  #edgeR package use inverted order compared to general convention in other DE tools
  # i.e: Rx vs. UnRx
  # edgeR: 1_UnRx vs. 2_Rx
  meta_df$condition <- ifelse(grepl('Primary', meta_df$condition),paste0("1_",'Primary'),
                              paste0("2_",'Metastasis'))
  print(meta_df$condition)
  meta_df$condition <- factor(meta_df$condition,levels = c("1_Primary","2_Metastasis"))
  # rownames(meta_df) <- meta_df$sample
  # counts_df <- load_counts_data(meta_df$sample, input_dir, datatag)
  # counts_df <- counts_df %>%
  #   dplyr::select(rownames(meta_df))
  
  print(colnames(counts_df))
  dge <- get_dgeList(counts_df, meta_df)
  return(dge)
  
}
get_dgeList <- function(counts_df, meta_df){
  # Creating a DGEList object for use in edgeR.
  groups_use <- meta_df$condition
  dge <- edgeR::DGEList(counts_df, group = groups_use)
  
  # dge <- edgeR::calcNormFactors(cts)
  # dge$genes <- meta_genes
  # normMat <- edgeR::cpm(cts, log=TRUE) 
  # dge <- edgeR::scaleOffset(dge, normMat)
  # filtering
  print(dim(dge))
  dge$samples$lib.size <- colSums(dge$counts)
  keep <- edgeR::filterByExpr(dge)
  dge <- dge[keep, ]
  print(dim(dge))
  print(class(dge))
  
  design <- model.matrix(~condition, data = meta_df)
  if(is.null(rownames(design))){
    rownames(design) <- meta_df$sample
  }
  # dge <- calcNormFactors(dge, method = "TMM")
  # sf_df <- data.table::fread(paste0(save_dir, datatag, '_sizefactors_Scran.csv.gz'))
  # sf_df <- sf_df %>%
  #   tibble::column_to_rownames('sample')
  # dge$samples$norm.factors <- sf_df[rownames(dge$samples),'size_factor']
  
  # dge$samples$norm.factors <- colData(sce)[rownames(dge$samples),'size_factor']
  dge$samples$norm.factors <- meta_df[rownames(dge$samples),'size_factor']
  dge <- edgeR::estimateDisp(dge, design = design)
  # dge <- edgeR::estimateCommonDisp(dge, design = design)
  # dge <- edgeR::estimateTagwiseDisp(dge)
  # dge <- edgeR::estimateGLMCommonDisp(dge, design = design)
  
  print(length(dge$tagwise.dispersion))
  return(dge)
}
viz_BCV_dispersion <- function(dge, meta_genes, cis_ens_genes, 
                               obs_ens_genes, save_dir, datatag){
  # cnv <- data.table::fread(paste0(save_dir, 'cnv_variance_CB.csv.gz'))
  # de_genes_df <- data.table::fread(paste0(input_dir, 'SA919_CMetMainMix_Bpri_DESeq2/DE_signif_genes.csv.gz'))
  # dim(de_genes_df)
  
  #Visualize the dispersion estimates with a BCV plot
  png(paste0(save_dir,datatag,"_glmQLF_plotBCV_check.png"), height = 2*350, width=2*400, res = 2*72)
  # jpeg("glmQLF_plotBCV.jpg")
  plotBCV(dge)
  dev.off()
  
  # plotMDS(dge)
  
  A <- dge$AveLogCPM
  if (is.null(A)) 
    A <- aveLogCPM(dge$counts, offset = getOffset(y))
  if (!is.null(dge$tagwise.dispersion)) {
    print('Hello')
    print(length(dge$tagwise.dispersion))
    df <- data.frame(ens_gene_id=rownames(dge), 
                     tagwise_dispersion=dge$tagwise.dispersion)
  }else{
    stop('Issue with dispersion calculation, double check!!!')
  }  
  # meta_genes <- rowData(sce) %>% as.data.frame() 
  meta_genes <- meta_genes %>%
    dplyr::inner_join(df, by='ens_gene_id')
  dim(meta_genes)
  
  meta_genes <- meta_genes %>%
    dplyr::filter(ens_gene_id %in% obs_ens_genes)%>%
    dplyr::mutate(gene_type=
                    case_when(ens_gene_id %in% cis_ens_genes ~ 'cis',
                              TRUE ~ 'trans'))
  meta_genes$tagwise_dispersion <- round(meta_genes$tagwise_dispersion, 3)
  data.table::fwrite(meta_genes, paste0(save_dir, datatag, '_dispersion_cis_trans_genes.csv'))
  
  dim(meta_genes)
  # View(head(meta_genes))
  cis_df <- meta_genes %>%
    dplyr::filter(gene_type=='cis')
  data.table::fwrite(cis_df, paste0(save_dir, datatag, '_dispersion_cis_genes.csv'))
  
  
  # library("ggplot2")
  p <- ggplot(meta_genes, aes(x=gene_type, y=tagwise_dispersion)) + #
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(size=0.4, aes(colour = gene_type), position=position_jitter(0.2))  + 
    theme_bw() + 
    labs(x='Gene type', y='Gene wise dispersion')
  
  png(paste0(save_dir,datatag,"_cis_trans_disperion.png"), height = 2*300, width=2*300, res = 2*72)
  print(p)
  dev.off()
  
  
  cis_disp <- meta_genes %>%
    dplyr::filter(gene_type=='cis') %>%
    dplyr::pull(tagwise_dispersion)
  trans_disp <- meta_genes %>%
    dplyr::filter(gene_type=='trans') %>%
    dplyr::pull(tagwise_dispersion)
  out_stat <- ks.test(cis_disp, trans_disp, alternative='less')
  
  res <- list(p_cistrans=p, out_stat=out_stat, meta_genes=meta_genes)
  return(res)
}
get_gene_wise_dispersion_edgeR <- function(sce, datatag, save_dir){
  
  
  ## B met vs. A pri
  ## C met vs. A pri
  library("edgeR")
  input_dir <- "/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/metastasis_material/materials/bulkRNAseq/preprocessed_09April2024/"
  save_dir <- "/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/results_bulkRNAseq/SA919_full/"
  input_data_dir <- "/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/bulk_SA919/"
  datatag <- 'SA919_full'
  
  ## B met vs. A pri
  ## C met vs. A pri
  sce <- readRDS(paste0(save_dir, datatag, '_sce.rds'))
  meta_genes <- rowData(sce) %>% as.data.frame() 
  meta_df <- as.data.frame(colData(sce))
  counts_df <- as.data.frame(counts(sce))
  dim(meta_df)
  meta_df$desc <- paste0(meta_df$main_clone, '_',meta_df$mainsite) 
  unique(meta_df$desc)
  meta_df$condition <- meta_df$mainsite
  meta_df$sample <- rownames(meta_df)
  
  
  meta_tmp <- meta_df %>%
    dplyr::filter(desc %in% c('B_Metastasis','A_Primary'))
  # meta_df <- meta_df %>%
  #   dplyr::filter(desc %in% c('A_Metastasis','A_Primary'))
  dim(meta_tmp)
  # View(meta_df)
  counts_tmp <- counts_df %>%
    dplyr::select(rownames(meta_tmp))
  dim(counts_tmp)
  head(counts_tmp)
  # sce$main_clone
  dge <- create_edgeR_object(counts_tmp, meta_tmp)
  
  subtag <- 'Bmet_Apri'
  save_fig_dir <- paste0(save_dir, subtag, '/')
  # dir.create(save_fig_dir)
  saveRDS(dge, paste0(save_fig_dir,'dge.rds'))
  dge <- readRDS(paste0(save_fig_dir,'dge.rds'))
  
  ## Loading cnv file
  obs_clones <- c('A','B')
  cnv <- get_cis_trans_gene_type(obs_clones)
  dim(cnv)
  cis_ens_genes <- cnv$ensembl_gene_id
  de_genes_df <- data.table::fread(paste0(save_fig_dir, subtag, '_DE_genes.csv.gz'))
  obs_ens_genes <- de_genes_df$ensembl_gene_id
  res <- viz_BCV_dispersion(dge, meta_genes, cis_ens_genes, obs_ens_genes, save_fig_dir, subtag)
  res$out_stat
  
  
  
  
  meta_tmp <- meta_df %>%
    dplyr::filter(desc %in% c('C_Metastasis','A_Primary'))
  # meta_df <- meta_df %>%
  #   dplyr::filter(desc %in% c('A_Metastasis','A_Primary'))
  dim(meta_tmp)
  # View(meta_df)
  counts_tmp <- counts_df %>%
    dplyr::select(rownames(meta_tmp))
  dim(counts_tmp)
  head(counts_tmp)
  # sce$main_clone
  dge <- create_edgeR_object(counts_tmp, meta_tmp)
  
  subtag <- 'Cmet_Apri'
  save_fig_dir <- paste0(save_dir, subtag, '/')
  # dir.create(save_fig_dir)
  saveRDS(dge, paste0(save_fig_dir,'dge.rds'))
  dge <- readRDS(paste0(save_fig_dir,'dge.rds'))
  
  ## Loading cnv file
  obs_clones <- c('A','C')
  cnv <- get_cis_trans_gene_type(obs_clones)
  dim(cnv)
  cis_ens_genes <- cnv$ensembl_gene_id
  length(cis_ens_genes)
  de_genes_df <- data.table::fread(paste0(save_fig_dir, subtag, '_DE_genes.csv.gz'))
  obs_ens_genes <- de_genes_df$ensembl_gene_id
  length(obs_ens_genes)
  sum(cis_ens_genes %in% obs_ens_genes)
  res <- viz_BCV_dispersion(dge, meta_genes, cis_ens_genes, obs_ens_genes, save_fig_dir, subtag)
  res$out_stat
  
  
  
}



get_cis_trans_gene_type <- function(obs_clones){
  ## Loading cnv file
  # script_dir <- "/home/htran/Projects/hakwoo_project/metastasis_material/materials/bulkRNAseq/"
  script_dir <- "/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/metastasis_material/materials/bulkRNAseq/"
  cnv <- data.table::fread(paste0(script_dir, 'SA919_DE_analysis_DESeq2_Hoa_09April2024/mapping_gene_cnv_SA919.csv.gz'))  
  dim(cnv)
  cnv <- cnv[!duplicated(cnv$ensembl_gene_id),]
  # head(cnv)
  rownames(cnv) <- cnv$ensembl_gene_id
  cnv_vals <- cnv %>%
    dplyr::select(all_of(obs_clones))
  rv <- rowVars(as.matrix(cnv_vals))  # median copy number profile of clone B, C 
  # length(rv)
  cnv <- as.data.frame(cnv)
  genes_used <- cnv$ensembl_gene_id[rv>0]
  # length(genes_used)
  cnv <- cnv %>%
    dplyr::filter(ensembl_gene_id %in% genes_used)
  dim(cnv)
  return(cnv)
}
