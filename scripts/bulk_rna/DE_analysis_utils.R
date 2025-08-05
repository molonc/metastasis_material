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
get_BCV_dispersion <- function(dge, meta_genes, cis_ens_genes, 
                               obs_ens_genes, save_dir, datatag){
  # cnv <- data.table::fread(paste0(save_dir, 'cnv_variance_CB.csv.gz'))
  # de_genes_df <- data.table::fread(paste0(input_dir, 'SA919_CMetMainMix_Bpri_DESeq2/DE_signif_genes.csv.gz'))
  # dim(de_genes_df)
  
  # #Visualize the dispersion estimates with a BCV plot
  # png(paste0(save_dir,datatag,"_glmQLF_plotBCV_check.png"), height = 2*350, width=2*400, res = 2*72)
  # # jpeg("glmQLF_plotBCV.jpg")
  # plotBCV(dge)
  # dev.off()
  
  # plotMDS(dge)
  
  # A <- dge$AveLogCPM
  # if (is.null(A)) 
  #   A <- aveLogCPM(dge$counts, offset = getOffset(y))
  if (!is.null(dge$tagwise.dispersion)) {
    print('Gene wise dispersion output: ')
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
  
  # meta_genes <- meta_genes %>%
  #   dplyr::filter(ens_gene_id %in% obs_ens_genes)%>%
  #   dplyr::mutate(gene_type=
  #                   case_when(ens_gene_id %in% cis_ens_genes ~ 'cis',
  #                             TRUE ~ 'trans'))
  meta_genes$tagwise_dispersion <- round(meta_genes$tagwise_dispersion, 3)
  data.table::fwrite(meta_genes, paste0(save_dir, datatag, '_dispersion_cis_trans_genes.csv'))
  print(summary(as.factor(meta_genes$gene_type)))  
  dim(meta_genes)
  # View(head(meta_genes))
  # cis_df <- meta_genes %>%
  #   dplyr::filter(gene_type=='cis')
  # data.table::fwrite(cis_df, paste0(save_dir, datatag, '_dispersion_cis_genes.csv'))
  
  
  # library("ggplot2")
  p <- ggplot(meta_genes, aes(x=gene_type, y=tagwise_dispersion)) + #
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(size=0.4, aes(colour = gene_type), position=position_jitter(0.2))  + 
    theme_bw() + 
    labs(x='Gene type', y='Gene wise dispersion')
  
  png(paste0(save_dir,datatag,"_cis_trans_disperion.png"), height = 2*300, width=2*300, res = 2*72)
  print(p)
  dev.off()
  saveRDS(p, paste0(save_dir,datatag,"_cis_trans_disperion_plot.rds"))
  
  res <- list(p_cistrans=p, meta_genes=meta_genes)
  return(res)
}
