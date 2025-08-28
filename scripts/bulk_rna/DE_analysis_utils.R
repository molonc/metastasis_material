



create_DESeq2_object <- function(sce){
  coldata <- colData(sce) %>% as.data.frame()
  dim(coldata)
  cts <- counts(sce)
  dim(cts)
  dds <- create_DESeq2_obj(coldata, cts, use_existing_size_factor=T)
  return(dds)
}


get_DE_results <- function(subtag, save_fig_dir, dds, meta_genes, obs_clones, exps=c('main_exp','mixing_exp')){
    
  filtered_conds <- c(paste0(obs_clones[2],'_Metastasis'),paste0(obs_clones[1],'_Primary'))
  print(filtered_conds)
  meta_df <- colData(dds) %>%
    as.data.frame() %>%
    dplyr::mutate(desc=paste0(main_clone, '_',mainsite)) %>%
    filter(desc %in% filtered_conds & experiment %in% exps)
  print('Number of observed samples:')
  print(dim(meta_df))
  print(meta_df)
  sids <- rownames(meta_df)
  
  ## Special case, only check within a mouse
  # if(obs_clones[1]=='A' & obs_clones[2]=='A'){
  #   input_dir <- "/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/preprocessed_09April2024/"
  #   meta_samples <- data.table::fread(paste0(input_dir, '../metadata/metadata_Hakwoo_bulkRNA_mixing_main_exp.csv'))
  #   dim(meta_samples)
  #   # View(meta_samples)
  #   meta_samples <- meta_samples %>%
  #     dplyr::filter(experiment=='main_exp' & main_clone=='A' & pdxid=='X08472164')
  #   sids <- meta_samples$sample_id
  # }
  dds_tmp <- dds[,sids]
  print(sizeFactors(dds_tmp)) ## using existing size factors calculated using Scran method
  res <- get_DE_genes_DESeq2(dds_tmp, DE_comp=c("Metastasis","Primary"),
                              filter_genes=F, min_total_exp_samples=10)
  dim(res)
  res <- res %>%
    as.data.frame() %>%
    filter(pvalue<0.05 & abs(log2FoldChange)>=1) # 
  print(dim(res))
  
  res <- res %>%
    dplyr::left_join(meta_genes, by=c('ensembl_gene_id'='ens_gene_id'))
  data.table::fwrite(res, paste0(save_fig_dir, subtag, '_DE_genes.csv.gz'))
  
  return(res)
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
