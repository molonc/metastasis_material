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
  library(edgeR)
})

## Loading utility functions
script_dir <- '/home/htran/Projects/hakwoo_project/metastasis_material/'
script_dir <- '/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/metastasis_material/'
script_dir <- '/Users/hoatran/Documents/BCCRC_projects/hakwoo_project/code/metastasis_material/'
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
  
  input_dir <- "/Users/hoatran/Documents/BCCRC_projects/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/preprocessed_09April2024/"
  save_dir <- "/Users/hoatran/Documents/BCCRC_projects/hakwoo_project/code/results_bulkRNAseq/SA919_full/"
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
  df_normalized <- normalize_by_size_factor_v2(df_counts_fn, datatag, save_dir)
  dim(df_normalized)
  sce <- readRDS(paste0(save_dir, datatag,'_sizefactor_normalized.rds'))
  dim(sce)
  colSums(normcounts(sce))
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
  
  # colSums(df_normalized)
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
get_statistical_test_gene_types <- function(meta_genes, alternative_theory="two.sided"){
  print('KS test for 2 groups: ')
  alternative_theory="two.sided"
  meta_genes <- meta_genes %>%
    dplyr::mutate(gene_type=
                    case_when(
                      gene_type %in% c('trans','trans_whole_genome') ~ 'trans_gene',
                      TRUE ~ 'cis_gene'
                    )) 
  cis_disp <- meta_genes %>%
    dplyr::filter(gene_type=='cis_gene') %>%
    dplyr::pull(tagwise_dispersion)
  trans_disp <- meta_genes %>%
    dplyr::filter(gene_type=='trans_gene') %>%
    dplyr::pull(tagwise_dispersion)
  out_stat <- ks.test(cis_disp, trans_disp, alternative=alternative_theory)
  print(out_stat)
  
  print('F-test variance testing for 2 groups: ')
  unique(meta_genes$gene_type)
  res.ftest <- var.test(tagwise_dispersion ~ gene_type, data = meta_genes)
  print(res.ftest)
  
  print('Bootstrap testing for 2 groups: ')
  bootstrap_stat <- get_bootstrap_stat_sampling(cis_disp, trans_disp, 
                              sampling_fraction=0.7, nsamples=1000, alternative_theory=alternative_theory)
  
  return(list(ks_stat=out_stat, ftest_stat=res.ftest, bootstrap_stat=bootstrap_stat))
}
get_cis_trans_gene_type <- function(obs_clones){
  ## Loading cnv file
  # script_dir <- "/home/htran/Projects/hakwoo_project/metastasis_material/materials/bulkRNAseq/"
  script_dir <- "/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/metastasis_material/materials/bulkRNAseq/"
  script_dir <- "/Users/hoatran/Documents/BCCRC_projects/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/"
  cnv <- data.table::fread(paste0(script_dir, 'SA919_DE_analysis_DESeq2_Hoa_09April2024/mapping_gene_cnv_SA919.csv.gz'))  
  dim(cnv)
  cnv <- cnv[!duplicated(cnv$ensembl_gene_id),]
  colnames(cnv)
  # head(cnv)
  cols_use <- c('ensembl_gene_id','gene_symbol','chr')
  cols_use <- c(cols_use, obs_clones)
  rownames(cnv) <- cnv$ensembl_gene_id
  cnv <- cnv %>%
    dplyr::select(all_of(cols_use))
  cnv_vals <- cnv %>%
    dplyr::select(all_of(obs_clones))
  rv <- rowVars(as.matrix(cnv_vals))  # median copy number profile of obs clones
  # length(rv)
  cnv <- as.data.frame(cnv)
  genes_used <- cnv$ensembl_gene_id[rv>0]
  # length(genes_used)
  cnv <- cnv %>%
    dplyr::filter(ensembl_gene_id %in% genes_used)
  dim(cnv)
  return(cnv)
}

# meta_genes <- res$meta_genes
# (sce_extract, res$meta_genes, obs_conds)
# save_dir <- save_fig_dir
get_avg_exp <- function(sce_extract, meta_genes, obs_conds, save_dir){
  
  tag <- paste(obs_conds, collapse='_vs_')
  # unique(sce_extract$experiment)
  # sce_extract <- sce_extract[,sce_extract$experiment=='main_exp']
  # colSums(normcounts(sce_extract))
  ## Check if we have clone B, C in sce meta samples? 
  dim(meta_genes)
  ## To Do: report for clone B, and clone C separately
  meta_genes_df <- meta_genes %>%
    dplyr::filter(ens_gene_id %in% rownames(sce_extract))
  cis_genes <- meta_genes %>%
    dplyr::filter(gene_type=='cis') %>%
    dplyr::pull(ens_gene_id)
  trans_genes <- meta_genes %>%
    dplyr::filter(gene_type=='trans') %>%
    dplyr::pull(ens_gene_id)
  trans_genes_total <- meta_genes %>%
    dplyr::filter(gene_type %in% c('trans','trans_whole_genome')) %>%
    dplyr::pull(ens_gene_id)
  summary(as.factor(meta_genes$gene_type))
  
  ## To do: get avg exp per clone here
  exp_df <- tibble::tibble()
  total_exp_df <- tibble::tibble()
  for(cds in obs_conds){
    print(cds)
    # cis_logCPM <- aveLogCPM(as.matrix(normcounts(sce_extract[cis_genes,sce_extract$desc==cds])),
    #                         normalized.lib.sizes=FALSE)
    cis_logCPM <- log2(as.matrix(normcounts(sce_extract[cis_genes,sce_extract$desc==cds]))+1)
    dim(cis_logCPM)
    cis_logCPM <- round(cis_logCPM, 2)
    # data.table::fwrite(cis_logCPM, paste0(save_dir, tag, '_cis_logCPM.csv.gz'))
    # trans_logCPM <- aveLogCPM(as.matrix(normcounts(sce_extract[trans_genes,sce_extract$desc==cds])),
    #                           normalized.lib.sizes=FALSE)
    trans_logCPM <- log2(as.matrix(normcounts(sce_extract[trans_genes,sce_extract$desc==cds]))+1)
    dim(trans_logCPM)
    trans_logCPM <- round(trans_logCPM, 2)
    # data.table::fwrite(trans_logCPM, paste0(save_dir, tag, '_trans_logCPM.csv.gz'))
    # trans_total_logCPM <- aveLogCPM(as.matrix(normcounts(sce_extract[trans_genes_total,sce_extract$desc==cds])),
    #                                 normalized.lib.sizes=FALSE)

    trans_total_logCPM <- log2(as.matrix(normcounts(sce_extract[trans_genes_total,sce_extract$desc==cds]))+1)
    # data.table::fwrite(trans_total_logCPM, paste0(save_dir, tag, '_trans_total_logCPM.csv.gz'))
    trans_total_logCPM <- round(trans_total_logCPM, 2)
    
    ## To Do: get trans - whole genome here
    df1 <- tibble::tibble(obs_condition=cds, gene_type='cis', mean_logCPM=round(mean(cis_logCPM),2), 
                          median_logCPM=round(median(cis_logCPM),2), sd_logCPM=round(sd(cis_logCPM),2))
    df2 <- tibble::tibble(obs_condition=cds, gene_type='trans',mean_logCPM=round(mean(trans_logCPM),2), 
                          median_logCPM=round(median(trans_logCPM),2), sd_logCPM=round(sd(trans_logCPM),2))
    df3 <- tibble::tibble(obs_condition=cds, gene_type='trans_whole_genome',mean_logCPM=round(mean(trans_total_logCPM),2), 
                          median_logCPM=round(median(trans_total_logCPM),2), sd_logCPM=round(sd(trans_total_logCPM),2))
    exp_df <- dplyr::bind_rows(exp_df, df1)
    exp_df <- dplyr::bind_rows(exp_df, df2)  
    exp_df <- dplyr::bind_rows(exp_df, df3)
    
    cis_logCPM <- cis_logCPM %>% 
      as_tibble()%>% 
      mutate(ens_gene_id=rownames(cis_logCPM))%>% 
      tidyr::pivot_longer(!ens_gene_id, names_to='sample_id', values_to='gene_exp')%>% 
      dplyr::mutate(obs_conds=cds, gene_type='cis')
    # trans_total_logCPM <- trans_total_logCPM %>% 
    #   as_tibble() %>%
    #   mutate(ens_gene_id=rownames(trans_total_logCPM))%>% 
    #   tidyr::pivot_longer(!ens_gene_id, names_to='sample_id', values_to='gene_exp')%>% 
    #   dplyr::mutate(obs_conds=cds, gene_type='trans_whole_genome')
    # 
    total_exp_df <- dplyr::bind_rows(total_exp_df, cis_logCPM)
    # total_exp_df <- dplyr::bind_rows(total_exp_df, trans_logCPM)
    # total_exp_df <- dplyr::bind_rows(total_exp_df, trans_total_logCPM)
  }
  # dim(total_exp_df)
  # View(head(total_exp_df))
  meta_genes$gene_type <- NULL
  total_exp_df <- total_exp_df %>%
    left_join(meta_genes, by='ens_gene_id')
  # View(exp_df)
  # for(cds in obs_conds){
  #   print(cds)
  #   cis_exp1 <- log2(as.matrix(normcounts(sce_extract[cis_genes,sce_extract$desc==obs_conds[1]]))+1)
  #   cis_exp2 <- log2(as.matrix(normcounts(sce_extract[cis_genes,sce_extract$desc==obs_conds[2]]))+1)
  #   # trans_logCPM <- log2(as.matrix(normcounts(sce_extract[trans_genes,sce_extract$desc==cds]))+1)
  #   # trans_total_exp <- log2(as.matrix(normcounts(sce_extract[trans_genes_total,sce_extract$desc==cds]))+1)
  #   alternative_theory <- 'two.sided'
  #   alternative_theory <- 'greater'
  #   out_stat <- ks.test(cis_exp1, cis_exp2, alternative=alternative_theory)
  #   print(out_stat)
  # }  
  data.table::fwrite(total_exp_df, paste0(save_dir, tag, '_total_exp.csv.gz'))
  data.table::fwrite(exp_df, paste0(save_dir, tag, '_avg_gene_exp.csv.gz'))
  return(list(exp_df=exp_df, total_exp_df=total_exp_df))

}
viz_gene_exp_comparison <- function(df, obs_clones, legend_pos='none'){
  # df <- total_exp_df
  
  ptittle <- paste0(obs_clones[2],' vs. ', obs_clones[1])
  # Convert the variable dose from a numeric to a factor variable
  df$chr <- as.factor(df$chr)
  my_font <- "Helvetica"
  p <- ggplot(df, aes(x=chr, y=gene_exp, color=obs_conds)) +
    # geom_jitter(aes(color=obs_conds), shape=16, position=position_jitter(0.2)) +
    geom_violin(position=position_dodge(0.8), alpha=0) +
    geom_jitter(position=position_jitterdodge(0.2), size=0.1) +
    # geom_boxplot(width=1)
    theme_bw() + 
    theme(legend.position = legend_pos,
        panel.grid = element_blank(), 
        text = element_text(color="black",hjust = 0.5, family=my_font),
        axis.text.x = element_text(color="black",size=12, hjust = 0.5, family=my_font),
        axis.title.x=element_text(color="black",size=10, hjust = 0.5, family=my_font),
        axis.text.y = element_text(color="black",size=6, hjust = 0.5, family=my_font),
        axis.title.y=element_text(color="black",size=10, hjust = 0.5, family=my_font)) +
    labs(title = ptittle, x='Chr position', y= 'Log2(exp + 1)') 
  # p
    # geom_jitter(shape=16, position=position_jitter(0.2)) #+
    # geom_dotplot(aes(y=len), dotsize=0.5)# + #binaxis='y', stackdir='center', 
    # geom_boxplot(aes()width=0.1, fill="white")#+
    # scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  
  return(p)
}
debug_program <- function(){
  save_dir <- '/Users/hoatran/Documents/BCCRC_projects/hakwoo_project/code/results_bulkRNAseq/SA919_full/'
  datatag <- 'SA919_full'
  
  ## B met vs. A pri
  ## C met vs. A pri
  sce <- readRDS(paste0(save_dir, datatag, '_sce.rds'))
  sce$desc <- paste0(sce$main_clone, '_', sce$mainsite)
  sce$condition <- sce$mainsite
  norm_df <- as.matrix(normcounts(sce))
  DelayedArray::colSums(norm_df)  
  DelayedArray::colSums(counts(sce))  
}
get_gene_wise_dispersion_edgeR <- function(sce, datatag, save_dir){
  
  
  ## B met vs. A pri
  ## C met vs. A pri
  # library("edgeR")
  input_dir <- "/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/metastasis_material/materials/bulkRNAseq/preprocessed_09April2024/"
  save_dir <- "/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/results_bulkRNAseq/SA919_full/"
  input_data_dir <- "/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/bulk_SA919/"
  
  save_dir <- '/Users/hoatran/Documents/BCCRC_projects/hakwoo_project/code/results_bulkRNAseq/SA919_full/'
  datatag <- 'SA919_full'
  
  ## B met vs. A pri
  ## C met vs. A pri
  sce <- readRDS(paste0(save_dir, datatag, '_sce.rds'))
  sce$desc <- paste0(sce$main_clone, '_', sce$mainsite)
  sce$condition <- sce$mainsite
  meta_genes <- rowData(sce) %>% as.data.frame() 
  counts_df <- as.data.frame(counts(sce))
  meta_df <- as.data.frame(colData(sce))
  # summary(as.factor(meta_df$main_clone))  
  dim(meta_df)
  meta_df$sample <- rownames(meta_df)
  
  
  obs_conds <- c('B_Metastasis','A_Primary')
  subtag <- 'Bmet_Apri'
  get_BMet_APri_dispersion(sce, meta_genes, meta_df, counts_df, save_dir, obs_conds, subtag)
  
  
  obs_conds <- c('C_Metastasis','A_Primary')
  subtag <- 'Cmet_Apri'
  get_CMet_APri_dispersion(sce, meta_genes, meta_df, counts_df, save_dir, obs_conds, subtag)
  
  obs_conds <- c('C_Metastasis','B_Primary')
  subtag <- 'Cmet_Bpri'
  get_CMet_BPri_dispersion(sce, meta_genes, meta_df, counts_df, save_dir, obs_conds, subtag)
}

#obs_conds <- c('C_Metastasis','B_Primary')

get_BMet_APri_dispersion <- function(sce, meta_genes, meta_df, counts_df, 
                                     save_dir, obs_conds, subtag='Bmet_Apri'){
  meta_tmp <- meta_df %>%
    dplyr::filter(desc %in% obs_conds)
  print(obs_conds)
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
  
  
  save_fig_dir <- paste0(save_dir, subtag, '/')
  # dir.create(save_fig_dir)
  saveRDS(dge, paste0(save_fig_dir,'dge.rds'))
  dge <- readRDS(paste0(save_fig_dir,'dge.rds'))
  
  ## Loading cnv file
  obs_clones <- c('A','B')
  cnv <- get_cis_trans_gene_type(obs_clones)
  dim(cnv)
  summary(as.factor(cnv$chr))
  # head(cnv)
  cis_ens_genes <- cnv$ensembl_gene_id
  
  ## To Do: filtering this list with logFC threshold
  de_genes_df <- data.table::fread(paste0(save_fig_dir, subtag, '_DE_genes.csv.gz'))
  summary(de_genes_df$log2FoldChange)
  sum(abs(de_genes_df$log2FoldChange)<0.5)
  obs_ens_genes <- de_genes_df$ensembl_gene_id
  length(obs_ens_genes)
  cis_ens_genes <- cis_ens_genes[cis_ens_genes %in% obs_ens_genes] # condition to be cis gene = cis overlapping regions, and DE genes
  trans_genes <- obs_ens_genes[!obs_ens_genes %in% cis_ens_genes]
  dim(meta_genes)
  
  ## Checking chr position of cis CNA region
  # t <- cnv %>%
  #   dplyr::filter(ensembl_gene_id %in% cis_ens_genes)
  # summary(as.factor(t$chr))
  
  
  meta_genes <- meta_genes %>%
    # dplyr::filter(ens_gene_id %in% obs_ens_genes)%>%
    dplyr::mutate(gene_type=
                    case_when(ens_gene_id %in% cis_ens_genes ~ 'cis',
                              ens_gene_id %in% trans_genes ~ 'trans',
                              TRUE ~ 'trans_whole_genome'))
  summary(as.factor(meta_genes$gene_type))
  res <- get_BCV_dispersion(dge, meta_genes, cis_ens_genes, obs_ens_genes, save_fig_dir, subtag)
  # meta_genes <- res$meta_genes
  sce_extract <- sce[,sce$desc %in% obs_conds]
  dim(sce_extract)
  dim(sce)
  res <- get_avg_exp(sce_extract, res$meta_genes, obs_conds)
  p <- viz_gene_exp_comparison(res$total_exp_df, obs_clones, legend_pos='none')
  saveRDS(p, paste0(save_fig_dir, subtag, '_gene_exp_plt.rds'))
  # View(exp_df)
  # exp_df1 <- exp_df
  # View(exp_df1)
  # View(exp_df)
}

get_CMet_APri_dispersion <- function(sce, meta_genes, meta_df, counts_df, save_dir, obs_conds, subtag='Cmet_Apri'){
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
  
  # subtag <- 'Cmet_Apri'
  save_fig_dir <- paste0(save_dir, subtag, '/')
  # dir.create(save_fig_dir)
  saveRDS(dge, paste0(save_fig_dir,'dge.rds'))
  # dge <- readRDS(paste0(save_fig_dir,'dge.rds'))
  
  ## Loading cnv file
  obs_clones <- c('A','C')
  cnv <- get_cis_trans_gene_type(obs_clones)
  dim(cnv)
  cis_ens_genes <- cnv$ensembl_gene_id
  length(cis_ens_genes)
  de_genes_df <- data.table::fread(paste0(save_fig_dir, subtag, '_DE_genes.csv.gz'))
  sum(abs(de_genes_df$log2FoldChange)<0.5)
  obs_ens_genes <- de_genes_df$ensembl_gene_id
  length(obs_ens_genes)
  sum(cis_ens_genes %in% obs_ens_genes)
  cis_ens_genes <- cis_ens_genes[cis_ens_genes %in% obs_ens_genes] # condition to be cis gene = cis overlapping regions, and DE genes
  trans_genes <- obs_ens_genes[!obs_ens_genes %in% cis_ens_genes]
  dim(meta_genes)
  meta_genes <- meta_genes %>%
    # dplyr::filter(ens_gene_id %in% obs_ens_genes)%>%
    dplyr::mutate(gene_type=
                    case_when(ens_gene_id %in% cis_ens_genes ~ 'cis',
                              ens_gene_id %in% trans_genes ~ 'trans',
                              TRUE ~ 'trans_whole_genome'))
  summary(as.factor(meta_genes$gene_type))
  res <- get_BCV_dispersion(dge, meta_genes, cis_ens_genes, obs_ens_genes, save_fig_dir, subtag)
  sce_extract <- sce[,sce$desc %in% obs_conds]
  dim(sce_extract)
  dim(sce)
  meta_genes <- res$meta_genes
  exp_df <- get_avg_exp(sce_extract, meta_genes, obs_conds)
  exp_df$datatag <- subtag
  data.table::fwrite(exp_df, paste0(save_fig_dir, subtag, '_avg_gene_exp.csv.gz'))
  # View(exp_df)
  
}
get_CMet_BPri_dispersion <- function(sce, meta_genes, meta_df, counts_df, save_dir, obs_conds, subtag='Cmet_Bpri'){
  meta_tmp <- meta_df %>%
    dplyr::filter(desc %in% obs_conds)
  print(obs_conds)
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
  
  
  save_fig_dir <- paste0(save_dir, subtag, '/')
  # dir.create(save_fig_dir)
  saveRDS(dge, paste0(save_fig_dir,'dge.rds'))
  # dge <- readRDS(paste0(save_fig_dir,'dge.rds'))
  
  ## Loading cnv file
  obs_clones <- c('B','C')
  cnv <- get_cis_trans_gene_type(obs_clones)
  dim(cnv)
  head(cnv)
  cis_ens_genes <- cnv$ensembl_gene_id
  
  ## To Do: filtering this list with logFC threshold
  # de_genes_df <- data.table::fread(paste0(save_fig_dir, subtag, '_DE_genes.csv.gz'))
  input_dir <- "/Users/hoatran/Documents/BCCRC_projects/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/preprocessed_09April2024/"
  de_genes_df <- data.table::fread(paste0(input_dir, 'SA919_CMetMainMix_Bpri_DESeq2/DE_signif_genes.csv.gz'))
  dim(de_genes_df)
  sum(abs(de_genes_df$log2FoldChange)<0.5)
  obs_ens_genes <- de_genes_df$ens_gene_id
  # obs_ens_genes <- de_genes_df$ensembl_gene_id
  length(obs_ens_genes)
  cis_ens_genes <- cis_ens_genes[cis_ens_genes %in% obs_ens_genes] # condition to be cis gene = cis overlapping regions, and DE genes
  trans_genes <- obs_ens_genes[!obs_ens_genes %in% cis_ens_genes]
  dim(meta_genes)
  meta_genes <- meta_genes %>%
    # dplyr::filter(ens_gene_id %in% obs_ens_genes)%>%
    dplyr::mutate(gene_type=
                    case_when(ens_gene_id %in% cis_ens_genes ~ 'cis',
                              ens_gene_id %in% trans_genes ~ 'trans',
                              TRUE ~ 'trans_whole_genome'))
  summary(as.factor(meta_genes$gene_type))
  res <- get_BCV_dispersion(dge, meta_genes, cis_ens_genes, obs_ens_genes, save_fig_dir, subtag)
  # meta_genes <- res$meta_genes
  sce_extract <- sce[,sce$desc %in% obs_conds]
  dim(sce_extract)
  dim(sce)
  exp_df <- get_avg_exp(sce_extract, res$meta_genes, obs_conds)
  exp_df$datatag <- subtag
  data.table::fwrite(exp_df, paste0(save_fig_dir, subtag, '_avg_gene_exp.csv.gz'))
  # View(exp_df)
}
viz_stat <- function(){
  stat_total_df <- meta_genes %>%
    dplyr::mutate(gene_type_desc=
                    case_when(
                      gene_type %in% c('trans','trans_whole_genome') ~ 'trans_gene',
                      TRUE ~ 'cis_gene'
                    )) #%>%
    # dplyr::group_by(gene_type_desc) #%>%
    # dplyr::summarise(mean_dispersion=round(mean(tagwise_dispersion),2),
    #                  median_dispersion=round(median(tagwise_dispersion),2),
    #                  sd_dispersion=round(sd(tagwise_dispersion),2), nb_genes=n())
  p <- ggplot(stat_total_df, aes(x=gene_type_desc, y=tagwise_dispersion)) + #
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(size=0.4, aes(colour = gene_type_desc), position=position_jitter(0.2))  + 
    theme_bw() + 
    theme(legend.position = 'none') + 
    labs(x='Gene type', y='Gene wise dispersion')
  summary(as.factor(stat_total_df$gene_type_desc))
  png(paste0(save_dir, subtag, '/',subtag,"_cis_trans_disperion.png"), height = 2*300, width=2*250, res = 2*72)
  print(p)
  dev.off()
  
  
  
}
get_statistical_test_gene_types <- function(exp_df, alternative_theory="two.sided"){
  colnames(exp_df)
  print('KS test for 2 groups: ')
  cds <- unique(exp_df$obs_condition)
  cd_met <- cds[grepl('Met',cds)]
  cd_pri <- cds[grepl('Pri',cds)]
  exp_df <- exp_df %>%
    dplyr::filter(gene_type=='cis')
  g1 <- exp_df %>%
    dplyr::filter(obs_condition==cd_met) %>%
    # dplyr::mutate(gene_type_desc=paste0(gene_type, '_', obs_condition))%>%
    dplyr::pull(median_logCPM)
  g2 <- exp_df %>%
    dplyr::filter(obs_condition==cd_pri) %>%
    # dplyr::mutate(gene_type_desc=paste0(gene_type, '_', obs_condition))%>%
    dplyr::select(median_logCPM)
  alternative_theory <- 'greater'
  out_stat <- ks.test(g1, g2, alternative=alternative_theory)
  print(out_stat)
  
  
  print('F-test variance testing for 2 groups: ')
  res.ftest <- var.test(median_logCPM ~ obs_condition, data = exp_df)
  print(res.ftest)
  
  # print('Bootstrap testing for 2 groups: ')
  # bootstrap_stat <- get_bootstrap_stat_sampling(cis_genes, trans_genes, 
  #                                               sampling_fraction=0.7, nsamples=1000, alternative_theory=alternative_theory)
  return(list(ks_stat=out_stat, ftest_stat=res.ftest, bootstrap_stat=bootstrap_stat))
}
summary_delta_results <- function(){
  save_dir <- '/Users/hoatran/Documents/BCCRC_projects/hakwoo_project/code/results_bulkRNAseq/SA919_full/'
  subtag <- 'Bmet_Apri'
  # subtag <- 'Cmet_Apri'
  # subtag='Cmet_Bpri'
  exp_df <- data.table::fread(paste0(save_dir, subtag, '/', subtag, '_avg_gene_exp.csv.gz'))
  # View(exp_df)
  colnames(exp_df)
  exp_df1 <- exp_df %>%
    # dplyr::mutate(gene_type_desc=paste0(gene_type, '_', obs_condition))%>%
    dplyr::select(gene_type, obs_condition, median_logCPM) %>%
    tidyr::pivot_wider(names_from = 'gene_type', values_from='median_logCPM')
  # ?tidyr::pivot_wider id_cols=c('median_logCPM','gene_type'),
  # View(exp_df1)
  round((exp_df1$cis[1] - exp_df1$cis[2])/exp_df1$cis[2],2) * 100
  
  meta_genes <- data.table::fread(paste0(save_dir, subtag, '/',subtag, '_dispersion_cis_trans_genes.csv'))
  stat_df <- meta_genes %>%
    dplyr::group_by(gene_type) %>%
    dplyr::summarise(mean_dispersion=mean(tagwise_dispersion), median_dispersion=median(tagwise_dispersion),
                     sd_dispersion=sd(tagwise_dispersion), nb_genes=n())
  
  stat_total_df <- meta_genes %>%
    dplyr::mutate(gene_type_desc=
                    case_when(
                      gene_type %in% c('trans','trans_whole_genome') ~ 'trans_gene',
                      TRUE ~ 'cis_gene'
                    )) %>%
    dplyr::group_by(gene_type_desc) %>%
    dplyr::summarise(mean_dispersion=round(mean(tagwise_dispersion),2),
                     median_dispersion=round(median(tagwise_dispersion),2),
                     sd_dispersion=round(sd(tagwise_dispersion),2), nb_genes=n())
  
  # View(stat_total_df)
  # View(stat_df)
  
  save_dir <- '/Users/hoatran/Documents/BCCRC_projects/hakwoo_project/code/results_bulkRNAseq/SA919_full/'
  subtag <- 'Cmet_Apri'
  meta_genes <- data.table::fread(paste0(save_dir, subtag, '/',subtag, '_dispersion_cis_trans_genes.csv'))
  exp_CA <- data.table::fread(paste0(save_dir, subtag, '/', subtag, '_avg_gene_exp.csv.gz'))
  exp_CA$datatag <- subtag
  df <- dplyr::bind_rows(exp_BA, exp_CA)
  
  
  save_dir <- '/Users/hoatran/Documents/BCCRC_projects/hakwoo_project/code/results_bulkRNAseq/SA919_full/'
  subtag='Cmet_Bpri'
  meta_genes <- data.table::fread(paste0(save_dir, subtag, '/',subtag, '_dispersion_cis_trans_genes.csv'))
  
  
}
