



suppressPackageStartupMessages({
  library("dplyr")
  library("ggplot2")
})


script_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/'
source(paste0(script_dir, 'scripts/bulk_rna/bulk_utils.R'))
source(paste0(script_dir, 'scripts/bulk_rna/DE_analysis_utils.R'))


get_stat_per_chr <- function(exp_df, obs_conds){
  datatag <- paste(obs_conds, collapse='_vs_')
  print(datatag)
  # unique(exp_df$gene_type)
  cds <- unique(exp_df$obs_conds)
  cd_met <- cds[grepl('Met',cds)]
  cd_pri <- cds[grepl('Pri',cds)]
  res <- tibble::tibble()
  
  for(obs_chr in unique(exp_df$chr)){
    print(obs_chr)
    exp_df1 <- exp_df %>%
      dplyr::filter(gene_type=='cis' & chr==obs_chr)
    tmp1 <- exp_df %>%
      dplyr::filter(gene_type=='cis' & obs_conds==cd_met & chr==obs_chr) #%>%
      # dplyr::mutate(gene_type_desc=paste0(gene_type, '_', obs_condition))%>%
      # dplyr::pull(gene_exp)
    tmp2 <- exp_df %>%
      dplyr::filter(gene_type=='cis' & obs_conds==cd_pri & chr==obs_chr) #%>%
      # dplyr::mutate(gene_type_desc=paste0(gene_type, '_', obs_condition))%>%
      # dplyr::pull(gene_exp)
    dim(tmp1)
    head(tmp1)
    dim(tmp2)
    exp_df$ens_gene_id
    tmp <- exp_df %>% 
      dplyr::group_by(ens_gene_id, chr, obs_conds) %>% 
      summarize(median=median(gene_exp)
                                                             # Lower=quantile(gene_exp, probs=0.25),
                                                             # Upper=quantile(gene_exp, probs=0.75)
    ) %>% 
      tidyr::pivot_longer(cols = median, values_to = "median_gene_exp")
    dim(tmp)
    head(tmp)
    colnames(tmp)
    unique(tmp$chr)
    tmp$name <- NULL
    
    # tmp <- tmp %>% 
    #   dplyr::mutate(desc = paste0(chr, '_', ens_gene_id)) 
    # tmp$ens_gene_id <- NULL
    tmp <- tmp %>% 
      # dplyr::select(desc, obs_conds, median_gene_exp) %>% 
      tidyr::pivot_wider(names_from = obs_conds, 
                         values_from = median_gene_exp)
    tmp <- tmp %>% 
      dplyr::rename(group2 = !!sym(cd_met), group1 = !!sym(cd_pri))
    100*sum(tmp$group2 - tmp$group1>0)/dim(tmp)[1]
    obs_chr <- '7'
    tmp <- tmp %>% 
      dplyr::filter(chr %in% obs_chr)
    data.table::fwrite(tmp, paste0(save_fig_dir, 'median_gene_exp.csv'))
    
    # , '_',obs_conds
    print('KS-test variance testing for 2 groups: ')
    # alternative_theory <- 'greater'
    # out_stat <- ks.test(g1, g2, alternative=alternative_theory)
    out_stat <- ks.test(g1, g2)
    print(out_stat)
    length(g2)
    length(g1)
    100*sum(g1-g2>0)/length(g1)
    dim(exp_df)
    res1 <- tibble::tibble(obs_chr=obs_chr, nb_genes=length(g1), 
                           p_value=out_stat$p.value,
                           D_value=as.numeric(out_stat$statistic),
                           hypothesis=out_stat$alternative)
    res <- dplyr::bind_rows(res, res1)
    # print('F-test variance testing for 2 groups: ')
    # res.ftest <- var.test(gene_exp ~ obs_conds, data = exp_df1)
    # print(res.ftest)
    # 
    # print('Bootstrap testing for 2 groups: ')
    # alternative_theory <- 'greater'
    # bootstrap_stat <- get_bootstrap_stat_sampling(g1, g2,
    #                                               sampling_fraction=0.7, nsamples=1000, alternative_theory=alternative_theory)
    
  }
  res$obs_condition <- datatag
  return(res)
}
permutation_test_cis_genes <- function(nsamples=2){
  obs_genes <- tmp$ens_gene_id
  resamples <- lapply(1:nsamples, function(i) sample(obs_genes, size=length(obs_genes), replace = T))
  length(resamples)
  resamples[[1]]
  our <- sum(our_obs %in% ref_set)
  # occurs <- sapply(resamples, count_occurance)
  occurs <- sapply(resamples, function(s) {return(sum(s %in% ref_set))})
  names(occurs) <- paste0('R',seq(1:length(occurs)))
  occurs['our'] <- our
  # length(occurs)
  r <- rank(occurs) #,ties.method = "max"
  pval <- (nsamples + 1 - r['our'])/nsamples
  return(list(CI=r['our'],pval=pval))
}
viz_Fig6_partB <- function(){
  save_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/results_bulkRNAseq/SA919_full/'
  datatag <- 'SA919_full'
  
  obs_conds <- c('B_Metastasis','A_Primary')
  subtag <- 'Bmet_Apri'
  save_fig_dir <- paste0(save_dir, subtag, '/')
  obs_clones <- c('A','B')
  df <- data.table::fread(paste0(save_fig_dir, paste(obs_conds, collapse='_vs_'), '_total_exp.csv.gz'))
  dim(df)
  obs_genes <- unique(df$ens_gene_id)
  length(obs_genes)
  unique(df$chr)
  # chr 7: 
  stat <- get_median_cnv_cis_genes(unique(df$ens_gene_id), obs_clones)
  stat
  p10 <- viz_cnv_profiles_specific_chrs_horizontal(stat)
  p11 <- viz_gene_exp_comparison(df, obs_clones, legend_pos='bottom')
  # saveRDS(p1, paste0(save_fig_dir, subtag, '_gene_exp_plt_with_legend.rds'))
  
  
  obs_conds <- c('C_Metastasis','A_Primary')
  subtag <- 'Cmet_Apri'
  save_fig_dir <- paste0(save_dir, subtag, '/')
  obs_clones <- c('A','C')
  df <- data.table::fread(paste0(save_fig_dir, paste(obs_conds, collapse='_vs_'), '_total_exp.csv.gz'))
  dim(df)
  stat <- get_median_cnv_cis_genes(unique(df$ens_gene_id), obs_clones)
  stat
  p20 <- viz_cnv_profiles_specific_chrs_horizontal(stat)
  p21 <- viz_gene_exp_comparison(df, obs_clones, legend_pos='bottom')
  # saveRDS(p2, paste0(save_fig_dir, subtag, '_gene_exp_plt_with_legend.rds'))
  # p2
  
  obs_conds <- c('C_Metastasis','B_Primary')
  subtag <- 'Cmet_Bpri'
  save_fig_dir <- paste0(save_dir, subtag, '/')
  obs_clones <- c('B','C')
  df <- data.table::fread(paste0(save_fig_dir, paste(obs_conds, collapse='_vs_'), '_total_exp.csv.gz'))
  dim(df)
  stat <- get_median_cnv_cis_genes(unique(df$ens_gene_id), obs_clones)
  stat
  p30 <- viz_cnv_profiles_specific_chrs_horizontal(stat)
  p31 <- viz_gene_exp_comparison(df, obs_clones, legend_pos='bottom')
  # saveRDS(p3, paste0(save_fig_dir, subtag, '_gene_exp_plt_with_legend.rds'))
  # p3
  
  # xpos <- 0.2
  # ypos <- 0.9
  # + theme(legend.position = c(xpos, ypos))
  p13 <- cowplot::plot_grid(p10 + theme(legend.position = 'top'),
                            p11,
                            ncol=1, rel_heights = c(0.7,3))
  p23 <- cowplot::plot_grid(p20 + theme(legend.position = 'top'),
                            p21,
                            ncol=1, rel_heights = c(0.7,3))
  
  p33 <- cowplot::plot_grid(p30 + theme(legend.position = 'top'),
                            p31,
                            ncol=1, rel_heights = c(0.7,3))
  
  # p13
  
  p_total_partB <- cowplot::plot_grid(p13, 
                                      p23,
                                      p33,
                                      nrow=1,
                                      rel_widths = c(1.8,3.1,1.8))
  p_total_partB
  saveRDS(p_total_partB, paste0(save_dir,"figs/","gene_exp_Fig6_partB.rds"))
  p_total_partB <- readRDS(paste0(save_dir,"figs/","gene_exp_Fig6_partB.rds"))
  png(paste0(save_dir,"figs/","gene_exp_Fig6_partB.png"), 
      height = 2*200, width=2*650, res = 2*72)
  print(p_total_partB)
  dev.off() 
  
  
}
get_stat_test <- function(){
  
  # Get list of bulk RNA-seq genes: 19K genes: genome genes
  script_dir <- "/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/"
  meta_genes <- data.table::fread(paste0(script_dir, 'SA919_DE_analysis_DESeq2_Hoa_09April2024/meta_genes_SA919.csv.gz'))  
  dim(meta_genes) #19588 genes
  head(meta_genes)
  genome_genes <- meta_genes$ens_gene_id
  
  # Get # cis, trans genes for each DE comparison
  save_dir <- "/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/results_bulkRNAseq/SA919_full/"
  # save_fig_dir <- paste0(save_dir,'revision/')
  save_fig_dir <- paste0(save_dir,'revision_v2/')
  dir.create(save_fig_dir)
  de_comps <- data.table::fread(paste0(save_dir,'list_DE_comparisons.csv'))
  de_comps <- de_comps %>%
    dplyr::filter(exp_cond1=='main_exp' & exp_cond2=='main_exp')
  de_comps
  
  obs_chrs <- c(as.character(1:22), "X")
  
  de <- de_comps$desc[1]
  
  obs_gene_type <- 'cis_gene'
  for(de in de_comps$desc){
    de_comps_tmp <- de_comps %>%
      dplyr::filter(desc==de)
    print('---------------------------')
    print(de)
    de_comps_tmp
    de_fn <- paste0(save_dir, de,'/',de,'_DE_genes_gene_type.csv.gz')
    if(file.exists(de_fn)){
      df_full <- data.table::fread(de_fn)
      df_full <- df_full %>%
        # dplyr::filter(chr %in% unique(df$chr))
        dplyr::filter(chr %in% obs_chrs)
      df <- df_full %>%
        dplyr::filter(gene_type==obs_gene_type)
      print(dim(df_full))
      print(dim(df))
      
      obs_clones <- c(de_comps_tmp$obs_clone1, de_comps_tmp$obs_clone2)
      if(obs_clones[1] != obs_clones[2]){
        cnv <- get_cis_trans_gene_type(obs_clones)
        dim(cnv)
        head(cnv)
        cnv <- cnv %>%
          dplyr::filter(chr %in% obs_chrs)
        df_full <- df_full %>%
          dplyr::filter(chr %in% unique(cnv$chr))
        sample_size <- dim(df_full)[1]
        # sample_size
        res <- get_bootstrap_stat(df$ensembl_gene_id, cnv$ensembl_gene_id, genome_genes, sample_size)  
        print(paste0('Confidence interval: ',res$CI, ' & p-value:',res$pval))      
        
        
      }
      
    }
  }
  # unique(cnv$chr)
  # dim(total_df)
  # cnv$total_cnv <- cnv$A + cnv$B
  sum(abs(df_full$log2FoldChange)<1)
  summary(as.factor(df_full$chr))
  summary(as.factor(df$chr))
  # Get bootstrap test for number of occurrence
  
}


get_bootstrap_stat <- function(our_obs, ref_set, genome_genes, sample_size, nsamples=10000){
  resamples <- lapply(1:nsamples, function(i) sample(genome_genes, size=sample_size, replace = T))
  
  our <- sum(our_obs %in% ref_set)
  # occurs <- sapply(resamples, count_occurance)
  occurs <- sapply(resamples, function(s) {return(sum(s %in% ref_set))})
  names(occurs) <- paste0('R',seq(1:length(occurs)))
  occurs['our'] <- our
  # length(occurs)
  r <- rank(occurs) #,ties.method = "max"
  pval <- (nsamples + 1 - r['our'])/nsamples
  return(list(CI=r['our'],pval=pval))
}
