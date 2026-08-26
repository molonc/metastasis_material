



suppressPackageStartupMessages({
  library("dplyr")
  library("ggplot2")
})


script_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/'
source(paste0(script_dir, 'scripts/bulk_rna/bulk_utils.R'))
source(paste0(script_dir, 'scripts/bulk_rna/DE_analysis_utils.R'))


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
