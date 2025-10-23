suppressPackageStartupMessages({
  library("dplyr")
  library("ggplot2")
  library("DESeq2")
  library("SingleCellExperiment")
})

script_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/'
source(paste0(script_dir, 'scripts/bulk_rna/bulk_utils.R'))
source(paste0(script_dir, 'scripts/bulk_rna/DE_analysis_utils.R'))
source(paste0(script_dir, 'scripts/bulk_rna/bulk_pathway_utils.R'))

get_pathways_Fig6_partD <- function(){
  # Loading DE genes files 
  # A met vs. A pri
  # B met vs. B pri
  # C met vs. B pri
  save_dir <- "/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/results_bulkRNAseq/SA919_full/"
  aa_pw <- data.table::fread(paste0(save_dir, 'Amet_Apri/pathways_Amet_Apri.csv.gz'))
  bb_pw <- data.table::fread(paste0(save_dir, 'Bmet_Bpri/pathways_Bmet_Bpri.csv.gz'))
  cb_pw <- data.table::fread(paste0(save_dir, 'Cmet_Bpri/pathways_Cmet_Bpri.csv.gz'))
  
  # aa_pw$comp <- 'A met vs. A pri'
  # bb_pw$comp <- 'B met vs. B pri'
  # cb_pw$comp <- 'C met vs. B pri'
  
  aa_pw$comp <- 'Met clone A vs. Pri clone A'
  bb_pw$comp <- 'Met clone B vs. Pri clone B'
  cb_pw$comp <- 'Met clone C vs. Pri clone B'
  
  df <- dplyr::bind_rows(aa_pw, bb_pw, cb_pw)
  # View(df)
  # colnames(df)
  p <- viz_pathway_sets(df)
  # p
  saveRDS(p, paste0(save_dir, 'figs/pathways_Fig6.rds'))
  ggsave(  
    filename = paste0(save_dir,"figs/","pathways_Fig6.svg"),  
    plot = p,  
    height = 3, width = 3, dpi = 150)
  
  
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
  
  obs_clones <- c('A','B')
  res1 <- get_DE_results(dds, obs_clones, save_dir)
  # subtag <- 'Bmet_Apri'
  # save_fig_dir <- paste0(save_dir, subtag, '/')
  # # dir.create(save_fig_dir)
  # meta_df <- colData(dds) %>%
  #   as.data.frame() %>%
  #   dplyr::mutate(desc=paste0(main_clone, '_',mainsite)) %>%
  #   filter(desc %in% c('B_Metastasis','A_Primary'))
  # dim(meta_df)
  # sids <- rownames(meta_df)
  # dds_tmp <- dds[,sids]
  # print(sizeFactors(dds_tmp)) ## using existing size factors calculated using Scran method
  # res1 <- get_DE_genes_DESeq2(dds_tmp, DE_comp=c("Metastasis","Primary"),
  #                             filter_genes=F, min_total_exp_samples=10)
  # dim(res1)
  # res1 <- res1 %>%
  #   as.data.frame() %>%
  #   filter(abs(log2FoldChange)>=1) #pvalue<0.05 & 
  # dim(res1)
  # head(res1)
  # data.table::fwrite(res1, paste0(save_fig_dir, subtag, '_DE_genes.csv.gz'))
  
  obs_clones <- c('A','C')
  res2 <- get_DE_results(dds, obs_clones, save_dir)
  # subtag <- 'Cmet_Apri'
  # save_fig_dir <- paste0(save_dir, subtag, '/')
  # dir.create(save_fig_dir)
  # meta_df <- colData(dds) %>%
  #   as.data.frame() %>%
  #   dplyr::mutate(desc=paste0(main_clone, '_',mainsite)) %>%
  #   filter(desc %in% c('C_Metastasis','A_Primary'))
  # dim(meta_df)
  # sids <- rownames(meta_df)
  # dds_tmp <- dds[,sids]
  # print(sizeFactors(dds_tmp)) ## using existing size factors calculated using Scran method
  # res2 <- get_DE_genes_DESeq2(dds_tmp, DE_comp=c("Metastasis","Primary"),
  #                             filter_genes=F, min_total_exp_samples=10)
  # dim(res2)
  # res2 <- res2 %>%
  #   as.data.frame() %>%
  #   filter(abs(log2FoldChange)>=1) #pvalue<0.05 & 
  # dim(res2)
  # head(res2)
  # data.table::fwrite(res2, paste0(save_fig_dir, subtag, '_DE_genes.csv.gz'))
  
  
}  
extract_DE_genes_DESeq2_cloneA <- function(){
  
  save_dir <- "/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/results_bulkRNAseq/SA919_full/"
  datatag <- 'SA919_full'
  sce <- readRDS(paste0(save_dir, datatag, '_sce.rds'))
  dds <- create_DESeq2_object(sce)
  obs_clones <- c('A','A')
  meta_genes <- rowData(sce) %>% as.data.frame()
  subtag <- paste0(obs_clones[2],'met','_',obs_clones[1],'pri')
  print(subtag)
  save_fig_dir <- paste0(save_dir, subtag, '/')
  dir.create(save_fig_dir)
  
  res <- get_DE_results(subtag, save_fig_dir, dds, meta_genes, obs_clones, exps=c('main_exp'))
  obs_genes_symb <- unique(res$symbol)
  length(obs_genes_symb)
  pw_res <- get_gprofiler_pathways_obsgenes_v2(obs_genes_symb, save_fig_dir, datatag, correction_method='gSCS',
                                            custom_id="gp__RYib_WEFZ_gEw", pathway_fn=NULL, save_data=T)
  sum(pw_res$stat$p_value<0.05)
  
  meta_cells <- colData(sce)
  summary(res$log2FoldChange)
  
  
  obs_clones <- c('B','B')
  meta_genes <- rowData(sce) %>% as.data.frame()
  subtag <- paste0(obs_clones[2],'met','_',obs_clones[1],'pri')
  print(subtag)
  save_fig_dir <- paste0(save_dir, subtag, '/')
  dir.create(save_fig_dir)
  
  res <- get_DE_results(subtag, save_fig_dir, dds, meta_genes, obs_clones, exps=c('main_exp'))
  obs_genes_symb <- unique(res$symbol)
  pw_res <- get_gprofiler_pathways_obsgenes_v2(obs_genes_symb, save_fig_dir, datatag, correction_method='gSCS',
                                            custom_id="gp__RYib_WEFZ_gEw", pathway_fn=NULL, save_data=T)
  sum(pw_res$stat$p_value<0.05)
  
  
  obs_clones <- c('B','C')
  meta_genes <- rowData(sce) %>% as.data.frame()
  subtag <- paste0(obs_clones[2],'met','_',obs_clones[1],'pri')
  print(subtag)
  save_fig_dir <- paste0(save_dir, subtag, '/')
  dir.create(save_fig_dir)
  
  res <- get_DE_results(subtag, save_fig_dir, dds, meta_genes, obs_clones, exps=c('main_exp'))
  obs_genes_symb <- unique(res$symbol)
  pw_res <- get_gprofiler_pathways_obsgenes(obs_genes_symb, save_fig_dir, datatag, correction_method='gSCS',
                                            custom_id="gp__RYib_WEFZ_gEw", pathway_fn=NULL, save_data=T)
  sum(pw_res$stat$p_value<0.05)
  
  
}
main(){
  ## Loading dds object
  extract_DE_genes_DESeq2(dds)
}

# main()