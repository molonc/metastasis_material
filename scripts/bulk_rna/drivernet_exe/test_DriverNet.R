# BiocManager::install("DriverNet")
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
  library(DriverNet)
})  

script_dir <- "/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/scripts/bulk_rna/drivernet_exe/"
source(paste0(script_dir, 'drivernet_utils.R'))
# save_dir <- save_output_dir

main <- function(){
  
  input_dir <- "/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/preprocessed_09April2024/"
  save_dir <- paste0(input_dir, 'drivernet_demo/')
  # obs_clones <- c('B','C')
  # obs_clones <- c('A','B')
  obs_clones <- c('A','C')
  
  
  obs_genes <- get_obs_genes(obs_clones)
  sampleGenes <- obs_genes
  length(sampleGenes)
  length(obs_genes)
  
  subtag <- paste0(obs_clones[2],'met_',obs_clones[1],'pri')
  save_output_dir <- paste0(save_dir, subtag,'/')
  # dir.create(save_output_dir)
  
  res <- get_expression_data(obs_clones, obs_genes, subtag, save_output_dir)
  exp_mtx <- res$exp_mtx
  obs_genes <- res$obs_genes
  length(obs_genes)

  # edges_ig1 <- get_influence_graph(obs_genes, pathway_fn=NULL) ## version 1
  # dim(edges_ig1)
  
  
  df <- tibble::tibble(sampleGenes=sampleGenes)
  data.table::fwrite(df, paste0(save_dir, 'sampleGenes.csv.gz'))
  
  res_gr <- get_influence_graph_from_InteractomeFI(obs_genes, save_output_dir, save_data = T) ## version 2
  influenceGraph <- res_gr$influenceGraph
  
  mutation_mtx <- get_mutation_data(obs_genes, obs_clones, save_output_dir)
  # mutation_mtx: creating a matrix, with var genes, samples have values 1
    
  # saveRDS(obs_genes, paste0(save_dir, 'obs_genes.rds'))
  # saveRDS(mutation_mtx, paste0(save_dir, 'mutation_mtx.rds'))
  # saveRDS(exp_mtx, paste0(save_dir, 'exp_mtx.rds'))
  # saveRDS(influenceGraph, paste0(save_dir, 'influenceGraph.rds'))
  
  # influenceGraph <- readRDS(paste0(save_dir, 'influenceGraph.rds'))
      # dim(mutation_mtx)
  # head(mutation_mtx[1:2,1:5])
  # head(exp_mtx[1:2,1:5])
  # head(influenceGraph[1:2,1:2])
  # dim(exp_mtx)
  # dim(influenceGraph)
  # sum(rownames(influenceGraph)==colnames(exp_mtx))
  # sum(rownames(influenceGraph)==colnames(mutation_mtx))
  # comparing cis versus trans in term of connection/ weight
  
  
  length(obs_genes)
  sum(obs_genes %in% sampleGenes)
  dim(mutation_mtx)
  dim(exp_mtx)
  dim(influenceGraph)
  length(sampleGenes)
  run_DriverNet(mutation_mtx, exp_mtx, influenceGraph, sampleGenes, save_output_dir)
    
}

check_associated_pathway <- function(total_df){
  ## get list of significant pathways
  ## check if driver genes are in pathway
  ## check if cover genes are in pathway
  gmt_fn <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/pathway_set/h.all.v7.0.symbols.gmt'
  ref_set <- fgsea::gmtPathways(gmt_fn)
  
  signif_pw_ls <- c('MYOGENESIS','EPITHELIAL_MESENCHYMAL_TRANSITION', 
                    'APICAL_JUNCTION','MITOTIC_SPINDLE')
  prefix <- 'HALLMARK_'
  
  class(ref_set)
  ref_set <- ref_set[paste0(prefix, signif_pw_ls)]
  for(g in total_df$gene_symbol){
    print('-------------------')
    print(g)
    tmp <- total_df %>%
      dplyr::filter(gene_symbol==g)
    cover_genes <- unlist(strsplit(total_df$cover_genes, split=','))
    for(rf in names(ref_set)){
      cover_included <- cover_genes[cover_genes %in% ref_set[[rf]]]
      st1 <- g %in% ref_set[[rf]]
      st2 <- length(cover_included)>0
      if(st1 || st2){
        print('         ')
        print(rf)  
        if(st1){
          print('Gene is included ')
        }else{
          print('Gene is NOT included')
        }
        
        if(st2){
          print('Covered gene is included')
          print(paste(cover_included, collapse = ', '))
        }
      }
      
      
    }  
  }
  
  
  
}

summary_driver_genes_output <- function(){
  input_dir <- "/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/preprocessed_09April2024/"
  save_dir <- paste0(input_dir, 'drivernet_output_metastasis_proj/')
  # obs_clones <- c('B','C')
  # obs_clones <- c('A','B')
  obs_clones <- c('A','C')
  subtag <- paste0(obs_clones[2],'met_',obs_clones[1],'pri')
  ac_df <- data.table::fread(paste0(save_dir, subtag,'/', 'significant_genes_DriverNet.csv'))
  dim(ac_df)
  # View(ac_df)
  ac_df <- get_cis_trans_gene_type(obs_clones, ac_df)
  dim(ac_df)
  data.table::fwrite(ac_df, paste0(save_dir, subtag,'/', subtag,'_significant_genes_DriverNet.csv'))
  # act <- ac_df %>%
  #   dplyr::filter(`p-value`<0.05)
  
  obs_clones <- c('A','B')
  subtag <- paste0(obs_clones[2],'met_',obs_clones[1],'pri')
  subtag
  ab_df <- data.table::fread(paste0(save_dir, subtag,'/', 'significant_genes_DriverNet.csv'))
  dim(ab_df)
  ab_df <- get_cis_trans_gene_type(obs_clones, ab_df)
  dim(ab_df)
  colnames(ab_df)
  ab_df$gene_type
  data.table::fwrite(ab_df, paste0(save_dir, subtag,'/', subtag,'_significant_genes_DriverNet.csv'))
  
  abt <- ab_df %>%
    dplyr::filter(`p-value`<0.05)
  dim(abt)
  sum(ac_df$gene_symbol %in% ab_df$gene_symbol)  
  similar_genes <- intersect(ac_df$gene_symbol, ab_df$gene_symbol)
  ac1 <- ac_df %>%
    dplyr::filter(gene_symbol %in% similar_genes)
  ab1 <- ab_df %>%
    dplyr::filter(gene_symbol %in% similar_genes)
  ac1$chr  
  ab1$chr  
  
  obs_clones <- c('B','C')
  subtag <- paste0(obs_clones[2],'met_',obs_clones[1],'pri')
  bc_df <- data.table::fread(paste0(save_dir, subtag,'/', 'significant_genes_DriverNet.csv'))
  dim(bc_df)
  bc_df <- get_cis_trans_gene_type(obs_clones, bc_df)
  dim(bc_df)
  colnames(bc_df)
  bc_df$gene_type
  data.table::fwrite(bc_df, paste0(save_dir, subtag,'/', subtag,'_significant_genes_DriverNet.csv'))
  
  
  
  bct <- bc_df %>%
    dplyr::filter(`p-value`<0.05)
  dim(bct)
  sum(ac_df$gene_symbol %in% bc_df$gene_symbol)  
  similar_genes <- intersect(ac_df$gene_symbol, bc_df$gene_symbol)
  ac2 <- ac_df %>%
    dplyr::filter(gene_symbol %in% similar_genes)
  bc1 <- bc_df %>%
    dplyr::filter(gene_symbol %in% similar_genes)
  ac2$chr
  bc1$chr  
}

z








# View(head(exp_df))
# exp_mtx <- matrix(FALSE, nrow=2, ncol=dim(exp_df)[1])
# rownames(exp_mtx) <- sids
# colnames(exp_mtx) <- exp_df$symbol
# head(exp_mtx[,1:4])
# exp_mtx['cloneC',] <- TRUE
# dim(exp_mtx)
# 
# 
# edges_ig <- get_influence_graph(exp_df$symbol, pathway_fn=NULL)
# dim(edges_ig)
# head(edges_ig)
# influenceGraph <- matrix(0, nrow=length(obs_genes), ncol=length(obs_genes))
# rownames(influenceGraph) <- obs_genes
# colnames(influenceGraph) <- obs_genes
# for(i in seq(1:dim(edges_ig)[1])){
#   influenceGraph[edges_ig[i,]$g1,edges_ig[i,]$g2] <- 1
# }
# 
# cnv_var_genes <- get_mutation_data(obs_genes)
# # mutation_mtx: creating a matrix, with var genes, samples have values 1
# length(cnv_var_genes)
# cnv_var_genes[1:10]
# 
# sum(cnv_var_genes %in% obs_genes)
# mutation_mtx <- matrix(0, nrow=2, ncol=dim(exp_df)[1])
# rownames(mutation_mtx) <- sids
# colnames(mutation_mtx) <- exp_df$symbol
# mutation_mtx['cloneC',cnv_var_genes] <- 1
# dim(mutation_mtx)
# 
# 
