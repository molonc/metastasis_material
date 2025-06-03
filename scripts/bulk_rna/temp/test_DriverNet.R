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


get_influence_graph_from_InteractomeFI <- function(obs_genes, save_dir, save_data=T){
  
  input_dir <- '/Users/hoatran/Documents/BCCRC_projects/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/preprocessed_09April2024/drivernet_demo/databases/'
  # graph_df1 <- data.table::fread(paste0(input_dir, 'FIsInGene_Predicted_041709.txt'))
  # graph_df2 <- data.table::fread(paste0(input_dir, 'FIsInGene_Pathway_041709.txt'), header=F)
  # graph_df3 <- data.table::fread(paste0(input_dir, 'FIsInGene_041709.txt'), header=F)
  # colnames(graph_df3)
  # graph_df1 <- graph_df1 %>%
  #   dplyr::mutate(desc=paste0(V1,V2))
  # graph_df2 <- graph_df2 %>%
  #   dplyr::mutate(desc=paste0(V1,V2))
  # graph_df3 <- graph_df3 %>%
  #   dplyr::mutate(desc=paste0(V1,V2))
  # dim(graph_df1)
  # dim(graph_df2)
  # dim(graph_df3)
  # sum(graph_df1$desc %in% graph_df3$desc)
  # sum(graph_df2$desc %in% graph_df3$desc)
  # head(graph_df3)
  
  # edge_df <- data.table::fread(paste0(input_dir, 'FIsInGene_041709.txt'), header=F)
  # edge_df <- edge_df %>%
  #   dplyr::rename(g1=V1, g2=V2) %>%
  #   dplyr::filter(g1 %in% obs_genes & g2 %in% obs_genes)
  # dim(edge_df)
  
  edge_df <- data.table::fread(paste0(input_dir, 'FIsInGene_070323_with_annotations_2022.txt'))
  # summary(edge_df$Score)
  # dim(edge_df)
  # View(head(edge_df))
  edge_df <- edge_df %>%
    dplyr::filter(Gene1 %in% obs_genes & Gene2 %in% obs_genes) %>%
    dplyr::rename(g1=Gene1, g2=Gene2)
  print(dim(edge_df))
  if(save_data){
    data.table::fwrite(edges_ig, paste0(save_dir, 'InteractomeFI_2022_wt_Score.csv.gz'))  
  }
  
  # head(edges_ig)
  # edges_ig <- test_df
  influenceGraph <- matrix(0, nrow=length(obs_genes), ncol=length(obs_genes))
  rownames(influenceGraph) <- obs_genes
  colnames(influenceGraph) <- obs_genes
  for(i in seq(1:dim(edges_ig)[1])){
    influenceGraph[edges_ig[i,]$g1,edges_ig[i,]$g2] <- 1
    influenceGraph[edges_ig[i,]$g2,edges_ig[i,]$g1] <- 1
  }
  # total_genes_ig <- union(edges_ig$g1, edges_ig$g2)
  for(i in seq(1:length(obs_genes))){
    influenceGraph[obs_genes[i],obs_genes[i]] <- 1
  }
  # head(influenceGraph[1:20,1:10])
  print(dim(influenceGraph))
  if(save_data){
    data.table::fwrite(influenceGraph, paste0(save_dir, 'influenceGraph.csv.gz'), row.names = T, col.names = T)
  }
  
  return(influenceGraph)
  
}


# influence graph: using pathway ref to create gene-gene connection or using existing inf graph
# binary matrix, 1: connected, 0: not connected. Same gene - connected 
# samplePatientMutationMatrix: using cnv gene mapping matrix, 
# test with only cis genes + genes in pathway first to see if algo works? 
# how many cis genes are in pathway sets?
# binary form: if change in CN between C versus B, true else false, check again binary format 
# samplePatientOutlierMatrix: 2 samples first, if up-reg true, else false

get_influence_graph <- function(obs_genes, pathway_fn=NULL){
  ## based on hallmark set or interactome database, testing hallmark first
  if(is.null(pathway_fn)){
    # pathway_fn = '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/pathway_set/c2.cp.kegg.v7.1.symbols.gmt'  
    pathway_fn = '/Users/hoatran/Documents/BCCRC_projects/hakwoo_project/testing/DriverNet_testing/h.all.v7.0.symbols.gmt'  
  }
  ref_set <- fgsea::gmtPathways(pathway_fn)
  dim(ref_set)
  total_genes <- c()
  rs_ls <- names(ref_set)
  for(rs in rs_ls){
    tmp_genes <- ref_set[[rs]]
    total_genes <- c(total_genes, tmp_genes)
  }
  total_genes <- intersect(unique(total_genes), obs_genes)
  tg <- length(total_genes)
  adj_mtx <- matrix(0, nrow=tg, ncol=tg)
  rownames(adj_mtx) <- total_genes
  colnames(adj_mtx) <- total_genes
  edge_df <- tibble::tibble()
  for(rs in rs_ls){
    tmp_genes <- ref_set[[rs]]
    tmp_genes <- intersect(tmp_genes, obs_genes)
    g1 <- tmp_genes
    g2 <- tmp_genes
    df <- tidyr::crossing(g1, g2)
    edge_df <- dplyr::bind_rows(edge_df, df)
  }
  # dim(edge_df)
  # library(igraph)
  # el=matrix(c('a','b','c','d','a','d','a','b','c','d'),ncol=2,byrow=TRUE) #a sample edgelist
  # g=graph.data.frame(el)
  # get.adjacency(g,sparse=FALSE)
  # grp=graph.data.frame(as.matrix(edge_df))
  # adj_mtx <- get.adjacency(g,sparse=FALSE)
  # 
  return(edge_df)
}

# patients x genes, values 0, 1
# with row names and col names are pt and genes
get_mutation_data <- function(obs_genes, obs_clones){
  script_dir <- "/Users/hoatran/Documents/BCCRC_projects/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/"
  cnv <- data.table::fread(paste0(script_dir, 'SA919_DE_analysis_DESeq2_Hoa_09April2024/mapping_gene_cnv_SA919.csv.gz'))  
  dim(cnv)
  cnv <- cnv[!duplicated(cnv$ensembl_gene_id),]
  # head(cnv)
  
  ## need to generate a cnv profile for each patient
  rownames(cnv) <- cnv$ensembl_gene_id
  
  # length(rv)
  
  cnv <- as.data.frame(cnv)
  cols_use <- c('ensembl_gene_id','gene_symbol','chr')
  cols_use <- c(cols_use, obs_clones)
  rownames(cnv) <- cnv$ensembl_gene_id
  cnv <- cnv %>%
    dplyr::select(all_of(cols_use))
  # genes_used <- cnv$ensembl_gene_id[rv>0]
  # length(genes_used)
  cnv <- cnv %>%
    dplyr::filter(gene_symbol %in% obs_genes)
  dim(cnv)
  cnv_vals <- cnv %>%
    dplyr::select(all_of(obs_clones))
  rv <- rowVars(as.matrix(cnv_vals))  # median copy number profile of obs clones
  # rv <- rowVars(as.matrix(cnv[,c(5, 6)]))  # median copy number profile of clone B, C 
  genes_used <- cnv$gene_symbol[rv>0]
  
  # sum(genes_used %in% obs_genes)
  # length(obs_genes)
  mutation_mtx <- matrix(0, nrow=2, ncol=length(obs_genes))
  dim(mutation_mtx)
  rownames(mutation_mtx) <- sids
  colnames(mutation_mtx) <- obs_genes
  # mutation_mtx['cloneC',cnv_var_genes] <- 1
  mutation_mtx[obs_clones[2],genes_used] <- 1
  print(dim(mutation_mtx))
  
  return(mutation_mtx)
}

get_expression_data <- function(){
  data_dir <- '/Users/hoatran/Documents/BCCRC_projects/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/preprocessed_09April2024/SA919_DE_analysis_DESeq2_Hoa_09April2024/SA919_CMetMainMix_Bpri_DESeq2/'
  de_df <- data.table::fread(paste0(data_dir, 'DE_signif_genes.csv.gz'))  
  dim(de_df)
  # View(head(de_df))
  ## To Do: get true, false for each sample, and gene id
  ## get list of common genes
  # exp_mtx <- de_df
  dim(exp_df)
  exp_df <- exp_df %>%
    dplyr::filter(symbol %in% obs_genes)
  dim(exp_df)
  obs_genes <- intersect(obs_genes, exp_df$symbol)
  length(unique(obs_genes))
  # View(head(exp_df))
  
  exp_mtx <- matrix(FALSE, nrow=2, ncol=length(obs_genes))
  dim(exp_mtx)
  rownames(exp_mtx) <- sids
  colnames(exp_mtx) <- obs_genes
  head(exp_mtx[,1:4])
  sum(exp_df$symbol %in% obs_genes)
  outlier_genes_B <- exp_df %>%
    dplyr::filter(log2FoldChange<0 & symbol %in% obs_genes) %>%
    dplyr::pull(symbol)
  
  outlier_genes_C <- exp_df %>%
    dplyr::filter(log2FoldChange>0 & symbol %in% obs_genes) %>%
    dplyr::pull(symbol)
  length(outlier_genes_C)
  exp_mtx['cloneC', outlier_genes_C] <- TRUE
  dim(exp_mtx)
  res <- list(exp_df=de_df, exp_mtx=exp_mtx, obs_genes=obs_genes)
  return(res)
}
get_obs_genes <- function(){
  input_dir <- "/Users/hoatran/Documents/BCCRC_projects/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/preprocessed_09April2024/"
  meta_genes <- data.table::fread(paste0(input_dir, 'dispersion_cis_trans_genes_cloneCB.csv'))
  dim(meta_genes)
  summary(as.factor(meta_genes$gene_type))
  
  script_dir <- "/Users/hoatran/Documents/BCCRC_projects/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/"
  cnv <- data.table::fread(paste0(script_dir, 'SA919_DE_analysis_DESeq2_Hoa_09April2024/mapping_gene_cnv_SA919.csv.gz'))  
  dim(cnv)
  cnv <- cnv[!duplicated(cnv$ensembl_gene_id),]
  obs_genes <- intersect(meta_genes$symbol, cnv$gene_symbol)
  length(unique(obs_genes))
  return(obs_genes)  
}
main <- function(){
  
  obs_genes <- get_obs_genes()
  # sids <- c('cloneB','cloneC')
  obs_clones <- c('B','C')
  sids <- paste0('clone', obs_clones)
  res <- get_expression_data(obs_genes)
  obs_genes <- res$obs_genes
  length(obs_genes)
  # edges_ig1 <- get_influence_graph(obs_genes, pathway_fn=NULL) ## version 1
  # dim(edges_ig1)
  save_dir <- paste0(input_dir, 'drivernet_demo/')
  # dir.create(save_dir)
  res_gr <- get_influence_graph_from_InteractomeFI(obs_genes, save_dir, save_data = T) ## version 2
  influenceGraph <- res_gr$influenceGraph
  
  
  mutation_mtx <- get_mutation_data(obs_genes, obs_clones)
  # mutation_mtx: creating a matrix, with var genes, samples have values 1
    
  # saveRDS(obs_genes, paste0(save_dir, 'obs_genes.rds'))
  # saveRDS(mutation_mtx, paste0(save_dir, 'mutation_mtx.rds'))
  # saveRDS(exp_mtx, paste0(save_dir, 'exp_mtx.rds'))
  # saveRDS(influenceGraph, paste0(save_dir, 'influenceGraph.rds'))
  
  # influenceGraph <- readRDS(paste0(save_dir, 'influenceGraph.rds'))
  dim(mutation_mtx)
  dim(exp_mtx)
  dim(influenceGraph)
  driversList = DriverNet::computeDrivers(mutation_mtx, 
                                          exp_mtx,
                                          influenceGraph,
                                          outputFolder=NULL, printToConsole=FALSE)
  
  
  class(driversList)
  driversList@drivers

  
  
    # dim(mutation_mtx)
  # head(mutation_mtx[1:2,1:5])
  # head(exp_mtx[1:2,1:5])
  # head(influenceGraph[1:2,1:2])
  # dim(exp_mtx)
  # dim(influenceGraph)
  # sum(rownames(influenceGraph)==colnames(exp_mtx))
  # sum(rownames(influenceGraph)==colnames(mutation_mtx))
  # comparing cis versus trans in term of connection/ weight
  
  sampleGenes <- unique(meta_genes$symbol)
  length(sampleGenes)
  length(obs_genes)
  sum(obs_genes %in% sampleGenes)
  dim(mutation_mtx)
  dim(exp_mtx)
  dim(influenceGraph)
  
  ## Main version from package
  # randomDriversResult = DriverNet::computeRandomizedResult(
  #   patMutMatrix=mutation_mtx,
  #   patOutMatrix=exp_mtx, 
  #   influenceGraph=influenceGraph,
  #   geneNameList=sampleGenes, outputFolder=NULL, 
  #   printToConsole=F,numberOfRandomTests=100, weight=FALSE, 
  #   purturbGraph=FALSE, purturbData=TRUE)
  length(sampleGenes)
  ## Correct errors version
  randomDriversResult = computeRandomizedResultV2(
    patMutMatrix=mutation_mtx,
    patOutMatrix=exp_mtx, 
    influenceGraph=influenceGraph,
    geneNameList=sampleGenes, outputFolder=NULL, 
    printToConsole=F,numberOfRandomTests=1050, weight=FALSE, 
    purturbGraph=FALSE, purturbData=TRUE)
  
  randomDriversResult <- randomDriversResult[!sapply(randomDriversResult,is.null)]
  
  
  # Summarize the results
  res = resultSummary(driversList, randomDriversResult, 
                      mutation_mtx, influenceGraph,
                      outputFolder=NULL, printToConsole=FALSE)
  res <- as.data.frame(res)
  
  total_df <- print_output(driversList, res, save_dir)
  View(total_df)
  
  total_df1 <- total_df %>%
    dplyr::filter(`p-value`<=0.05)
  total_df2 <- total_df %>%
    dplyr::filter(`p-value`>0.05)
  # View(total_df1)
  # View(total_df2)
  paste(total_df1$gene_symbol, collapse = ', ')
  saveRDS(sampleGenes, paste0(save_dir, 'sampleGenes.rds'))
  saveRDS(randomDriversResult, paste0(save_dir, 'randomDriversResult_1050iterations.rds'))
  saveRDS(res, paste0(save_dir, 'res.rds'))
  saveRDS(driversList, paste0(save_dir, 'driversList.rds'))
  
  save.image("~/Documents/BCCRC_projects/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/preprocessed_09April2024/drivernet_demo/results_5May_50iterations.RData")
  
}

check_associated_pathway <- function(total_df){
  ## get list of significant pathways
  ## check if driver genes are in pathway
  ## check if cover genes are in pathway
  gmt_fn <- '/Users/hoatran/Documents/BCCRC_projects/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/pathway_set/h.all.v7.0.symbols.gmt'
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
          print('Gene is INCLUDED ')
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

print_output <- function(driversList, res, save_dir){
  # res$`p-value`
  # signif_genes <- res %>%
  #   dplyr::filter(`p-value`<=0.05) %>%
  #   dplyr::pull(gene)
  
  total_df <- tibble::tibble()
  
  for(g in res$gene){
    print(g)
    tmp <- as.data.frame(driversList@actualEvents[[g]])
    connect_genes <- tmp$actual_events_gen[tmp$actual_events_gen != g]
    temp <- res %>%
      dplyr::filter(gene==g) %>%
      dplyr::mutate(cover_genes=paste(connect_genes, collapse=','))
    # View(temp)
    
    # 
    # df <- tibble::tibble(gene_symbol=g, rank=rank, nb_connected_genes=dim(tmp)[1],
    #                      connected_genes=paste(tmp$actual_events_gen, collapse=','))
    total_df <- dplyr::bind_rows(total_df, temp)
    # $actual_events_gen
  }
  total_df <- total_df %>%
    dplyr::mutate(`p-value`=round(as.numeric(`p-value`),3)) %>%
    dplyr::select(gene, rank, `p-value`, covered_events, cover_genes, everything()) %>%
    dplyr::rename(gene_symbol=gene)
  # View(total_df)
  
  obs_chrs <- c(as.character(1:22),'X')
  ref_df1 <- annotables::grch38 %>%
    dplyr::select(ensgene, symbol, chr, description) %>%
    dplyr::rename(gene_symbol=symbol, ens_gene_id=ensgene) %>%
    dplyr::filter(gene_symbol %in% total_df$gene_symbol & chr %in% obs_chrs) 
  # as.data.frame()
  dim(ref_df1)
  
  # View(total_df)
  total_df <- total_df %>%
    dplyr::left_join(ref_df1, by='gene_symbol') %>%
    dplyr::select(gene_symbol, rank, `p-value`, covered_events, cover_genes, chr, description, everything()) 
  
  data.table::fwrite(total_df, paste0(save_dir, 'significant_genes_DriverNet_cloneCB_v2.csv'))
  
  return(total_df)
}
get_weight_per_gene <- function(){
  cis_genes <- cnv_var_genes
  trans_genes <- obs_genes[!obs_genes %in% cis_genes]
  length(cis_genes)
  length(trans_genes)
  length(obs_genes)
  head(edges_ig)
  gene_weight_ls <- list()
  for(g in obs_genes){
    edges_tmp <- edges_ig %>%
      dplyr::filter(g1==g)
    dim(edges_tmp)
    gene_weight_ls[[g]] <- c(gene_symbol=g, weight=dim(edges_tmp)[1])
  }
  gene_weight_df <- as.data.frame(gene_weight_ls)
  dim(gene_weight_df)
  gene_weight_df <- as.data.frame(t(gene_weight_df))
  head(gene_weight_df)
  
  # cis <- gene_weight_df %>%
  #   dplyr::filter(gene_symbol %in% cis_genes) %>%
  #   dplyr::mutate(weight=as.numeric(weight))%>%
  #   dplyr::summarise(avg_weight=mean(weight), sd_weight=sd(weight))
  # trans_weight <- gene_weight_df %>%
  #   dplyr::filter(gene_symbol %in% trans_genes) %>%
  #   dplyr::mutate(weight=as.numeric(weight))%>%
  #   dplyr::summarise(avg_weight=mean(weight), sd_weight=sd(weight))
  # 
  gene_weight_df <- gene_weight_df %>%
    dplyr::mutate(gene_type=
                    case_when(gene_symbol %in% cis_genes ~ 'cis',
                              TRUE ~ 'trans'),
                  weight=as.numeric(weight))
  stat_df <- gene_weight_df %>%
    group_by(gene_type) %>%
    dplyr::summarise(avg_wd=mean(weight), #median_wd=median(weight), 
                     sd_wd = sd(weight))
  stat_df
  
  stat_df1 <- gene_weight_df %>%
    dplyr::filter(weight>0) %>%
    group_by(gene_type) %>%
    dplyr::summarise(avg_wd=mean(weight), #median_wd=median(weight), 
                     sd_wd = sd(weight))
  stat_df1
  library("ggplot2")
  p <- ggplot(gene_weight_df, aes(x=gene_type, y=weight)) + #
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(size=0.2, aes(colour = gene_type), position=position_jitter(0.2))  + 
    theme_bw() + 
    labs(x='Gene type', y='Gene connected weight')
  p
  png(paste0(save_dir,"cis_trans_connected_weight_cloneBC.png"), height = 2*300, width=2*350, res = 2*72)
  print(p)
  dev.off()
  
  cis_genes <- gene_weight_df %>%
    dplyr::filter(gene_type=='cis' & weight>0) %>%
    dplyr::pull(weight)
  trans_genes <- gene_weight_df %>%
    dplyr::filter(gene_type=='trans' & weight>0) %>%
    dplyr::pull(weight)
  summary_stat_df <- get_bootstrap_stat_sampling(cis_genes, trans_genes, 
                                                 sampling_fraction=0.9, nsamples=100, alternative_theory='greater')
  
}





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
