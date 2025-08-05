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

## input, output? 

viz_graph_results <- function(datatag, save_output_dir){
  library(ggraph)
  library(igraph)
  library(tidygraph)
  
  input_dir <- "/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/preprocessed_09April2024/"
  save_dir <- paste0(input_dir, 'drivernet_demo/')
  obs_clones <- c('B','C') ## check % cis genes, a trackplot
  subtag <- paste0(obs_clones[2],'met_',obs_clones[1],'pri')
  save_output_dir <- paste0(save_dir, subtag,'/')
  pw_df <- data.table::fread(paste0(save_output_dir, 'pathways_',subtag,'.csv.gz'))
  # df <- data.table::fread(paste0(save_dir, 'significant_genes_DriverNet_cloneCB_v2.csv'))
  df <- data.table::fread(paste0(save_output_dir, 'significant_genes_DriverNet.csv'))
  # edges_df <- data.table::fread(paste0(save_output_dir, 'InteractomeFI_2022_wt_Score.csv.gz'))  
  
  data_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/results_bulkRNAseq/SA919_full/'
  exp_df <- data.table::fread(paste0(data_dir, subtag, '/', subtag, '_DE_genes.csv.gz'))
  print(subtag)
  print(dim(exp_df))
  head(exp_df)
  # 
  # dim(edges_df)  
  # head(edges_df)
  # dim(df)
  # df[1:5,1:5]
  # names(df)
  # unique(edges_df$Direction)
  
  # dolphin %>% 
  #   activate(nodes) %>%
  #   mutate(community = as.factor(group_louvain())) %>% 
  #   ggraph() + 
  #   geom_edge_link() + 
  #   geom_node_point(aes(colour = community), size = 5)
  
  
  driver_df <- df %>%
    dplyr::filter(`p-value`<0.05) %>%
    # dplyr::select(gene_symbol, ens_gene_id) %>%
    # dplyr::rename(name=gene_symbol) %>%
    dplyr::mutate(gene_type='driver gene')
  nodes_df <- df %>%
    # dplyr::filter(`p-value`<0.05) %>%
    dplyr::select(gene_symbol, ens_gene_id) %>%
    dplyr::rename(name=gene_symbol) %>%
    dplyr::mutate(gene_type='driver gene', label=name)
  
  cover_genes <- c()
  edges_df <- tibble::tibble()
  for(g in df$gene_symbol){
    gls <- df %>%
      dplyr::filter(gene_symbol==g) %>%
      dplyr::pull(cover_genes)
    gls <- unlist(strsplit(gls,','))
    
    if(length(gls)>0){
      print('______________')
      print(g)
      print(gls)
      cover_genes <- c(cover_genes, gls)
      ed_tmp <- expand.grid(c(g), gls) %>%
        dplyr::rename(x=Var1, y=Var2) %>%
        dplyr::mutate(edge_type='driver_cover_connected')
      edges_df <- dplyr::bind_rows(edges_df, ed_tmp)
    }
  }
  cover_genes <- unique(cover_genes)
  cover_genes <- cover_genes[!cover_genes %in% nodes_df$name]
  nodes_df2 <- tibble::tibble(name=cover_genes, gene_type='cover gene',label='') # empty label for covered genes
  nodes_df <- dplyr::bind_rows(nodes_df, nodes_df2)
  dim(nodes_df)
  dim(edges_df)
  
  # nodes_df <- nodes_df %>%
  #   dplyr::mutate(size=
  #                   case_when(gene_type=='driver gene' ~ 5,
  #                             TRUE ~ 4))
  # dim(nodes_df)
  # rs <- runif(dim(nodes_df)[1], 1, 4)
  # nodes_df$size <- rs
  
  exp_df <- exp_df %>%
    dplyr::filter(symbol %in% nodes_df$name) %>%
    dplyr::select(log2FoldChange, symbol) %>%
    dplyr::mutate(log2FoldChange=abs(log2FoldChange))%>%
    dplyr::rename(size=log2FoldChange)
  nodes_df <- nodes_df %>%
    dplyr::left_join(exp_df, by=c('name'='symbol'))
  # sum(nodes_df$name %in% exp_df$symbol)
  # nodes_df$size # check if cover genes with smaller size to emphasize only driver genes
  nodes_df <- nodes_df %>%
    dplyr::mutate(
      size=case_when(
        size > 3 ~ 3, 
        size < -3 ~ -3, 
        TRUE ~ size
      )
    )
  nodes_df <- nodes_df %>%
    dplyr::mutate(
      size=case_when(
        gene_type== "cover gene"  ~ 0.5, 
        TRUE ~ size
      )
    )
    
  ## Adding pathways labels here
  ## If nodes are in pathways, coloring nodes by different colors, and edge color between nodes as well
  
  pw_total <- tibble::tibble()
  for(rf in pw_df$reference_set){
    pw_tmp <- pw_df %>%
      dplyr::filter(reference_set==rf)
    genes_used <- as.character(unlist(strsplit(pw_tmp$signif_genes, split = ',')))
    stat_tmp <- tibble::tibble(gene_name=genes_used, pathway=rf)
    pw_total <- dplyr::bind_rows(pw_total, stat_tmp)
  }
  dim(pw_total)
  # View(pw_total)
  pw_total <- pw_total[!duplicated(pw_total$gene_name),]
  # nodes_df$label
  
  nodes_df <- nodes_df %>%
    dplyr::left_join(pw_total, by=c('name'='gene_name')) 
  dim(nodes_df )
  summary(as.factor(nodes_df$pathway))
  nodes_df <- nodes_df %>%
    dplyr::mutate(
      pathway=case_when(
        !is.na(pathway) ~ pathway,
        TRUE ~ 'no_pathway')
    )
  ## if genes are in pathways --> show genes name, otherwise, empty label text
  nodes_df <- nodes_df %>%
    dplyr::mutate(
      label=case_when(
        pathway!='no_pathway' ~ name,
        TRUE ~ label
      )
    )
  dim(nodes_df)
  # edges_df <- edges_main
  # datatag <- subtag
  pg <- viz_graph(nodes_df, edges_df, subtag, save_output_dir)
  pg
  ggsave(paste0(save_output_dir,subtag,"_network_2.svg"),
         plot = pg,
         height = 7,
         width = 10,
         # useDingbats=F
  )
  
  
  
  return(pg)
  
}
viz_graph <- function(nodes_df, edges_df, datatag, save_output_dir){
  nodes_df$pathway <- as.factor(nodes_df$pathway)
  ls_pws <- unique(nodes_df$pathway)
  ls_pws <- ls_pws[ls_pws!='no_pathway']
  for(pw in ls_pws){
    nodes_pw2 <- nodes_df %>%
      dplyr::filter(pathway==pw)  
    edges_df <- edges_df %>%
      dplyr::mutate(
        edge_type=case_when(
          x %in% nodes_pw2$name & y %in% nodes_pw2$name ~ pw,
          TRUE ~ edge_type)
      )
    
  }
  # library(RColorBrewer)
  
  cols_ref <- brewer.pal(8, "Set2")
  cols_use <- c('#D3D3D3', cols_ref[1:length(ls_pws)])
  names(cols_use) <- c('no_pathway',as.character(ls_pws))
  
  cols_use <- c('#b2d8b2','#800080')
  names(cols_use) <- c("cover gene", "driver gene")
  # cols_use
  # summary(as.factor(edges_df$edge_type))
  # summary(as.factor(nodes_df$gene_type))
  # nodes_df$desc <- paste0(nodes_df$gene_type, '_',nodes_df$pathway)
  # summary(as.factor(nodes_df$pathway))
  # edge_cols_use <- c('grey', 'darkgreen','darkred')
  # names(edge_cols_use) <- c('driver_cover_connected','pathway2','pathway3')
  # 
  # edge_cols_use <- c('grey', 'darkgreen','darkred')
  # names(edge_cols_use) <- c('driver_cover_connected','pathway2','pathway3')
  # 
  # edges_wd <- c(0.4, 0.8, 0.9)
  # names(edges_wd) <- c('driver_cover_connected','pathway2','pathway3')
  # 
  gt_cols <- c(0.9, 0.6)
  names(gt_cols) <- c("driver gene","cover gene")
  sz <- c(0.9, 0.6)
  names(sz) <- c("driver gene","cover gene")
  # fc1 <- c('bold', 'italic')
  # names(fc1) <- c("driver gene","cover gene")
  ## get pathway labels here
  # make a tidy graph
  # nodes = NULL, edges = NULL, directed = TRUE, node_key = "name"
  ig <- tidygraph::tbl_graph(nodes = nodes_df, edges = edges_df, directed = FALSE)
  # ig
  # setting theme_graph 
  set_graph_style()
  
  ## To Do: need to check edge directions here 
  # basic plot
  
  ## To do: alpha scale, gene type=signf, not signif, cover genes
  pg <- ggraph(ig, layout='kk') + 
    geom_edge_link(alpha=.3, colour='darkgrey') + #, colour='darkgrey' aes(color = edge_type, width=edge_type)
    geom_node_point(aes(size=size, colour = gene_type), alpha=0.7) + #shape = gene_type, 
    geom_node_point(aes(size=size, alpha=gene_type), shape = 1) + #shape = gene_type, , colour='black'
    geom_node_text(aes(label = label), repel=TRUE, max.overlaps = Inf, size=5.5) + 
    # scale_edge_colour_manual(values=edge_cols_use) + 
    # scale_edge_width_manual(values=edges_wd) + 
    # scale_size_manual(values = sz) + 
    scale_color_manual(values = cols_use) + 
    scale_alpha_manual(values = gt_cols) +
    theme(legend.position = 'bottom',
          legend.text  = element_text(size=11, color='black'))
  # pg
  # scale_size_continuous() + 
  # scale_size_area(values=sz)
  # pg
  # geom_edge_link(aes(color = factor(to), width = log(weight)), alpha = 0.5, 
  #                start_cap = circle(2, 'mm'), end_cap = circle(2, 'mm'))
  # Define colors of locations and characters 
  # V(pg)$color <- "grey20"
  # V(pg)$color[V(pg)$name %in% driver_df$gene_symbol] <- "red"
  
  # BiocManager::install('svglite', ask=F)
  ggsave(  
    filename = paste0(save_output_dir,datatag,"_graph.svg"),  
    plot = pg,  
    height = 8, width = 12, dpi = 150)
  
  saveRDS(pg, paste0(save_output_dir, datatag, '_graph.rds'))
  ## to do: significant genes with red color node
  ## size = # connection, with log(size)
  ## pathway genes: including cover genes, with label text
  return(pg)
  
}
run_DriverNet <- function(mutation_mtx, exp_mtx, influenceGraph, sampleGenes, save_output_dir){
  
  dim(mutation_mtx)
  dim(exp_mtx)
  dim(influenceGraph)
  driversList = DriverNet::computeDrivers(mutation_mtx, 
                                          exp_mtx,
                                          influenceGraph,
                                          outputFolder=NULL, printToConsole=FALSE)
  
  
  class(driversList)
  # driversList@drivers
  
  ## NOTE: Main version from package
  # randomDriversResult = DriverNet::computeRandomizedResult(
  #   patMutMatrix=mutation_mtx,
  #   patOutMatrix=exp_mtx, 
  #   influenceGraph=influenceGraph,
  #   geneNameList=sampleGenes, outputFolder=NULL, 
  #   printToConsole=F,numberOfRandomTests=100, weight=FALSE, 
  #   purturbGraph=FALSE, purturbData=TRUE)
  
  ## corrected errors version: computeRandomizedResultV2
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
  
  total_df <- print_output(driversList, res, save_output_dir)
  # View(total_df)
  
  total_df1 <- total_df %>%
    dplyr::filter(`p-value`<=0.05)
  total_df2 <- total_df %>%
    dplyr::filter(`p-value`>0.05)
  print(total_df1)
  # View(total_df1)
  # View(total_df2)
  # paste(total_df1$gene_symbol, collapse = ', ')
  saveRDS(sampleGenes, paste0(save_output_dir, 'sampleGenes.rds'))
  saveRDS(randomDriversResult, paste0(save_output_dir, 'randomDriversResult_1050iterations.rds'))
  saveRDS(res, paste0(save_output_dir, 'res.rds'))
  saveRDS(driversList, paste0(save_output_dir, 'driversList.rds'))
  # save.image(paste0(save_dir, "resuls_DriverNet.RData")
  
}

get_influence_graph_from_InteractomeFI <- function(obs_genes, save_dir, save_data=T){
  
  input_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/preprocessed_09April2024/drivernet_demo/databases/'
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
    data.table::fwrite(edge_df, paste0(save_dir, 'InteractomeFI_2022_wt_Score.csv.gz'))  
  }
  
  influenceGraph <- matrix(0, nrow=length(obs_genes), ncol=length(obs_genes))
  rownames(influenceGraph) <- obs_genes
  colnames(influenceGraph) <- obs_genes
  for(i in seq(1:dim(edge_df)[1])){
    influenceGraph[edge_df[i,]$g1,edge_df[i,]$g2] <- 1
    influenceGraph[edge_df[i,]$g2,edge_df[i,]$g1] <- 1
  }
  # total_genes_ig <- union(edges_ig$g1, edges_ig$g2)
  for(i in seq(1:length(obs_genes))){
    influenceGraph[obs_genes[i],obs_genes[i]] <- 1
  }
  # head(influenceGraph[1:20,1:10])
  print(dim(influenceGraph))
  if(save_data){
    data.table::fwrite(as.data.frame(influenceGraph), paste0(save_dir, 'influenceGraph.csv.gz'), row.names = T, col.names = T)
  }
  res <- list(influenceGraph=influenceGraph, edge_df=edge_df)
  return(res)
  
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
    pathway_fn = '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/testing/DriverNet_testing/h.all.v7.0.symbols.gmt'  
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
get_mutation_data <- function(obs_genes, obs_clones, save_dir){
  script_dir <- "/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/"
  cnv <- data.table::fread(paste0(script_dir, 'SA919_DE_analysis_DESeq2_Hoa_09April2024/mapping_gene_cnv_SA919.csv.gz'))  
  cnv <- cnv[!duplicated(cnv$ensembl_gene_id),]
  dim(cnv)
  # head(cnv)
  ## need to generate a cnv profile for each patient
  rownames(cnv) <- cnv$ensembl_gene_id
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
  sids <- paste0('clone', obs_clones)
  rownames(mutation_mtx) <- sids
  colnames(mutation_mtx) <- obs_genes
  # mutation_mtx['cloneC',cnv_var_genes] <- 1
  mutation_mtx[paste0('clone',obs_clones[2]),genes_used] <- 1
  print(dim(mutation_mtx))
  data.table::fwrite(as.data.frame(mutation_mtx), paste0(save_dir, 'mutation_mtx.csv.gz'), row.names=T)
  return(mutation_mtx)
}

get_expression_data <- function(obs_clones, obs_genes, subtag, save_dir){
  # data_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/preprocessed_09April2024/SA919_DE_analysis_DESeq2_Hoa_09April2024/SA919_CMetMainMix_Bpri_DESeq2/'
  # de_df <- data.table::fread(paste0(data_dir, 'DE_signif_genes.csv.gz'))  
  # dim(de_df)
  sids <- paste0('clone', obs_clones)
  data_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/results_bulkRNAseq/SA919_full/'
  save_fig_dir <- paste0(data_dir, subtag, '/')
  exp_df <- data.table::fread(paste0(save_fig_dir, subtag, '_DE_genes.csv.gz'))
  print(subtag)
  print(dim(exp_df))
  
  if(!'symbol' %in% colnames(exp_df)){
    ref <- annotables::grch38 %>%
      dplyr::select(ensgene,symbol) %>%
      dplyr::rename(ensembl_gene_id=ensgene)
    ref <- ref[!duplicated(ref$ensembl_gene_id),]
    exp_df <- exp_df %>%
      dplyr::left_join(ref, by='ensembl_gene_id')
  }
  exp_df <- exp_df %>%
    dplyr::filter(symbol %in% obs_genes)
  print(dim(exp_df))
  obs_genes <- intersect(obs_genes, exp_df$symbol)
  length(unique(obs_genes))
  # View(head(exp_df))
  
  exp_mtx <- matrix(FALSE, nrow=2, ncol=length(obs_genes))
  dim(exp_mtx)
  rownames(exp_mtx) <- sids
  colnames(exp_mtx) <- obs_genes
  head(exp_mtx[,1:4])
  sum(exp_df$symbol %in% obs_genes)
  # outlier_genes_B <- exp_df %>%
  #   dplyr::filter(log2FoldChange<0 & symbol %in% obs_genes) %>%
  #   dplyr::pull(symbol)
  
  # outlier_genes_C <- exp_df %>%
  #   dplyr::filter(log2FoldChange>0 & symbol %in% obs_genes) %>%
  #   dplyr::pull(symbol)
  # length(outlier_genes_C)
  outlier_genes <- exp_df %>%
    dplyr::filter(log2FoldChange>0 & symbol %in% obs_genes) %>%
    dplyr::pull(symbol)
  length(outlier_genes)
  exp_mtx[paste0('clone',obs_clones[2]), outlier_genes] <- TRUE
  dim(exp_mtx)
  data.table::fwrite(as.data.frame(exp_mtx), paste0(save_dir, 'exp_mtx.csv.gz'), row.names=T)
  res <- list(subtag=subtag, exp_df=exp_df, exp_mtx=exp_mtx, obs_genes=obs_genes)
  return(res)
}
get_obs_genes <- function(obs_clones){
  
  ## To Do: check negative, positive tendency, may need to redo DE analysis
  
  input_dir <- "/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/preprocessed_09April2024/"
  input_dir <- "/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/preprocessed_09April2024/"
  
  meta_genes <- data.table::fread(paste0(input_dir, 'dispersion_cis_trans_genes_cloneCB.csv'))
  dim(meta_genes)
  summary(as.factor(meta_genes$gene_type))
  
  script_dir <- "/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/"
  cnv <- data.table::fread(paste0(script_dir, 'SA919_DE_analysis_DESeq2_Hoa_09April2024/mapping_gene_cnv_SA919.csv.gz'))  
  dim(cnv)
  cnv <- cnv[!duplicated(cnv$ensembl_gene_id),]
  obs_genes <- intersect(meta_genes$symbol, cnv$gene_symbol)
  length(unique(obs_genes))
  # res <- list(obs_genes=obs_genes, total_genes=meta_genes$symbol)
  return(obs_genes)  
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
  
  data.table::fwrite(total_df, paste0(save_dir, 'significant_genes_DriverNet.csv'))
  
  return(total_df)
}
# patMutMatrix <- mutation_mtx
# patOutMatrix <- exp_mtx
# # influenceGraph influenceGraph
# geneNameList <- sampleGenes
# outputFolder = NULL
# printToConsole = FALSE
# numberOfRandomTests = 2
# weighted = FALSE
# purturbGraph = FALSE
# purturbData = TRUE


computeRandomizedResultV2 <- function(patMutMatrix, patOutMatrix, influenceGraph, geneNameList, outputFolder = NULL, printToConsole = FALSE, numberOfRandomTests = 500,  weighted = FALSE, purturbGraph = FALSE, purturbData = TRUE) {
  if (weighted) {
    stop("Weighted algorithm is not implemented.")
  }
  if (purturbGraph) {
    stop("Purturb graph option is not implemented.")
  }
  if (purturbData == FALSE) {
    stop("Purturb data option must be TRUE at this moment.")
  }  
  ## set maxNumOfDrivers to the number of mutations
  maxNumOfDrivers <- length(which(colSums(patMutMatrix)>0))
  out_fname <- c()
  if (!identical(outputFolder, NULL)) {
    out_fname <- c(paste(outputFolder, "randomized_result_", format(Sys.time(), "%Y_%m_%d_%H%M%S"), ".txt", sep=""))	
  }
  if (printToConsole) {
    out_fname <- c(out_fname, "")
  }  
  
  if (purturbGraph == FALSE) {
    nG_out <- DriverNet:::.neighborGraph(influenceGraph[intersect(colnames(patMutMatrix), rownames(influenceGraph)), intersect(colnames(influenceGraph), colnames(patOutMatrix))])   ## genes whose expression is affected by dna mutation in g
  }
  
  coverageResults <- vector(mode="list", length=numberOfRandomTests)
  i <- 1
  for (i in 1:numberOfRandomTests) {
    # print('---------------------------')
    if(i %% 50 ==0){
      print(paste0('i: ', i))  
    }
    
    randomPatMutMatrix <- patMutMatrix
    randomPatOutMatrix <- patOutMatrix
    if (purturbData) {
      ## purturb the data by randomizing the gene names
      randomizedOutlierNames <- geneNameList[sample(1:length(geneNameList))[1:ncol(patOutMatrix)]] 
      randomizedMutationNames <- geneNameList[sample(1:length(geneNameList))[1:ncol(patMutMatrix)]]
      colnames(randomPatOutMatrix) <- randomizedOutlierNames
      colnames(randomPatMutMatrix) <- randomizedMutationNames
    }
    # print(dim(randomPatOutMatrix))
    
    
    if (purturbGraph) {
      numColInfluenceGraph=ncol(influenceGraph) 
      ## purturb influenceGraph as well
      infGraphNewGeneNames <- geneNameList[sample(1:length(geneNameList))[1:numColInfluenceGraph]]
      colnames(influenceGraph) <- infGraphNewGeneNames
      rownames(influenceGraph) <- infGraphNewGeneNames
      nG_out <- DriverNet:::.neighborGraph(influenceGraph[intersect(colnames(randomPatMutMatrix), rownames(influenceGraph)), intersect(colnames(influenceGraph), colnames(randomPatOutMatrix))])   ## genes whose expression is affected by dna mutation in g
    }
    ## reuse these intersections
    affectedGenesIntersection <- intersect(colnames(influenceGraph), colnames(randomPatOutMatrix))
    patientIntersection <- intersect(rownames(randomPatMutMatrix), rownames(randomPatOutMatrix))
    mutatedGenesIntersection <- intersect(rownames(influenceGraph), colnames(randomPatMutMatrix))
    ## pre-process the new matrices
    influenceGraph2 <- influenceGraph[mutatedGenesIntersection, affectedGenesIntersection]
    randomPatOutMatrix2 <- randomPatOutMatrix[patientIntersection, affectedGenesIntersection]
    randomPatMutMatrix2 <- randomPatMutMatrix[patientIntersection, mutatedGenesIntersection] 
    
    # print(dim(influenceGraph2))
    # print(dim(randomPatOutMatrix2))
    # print(dim(randomPatMutMatrix2))
    if (weighted) {
      ## Not implemented
      # drivers[[i]] <- .greedyGeneDriverSelection_weighted(out_fname, randomPatOutMatrix2, randomPatMutMatrix2, influenceGraph2, maxNumOfDrivers)[[1]]   
    } else {
      # print('Debug 1')
      # runResult <- DriverNet:::.greedyGeneDriverSelection(out_fname, randomPatOutMatrix2, randomPatMutMatrix2,
      #                                                     influenceGraph2, nG_out, NULL, maxNumOfDrivers)
      runResult <- NULL
      runResult <- tryCatch({
        DriverNet:::.greedyGeneDriverSelection(out_fname, randomPatOutMatrix2, randomPatMutMatrix2,
                                               influenceGraph2, nG_out, NULL, maxNumOfDrivers)
      }, error = function(e) {
        # message("Error: ", e$message)
        NULL  # Returning NA in case of an error
      })
      
    # }
    
      # length(intersect(colnames(randomPatMutMatrix2), colnames(randomPatOutMatrix2)))
      # length(intersect(colnames(randomPatMutMatrix2), names(nG_out)))
      # length(intersect(colnames(randomPatOutMatrix2), names(nG_out)))
      # length(intersect(colnames(influenceGraph2), names(nG_out)))
      # length(intersect(rownames(influenceGraph2), names(nG_out)))
      # print(runResult)
      # print('Debug 2')
      if(!is.null(runResult)){
      if (!is.null(runResult$drivers)) {
        len <- length(runResult$drivers)
        # print(len)
        # coverageVector <- vector(mode="integer", length=len)
        # names(coverageVector) <- runResult$drivers
        if (len>0) {
          coverageVector <- vector(mode="integer", length=len)
          names(coverageVector) <- runResult$drivers
          # print(runResult$actualEvents)
          
          for (k in 1:len) {
            ## nrow(actualEvents[[i]][[k]]) equals the number of events covered by the k-th driver found in the i-th run
            coverageVector[[k]] <- nrow(runResult$actualEvents[[k]])
          }
          # print('coverageVector: ')
          # print(coverageVector)
          coverageResults[[i]] <- coverageVector
          
        }
        # print('coverageVector: ')
        # print(coverageVector)
        # coverageResults[[i]] <- coverageVector
        # print('---------------------------')
        
      }
      }  
    }
  }
  
  if (!identical(outputFolder, NULL)) {
    message(paste("Log file written to: ", out_fname[[1]], "\n", sep=""))
  }  
  
  coverageResults
}


# example_function <- function(out_fname, randomPatOutMatrix2, randomPatMutMatrix2, 
#                              influenceGraph2, nG_out, NULL, maxNumOfDrivers) {
#   result <- tryCatch({
#     DriverNet:::.greedyGeneDriverSelection(out_fname, randomPatOutMatrix2, randomPatMutMatrix2, 
#                                            influenceGraph2, nG_out, NULL, maxNumOfDrivers)
#   }, error = function(e) {
#     message("Error: ", e$message)
#     NA  # Returning NA in case of an error
#   })
#   return(result)
# }



check_cis_genes_tendency <- function(){
  obs_clones <- c('B','C') ## check % cis genes, a trackplot
  # obs_clones <- c('A','B')
  # obs_clones <- c('A','C')
  sids <- paste0('clone', obs_clones)
  subtag <- paste0(obs_clones[2],'met_',obs_clones[1],'pri')
  data_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/results_bulkRNAseq/SA919_full/'
  save_fig_dir <- paste0(data_dir, subtag, '/')
  exp_df <- data.table::fread(paste0(save_fig_dir, subtag, '_DE_genes.csv.gz'))
  print(subtag)
  print(dim(exp_df))
  sum(abs(exp_df$log2FoldChange)>1)
  
  if(!'symbol' %in% colnames(exp_df)){
    ref <- annotables::grch38 %>%
      dplyr::select(ensgene,symbol) %>%
      dplyr::rename(ensembl_gene_id=ensgene)
    ref <- ref[!duplicated(ref$ensembl_gene_id),]
    exp_df <- exp_df %>%
      dplyr::left_join(ref, by='ensembl_gene_id')
  }
  obs_genes_symb <- exp_df$symbol
  length(obs_genes_symb)
  # paste(obs_genes_symb, collapse = '')
  pw_df <- get_gprofiler_pathways_obsgenes(obs_genes_symb, save_fig_dir, subtag, 
                                           custom_id=NULL, 
                                           pathway_fn=NULL, save_data=T, correction_func = 'gSCS')
  
  
  dim(pw_df)
  # pw_df
  
  # t <- data.table::fread('/Users/hoatran/Downloads/gProfiler_All_C_Up_Th1.txt')
  # t <- t %>%
  #   dplyr::filter(p_value<0.05)
  # t$p_value
  # write.table(obs_genes_symb,'/Users/hoatran/Downloads/gProfiler_testBC.txt',sep="\t",quote=F,row.names = F)
  # write_file(obs_genes_symb, file = '/Users/hoatran/Downloads/gProfiler_testBC.txt') #, append = '\n'
  
  cnv <- load_cis_genes(obs_clones)
  dim(cnv)
  exp_df <- exp_df %>%
    dplyr::inner_join(cnv, by=c('symbol'='gene_symbol'))
  print(dim(exp_df))
  # View((head(exp_df)))
  pp <- sum(exp_df$cn_change>0 & exp_df$log2FoldChange>0)
  nn <- sum(exp_df$cn_change<0 & exp_df$log2FoldChange<0)
  # 100*(pp+nn)/dim(exp_df)[1]
  
  
}

load_cis_genes <- function(obs_clones){
  script_dir <- "/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/"
  cnv <- data.table::fread(paste0(script_dir, 'SA919_DE_analysis_DESeq2_Hoa_09April2024/mapping_gene_cnv_SA919.csv.gz'))  
  cnv <- cnv[!duplicated(cnv$ensembl_gene_id),]
  dim(cnv)
  # head(cnv)
  ## need to generate a cnv profile for each patient
  rownames(cnv) <- cnv$ensembl_gene_id
  cnv <- as.data.frame(cnv)
  cols_use <- c('gene_symbol','chr') #'ensembl_gene_id',
  cols_use <- c(cols_use, obs_clones)
  rownames(cnv) <- cnv$ensembl_gene_id
  cnv <- cnv %>%
    dplyr::select(all_of(cols_use))
  # genes_used <- cnv$ensembl_gene_id[rv>0]
  # length(genes_used)
  # cnv <- cnv %>%
  #   dplyr::filter(gene_symbol %in% obs_genes)
  dim(cnv)
  cnv_vals <- cnv %>%
    dplyr::select(all_of(obs_clones))
  rv <- rowVars(as.matrix(cnv_vals))  # median copy number profile of obs clones
  # rv <- rowVars(as.matrix(cnv[,c(5, 6)]))  # median copy number profile of clone B, C 
  genes_used <- cnv$gene_symbol[rv>0]
  
  cnv <- cnv %>%
    dplyr::filter(gene_symbol %in% genes_used) %>%
    dplyr::mutate(cn_change=as.integer(!!sym(obs_clones[2])-!!sym(obs_clones[1])))
  
  
  return(cnv)
}

get_cis_trans_gene_type <- function(obs_clones, drivernet_df){
  obs_genes <- drivernet_df$gene_symbol
  script_dir <- "/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/"
  cnv <- data.table::fread(paste0(script_dir, 'SA919_DE_analysis_DESeq2_Hoa_09April2024/mapping_gene_cnv_SA919.csv.gz'))  
  cnv <- cnv[!duplicated(cnv$ensembl_gene_id),]
  dim(cnv)
  # head(cnv)
  ## need to generate a cnv profile for each patient
  rownames(cnv) <- cnv$ensembl_gene_id
  cnv <- as.data.frame(cnv)
  cols_use <- c('gene_symbol') #'ensembl_gene_id',
  cols_use <- c(cols_use, obs_clones)
  rownames(cnv) <- cnv$ensembl_gene_id
  cnv <- cnv  %>%
    dplyr::filter(gene_symbol %in% obs_genes)%>%
    dplyr::select(all_of(cols_use))
  # genes_used <- cnv$ensembl_gene_id[rv>0]
  # length(genes_used)
  # cnv <- cnv %>%
  #   dplyr::filter(gene_symbol %in% obs_genes)
  dim(cnv)
  # cnv_vals <- cnv %>%
  #   dplyr::select(all_of(obs_clones))
  # rv <- rowVars(as.matrix(cnv_vals))  # median copy number profile of obs clones
  # # rv <- rowVars(as.matrix(cnv[,c(5, 6)]))  # median copy number profile of clone B, C 
  # genes_used <- cnv$gene_symbol[rv>0]
  
  cnv <- cnv %>%
    # dplyr::filter(gene_symbol %in% genes_used) %>%
    dplyr::mutate(cn_change=as.integer(!!sym(obs_clones[2])-!!sym(obs_clones[1])))
  
  drivernet_df <- drivernet_df %>%
    dplyr::left_join(cnv, by='gene_symbol') %>%
    dplyr::mutate(gene_type=
                    case_when(
                      cn_change != 0 ~ 'cis',
                      TRUE ~ 'trans'
                    )) %>%
    dplyr::select(-cn_change)
  cols_use <- c('gene_symbol','rank','p-value','cover_genes','chr','gene_type', obs_clones,'description')
  
  drivernet_df <- drivernet_df %>%
    dplyr::select(all_of(cols_use), everything())
  
  return(drivernet_df)
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

viz_heatmap <- function(){
  ## common driver genes between clones
  ## plot for CNA, and for gene exp
}
