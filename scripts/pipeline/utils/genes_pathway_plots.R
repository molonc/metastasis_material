suppressPackageStartupMessages({
  require(RColorBrewer)
  # require(glue)
  # require(ggrepel)
  # require(tidyr)
  # require(fgsea)
  # require(scales)
  # require(ComplexHeatmap)
  # require(DOSE)
  # require(enrichplot)
  require(stringr)
  require(ggraph)
  require(ggplot2)
  require(igraph)
  require(tidyverse)
})
# generate_genes_network(input_dir, datatag='SA609_R-UTTTT_vs_E-UUUUUUUU', 
#                        obs_pw='top_up_pathway_cls_',
#                        de_desc='SA609_UTTTT_UUUUUUUU_R_E',
#                        genes_set='hallmark',
#                        pattern_use='^hallmark_'
# )
cnetplot_v2 <- function(geneSets, datatag,
                        tag, save_dir,
                        foldChange=NULL,
                        layout = "kk",
                        colorEdge = FALSE,
                        circular = FALSE,
                        node_label = "all",
                        node_scale = 1) 
{
  
  node_label <- match.arg(node_label, c("category", "gene", "all", "none"))
  
  # if (circular) {
  #   layout <- "linear"
  #   geom_edge <- geom_edge_arc
  # } else {
  #   geom_edge <- geom_edge_link
  # }
  # foldChange <- enrichplot:::fc_readable(x, foldChange)
  # geneSets <- enrichplot:::extract_geneSets(x, showCategory) #showCategory = 5,
  g <- enrichplot:::list2graph(geneSets)
  size <- sapply(geneSets, length)
  V(g)$size <- min(size)/2
  n <- length(geneSets)
  V(g)$size[1:n] <- size
  
  # nb_genes <- sapply(geneSets, length)
  # V(g)$nb_genes <- min(nb_genes)/2
  # n <- length(geneSets)
  # V(g)$nb_genes[1:n] <- nb_genes 
  
  
  if (colorEdge) {
    E(g)$category <- rep(names(geneSets), sapply(geneSets, length))
    edge_layer <- ggraph::geom_edge_link(aes_(color = ~category), alpha=.9)
  } else {
    # edge_layer <- geom_edge(alpha=.8, colour='darkgrey')
    edge_layer <- ggraph::geom_edge_link(alpha=.2, colour='darkgrey')
    
  }
  
  if (!is.null(foldChange)) {
    fc <- foldChange[V(g)$name[(n+1):length(V(g))]]
    V(g)$color <- NA
    V(g)$color[(n+1):length(V(g))] <- fc
    #palette <- fc_palette(fc)
    p <- ggraph(g, layout=layout, circular = circular) +
      edge_layer +
      geom_node_point(aes_(color=~as.numeric(as.character(color)),
                           size=~size)) +
      scale_colour_gradient2(name = "logFC", low = "green",
                             mid = "blue", high = "red")
  } else {
    V(g)$color <- "#B3B3B3"
    V(g)$color[1:n] <- "#E5C494"
    p <- ggraph::ggraph(g, layout=layout, circular=circular) +
      edge_layer +
      geom_node_point(aes_(color=~I(color), size=~size))
    # p
  }
  
  p <- p + scale_size(range=c(3, 15) * node_scale,
                      breaks=unique(round(seq(min(size), max(size), length.out=4)))) +
    theme_void()
  if (node_label == "category") {
    p <- p + geom_node_text(aes_(label=~name), data = p$data[1:n,]) #, fontface ="bold"
  } else if (node_label == "gene") {
    p <- p + geom_node_text(aes_(label=~name), data = p$data[-c(1:n)], repel=TRUE)
  } else if (node_label == "all") {
    if (utils::packageVersion("ggrepel") >= "0.9.0") {
      p <- p + geom_node_text(aes_(label=~name), repel=TRUE, bg.color = "white", size = 7)
    } else {
      p <- p + geom_node_text(aes_(label=~name), repel=TRUE, size=5.7)
    }
    
  }
  
  p <- p + ggtitle(paste0(datatag," \n",tag)) + 
            theme(
              plot.title=element_text(color="black", size=15, hjust = 0.5, face='bold')
            )
  if(length(geneSets)>=3 & length(geneSets)<=5){
    ht <- 750
    wd <- 1000
  } else if(length(geneSets)<3){
    ht <- 350
    wd <- 600
  }else{
    ht <- 1000
    wd <- 1300
  }
  datatag <- gsub(' ','_',datatag)
  tag <- gsub(' ','_',tag)
  png(paste0(save_dir, tag,'_genes_network.png'), height = 2*ht, width=2*wd,res = 2*72)
  print(p)
  dev.off()
  
  # return(p)
}

generate_genes_network_plot <- function(pathway_df, datatag, tag, save_dir, foldChange, selected_pathways=FALSE, obsvered_pws){
  rownames(pathway_df) <- pathway_df$Term
  
  obsvered_pws <- paste0('HALLMARK_',obsvered_pws)
  if(selected_pathways){
    pathway_df <- pathway_df %>%
      dplyr::filter(Term %in% obsvered_pws)
  }
  print(dim(pathway_df))
  pattern_use <- '^hallmark_'
  geneSets <- list()
  for(p in rownames(pathway_df)){
    gene_ls <- strsplit(as.character(pathway_df[p,'ledge_genes']),';')
    gene_ls <- unlist(gene_ls, use.names=FALSE)
    gene_ls <- gene_ls[gene_ls %in% names(foldChange)]
    p <- gsub(pattern_use,'',tolower(p))
    geneSets[[p]] <- gene_ls
    # print(length(gene_ls))
  }
  cnetplot_v2(geneSets, datatag, tag, save_dir, foldChange)
}


viz_genes_network <- function(save_dir, 
                                   pw_fn,
                                   logFC_fn,
                                   datatag='SA',
                                   tag=''){
  if(!dir.exists(save_dir)){
    dir.create(save_dir, showWarnings = F)
  }
  stat_genes <- read.csv(logFC_fn, check.names=F, stringsAsFactors=F) #total_markers.txt
  # View(head(stat_genes))
  foldChange <- stat_genes$logFC
  names(foldChange) <- stat_genes$gene_symbol
  
  # pathway_df <- read.csv(paste0(input_dir, genes_set,'/SA609_UTTTT_UUUUUUUU_up_pathway.csv'), check.names=F, row.names=1)
  
  
  if(file.exists(pw_fn)){
    pathway_df <- read.csv(pw_fn, check.names=F, stringsAsFactors=F)
    # View(head(pathway_df))
    # Get up_pathways:
    up_pathway_df <- pathway_df %>%
      dplyr::filter((es > 0) & (fdr < 0.05))
    print(dim(up_pathway_df))
    
    down_pathway_df <- pathway_df %>%
      dplyr::filter((es < 0) & (fdr < 0.05))
    print(dim(down_pathway_df))
    
    generate_genes_network_plot(up_pathway_df, datatag, paste0('Up-regulated pathways ',tag), 
                                save_dir, foldChange, selected_pathways=FALSE, obsvered_pws=NULL)
    
    generate_genes_network_plot(down_pathway_df, datatag, paste0('Down-regulated pathways ',tag), 
                                save_dir, foldChange, selected_pathways=FALSE, obsvered_pws=NULL)
    # obsvered_pws <- c('EPITHELIAL_MESENCHYMAL_TRANSITION','OXIDATIVE_PHOSPHORYLATION','TGF_BETA_SIGNALING',
    #                   'HYPOXIA','MTORC1_SIGNALING','MYC_TARGETS_V1','PEROXISOME')
    # generate_genes_network_plot(pathway_df, tag, datatag, save_dir, foldChange, selected_pathways=TRUE, obsvered_pws)
    # 
    
  }else{
    print("Double check the file name ")
    print(pw_fn)
  }
  
}

# SA1035: UTTTT-H vs UUUUU-E
# SA609: UTTT-R vs UUUU-H (edited) 
# SA535 cisplatin: UUTTT S-T. vs. UUUUU - J
# SA535 CX: UXXXX-U vs UUUUU - J

input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/pathways/'
save_dir <- paste0(input_dir, 'pathway_plots/')
pw_fn <- paste0(input_dir, 'pathway12-SA535_CX5461-v6-HALLMARK-FDR-0.001-pathway.csv')
logFC_fn <- paste0(input_dir,'pathway12_SA535_UXXXX_U_UUUUU_Q_logfc_results.csv')
datatag <- 'SA535'
tag <- 'SA535_CX5461_UXXXX_U_UUUUU_Q'


input_dir <- '/home/htran/Projects/farhia_project/drug_resistance/rnaseq/differential_expression/results/'
save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/pathways/pathway_plots/'
pw_fn <- paste0(input_dir, 'pathway12-SA535_CX5461-v6-HALLMARK-FDR-0.001-pathway.csv')
logFC_fn <- paste0(input_dir,'pathway12_SA535_UXXXX_U_UUUUU_Q_logfc_results.csv')
datatag <- 'SA535'
tag <- 'SA535_CX5461_UXXXX_U_UUUUU_Q'

input_dir <- '/home/htran/Projects/farhia_project/drug_resistance/rnaseq/differential_expression/results/'
logFC_ls <- c('SA609-v6/comps/pathway9_SA609_UTTTT_R_UUUUU_H','SA1035-v6/comps/ressens8_SA1035_UTTTT_H_UUUUU_E',
           'SA535-v6/comps/pathway4Tvs4U_SA535_UUTTTT_S_T_UUUUUU_J_Q','SA535-v6/comps/ressens18_SA535_UXXXX_U_UUUUU_J')

pw_ls <- c('SA609-v6/comps/pathway9-SA609-v6-HALLMARK-FDR-0.001','SA1035-v6/comps/ressens8-SA1035-v6-HALLMARK-FDR-0.001',
           'SA535-v6/comps/pathway4Tvs4U_2-SA535_cisplatin-v6-HALLMARK-FDR-0.001','SA535-v6/comps/ressens18-SA535_CX5461-v6-HALLMARK-FDR-0.001')
tag_ls <- c('UTTTT-R versus UUUUU-H', 'UTTTT-H versus UUUUU-E', 'UUTTTT-S_T versus UUUUUU-J', 'UXXXX-U versus UUUUU-J')
datatag_ls <- c('SA609','SA1035','SA535 Cisplatin','SA535 CX5461')

for(i in seq(length(logFC_ls))){
  # i <- 3
  viz_genes_network(save_dir,
                    paste0(input_dir,pw_ls[i],'-pathway.csv'),
                    paste0(input_dir,logFC_ls[i],'_logfc_results.csv'),
                    datatag_ls[i],
                    tag_ls[i]
  )
}


