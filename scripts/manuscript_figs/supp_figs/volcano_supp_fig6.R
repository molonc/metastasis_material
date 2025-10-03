suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(ggrepel)
  options(dplyr.summarise.inform = FALSE)
  options(tidyverse.quiet = TRUE)
})




plot_DE_genes_ggplot <- function(df, topGenes=NULL, capstr='', 
                                          # FDRcutoff=0.01, logFCcutoff=0.25, pValuecutoff=0.05,
                                          plttitle="A versus B", save_dir="",legendVisible=F,
                                          iscaption=TRUE, legend_verbose='none', save_plot=TRUE, 
                                          xl=NULL, yl=NULL){  
  yl <- c(0, 40)
  anno_text_sz <- 4
  nb_top_genes <- 20
  df$logFC <- df$log2FoldChange
  
  if(!'symbol' %in% colnames(df)){
    ref <- annotables::grch38 %>%
      dplyr::select(ensgene,symbol) %>%
      dplyr::rename(ensembl_gene_id=ensgene, gene_symbol=symbol)
    ref <- ref[!duplicated(ref$ensembl_gene_id),]
    df <- df %>%
      dplyr::left_join(ref, by='ensembl_gene_id')
  }else{
    df <- df %>%
      dplyr::rename(gene_symbol=symbol)
  }
  if(length(topGenes)>nb_top_genes){
    tmp <- df %>%
      dplyr::filter(gene_symbol %in% topGenes) %>%
      dplyr::arrange(-abs(logFC))
    # print(tmp)
    topGenes <- tmp$gene_symbol[1:nb_top_genes]
  }
  print(topGenes)
  # df <- df[abs(df$logFC)>logFCcutoff & df$FDR<FDRcutoff  & df$PValue<pValuecutoff,]
  maxLogFC <- 10
  # df$logFC <- ifelse(df$logFC>maxLogFC,maxLogFC,
  #                    ifelse(df$logFC<-maxLogFC,-maxLogFC,df$logFC))
  df$logFC <- sapply(df$logFC, function(x) replace(x, x > maxLogFC, maxLogFC))
  df$logFC <- sapply(df$logFC, function(x) replace(x, x < (-maxLogFC), -maxLogFC))
  
  thrs <- 40
  df$log10Adj <- sapply(-log10(df$pvalue), function(x) replace(x, x>thrs, thrs))
  ylabel <- bquote(~-Log[10] ~ ' Pvalue')
  # if(sum(is.na(df$padj))>0){
  #   df$log10Adj <- sapply(-log10(df$pvalue), function(x) replace(x, x>thrs, thrs))
  #   ylabel <- bquote(~-Log[10] ~ ' P value')
  # }else{
  #   df$log10Adj <- sapply(-log10(df$padj), function(x) replace(x, x>thrs, thrs))  
  #   ylabel <- bquote(~-Log[10] ~ ' P adj')
  # }
  
  my_font <- "Helvetica"
  thesis_theme <- theme(text=element_text(size = 8,family=my_font),
                          # plot.title = element_blank(),
                          plot.title = element_text(color="black", size=13, hjust = 0, face = "bold", family=my_font),
                          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                          axis.text.y = element_text(size=8, hjust = 0.5, family=my_font),
                          axis.text.x = element_text(size=10, hjust = 0.5, family=my_font),
                          axis.title = element_text(size=11, hjust = 0.5, family=my_font),
                          plot.caption = element_text(size=8, hjust = 1, family=my_font),
                          legend.title = element_text(size=9, hjust = 0.5, family=my_font), #, angle=90
                          legend.text = element_text(size=9, hjust = 0, family=my_font),
                          strip.text.x = element_text(color="black",size=9, family=my_font),
                          strip.text.y = element_text(color="black",size=9, family=my_font),
                          legend.spacing.x = unit(0.1, 'mm'),
                          legend.spacing.y = unit(0.1, 'mm'),
                          legend.key.height=unit(1,"line"),
                          # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          legend.position = legend_verbose)
  
  
  
  
  # topGenes <- get_top_genes_v2(df, minLogFC=0.25, nbtopup=10, nbtopdown=10)
  p <- ggplot(df, aes(x = logFC, y = log10Adj)) +
    geom_point(size=1.5, alpha=1, color='#088F8F') + #aes(color = gene_type), shape=21
    # scale_color_manual(values = cols_use_cis, name = NULL) + #name = "Gene Type"
    thesis_theme + 
    geom_text_repel(family = my_font,
                    data = df[df$gene_symbol %in% topGenes, ],
                    aes(label = gene_symbol), size = anno_text_sz,
                    box.padding = unit(0.35, "lines"),
                    point.padding = unit(0.3, "lines"),
                    max.overlaps = Inf,
                    min.segment.length = 0,  # draw segment lines, not matter how short they are)
                    color='black', segment.alpha = 0.2)
  p <- p + labs(x= bquote(~Log[2] ~ ' Fold Change'), 
                        y=ylabel, title = plttitle) #+ #paste0(plttitle,', ','in cis DE genes'). bquote(~-Log[10]~italic(FDR))
    # guides(color = guide_legend(ncol=2, override.aes = list(size=3, shape=16, alpha=1, title='In cis genes')))
  
  if(!is.null(xl)){
    p <- p + xlim(xl[1],xl[2])
  }
  if(!is.null(yl)){
    p <- p + ylim(yl[1],yl[2])
  }
  return(p)
}



viz_volcano_plts <- function(){
  # Loading DE genes files
  save_fig_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/figures/supp_figs/volcano_plts/'
  dir.create(save_fig_dir)
  input_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/DE_results/'
  # obs_clones <- c('B','C')
  # obs_clones <- c('A','B')
  obs_clones <- c('A','B')
  subtag <- paste0(obs_clones[2],'met_',obs_clones[1],'pri')
  plt_title <- paste0('Clone ',obs_clones[2],' metastasis vs. clone ',obs_clones[1],' primary')
  ab_df <- data.table::fread(paste0(input_dir, subtag,'_drivernet/', 'significant_genes_DriverNet.csv'))
  dim(ab_df)
  ab_df <- ab_df %>%
    dplyr::filter(`p-value`<0.05)
  annot_genes <- unique(ab_df$gene_symbol)
  ab_de <- data.table::fread(paste0(input_dir, subtag,'/',subtag, '_DE_genes.csv.gz'))
  dim(ab_de)
  ab_de <- ab_de %>%
    dplyr::filter(pvalue<0.05 & abs(log2FoldChange)>=0.5)
  
  pab <- plot_DE_genes_ggplot(ab_de, topGenes=annot_genes, capstr='', 
                       # FDRcutoff=0.01, logFCcutoff=0.25, pValuecutoff=0.05,
                       plttitle=plt_title, save_fig_dir)
  # pab
  
  obs_clones <- c('B','C')
  subtag <- paste0(obs_clones[2],'met_',obs_clones[1],'pri')
  plt_title <- paste0('Clone ',obs_clones[2],' metastasis vs. clone ',obs_clones[1],' primary')
  drivernet_df <- data.table::fread(paste0(input_dir, subtag,'_drivernet/', 'significant_genes_DriverNet.csv'))
  drivernet_df <- drivernet_df %>%
    dplyr::filter(`p-value`<0.05)
  dim(drivernet_df)
  # View(drivernet_df)
  annot_genes <- unique(drivernet_df$gene_symbol)
  
  pathway_df <- data.table::fread(paste0(input_dir, subtag,'/','pathways_',subtag, '.csv.gz'))
  pathway_df <- pathway_df %>%
    dplyr::filter(p_value<0.05)
  dim(pathway_df)
  
  topGenes <- c()
  for(r in unique(pathway_df$reference_set)){
    pw_genes <- pathway_df %>%
      dplyr::filter(reference_set==r) %>%
      dplyr::pull(signif_genes)
    pw_genes <- as.character(unlist(strsplit(pw_genes, ',')))
    topGenes <- c(topGenes, pw_genes)
  }
  
 
  
  length(annot_genes)
  
  de_df <- data.table::fread(paste0(input_dir, subtag,'/',subtag, '_DE_genes.csv.gz'))
  de_df <- de_df %>%
    dplyr::filter(pvalue<0.05 & abs(log2FoldChange)>=0.5)
  dim(de_df)
  
  tmp <- de_df %>%
    dplyr::filter(symbol %in% topGenes) %>%
    dplyr::arrange(-abs(log2FoldChange))
  # print(tmp)
  topGenes <- tmp$symbol[1:9]
  annot_genes <- c(annot_genes, topGenes)
  sum(annot_genes %in% de_df$symbol)
  pbc <- plot_DE_genes_ggplot(de_df, topGenes=annot_genes, capstr='', 
                              # FDRcutoff=0.01, logFCcutoff=0.25, pValuecutoff=0.05,
                              plttitle=plt_title, save_fig_dir)
  
  # pbc
  
  obs_clones <- c('A','C')
  subtag <- paste0(obs_clones[2],'met_',obs_clones[1],'pri')
  plt_title <- paste0('Clone ',obs_clones[2],' metastasis vs. clone ',obs_clones[1],' primary')
  drivernet_df <- data.table::fread(paste0(input_dir, subtag,'_drivernet/', 'significant_genes_DriverNet.csv'))
  drivernet_df <- drivernet_df %>%
    dplyr::filter(`p-value`<0.05)
  dim(drivernet_df)
  annot_genes <- unique(drivernet_df$gene_symbol)
  de_df <- data.table::fread(paste0(input_dir, subtag,'/',subtag, '_DE_genes.csv.gz'))
  
  de_df <- de_df %>%
    dplyr::filter(pvalue<0.05 & abs(log2FoldChange)>=0.5)
  dim(de_df)
  
  pac <- plot_DE_genes_ggplot(de_df, topGenes=annot_genes, capstr='', 
                              # FDRcutoff=0.01, logFCcutoff=0.25, pValuecutoff=0.05,
                              plttitle=plt_title, save_fig_dir)
  
  # pac
  
  
  # Global change comparison
  
  obs_clones <- c('A','A')
  subtag <- paste0(obs_clones[2],'met_',obs_clones[1],'pri')
  plt_title <- paste0('Clone ',obs_clones[2],' metastasis vs. clone ',obs_clones[1],' primary')
  pathway_df <- data.table::fread(paste0(input_dir, subtag,'/','pathways_',subtag, '.csv.gz'))
  pathway_df <- pathway_df %>%
    dplyr::filter(p_value<0.05)
  dim(pathway_df)
  annot_genes <- c()
  for(r in unique(pathway_df$reference_set)){
    pw_genes <- pathway_df %>%
      dplyr::filter(reference_set==r) %>%
      dplyr::pull(signif_genes)
    pw_genes <- as.character(unlist(strsplit(pw_genes, ',')))
    annot_genes <- c(annot_genes, pw_genes)
  }
  length(annot_genes)
  de_df <- data.table::fread(paste0(input_dir, subtag,'/',subtag, '_DE_genes.csv.gz'))
  # de_df <- de_df %>%
  #   dplyr::filter(pvalue<0.05 & abs(log2FoldChange)>=1)
  
  de_df <- de_df %>%
    dplyr::filter(pvalue<0.05 & abs(log2FoldChange)>=0.5)
  dim(de_df)
  
  paa <- plot_DE_genes_ggplot(de_df, topGenes=annot_genes, capstr='', 
                              # FDRcutoff=0.01, logFCcutoff=0.25, pValuecutoff=0.05,
                              plttitle=plt_title, save_fig_dir)
  
  # paa
  
  
  obs_clones <- c('B','B')
  subtag <- paste0(obs_clones[2],'met_',obs_clones[1],'pri')
  plt_title <- paste0('Clone ',obs_clones[2],' metastasis vs. clone ',obs_clones[1],' primary')
  pathway_df <- data.table::fread(paste0(input_dir, subtag,'/','pathways_',subtag, '.csv.gz'))
  pathway_df <- pathway_df %>%
    dplyr::filter(p_value<0.05)
  dim(pathway_df)
  annot_genes <- c()
  for(r in unique(pathway_df$reference_set)){
    pw_genes <- pathway_df %>%
      dplyr::filter(reference_set==r) %>%
      dplyr::pull(signif_genes)
    pw_genes <- as.character(unlist(strsplit(pw_genes, ',')))
    annot_genes <- c(annot_genes, pw_genes)
  }
  length(annot_genes)
  de_df <- data.table::fread(paste0(input_dir, subtag,'/',subtag, '_DE_genes.csv.gz'))
  # de_df <- de_df %>%
  #   dplyr::filter(pvalue<0.05 & abs(log2FoldChange)>=1)
  
  de_df <- de_df %>%
    dplyr::filter(pvalue<0.05 & abs(log2FoldChange)>=0.5)
  dim(de_df)
  
  pbb <- plot_DE_genes_ggplot(de_df, topGenes=annot_genes, capstr='', 
                              # FDRcutoff=0.01, logFCcutoff=0.25, pValuecutoff=0.05,
                              plttitle=plt_title, save_fig_dir)
  
  # pbb
  
  
  p_total <- cowplot::plot_grid(pab, pac, pbc, 
                     paa, pbb, NULL, ncol=2,
                     labels = c('A.','B.', 'C.', 'D.', 'E.',' ')) + 
              theme(plot.background = element_rect(fill = "white", colour = "white")) 
  p_total
  
  ggsave(  
    filename = paste0(save_fig_dir,"Supp_Fig6_volcano_plts.svg"),  
    plot = p_total,  
    height = 12, width = 9, dpi = 150)
}


get_annotated_genes <- function(de_df, nb_genes=15){
  # loading genes from drivernet output
  # loading met genes from HW annotation file 
  # loading genes from pathways
  
}

