suppressPackageStartupMessages({
  require(RColorBrewer)
  require(glue)
  require(ggrepel)
  require(tidyr)
  # require(fgsea)
  require(scales)
  require(dplyr)
  require(stringr)
  require(ggplot2)
})


# SA535
results_10x_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_total/deg_analysis/'
tag <- 'SA535_GF_D'
input_deg_fn <- paste0(results_10x_dir,'SA535_GF_treated_D/','total_markers.txt')
View(head(deg_df))
deg_df$GENEID
colnames(deg_df)[which(colnames(deg_df)=='GENEID')] <- 'ensembl_gene_id'
clone_str <- 'G,F vs D'
sample_name <- 'SA535_DE_resistance_vs_sensitive'


#SA609
input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA609_rna/deg_analysis/'

tag <- 'SA609_R_UTTTT_E'
tag1 <- 'SA609_R_E'
top_pathway_ext <- read.table(paste0(input_dir,'SA609_R_UTTTT_E/',pt_use,'top_up_pathway_cls_',tag1,'.txt'),
                              sep = '\t', header=T, check.names=F)

input_deg_fn <- paste0(input_dir,'SA609_R_UTTTT_E/','total_markers.txt')
clone_str <- 'R_UTTTT vs E_UUUUUUUU'
sample_name <- 'SA609_DE_resistance_vs_sensitive, R_UTTTT vs E_UUUUUUUU'


get_genes_common_pathways<- function(){
  obs_pathways <- c('tnfa_signaling_via_nfkb','tgf_beta_signaling',
                 'epithelial_mesenchymal_transition','uv_response_dn',
                 'mitotic_spindle')
  input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_total/deg_analysis/'
  tag <- 'SA535_GF_D'
  pt_use <- 'hallmark/'
  pattern_use <- 'HALLMARK_'
  # up_SA535_GF_D <- read.csv(paste0(input_dir,'SA535_GF_treated_D/',pt_use,'SA535_up_pathway.csv'),
  #                           row.names=1, check.names=F, stringsAsFactors=F)
  # View(up_SA535_GF_D)
  # pt_ls <- rownames(up_SA535_GF_D)
  # pt_ls <- str_replace_all(tolower(pt_ls), "^hallmark_", "")
  # rownames(up_SA535_GF_D)  <- pt_ls
  # up_SA535_GF_D <- up_SA535_GF_D[rownames(up_SA535_GF_D) %in% obs_pathways,]
  # View(head(up_SA535_GF_D))
  
  top_pathway_ext <- read.table(paste0(input_dir,'SA535_GF_treated_D/',pt_use,'top_up_pathway_cls_',tag,'.txt'),
                                sep = '\t', header=T, check.names=F)
  View(head(top_pathway_ext))
  top_pathway_ext$pathway <- str_replace_all(tolower(top_pathway_ext$pathway), "^hallmark_", "")
  top_pathway_ext <- top_pathway_ext[top_pathway_ext$pathway %in% obs_pathways,]
  # View(top_pathway_ext)
  genes_use <- c()
  for(idx in rep(1:nrow(top_pathway_ext),1)){
    
    gene_ls <- strsplit(as.character(top_pathway_ext$leadingEdge[[idx]]),' ')
    gene_ls <- unlist(gene_ls, use.names=FALSE)
    genes_use <- c(genes_use, gene_ls)
  }
  
  genes_use <- unique(genes_use)
  length(genes_use)
}
plot_DE_genes_chr <- function(input_deg_fn, clones=c('E','H'), sample_name='',
                              additional_genes = NULL, n_genes_to_annotate = 25){
  
  # input_deg_fn <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_cys/clonealign/SA535X9XB03616/deg/SA535X9XB03616_EH/seurat_DEG_total_markers.txt'
  deg_df <- read.table(input_deg_fn, check.names=F, stringsAsFactors=F)
  dim(deg_df)
  deg_df <- deg_df[deg_df$p_val_adj<0.05,]
  dim(deg_df)
  colnames(deg_df)[which(colnames(deg_df)=='GENEID')] <- 'ensembl_gene_id'
  # deg_df <- deg_df %>% rownames_to_column("ensembl_gene_id")
  colnames(deg_df)[which(names(deg_df) == "gene_symb")] <- "gene_symbol"
  # clones <- c('E','H')
  clone_str <- paste(clones, collapse = " vs ")
  print(clone_str)
  save_dir <- paste0(dirname(input_deg_fn),'/')
  if (!file.exists(save_dir)){
    dir.create(save_dir)
  }
  chrs <- c(as.character(1:23), "X", "Y")
  df_track <- annotables::grch37 %>% 
    dplyr::select(gene_symbol = symbol, chr, start, end) %>% 
    inner_join(deg_df) %>% 
    dplyr::filter(chr %in% c(as.character(1:23), "X", "Y")) %>% 
    dplyr::mutate(chr = factor(chr, levels = chrs),
                  position = (start + end) / 2)
  
  # MA: sometimes a genes appears more than once, keep only the unique ones
  df_track <- distinct(df_track)
  
  cols <- rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu") )
  #-log10(p_val_adj)
  track_plot <- ggplot(df_track, aes(x = position, y = avg_logFC)) +
    geom_point(aes(colour = p_val_adj), size = 2) +   # MA: size was 1
    facet_wrap(~ chr, scales = "free_x", 
               strip.position = "bottom",
               # switch = "x", 
               nrow = 1) +
    scale_colour_gradientn(name = glue("pval_adj, {clone_str}"),
                           colours = cols, 
                           values = rescale(c(0, 0.03, 0.05)),
                           limits = c(0, 0.05)) +
    theme(strip.background = element_rect(fill = 'white', colour = 'white'),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size=14),        
          strip.placement = "outside",
          legend.position = "top",
          panel.spacing = unit(0.1, 'cm')) +   ## MA: was 0.2
    theme(text = element_text(size = 12)) +
    labs(x = "Chromosome",
         subtitle = sample_name) +
    scale_x_continuous(expand = c(0,0))
  
  # df_annot <- top_n(df_track, n_genes_to_annotate, -log10(p_val_adj)) 
  additional_genes <- genes_use
  length(additional_genes)
  df_annot <- df_track[df_track$gene_symbol %in% additional_genes,]
  dim(df_annot)
  # ext_genes <- setdiff(genes_use,df_annot$gene_symbol)
  # sum(!as.character(genes_use) %in% as.character(df_annot$gene_symbol))
  # genes_add <- setdiff(genes, colnames(df1))
  # sum(genes_use %in% genes_use)
  # length(unique(additional_genes))
  # length(ext_genes)
  # View(df_annot)
  # df_annot1 <- df_annot[df_annot$gene_symbol=='AMTN',]
  # View(df_annot1)
  print("The labeled genes")
  # print(df_annot)  
  # View(df_track[df_track$gene_symbol=='TPM2',])
  if(!is.null(additional_genes)) {
    df_annot <- df_annot  %>% 
      bind_rows(dplyr::filter(df_track, gene_symbol %in% additional_genes))
  }
  
  ## thhis is where the genes are added
  track_plot2 <- track_plot +
    geom_text_repel(data = df_annot, aes(label = gene_symbol), size = 2.3)
  
  # base_name <- paste(clones, collapse = "_vs_")
  base_name <- tag
  png(paste0(save_dir,"DE_chr_",base_name,".png"), height = 2*450, width=2*1200,res = 2*72)
  print(track_plot2)
  dev.off()
}

