suppressPackageStartupMessages({
  require(RColorBrewer)
  require(glue)
  require(ggrepel)
  require(tidyr)
  require(fgsea)
  # require(scales)
  require(ComplexHeatmap)
  # require(DOSE)
  # require(enrichplot)
  # require(ggraph)
  require(ggplot2)
  # require(igraph)
  require(tidyverse)
  require(data.table)
  require(Seurat)
  require(dplyr)
})


plot_DE_genes_chr <- function(input_deg_fn, clones=c('E','H'), sample_name='',
                              additional_genes = NULL, n_genes_to_annotate = 25){
 
  # input_deg_fn <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_cys/clonealign/SA535X9XB03616/deg/SA535X9XB03616_EH/seurat_DEG_total_markers.txt'
  deg_df <- read.table(input_deg_fn, check.names=F, stringsAsFactors=F)
  dim(deg_df)
  deg_df <- deg_df %>% rownames_to_column("ensembl_gene_id")
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
  
  cols <- rev(brewer.pal(n = 5, name = "RdBu") )
  
  track_plot <- ggplot(df_track, aes(x = position, y = -log10(p_val_adj))) +
    geom_point(aes(colour = avg_log2FC), size = 2) +   # MA: size was 1
    facet_grid(~ chr, scales = "free_x", space='free',  #facet_wrap
               strip.position = "bottom",
               # switch = "x", 
               nrow = 1) +
    scale_colour_gradientn(name = glue("log2 fold change, {clone_str}"),
                           colours = cols, 
                           values = scales::rescale(c(-2, -0.3, 0, 0.3, 2)),
                           limits = c(-3.61, 3.61)) +
    theme(strip.background = element_rect(fill = 'white', colour = 'white'),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size=14),        
          strip.placement = "outside",
          legend.position = "top",
          panel.spacing = unit(0.1, 'cm')) +   ## MA: was 0.2
    theme(text = element_text(size = 18)) +
    labs(x = "Chromosome",
         subtitle = sample_name) +
    scale_x_continuous(expand = c(0,0))
  
  df_annot <- top_n(df_track, n_genes_to_annotate, -log10(p_val_adj)) 
  
  
  print("The labeled genes")
  print(df_annot)  
  
  if(!is.null(additional_genes)) {
    df_annot <- df_annot  %>% 
      bind_rows(dplyr::filter(df_track, gene_symbol %in% additional_genes))
  }
  
  ## this is where the genes are added
  track_plot2 <- track_plot +
    geom_text_repel(data = df_annot, aes(label = gene_symbol), size = 3.0)
  
  base_name <- paste(clones, collapse = "_vs_")
  png(paste0(save_dir,"DE_chr_",base_name,".png"), height = 2*450, width=2*950,res = 2*72)
  print(track_plot2)
  dev.off()
}


# input_cnv_fn <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna_cys/clonealign/SA535X9XB03616/SA535X9XB03616.txt'
# df_cnv <- read.table(input_cnv_fn, check.names=F, stringsAsFactors=F, header=T)
# dim(df_cnv)
# colnames(df_cnv) should contain "ensembl_gene_id", 
# cluster is clone label, cnv is median copy number
# View(head(df_cnv))
plot_CNV <- function(df_cnv, clones, sample_name='SA535X9XB03616',
                     additional_genes = NULL, n_genes_to_annotate = 35){
  
  
  annots <- annotables::grch37
  
  annots <- dplyr::select(annots, ensembl_gene_id = ensgene, chr, start, end) 
  
  df_cnv <- inner_join(df_cnv, annots)
  
  df_cnv <- dplyr::select(df_cnv, ensembl_gene_id, chr, start) %>% 
    group_by(chr) %>% 
    dplyr::mutate(start_order = rank(start)) %>% 
    ungroup() %>% 
    inner_join(df_cnv)
  
  chr_levels <- c(as.character(1:23), "X", "Y")
  
  df_cnv$chr <- factor(df_cnv$chr, levels = chr_levels)
  
  df_cnv <- drop_na(df_cnv)
  
  # MA: 11 Apr 2020, setting the copy number colors that were used in the heatmap
  # TO COME BACK
  #print("data frame")
  #cnv_colors <- c('#4880B8', '#A7C9DF','#CCCCCC','#F5CE93','#ED9364','#D2553E','#A42116','#8B1A43','#CB3576','#D06CAD','#C196C4','#D0BAD8')
  #cnv_cols <- data.frame(cn=0:11,color<-cnv_cols)
  #colnames(cnv_cols) <- c("cn","color")
  #levels(cnv_cols$color) <- cnv_colors
  
  cnv_colors <- c('#4880B8', '#A7C9DF','#CCCCCC','#F5CE93','#ED9364','#D2553E','#A42116','#8B1A43','#CB3576','#D06CAD','#C196C4','#D0BAD8')
  
  cnv_cols <- c('0'='#4880B8', '1'='#A7C9DF','2'='#CCCCCC','3'='#F5CE93','4'='#ED9364',
                '5'='#D2553E','6'='#A42116','7'='#8B1A43','8'='#CB3576','9'='#D06CAD','10'='#C196C4','11'='#D0BAD8')
  #levels(cnv_cols) <- 0:11
  #  cnv_cols <- c("0" = "#2166ac",
  #                "1" = "#92c5de", 
  #                "2" = "grey80", 
  #                "3" = "#f4a582", 
  #                "4" = "#d6604d",
  #                "5" = "#b2182b",
  #                "6+" = "#67001f")
  
  df_cnv$cnv <- as.character(round(df_cnv$median_cnmode))
  levels(df_cnv$cnv) <- 0:(length(cnv_colors)-1)
  #print(levels(df_cnv$cnv))
  
  # MA: removing the 6+ restriction (clonealign still has that restriction though)
  # df_cnv$cnv[round(df_cnv$median_cnmode) >= 6] <- "6+"
  #!chr %in% c("X","Y")
  cnv_plot <- dplyr::filter(df_cnv, cluster %in% clones, !chr %in% c("X","Y")) %>% 
    ggplot(aes(x = start_order, y = cluster, fill = cnv)) +
    geom_raster() +
    # facet_wrap(~ chr, scales = "free_x", nrow = 1, switch = "x") +
    facet_grid(~ chr, scales = "free_x", switch = "x", space='free') +
    #theme(legend.position = "top", axis.text.x = element_blank()) +
    scale_y_discrete(expand = c(0, 0)) +
    #scale_fill_manual(values=cnv_colors, name = "Copy number", guide = 'legend',labels = 0:(length(cnv_colors)-1),drop=FALSE)  +
    #          theme(legend.position = "bottom") +   
    scale_fill_manual(values = cnv_cols, name = "Copy number") +  #, labels = 0:(length(cnv_cols)-1),drop=FALSE) +
    labs(x = "chromosome") +
    theme(strip.background = element_rect(fill = 'white', colour = 'white'),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size=14), 
          strip.placement = "outside",
          legend.position = "none",
          panel.spacing = unit(c(0.1), 'cm')) +   # MA: was 0.2
    theme(text = element_text(size = 18)) +
    labs(x = "Chromosome", y = "Clone")+
    guides(fill = guide_legend(nrow = 1)) +
    scale_x_continuous(expand = c(0,0))
  
  # cnv_plot
  main_plot <- cowplot::plot_grid(
    track_plot2 + theme(axis.title.x = element_blank(),
                        strip.text.x = element_blank()),
    cnv_plot,
    ncol = 1,
    rel_heights = c(2.5,1),
    align = 'v'
  )
  
  main_plot
}


calculate_DE_analysis_v2 <- function(base_name, meta_data, cells_use_g1, cells_use_g2, 
                                     feature_use="Grouping",
                                  srt, genes_map, groups_use, test_use='wilcox', 
                                  save_dir="", input_dir="",
                                  pAdjustThrs=0.05, minLogFC=0.25,
                                  nbtopup = 30, nbtopdown = 30, save_data=T, viz=F){
  
  # print(head(meta_data))
  print(dim(meta_data))
  # sub_dir <- paste0(base_name,"_",groups_use[1],"_",groups_use[2])
  sub_dir <- base_name
  save_dir_pw <- paste0(save_dir,sub_dir,"/")
  if (!file.exists(save_dir_pw)){
    dir.create(save_dir_pw)
  }
  cells_use <- c(cells_use_g1, cells_use_g2)
  # cells_use <- rownames(meta_data)
  cat("DE analysis ", file = paste0(save_dir_pw,"de_analysis_log.txt"))
  print(length(cells_use))
  cat(paste0("\n Observed clones: ",groups_use[1],' vs ',groups_use[2]," Nb obs cells: ",length(cells_use), 
             "  nb cells in ",groups_use[1], ":",length(cells_use_g1),
             "  nb cells in ",groups_use[2], ":",length(cells_use_g2)), file = paste0(save_dir_pw,"de_analysis_log.txt"), append = TRUE)
  
 
  # length(cells_use)
  srt2 <- srt[,as.character(cells_use)]  # subset cells
  dim(srt2)
  Idents(object = srt2) <- feature_use
  unique(srt2$library_id)
  # print(Idents(object = srt2))
  # rownames(genes_map) <- genes_map$gene_ens
  # Differential Expression Analysis 
  
  markers_ls <- Seurat::FindMarkers(srt2, slot="data",  # using normalized data   
                            ident.1=groups_use[1], ident.2=groups_use[2],    #Vector of cell names belonging to group 1, group 2
                            logfc.threshold=minLogFC, test.use=test_use)
  
  print(paste0("Total DE marker genes: ",nrow(markers_ls)))
  cat(paste0("\n Total DE marker genes: ",nrow(markers_ls)), file = paste0(save_dir_pw,"de_analysis_log.txt"), append = TRUE)
  markers_ls$gene_symb <- as.character(genes_map[rownames(markers_ls),"gene_symb"])
  markers_ls$GENEID <- rownames(markers_ls)
  # View(head(markers_ls))
  # markers_ls$avg_log2FC
  # Get top marker genes
  markers_ls_output <- markers_ls[markers_ls$p_val_adj<pAdjustThrs & abs(markers_ls$avg_log2FC)>minLogFC,]
  markers_ls_output <- markers_ls_output[order(markers_ls_output$avg_log2FC,decreasing = T),]  
  print(paste0("Nb selected marker genes: ",nrow(markers_ls_output)))
  cat(paste0("\n Nb selected marker genes: ",nrow(markers_ls_output)), file = paste0(save_dir_pw,"de_analysis_log.txt"), append = TRUE)
  # Get top up-regulated genes
  markers_ls_upreg <- markers_ls_output[markers_ls_output$avg_log2FC>minLogFC,]
  markers_ls_upreg <- markers_ls_upreg[order(markers_ls_upreg$avg_log2FC,decreasing = T),] 
  
  
  # dim(markers_ls_upreg)
  if(nrow(markers_ls_upreg) < nbtopup){
    nbtopup <- nrow(markers_ls_upreg)
  }
  markers_ls_upreg_top <- markers_ls_upreg[1:nbtopup,]
  
  # Get top down-regulated genes
  markers_ls_downreg <- markers_ls_output[markers_ls_output$avg_log2FC<(-minLogFC),]
  markers_ls_downreg <- markers_ls_downreg[order(markers_ls_downreg$avg_log2FC,decreasing = F),] 
  # dim(markers_ls_downreg)
  if(nrow(markers_ls_downreg) < nbtopdown){
    nbtopdown <- nrow(markers_ls_downreg)
  }
  markers_ls_downreg_top <- markers_ls_downreg[1:nbtopdown,]
  
  markers_ls_upreg_top$genes_deg <- 'up-regulated'
  markers_ls_downreg_top$genes_deg <- 'down-regulated'
  genes_deg_df <- rbind(markers_ls_upreg_top, markers_ls_downreg_top)
  
  if(save_data){
    # write.table(markers_ls_output,paste0(save_dir_pw,"de_significant_genes.txt"),sep="\t",col.names=NA,quote=F)
    data.table::fwrite(markers_ls_output,paste0(save_dir_pw,"de_significant_genes.csv.gz"))
    data.table::fwrite(markers_ls,paste0(save_dir_pw,"total_markers.csv.gz"))
    # write.table(markers_ls,paste0(save_dir_pw,"total_markers.txt"),sep="\t",col.names=NA,quote=F)
    # ,paste0(save_dir_pw,"de_significant_genes.txt")convert_de_csv(sub_dir, save_dir_pw)
    
    # write.table(genes_deg_df, paste0(save_dir_pw,"de_topup_",nbtopup,"_topdown_",nbtopdown,".txt"),sep="\t",col.names=NA,quote=F)
  }
    
  print("Plot DE genes")
  # plttitle <- paste0(base_name,"; ",groups_use[1]," vs. ",groups_use[2])
  plttitle <- base_name
  topGenes <- genes_deg_df$gene_symb
  # markers_ls_tmp <- markers_ls
  # rownames(markers_ls_tmp) <- markers_ls_tmp$gene_symb
  # plot_DE_genes(markers_ls_tmp, topGenes, nrow(markers_ls_output), 
  #               FCcutoff=minLogFC, pCutoff=pAdjustThrs, plttitle, save_dir_pw)
  # 
  # summary_genes_ls <- list(markers_ls=markers_ls,
                           # markers_ls_output=markers_ls_output,
                           # top_genes=genes_deg_df,
                           # markers_ls_upreg=markers_ls_upreg,
                           # markers_ls_downreg=markers_ls_downreg)
  # if(save_data){
  #   saveRDS(summary_genes_ls,file = paste0(save_dir_pw,'summary_genes_ls.rds'))
  # }
  
  p <- plot_DE_genes_ggplot(markers_ls, topGenes, capstr="", 
                       FDRcutoff=0.01, logFCcutoff=0.25, pValuecutoff=0.05,
                       plttitle, save_dir_pw, legendVisible=F,
                       iscaption=TRUE, legend_verbose='none', save_plot=TRUE, xl=NULL)
  # p  
  # Plot heatmap
  # norm_df <- GetAssayData(object = srt2, slot = "data")
  # ext <- rownames(genes_deg_df)[rownames(genes_deg_df) %in% rownames(norm_df)]
  # norm_df <- norm_df[ext,] # get correct order
  # rownames(norm_df) <- genes_deg_df[rownames(norm_df),'gene_symb']
  # clusters <- srt2@meta.data[colnames(norm_df), feature_use]
  # genes_cls <- genes_deg_df[ext,]$genes_deg
  # print("Plot heatmap")
  # p <- Heatmap(as.matrix(norm_df),show_column_names=FALSE, 
  #              cluster_rows=F,cluster_columns=F,
  #              column_split=clusters, row_split = genes_cls) 
  # plttitle_hm <- paste0(base_name,"_",groups_use[1],"_vs_",groups_use[2])
  # png(paste0(save_dir_pw,"heatmap_",plttitle_hm,".png"), height = 2*500, width=2*700,res = 2*72)
  # print(p)
  # dev.off()
  
  # Plot the top 20 reference genes cisplatin
  # gmt_dir <- paste0(input_dir,"biodatabase/")
  # cis_genes <- read.csv(paste0(gmt_dir,'cisplatin_20_genes.csv'), header = T, sep=',', check.names=F, stringsAsFactors=F)
  
  # markers_ls_output <- read.table(paste0(save_dir,'SA1035_E_D/de_significant_genes.txt'), header=T,sep='\t',row.names = 1)
  # head(markers_ls_output)
  
  
  # if(length(unique(markers_ls_output$gene_symb))!=nrow(markers_ls_output)){
  #   print("Duplicate gene symbol, remove duplicated genes from the mtx")
  #   deg_mtx <- markers_ls_output[!duplicated(markers_ls_output$gene_symb),]
  # } else{
  #   deg_mtx <- markers_ls_output
  # }
  # c_use <- markers_ls$gene_symb %in% cis_genes$gene
  # if(sum(c_use)>0){
  #   cat(paste0("\n Nb cells exist in top20 ref genes: ",sum(c_use)), file = paste0(save_dir_pw,"de_analysis_log.txt"), append = TRUE)
  #   norm_df <- GetAssayData(object = srt2, slot = "data")
  #   # dim(norm_df)
  #   norm_df <- norm_df[rownames(markers_ls),]
  #   norm_df <- norm_df[c_use,]
  #   rownames(norm_df) <- markers_ls[rownames(norm_df),'gene_symb']
  # 
  #   # norm_df <- norm_df[rownames(norm_df) %in% cis_genes$gene,]
  #   print(dim(norm_df))
  #   cell_clusters <- srt2@meta.data[colnames(norm_df), feature_use]
  # # # genes_cls <- rownames(norm_df) %in% cis_genes$gene
  # # # names(genes_cls) <- rownames(norm_df)
  #   # markers_ls$is_included_top20_ref <- markers_ls_output$gene_symb %in% cis_genes$gene
  #   
  #   # print("Plot heatmap")
  #   # p <- Heatmap(as.matrix(norm_df),show_column_names=F, show_row_names =T,
  #   #              cluster_rows=F,cluster_columns=F, column_split = cell_clusters)  #,row_split = genes_cls
  #   # ht <- nrow(norm_df)
  #   # png(paste0(save_dir_pw,"hm_ref_cisplatin_genes.png"), height = 2*40*ht+40, width=2*700,res = 2*72)
  #   # print(p)
  #   # dev.off()
  # } else{
  #   print("No gene match the reference 20 gene list")
  # }
  
  
  
  
  # Pathway analysis
  print("Pathway analysis....")
  # gmt_dir <- paste0(input_dir,"biodatabase/")
  gmt_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/pathway_set/'
  deg_stat <- markers_ls_output$avg_log2FC
  names(deg_stat) <- markers_ls_output$gene_symb
  # # print(head(deg_stat))
  
  # comment back
  method_use <- paste0("seurat_DEG_",test_use)
  gmt_ls <- c("h.all.v7.0.symbols.gmt") #,"c2.cp.kegg.v7.1.symbols.gmt","GO_c5.all.v7.1.symbols.gmt"
  # pathway_names <- c("hallmark","kegg") #,"go"
  pathway_names <- c("hallmark")
  
  gmtfile <- paste0(gmt_dir, gmt_ls[1])
  
  # de_genes <- markers_ls_upreg
  # if(viz){
  #   viz_pathway(gmt_dir, markers_ls_upreg, save_dir_pw, base_name, tag='up_pathway')
  #   viz_pathway(gmt_dir, markers_ls_downreg, save_dir_pw, base_name, tag='down_pathway')
  # }
    
  print("Get pathway")
  # if(!viz){
    # Deal with overlapping and NA gene symbols
    markers_ls <- markers_ls %>%
      dplyr::select(gene_symb, avg_log2FC) %>%
      na.omit() %>%
      distinct() %>%
      group_by(gene_symb) %>%
      summarize(avg_log2FC=mean(avg_log2FC))

    for(i in rep(1:length(gmt_ls),1)){
      # gmt_fn <- paste0(gmt_dir, gmt_ls[i])
      print(gmt_ls[i])
      cat(paste0("\n Pathway analysis, and using gene sets: ",gmt_ls[i]," pathway is: ",pathway_names[i]), file = paste0(save_dir_pw,"de_analysis_log.txt"), append = TRUE)
      print(pathway_names[i])
      markers_ls$avg_logFC <- markers_ls$avg_log2FC
      pathway_ls <- get_pathway_results(markers_ls, method_use,
                                        base_name, paste0(gmt_dir, gmt_ls[i]), 
                                        pathway_names[i],
                                        groups_use,
                                        save_dir_pw, 30)
    }

  # # return(pathway_ls)
  # }
  # print("Generate genes network")
  # generate_genes_network(save_dir_pw,
  #                        datatag=sub_dir,
  #                        de_desc=sub_dir,
  #                        genes_set='hallmark',
  #                        pattern_use='^hallmark_'
  # )
}



get_pathway_results <- function(markers_ls_output,                      # named vector of statistical significance 
                                method_use = "seurat_DEG_Wilcox",
                                base_name = '',    
                                gmt_fn="biodatabase/h.all.v7.0.symbols.gmt",
                                pathway_name='hallmark',
                                groups_use=c("UT","UU"),    # vector of 2 elements, 2 group name used for DE analysis
                                save_dir = "/home/htran/",
                                n_top=20){
  
  
  # library(fgsea)
  save_dir_log <- save_dir
  save_dir <- paste0(save_dir,pathway_name,'/')
  if(!dir.exists(save_dir)){
    dir.create(save_dir)
  }
  print(save_dir)
  # Can use deframe: first col as name, second col as value, deg_stat <- deframe(markers_ls_output)
  deg_stat <- markers_ls_output$avg_log2FC
  names(deg_stat) <- markers_ls_output$gene_symb   # TO DO: check duplicated genes, remove duplicated columns cause gene symbol are not unique in some case
  Hs.H <- gmtPathways(gmt_fn) 
 
  # Hs.c5 <- gmtPathways(paste0(gmt_dir, gmt_ls[2])) 
  # Hs.c6 <- gmtPathways(paste0(gmt_dir, gmt_ls[3])) 

  # pathways: List of gene sets to check, ex in data(examplePathways)
  # stats: Named vector of gene-level stats. Names should be the same as in ’pathways’
  # nperm: Number of permutations to do. Minimial possible nominal p-value is about 1/nperm
  gsea_hallmark <- fgsea(pathways=Hs.H, stats=deg_stat, nperm=10000)
  # gsea_go <- fgsea(pathways=Hs.c5, stats=deg_stat, nperm=10000)
  # gsea_oncogenic <- fgsea(pathways=Hs.c6, stats=deg_stat, nperm=10000)

  plottitle <- paste0(base_name,"_",groups_use[1],"_vs_",groups_use[2])
  
  
  topPathways5Up <- gsea_hallmark[(ES > 0) & (padj < 0.05)][head(order(padj), n=n_top), pathway]
  topPathways5Down <- gsea_hallmark[ES < 0 & (padj < 0.05)][head(order(padj), n=n_top), pathway]
  topPathways <- c(topPathways5Up, rev(topPathways5Down))
  save_fn = paste0("summary_pathway_cls_",base_name,"_",groups_use[1],"_",groups_use[2],".rds")
  save_fn_txt = paste0("pathway_cls_",base_name,"_",groups_use[1],"_",groups_use[2],".txt")
  
  cat(paste0("\n Number of up pathways: ",length(topPathways5Up)), file = paste0(save_dir_log,"de_analysis_log.txt"), append = TRUE)
  
  if(length(topPathways5Up)>0){
    
    # png(paste0(save_dir,"up_pathway_cls_",base_name,"_",groups_use[1],"_",groups_use[2],".png"), height = 2*500, width=2*900,res = 2*72)
    #   plotGseaTable(Hs.H[topPathways5Up], deg_stat, gsea_hallmark, gseaParam = 0.5)
    # dev.off()
    
    topPathwaysUp_df <- gsea_hallmark[match(topPathways5Up, pathway)]
    pathway_barplot(topPathwaysUp_df, plottitle, xtitle="padj", save_dir, typePth='upreg', x="Count", color='pval', showCategory=n_top)
    fwrite(gsea_hallmark[match(topPathways5Up, pathway)], file=paste0(save_dir, "top_up_",save_fn_txt), sep="\t", sep2=c("", " ", ""))
    
    top_pathway_ext <- topPathwaysUp_df[topPathwaysUp_df$padj < 0.05,]
    # pw_df_up <- get_hm_data(top_pathway_ext, markers_ls_output, save_dir, base_name, 'up_pathway')
    
  }
  
  
  cat(paste0("\n Number of down pathways: ",length(topPathways5Down)), file = paste0(save_dir_log,"de_analysis_log.txt"), append = TRUE)
  
  if(length(topPathways5Down)>0){
    # print("Contain down-pathway")
    # print(length(topPathways5Down))
    # png(paste0(save_dir,"down_pathway_cls_",base_name,"_",groups_use[1],"_",groups_use[2],".png"), height = 2*500, width=2*900,res = 2*72)
    #   plotGseaTable(Hs.H[topPathways5Down], deg_stat, gsea_hallmark, gseaParam = 0.5)
    # dev.off()
    topPathwaysDown_df <- gsea_hallmark[match(topPathways5Down, pathway)]
    pathway_barplot(topPathwaysDown_df, plottitle, xtitle="padj", save_dir, typePth='downreg', x="Count", color='pval', showCategory=n_top)
    fwrite(gsea_hallmark[match(topPathways5Down, pathway)], file=paste0(save_dir, "top_down_",save_fn_txt), sep="\t", sep2=c("", " ", ""))  
   
    top_pathway_ext <- topPathwaysDown_df[topPathwaysDown_df$padj < 0.05,]
    print(dim(top_pathway_ext))
    # pw_df_down <- get_hm_data(top_pathway_ext, markers_ls_output, save_dir, base_name, 'down_pathway')
    
  }  
  if(length(topPathways)>0){
    # png(paste0(save_dir,"summary_pathway_cls_",base_name,"_",groups_use[1],"_",groups_use[2],".png"), height = 2*800, width=2*900,res = 2*72)
    #   plotGseaTable(Hs.H[topPathways], deg_stat, gsea_hallmark, gseaParam = 0.5)
    # dev.off()
    
    # collapsedPathways <- collapsePathways(gsea_hallmark[order(pval)][padj < 1],   #padj < 0.01
    #                                       Hs.H, deg_stat)
    # mainPathways <- gsea_hallmark[pathway %in% collapsedPathways$mainPathways][
    #   order(-NES), pathway]
    # if(length(mainPathways)>0){
    #   png(paste0(save_dir,"summary_collapse_pathway_cls_",base_name,"_",groups_use[1],"_",groups_use[2],".png"), height = 2*800, width=2*900,res = 2*72)
    #   plotGseaTable(Hs.H[mainPathways], deg_stat, gsea_hallmark, 
    #                 gseaParam = 0.5)
    #   dev.off()  
    #   fgseaResMain <- gsea_hallmark[match(mainPathways, pathway)]
    #   # fgseaResMain[, leadingEdge := lapply(leadingEdge, mapIds, x=org.Hs.eg.db, keytype="ENTREZID", column="SYMBOL")]
    #   fwrite(fgseaResMain, file=paste0(save_dir, "main_",save_fn_txt), sep="\t", sep2=c("", " ", ""))
    # }
    
    topPathways_df <- gsea_hallmark[match(topPathways, pathway)]
    
    
    # pathway_barplot(topPathways_df, plottitle, xtitle="padj", save_dir, typePth='total', x="Count", color='pval', showCategory=n_top)
    
    # pathway_ls <- list(
    #   stat = deg_stat,
    #   topPathwaysUp = topPathways5Up,
    #   topPathwaysDown = topPathways5Down,
    #   gsea_hallmark = gsea_hallmark,
    #   method = method_use,
    #   base_name = base_name, 
    #   groups = groups_use
    # )
    
    # library(data.table)
    fwrite(gsea_hallmark, file=paste0(save_dir, save_fn_txt), sep="\t", sep2=c("", " ", ""))
    # Readable version
    # library(org.Hs.eg.db)
    fwrite(gsea_hallmark[match(topPathways, pathway)], file=paste0(save_dir, "summary_topPathways_",save_fn_txt), sep="\t", sep2=c("", " ", ""))
    
    # saveRDS(pathway_ls, file = paste0(save_dir, save_fn))
    # print(paste0("Save pathway ",save_fn," in the directory ",save_dir))
    # return(pathway_ls)
    
  } else{
    print("Do not exist any pathway!!!")
    return(FALSE)
  }
  
}


convert_de_csv <- function(tag, save_dir){
  marker_fn <- paste0(save_dir,'total_markers.txt')
  print(marker_fn)
  if(file.exists(marker_fn)){
    print('Converting..')
    markers_df <- read.table(marker_fn, 
                             row.names=1, stringsAsFactors=F, check.names=F)
    print(dim(markers_df))
    # View(head(markers_df))
    # colnames(markers_df)
    markers_df <- markers_df[,colnames(markers_df) %in% c('p_val','avg_log2FC','p_val_adj','gene_symb','GENEID')]
    colnames(markers_df)[which(colnames(markers_df) == "gene_symb")] <- "gene_symbol"
    colnames(markers_df)[which(colnames(markers_df) == "avg_log2FC")] <- "logFC"
    colnames(markers_df)[which(colnames(markers_df) == "p_val_adj")] <- "FDR"
    markers_df <- markers_df %>% 
      dplyr::select(gene_symbol, logFC) %>% 
      na.omit() %>% 
      distinct() %>% 
      group_by(gene_symbol) %>% 
      summarize(logFC=mean(logFC))
    
    
    write.csv(markers_df, file=paste0(save_dir,'total_markers_',tag,'.csv'),row.names=F, quote=F)
    
  }
  
}
viz_pathway <- function(gmt_dir, de_genes, save_dir, base_name, tag='up_pathway'){
  library(DOSE)
  library(gdata)
  
  
  # Read pathway reference set
  # gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
  h <- read.gmt(paste0(gmt_dir,'h.all.v7.0.entrez.gmt'))
  
  
  # Input of enricher is gene id
  # Read genes and avg_log2FC from file
  # convert it to entrezid genes
  # de_genes$GENEID <- rownames(de_genes)
  c <- 'GENEID'
  if(!c %in% colnames(de_genes)){
    de_genes$GENEID <- rownames(de_genes)
  }
  
  genes_map <- convert_genes_set(rownames(de_genes), fromSet="GENEID", toSet=c("ENTREZID")) 
  print(dim(genes_map))
  genes_map <- genes_map[!is.na(genes_map$ENTREZID),]
  de_genes <- de_genes %>% inner_join(genes_map, by='GENEID')
  rownames(de_genes) <- as.character(de_genes$ENTREZID)
  
  
  
  # dup <- duplicated2(de_genes_ext$ENTREZID)
  de_genes_ext <- de_genes
  # de_genes_ext <- de_genes_ext[!is.na(de_genes_ext$ENTREZID),]
  # rownames(de_genes_ext) <- de_genes_ext$ENTREZID
  # de_genes_ext <- de_genes_ext[de_genes_ext$avg_log2FC>=0.5,]
  de_genes_ext <- de_genes_ext[,c('avg_log2FC'), drop=F]
  # View(head(de_genes_ext))
  
  # Enrichment
  genes_ls <- de_genes_ext$avg_log2FC
  names(genes_ls) <- rownames(de_genes_ext)
  egmt <- clusterProfiler:::enricher(rownames(de_genes), TERM2GENE=h)
  
  # gseaplot(egmt, geneSetID=1)
  # cnetplot(egmt)
  # barplot(egmt, showCategory=20)
  # d <- dotplot(egmt)
  edox <- setReadable(egmt, 'org.Hs.eg.db', 'ENTREZID')
  
  # Gene Concept Network
  # p1 <- cnetplot(edox, foldChange=genes_ls)
  # p1
  
  ## categorySize can be scaled by 'pvalue' or 'geneNum'
  p_gc1 <- cnetplot(edox, categorySize="pvalue", foldChange=genes_ls)
  p_gc2 <- cnetplot(edox, foldChange=genes_ls, circular = TRUE, colorEdge = TRUE)
  png(paste0(save_dir, base_name,'_', tag,'_network_concept1_plt.png'), height = 2*700, width=2*1200,res = 2*72)
  print(p_gc1)
  dev.off()
  
  png(paste0(save_dir, base_name,'_', tag,'_network_concept2_plt.png'), height = 2*650, width=2*1000,res = 2*72)
  print(p_gc2)
  dev.off()
  # Heatmap plot
  p_hm <- heatplot(edox, foldChange=genes_ls)
  png(paste0(save_dir, base_name,'_', tag,'_hm_pathway_plt.png'), height = 2*550, width=2*1.45*length(genes_ls),res = 2*72)
  print(p_hm)
  dev.off()
  print("Enrichment map plot")
  p_enr <- emapplot(edox, pie_scale=1.5) #,layout="kk"
  png(paste0(save_dir, base_name,'_', tag,'_enrichment_map_plt.png'), height = 2*400, width=2*500,res = 2*72)
  print(p_enr)
  dev.off()
}


get_hm_data <- function(top_pathway_ext, logFC_df, save_dir,base_name, tag='pathway'){
  top_pathway_ext <- setDF(top_pathway_ext)
  logFC_df <- logFC_df %>% 
    dplyr::select(gene_symb, avg_log2FC) %>% 
    na.omit() %>% 
    distinct() %>% 
    group_by(gene_symb) %>% 
    summarize(avg_log2FC=mean(avg_log2FC))
  logFC_df <- as.data.frame(logFC_df)
  logFC_df <- column_to_rownames(logFC_df, var = "gene_symb")
  # rownames(logFC_df) <- logFC_df$gene_symb
  rownames(top_pathway_ext) <- as.character(top_pathway_ext$pathway)
  pw_ls <- list()
  genes_use <- c()
  pathway <- c()
  for(idx in rep(1:nrow(top_pathway_ext),1)){
    
    p <- as.character(top_pathway_ext$pathway[idx])
    gene_ls <- strsplit(as.character(top_pathway_ext$leadingEdge[[idx]]),' ')
    # gene_ls <- top_pathway_ext$leadingEdge[[idx]]
    gene_ls <- unlist(gene_ls, use.names=FALSE)
    genes_use <- c(genes_use, gene_ls)
    # p  <- str_replace_all(tolower(p), "^hallmark_", "")
    # p  <- str_replace_all(p, "_pathway$", "")
    # p <- paste0(base_name,'_',p)
    pathway <- c(pathway, p)
    pw_ls[[p]] <- gene_ls
  }
  
  genes_use <- unique(genes_use)
  length(genes_use)
  
  # pathway <- str_replace_all(tolower(top_uppw$pathway), "^hallmark_", "")
  pw_df <- as.data.frame(matrix(NA, nrow = length(top_pathway_ext$pathway), 
                                ncol =length(genes_use),
                                dimnames = list(pathway, genes_use)))
  
  for(p in names(pw_ls)){
    gs <- pw_ls[[p]]
    pw_df[p,gs] <- logFC_df[gs,'avg_log2FC']
  }
  
  print(dim(pw_df))
  write.csv(pw_df, file = paste0(save_dir,base_name,'_',tag,'.csv'), row.names=T, quote=F)
  hm <- ComplexHeatmap::Heatmap(as.matrix(pw_df), na_col = "white",
                                   show_column_names=T, 
                                   show_row_names = T,
                                   cluster_rows=F,cluster_columns=F,
                                   name = tag, 
                                   column_names_gp = grid::gpar(fontsize = 7), 
                                   row_names_gp = grid::gpar(fontsize = 4)) 
  sh <- dim(pw_df)[1]
  png(paste0(save_dir, base_name,'_',tag,'_HM.png'), height = 2*30*sh+100, width=2*1600,res = 2*72)
  print(hm)
  dev.off()
  
  return(pw_df)
}

plot_DE_genes <- function(df, topGenes, nbgenessig=0, FCcutoff=0.25, pCutoff=0.01, 
                          plttitle="A versus B", save_dir="",legendVisible=F){
  
  # library(EnhancedVolcano)
 
  # colnames(df)[which(names(df) == "avg_log2FC")] <- "log2FoldChange"
  # colnames(df)[which(names(df) == "p_val_adj")] <- "padj"
  capstr <- paste0("FC cutoff, ",FCcutoff,"; p-value cutoff, ",pCutoff, "; nb genes signif, ",nbgenessig)
  p <- EnhancedVolcano::EnhancedVolcano(df,
                        lab = df$gene_symb,   #df$symbol or rownames(df)
                        x = 'avg_log2FC',
                        y = 'p_val_adj',
                        # selectLab = as.character(topGenes),
                        selectLab = topGenes,
                        xlim=c(-2.5,2.5),
                        xlab = bquote(~Log ~ 'fold change'),
                        ylab = bquote(~-Log[10]~adjusted~italic(P)),
                        pCutoff = pCutoff,
                        FCcutoff = FCcutoff,
                        labSize = 4.0,
                        colAlpha = 0.9,
                        # cutoffLineType = 'blank',
                        # cutoffLineCol = 'black',
                        # cutoffLineWidth = 0.5,
                        # hline = c(0.01, 0.001),
                        # hlineCol = c('grey0', 'grey25'),
                        # hlineType = 'longdash',
                        legend=c('NS','Log2 FC','Adjusted p-value',
                                 'Adjusted p-value & Log2 FC'),
                        legendPosition = 'bottom',
                        legendLabSize = 10,
                        legendIconSize = 3.0,
                        legendVisible = legendVisible,
                        col=c('black', 'black', 'black', 'red3'),
                        title = plttitle,
                        subtitle = " ",
                        caption = capstr,
                        drawConnectors = TRUE,
                        widthConnectors = 0.23,
                        colConnectors = 'grey30')
  
  png(paste0(save_dir,"DE_",plttitle,".png"), height = 2*500, width=2*600,res = 2*72)
    print(p)
  dev.off()
  
  # return(p)
}


plot_DE_genes_edgeR <- function(df, topGenes, capstr='', FDRcutoff=0.01, logFCcutoff=0.25, 
                          plttitle="A versus B", save_dir="",legendVisible=F,
                          iscaption=TRUE){
  
  # library(EnhancedVolcano)
  
  # colnames(df)[which(names(df) == "avg_log2FC")] <- "log2FoldChange"
  # colnames(df)[which(names(df) == "p_val_adj")] <- "padj"
  # capstr <- paste0("FDR cutoff, ",FDRcutoff,"; logFC cutoff, ",logFCcutoff, "; nb genes signif, ",nbgenessig)
  # summary(as.factor(df$gene_type))
  Gene_Type=c('In_cis_Decrease_DownRegulated','In_cis_Decrease_UpRegulated',
              'In_cis_Increase_DownRegulated','In_cis_Increase_UpRegulated',
              'In_trans_DownRegulated','In_trans_UpRegulated'
              )
  gt <- data.frame(Gene_Type=Gene_Type,
                   gene_type=c('In-cis_DD','In-cis_DU',
                               'In-cis_ID','In-cis_IU',
                               'In-trans_D','In-trans_U'),
                   col_gt=c('#E61367','#D35B53',
                            '#EF816E','#E82503',
                            '#7BEADF','#093CB4'))
  
  df <- df %>% left_join(gt, by='Gene_Type')
  # unique(df$col_gt)
  keyvals_colour <- as.character(df$col_gt)
  names(keyvals_colour) <- df$gene_type
  # keyvals_colour <- factor(keyvals_colour, levels = unique(keyvals_colour))
  # unique(keyvals_colour)
  # names(keyvals_colour[1:3])
  st <- summary(df$gene_type)
  keyvals_shape <- ifelse(df$is_fitness_gene==T, 3, 16)
  names(keyvals_shape) <- ifelse(df$is_fitness_gene==T,'Fitness genes','ours')
  
  # df$gene_type <- factor(df$gene_type, levels = unique(df$gene_type))
  if(capstr=='' & iscaption){
    capstr <- paste0(capstr,names(st[1]),':',as.numeric(st[1]), ', ')
    capstr <- paste0(capstr,names(st[2]),':',as.numeric(st[2]), ', ')
    capstr <- paste0(capstr,names(st[3]),':',as.numeric(st[3]), ', ')
    capstr <- paste0(capstr,names(st[4]),':',as.numeric(st[4]), ' \n')
    capstr <- paste0(capstr,names(st[5]),':',as.numeric(st[5]), ', ')
    capstr <- paste0(capstr,names(st[6]),':',as.numeric(st[6]), ' ')
    # for(i in rep(1:length(st),1)){
    #   # print(st[i])
    #   capstr <- paste0(capstr,names(st[i]),':',as.numeric(st[i]), ' ')
    # }
    
  }
  df$mlog10FDR <- -log10(df$FDR)
  
  df$mlog10FDR <- sapply(df$mlog10FDR, function(x) replace(x, is.infinite(x), 300))
  if(legendVisible){
    lpos <- "right"
  }else{
    lpos <- "none"
  }
  p <- EnhancedVolcano::EnhancedVolcano(df,
                                        lab = df$gene_symbol,   #df$symbol or rownames(df)
                                        x = 'logFC',
                                        y = 'FDR',
                                        # selectLab = as.character(topGenes),
                                        selectLab = topGenes,
                                        xlim=c(-4.3,4.3),
                                        ylim = c(0, max(df$mlog10FDR, na.rm = TRUE)),
                                        xlab = bquote(~Log[2] ~ ' fold change'),
                                        ylab = bquote(~-Log[10]~italic(FDR)),
                                        pCutoff = logFCcutoff,
                                        FCcutoff = FDRcutoff,
                                        # pointSize = 3.2,
                                        pointSize = c(ifelse(df$is_fitness_gene==T, 4, 3)),
                                        labSize = 3,
                                        colAlpha = 0.7,
                                        gridlines.major = FALSE,
                                        gridlines.minor = FALSE,
                                        colCustom = keyvals_colour,
                                        shapeCustom = keyvals_shape,
                                        # cutoffLineType = 'blank',
                                        # cutoffLineCol = 'black',
                                        # cutoffLineWidth = 0.5,
                                        # hline = c(0.01, 0.001),
                                        # hlineCol = c('grey0', 'grey25'),
                                        # hlineType = 'longdash',
                                        legend=NULL,
                                        # legend=c('NS','Log2 FC','Adjusted p-value',
                                        #          'Adjusted p-value & Log2 FC'),
                                        legendPosition = lpos,
                                        legendLabSize = 8,
                                        legendIconSize = 6.0,
                                        legendVisible = legendVisible,
                                        # col=c('black', 'black', 'black', 'red3'),
                                        col=unique(keyvals_colour),
                                        title = plttitle,
                                        subtitle = NULL,
                                        caption = capstr,
                                        drawConnectors = TRUE,
                                        widthConnectors = 0.15,
                                        colConnectors = 'grey30')
  # p
  # p <- p + theme(legend.position="none")
  
  
  png(paste0(save_dir,"DE_",plttitle,".png"), height = 2*550, width=2*650,res = 2*72)
  print(p)
  dev.off()
  
  return(p)
}



pathway_barplot <- function(df, plottitle="Pathway", xtitle="pval", save_dir="",
                            typePth='upreg', x="Count", color='padj', showCategory=10) {
  # library(enrichplot)
  # library(DOSE)
  ## source code based on barplot function in enrichplot package 
  ## See more at https://bioconductor.org/packages/release/bioc/html/enrichplot.html
  if(nrow(df)>showCategory){
    df <- df[1:showCategory,]  
  }
  
  
  colorBy <- match.arg(color, c("pval", "p.adjust", "qvalue","padj"))
  if (x == "geneRatio" || x == "GeneRatio") {
    x <- "GeneRatio"
  } else if (x == "count" || x == "Count") {
    x <- "Count"
  } else{
    x <- "Count"
  }
  
  # change x to log fold change here 
  if(typePth=="upreg"){
    lowC <- "#006400"
    highC <- "#90EE90"
    plottitle <- paste0(plottitle,"_up-regulated")
  } else if(typePth=="downreg"){
    lowC <- "#4B0082"
    highC <- "#EE82EE"
    plottitle <- paste0(plottitle,"_down-regulated")
  } else{
    lowC <- "red"
    highC <- "blue"
    plottitle <- paste0(plottitle,"_total")
  }
  # df <- fortify(object, showCategory=showCategory, by=x, ...)
  if(!is.data.frame(df)){
    df <- as.data.frame(df)
  }
  # remove pathway_ str at the beginning of each description, change the text from capital to lower character
  # pathway_ls <- df$pathway
  # pathway_ls2 <- c()
  # for(p in pathway_ls){
  #   p <- str_sub(p, 10, str_length(p))
  #   pathway_ls2 <- c(pathway_ls2, p)
  # }
  # df$pathways <- pathway_ls2
  
  
  if(colorBy %in% colnames(df)) {
    p <- ggplot(df, aes(x = reorder(pathway,-padj), y = padj, fill = pval)) +
      DOSE::theme_dose(12) +
      # scale_fill_continuous(type = "gradient", name = color, guide=guide_colorbar(reverse=TRUE))
      scale_fill_continuous(low=lowC, high=highC, name = color, guide=guide_colorbar(reverse=TRUE)) 
  } else {
    p <- ggplot(df, aes_string(x = "pathway", y = xtitle)) + theme_dose(12) 
  }
  p <- p + geom_bar(stat = "identity", width=0.4) + coord_flip() 
    # ggtitle(title) + xlab(NULL) + ylab(xtitle)
  p <- p + labs(x="", y=xtitle, title=plottitle)
  p <- p + theme(plot.title = element_text(color="black", size=11, hjust = 0.5),
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                 axis.text.y = element_text(color="black", size=9),
                 axis.text.x = element_text(color="black", size=7, angle=90),
                 axis.title = element_text(color="black", size=15))
  
  
  png(paste0(save_dir,"pathway_barplot_",plottitle,".png"), height = 2*25*nrow(df)+100, width=2*600,res = 2*72)
  print(p)
  dev.off()
  
  return(p)
}



remove_bad_genes <- function(ens_genes){
  genes_map <- convert_genes_set(ens_genes, fromSet="GENEID", toSet=c("SYMBOL")) 

  bad_genes <- grepl("^MT-|^RP[L|S]", genes_map$SYMBOL)
  rv_genes <- rownames(genes_map)[bad_genes]
  genes_map <- genes_map[rv_genes,]    
  # genes_symb <- genes_map$SYMBOL[bad_genes]
  # genes_map_2 <- convert_genes_set(genes_symb, fromSet="SYMBOL", toSet=c("GENEID")) 
  print(paste0("Nb bad ens genes: ",length(rv_genes)))
  return(genes_map)   # return the bad genes

}


convert_genes_set <- function(genes_ls, fromSet="GENEID", toSet=c("ENTREZID", "SYMBOL")){

  library(EnsDb.Hsapiens.v79)
  # library(dplyr)
  # ls <- c("ENSG00000148908", "ENSG00000072506", "ENSG00000117394", "ENSG00000150753")
  genes_map <- ensembldb::select(EnsDb.Hsapiens.v79, key=genes_ls, columns=toSet, keytype=fromSet)
  
  # genes_map <- select(EnsDb.Hsapiens.v79, key=ls, columns="SYMBOL", keytype="GENEID")
  genes_map <- as.data.frame(genes_map) 
  genes_map <- dplyr::distinct(genes_map, GENEID, .keep_all= TRUE)   # remove duplicate
  rownames(genes_map) <- genes_map$GENEID

  map_gene <- length(genes_ls) - nrow(genes_map)
  if(map_gene==0){   # mapping 1:1
    print("1:1 mapping genes set")
  } else if(map_gene<0){  # 1 ensemble gene id correspond to multiple gene symbol, get the first appearance only
    print("Check the mapping results")
  } else{
    print(paste0("Could not find the mapping for ",map_gene," genes"))
  }
  
  return(genes_map)  

  # Second method
  # gene_symb <- bitr(rownames(sce), fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db", drop=FALSE)
  # sum(isUnique(gene_symb$ENSEMBL)==FALSE)
  # Third method
  # annots <- select(org.Hs.eg.db, keys=as.vector(rownames(sce)), columns="SYMBOL", keytype="ENSEMBL")
}
## get_metadata
# meta_data: dataframe with cells name as row name 
# save_dir: directory where to save output

get_metadata <- function(meta_data){
  
  print(colnames(meta_data))
  if(!is.data.frame(meta_data)){
    meta_data <- as.data.frame(meta_data)
  }
  # Reference, our definition based on biological experiment
  # Can detect pattern and automatically set
  # Can load from a sample file, do later
  treat_timepoint_map <- data.frame(ts=c("U","UT","UU","UUU","UTT","UTU","UUUU","UTTT","UTTU","UUUUU","UTTTT","UTTTU"),  # treatment status
                                    libid=c("SCRNA10X_SA_CHIP0142_002","TENX071","TENX076","SCRNA10X_SA_CHIP0079_001","SCRNA10X_SA_CHIP0079_002",
                                            "SCRNA10X_SA_CHIP0142_004","SCRNA10X_SA_CHIP0162_002","SCRNA10X_SA_CHIP0145_001","SCRNA10X_SA_CHIP0145_002",
                                            "SCRNA10X_SA_CHIP0175_002","SCRNA10X_SA_CHIP0151_001","SCRNA10X_SA_CHIP0162_001"),
                                    tp=c("T1","T2","T2","T3","T3","T3","T4","T4","T4","T5","T5","T5"),                   # correspongding time point 
                                    sm=c("untreated","treated","untreated","untreated","treated","previously_treated",
                                         "untreated","treated","previously_treated","untreated","treated","previously_treated"))
  rownames(treat_timepoint_map) <- treat_timepoint_map$libid
  
  # View(treat_timepoint_map)
  meta_data$treatment_summary <- treat_timepoint_map[meta_data$libid, 'sm']
  meta_data$timepoint <- treat_timepoint_map[meta_data$libid, 'tp']
  print(dim(meta_data))
  # print(colnames(meta_data))
  return(meta_data)
}


plot_clusters <- function(meta_data, save_dir){
  print(colnames(meta_data))
  meta_data <- get_metadata(meta_data)
  print(dim(meta_data))
  print(head(meta_data))
  meta_data <- 
  library(RColorBrewer)
  getPalette = colorRampPalette(brewer.pal(8, "Set1"))
  color_ls = getPalette(length(unique(meta_data$seurat_clusters)))
  cls_df <- get_pct_clusters(meta_data, color_ls, save_dir)
  meta_data$cls_freq <- cls_df[paste0(meta_data$seurat_clusters,"_",meta_data$treatment_status),'Freq']

  plot_passage_clusters(meta_data, color_ls, save_dir)

  write.table(meta_data,paste0(save_dir,'meta_data.txt'), sep="\t", quote=F, col.names=NA)

  return(meta_data)
}

plot_heatmap <- function(srt){
  
  # top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  sce <- readRDS(paste0(this.dir,"result_mt30/normv2/","SA0135_qc_filtered_overall_norm_normtotal_clusters_sce.rds"))
  genes_ens <- rownames(sce) 
  length(genes_ens)  
  sce <- sce[genes_map$GENEID,]
  sce
  rownames(sce) <- genes_map$SYMBOL
  length(genes_map$SYMBOL)
  srt2 <- as.Seurat(sce, counts = "normcounts", data = "logcounts")
  genes_ens <- rownames(x = srt2)
  genes_ens[1:3]
  genes_map <- convert_genes_set(genes_ens, fromSet="GENEID", toSet=c("SYMBOL"))
  # my_data %>% distinct(Sepal.Length, .keep_all = TRUE)
  genes_map <- dplyr::distinct(genes_map, SYMBOL, .keep_all= TRUE)   # remove duplicate  
  dim(genes_map)
  rownames(x = srt2) <-  genes_map$SYMBOL
  
  
  fn <- paste0(save_dir,"logfoldchange_cls_0_UTT_vs_UUU/","top_up_pathway_cls_0_UTT_UUU.txt")
  uppathway_df <- read.table(fn, sep="\t",header=T,row.names=NULL)
  print(head(uppathway_df))
  print(colnames(uppathway_df))
  rownames(uppathway_df) <- uppathway_df$pathway
  print(uppathway_df$leadingEdge)
  fn <- paste0(save_dir,"logfoldchange_cls_0_UTT_vs_UUU/","top_down_pathway_cls_0_UTT_UUU.txt")
  downpathway_df <- read.table(fn, sep="\t",header=T,row.names=NULL)
  
  downpathway_df$pathway
  
  marker_genes <- c()
  for(gl in uppathway_df$leadingEdge[1:10]){
    gs <- strsplit(as.character(gl), " ")
    gs <- unlist(gs, use.names=FALSE)
    marker_genes <- c(marker_genes, gs)
  }
  
  for(gl in downpathway_df$leadingEdge[1:10]){
    gs <- strsplit(as.character(gl), " ")
    gs <- unlist(gs, use.names=FALSE)
    marker_genes <- c(marker_genes, gs)
  }
  length(marker_genes)
  marker_genes[1:5]
  
  class(t)
  groups_use <- c("UTT","UUU")
  cells_use <- rownames(srt2@meta.data)[srt2@meta.data$treatment_status %in% groups_use]
  length(cells_use)
  srt3 <- srt2[,cells_use]
  srt3 <- ScaleData(object = srt3)
  df <- GetAssayData(object = srt3, slot = "scale.data")
  rownames(df)[1:3]
  p <- DoHeatmap(srt3, slot = "data", group.by = "treatment_status", features = marker_genes) + NoLegend()
  png(paste0(save_dir,"heatmap_up_downgenes.png"), height = 2*1000, width=2*1000,res = 2*72)
  print(p)
  dev.off()
  
}

plot_heatmap <- function(sce2, data_dir, save_dir, obs_cluster = 0, groups_use = c("UT","UU")){
  
  pAdjustThrs = 0.0125
  minLogFC = 0.25
  test_use = 'wilcox'
  save_dir_pw <- paste0(save_dir,"logfoldchange_","cls_",obs_cluster,"_",groups_use[1],"_vs_",groups_use[2],"/")
  markers_ls <- read.table(paste0(save_dir_pw,"seurat_wilcox_DEG_total_markers.txt"),sep="\t", header=T, row.names=1, check.names=F)
  
  
  markers_ls_upreg <- markers_ls[markers_ls$p_val_adj<pAdjustThrs & markers_ls$avg_log2FC>minLogFC,]
  markers_ls_upreg <- markers_ls_upreg[order(markers_ls_upreg$avg_log2FC,decreasing = T),] 
  markers_ls_downreg <- markers_ls[markers_ls$p_val_adj<pAdjustThrs & markers_ls$avg_log2FC<(-minLogFC),]
  markers_ls_downreg <- markers_ls_downreg[order(markers_ls_downreg$avg_log2FC,decreasing = F),] 
  
  # Can read from output file
  markers_ls_output <- markers_ls[markers_ls$p_val_adj<pAdjustThrs & abs(markers_ls$avg_log2FC)>minLogFC,]
  markers_ls_output <- markers_ls_output[order(markers_ls_output$avg_log2FC,decreasing = T),]  
  
  plttitle <- paste0("Cluster ",obs_cluster,"; ",groups_use[1]," versus ",groups_use[2])
  plot_DE_genes(markers_ls, nrow(markers_ls_output), FCcutoff=minLogFC, 
                pCutoff=pAdjustThrs, plttitle, save_dir_pw,legendVisible=F)
  # df, topGenes, nbgenessig=0, FCcutoff=0.25, pCutoff=0.01, 
  # plttitle="A versus B", save_dir="",legendVisible=F
  
  # Can read from output file
  nbtopup <- 30
  nbtopdown <- 30
  if(nrow(markers_ls_down)<nbtopdown){
    nbtopdown <- nrow(markers_ls_down)
  }
  if(nrow(markers_ls_upreg)<nbtopup){
    nbtopup <- nrow(markers_ls_upreg)
  }
  markers_ls_down <- markers_ls_down[1:nbtopdown,]
  markers_ls_upreg <- markers_ls_upreg[1:nbtopup,]
  markers_ls_upreg$genes_deg <- 'up-regulated'
  markers_ls_down$genes_deg <- 'down-regulated'
  genes_deg_df <- rbind(markers_ls_upreg, markers_ls_down)
  
  
  genes_ext <- rownames(sce2) %in% rownames(genes_deg_df)
  cells_use <- (sce2$seurat_clusters==obs_cluster) & (sce2$treatment_status %in% groups_use)
  if(length(cells_use)>0 && length(genes_ext)>0){
    sce_topgenes <- sce2[genes_ext,cells_use]
    
    norm_df <- logcounts(sce_topgenes)
    ext <- rownames(genes_deg_df)[rownames(genes_deg_df) %in% rownames(norm_df)]
    norm_df <- norm_df[ext,] # get correct order
    rownames(norm_df) <- genes_deg_df[rownames(norm_df),'symbol']
    clusters <- colData(sce_topgenes)[colnames(norm_df),'treatment_status']
    
    genes_cls <- genes_deg_df[ext,]$genes_deg
    print("Plot heatmap")
    p <- Heatmap(as.matrix(norm_df),show_column_names=FALSE, 
                 cluster_rows=F,cluster_columns=F,
                 column_split=clusters, row_split = genes_cls) 
    plttitle_hm <- paste0("cls_",obs_cluster,"_",groups_use[1],"_vs_",groups_use[2])
    png(paste0(save_dir_pw,"Heatmap_",plttitle_hm,".png"), height = 2*800, width=2*900,res = 2*72)
    print(p)
    dev.off()
    # # Pathway analysis
    # print("Pathway analysis....")
    # 
    # deg_stat <- markers_ls_output$avg_log2FC
    # names(deg_stat) <- markers_ls_output$symbol
    # print(head(deg_stat))
    # gmt_dir <- paste0(data_dir,"biodatabase/")
    # method_use <- paste0("seurat_DEG_",test_use)
    # gmt_ls = c("h.all.v7.0.symbols.gmt")
    # print("Get pathway")
    # get_pathway_results(deg_stat, method_use, obs_cluster, gmt_dir, groups_use, save_dir_pw, gmt_ls)
    return(TRUE)
  }
  
  return(FALSE)
}

cnetplot_v2 <- function(geneSets, 
                        datatag, save_dir,
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
  
  
  if (colorEdge) {
    E(g)$category <- rep(names(geneSets), sapply(geneSets, length))
    edge_layer <- ggraph::geom_edge_link(aes_(color = ~category), alpha=.8)
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
      scale_colour_gradient2(name = "fold change", low = "green",
                             mid = "blue", high = "red")
  } else {
    V(g)$color <- "#B3B3B3"
    V(g)$color[1:n] <- "#E5C494"
    p <- ggraph::ggraph(g, layout=layout, circular=circular) +
      edge_layer +
      geom_node_point(aes_(color=~I(color), size=~size))
    p
  }
  
  p <- p + scale_size(range=c(3, 15) * node_scale,
                      breaks=unique(round(seq(min(size), max(size), length.out=4)))) +
    theme_void()
  if (node_label == "category") {
    p <- p + geom_node_text(aes_(label=~name), data = p$data[1:n,])
  } else if (node_label == "gene") {
    p <- p + geom_node_text(aes_(label=~name), data = p$data[-c(1:n),], repel=TRUE)
  } else if (node_label == "all") {
    if (utils::packageVersion("ggrepel") >= "0.9.0") {
      p <- p + geom_node_text(aes_(label=~name), repel=TRUE, bg.color = "white", size = 5)
    } else {
      p <- p + geom_node_text(aes_(label=~name), repel=TRUE, size=5)
    }
    
  }
  ht <- 1000
  if(length(geneSets)<5){
    ht <- 750
  }
  png(paste0(save_dir, datatag,'_genes_network.png'), height = 2*ht, width=2*1300,res = 2*72)
  print(p)
  dev.off()
  
  # return(p)
}

# generate_genes_network(input_dir, datatag='SA609_R-UTTTT_vs_E-UUUUUUUU', 
#                        obs_pw='top_up_pathway_cls_',
#                        de_desc='SA609_UTTTT_UUUUUUUU_R_E',
#                        genes_set='hallmark',
#                        pattern_use='^hallmark_'
# )
generate_genes_network <- function(input_dir, datatag='SA609',
                                   de_desc='SA609_UTTTT_UUUUUUUU_R_E',
                                   genes_set='hallmark',
                                   pattern_use='^hallmark_'
){
  obs_pw <- c('top_up_pathway_cls_','top_down_pathway_cls_')
  pw_desc <- c('Up','Down')
  save_dir <- paste0(input_dir,genes_set,'/')
  stat_genes <- read.table(paste0(input_dir,'de_significant_genes.txt'), sep='\t', header=T, check.names=F, row.names = 1) #total_markers.txt
  # View(head(stat_genes))
  foldChange <- stat_genes$avg_log2FC
  names(foldChange) <- stat_genes$gene_symb
  
  # pathway_df <- read.csv(paste0(input_dir, genes_set,'/SA609_UTTTT_UUUUUUUU_up_pathway.csv'), check.names=F, row.names=1)
  
  for(w in rep(1:length(obs_pw),1)){
    pw_fn <- paste0(input_dir, genes_set,'/',obs_pw[w], de_desc,'.txt')
    if(file.exists(pw_fn)){
      pathway_df <- read.table(pw_fn, sep='\t',check.names=F, 
                               row.names=1, header=T)
      
      geneSets <- list()
      for(p in rownames(pathway_df)){
        gene_ls <- strsplit(as.character(pathway_df[p,'leadingEdge']),' ')
        gene_ls <- unlist(gene_ls, use.names=FALSE)
        gene_ls <- gene_ls[gene_ls %in% stat_genes$gene_symb]
        p <- gsub(pattern_use,'',tolower(p))
        geneSets[[p]] <- gene_ls
      }
      cnetplot_v2(geneSets, paste0(pw_desc[w],'_',datatag), save_dir, foldChange)
    }
  }
  
}

edgeR_de <- function(sce_de, save_dir){
  
  print("Filtering data...")
  sce_de <- sce_de[, sce_de$total_features_by_counts > 1500]
  # rs <- rowSums(as.matrix(counts(sce_de)))
  # qplot(rs, log='x') + geom_vline(xintercept = 100)
  
  sce_de <- sce_de[rowSums(as.matrix(counts(sce_de))) > 100, ]
  print(dim(sce_de))
  mycounts <- as.matrix(counts(sce_de))    # zeros are good
  print("Create DGE edgeR object...")
  # dge <- DGEList(counts=mycounts, group=sce_de$clone)
  dge <- DGEList(counts=mycounts, group=sce_de$treatmentSt)
  
  print("DE Analysis...")
  # which vs which ???
  # design <- model.matrix(~ clone, data = colData(sce_de))
  design <- model.matrix(~ treatmentSt, data = colData(sce_de))
  # This describes the edgeR user manual
  # http://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
  
  dge <- edgeR::estimateDisp(dge, design = design)
  fit <- edgeR::glmQLFit(dge, design = design)
  qlf <- edgeR::glmQLFTest(fit)
  tt <- edgeR::topTags(qlf, n = Inf)
  
  print("Generate output...")
  tt <- as.data.frame(tt) %>% 
    rownames_to_column("gene_id")
  # tt$gene_symbol <- rowData(sce_de[tt$gene_id,])$ID
  tt$gene_symbol <- rowData(sce_de[tt$gene_id,])$Symbol
  
  # Saving the logFC file
  # tt <- tt[tt$FDR<0.01 & tt$PValue<0.05 & abs(tt$logFC)>0.25,] # 
  write.csv(tt, file=paste0(save_dir, 'edgeR_significant_genes.csv'), row.names=FALSE, quote=FALSE)
  print("With FDR<0.01 and PValue<0.05 and abs(tt$logFC)>0.25, number of significant genes is: ")
  print(dim(tt))
  # print("With threshold logFC>0.25, number of significant genes is: ")
  # print(nrow(tt[abs(tt$logFC)>0.25,]))
  return(tt)
}

edgeR_DE_analysis_by_clones_treatment_status <- function( input_dir, input_file, output_file,
                               datatag='SA', fdr_thrs=0.01){
  input_dir=paste0(input_dir,'/')
  save_dir <- paste0(dirname(output_file),'/')
  if (!file.exists(save_dir)){
    dir.create(save_dir)
  }
  
  sce <- readRDS(input_file)
  print("Initialized sce")
  print(dim(sce))
  print(summary(as.factor(sce$clone)))
  
  # mito_genes <- str_detect(rowData(sce)$Symbol, "^MT\\-")
  # sum(mito_genes==TRUE)
  # 
  # ribo_genes <- str_detect(rowData(sce)$Symbol, "^RP(L|S)")  # or ^RP[L|S]?
  # sum(ribo_genes==TRUE)
  # sce <- sce[(!mito_genes) & (!ribo_genes), ]
  print("Observed sce: ")
  print(paste0('Removing mito and ribo genes: ',dim(sce)[1],'_',dim(sce)[2]))
  
  groups_use_ls <- list(T1=c('R','C'))
  treatment_st <- list(R=c("UTTT"),C=c("UUUU"))
  # 
  # meta_data <- meta_data[meta_data$clone %in% c('R','C'),]
  # table(meta_data$clone,meta_data$treatmentSt)
  meta_data <- as.data.frame(colData(sce))
  # summary(as.factor(meta_data$clone))
  # rownames(meta_data)[1:3]
  # View(head(meta_data))
  dim(meta_data)
  for(groups_use in groups_use_ls){
    cells_use_g1 <- rownames(meta_data)[meta_data$clone==groups_use[1] & meta_data$treatmentSt %in% treatment_st[[groups_use[1]]]]
    cells_use_g2 <- rownames(meta_data)[meta_data$clone==groups_use[2] & meta_data$treatmentSt %in% treatment_st[[groups_use[2]]]]
    if(length(cells_use_g1) < 100 || length(cells_use_g2) < 100){
      print(groups_use)
      print("There are no cells or small nb cells 
            only which satisfy the input condition ")
      print(paste0("\n Observed clones: ",groups_use[1],' vs ',groups_use[2], 
                   "  nb cells in ",groups_use[1], ":",length(cells_use_g1),
                   "  nb cells in ",groups_use[2], ":",length(cells_use_g2)))
      
    } else{
      sub_dir <- paste0(datatag,"_",groups_use[1],"-",
                        paste(treatment_st[groups_use[1]], collapse = "-"),'_',
                        groups_use[2],"-",paste(treatment_st[groups_use[1]], collapse = "-"))
      save_dir_pw <- paste0(save_dir,sub_dir,"/")
      if (!file.exists(save_dir_pw)){
        dir.create(save_dir_pw)
      }
      cells_use <- c(cells_use_g1, cells_use_g2)
      # cells_use <- rownames(meta_data)
      cat("DE analysis ", file = paste0(save_dir_pw,"de_analysis_log.txt"))
      print(length(cells_use))
      cat(paste0("\n Observed clones: ",groups_use[1],' vs ',groups_use[2]," Nb obs cells: ",length(cells_use), 
                 "  nb cells in ",groups_use[1], ":",length(cells_use_g1),
                 "  nb cells in ",groups_use[2], ":",length(cells_use_g2)), file = paste0(save_dir_pw,"de_analysis_log.txt"), append = TRUE)
      
      sce_tmp <- sce[,cells_use]
      print(dim(sce_tmp))
      
      # Similar to findMarkers func in Seurat, only test genes that are detected in a minimum fraction 
      # of min.pct cells in either of the two populations. Meant to speed up the function
      # by not testing genes that are very infrequently expressed. Default is 0.1
      min_pct = 0.1 
      zero_g1 <- DelayedArray::rowMeans(assay(sce_tmp[,sce_tmp$clone==groups_use[1]], "counts") == 0)
      obs_genes1 <- names(zero_g1[zero_g1 <= (1 - min_pct)])
      print(length(obs_genes1))
      zero_g2 <- DelayedArray::rowMeans(assay(sce_tmp[,sce_tmp$clone==groups_use[2]], "counts") == 0)
      obs_genes2 <- names(zero_g2[zero_g2 <= (1 - min_pct)])
      print(length(obs_genes2))
      genes_use <- intersect(obs_genes1, obs_genes2)
      sce_tmp <- sce_tmp[genes_use,]
      print("Filtering by pct minimum fraction genes:")
      print(dim(sce_tmp))
      sce_tmp <- sce_tmp[, sce_tmp$total_features_by_counts > 1500]
      print("Filtering by total_features_by_counts greater than 1500:")
      print(dim(sce_tmp))
      sce_tmp$clone <- ifelse(sce_tmp$clone==groups_use[1],paste0("2_",sce_tmp$clone),paste0("1_",sce_tmp$clone))
      print(summary(as.factor(sce_tmp$clone)))
      
      edgeR_de(sce_tmp, datatag, save_dir_pw)
      
    }  
  }
  
}


plot_DE_genes_ggplot <- function(df, topGenes, capstr='', FDRcutoff=0.01, logFCcutoff=0.25, pValuecutoff=0.05,
                                 plttitle="A versus B", save_dir="",legendVisible=F,
                                 iscaption=TRUE, legend_verbose='none', save_plot=TRUE, 
                                 xl=NULL, yl=NULL){  
  # library(ggplot2)
  # library(ggrepel)
  #legend_verbose is none or 'right', 'left','bottom'
  # df <- de_genes
  # library(EnhancedVolcano)
  
  # colnames(df)[which(names(df) == "avg_log2FC")] <- "logFC"
  # colnames(df)[which(names(df) == "p_val_adj")] <- "padj"
  df <- df %>%
    dplyr::rename(logFC=avg_log2FC, padj=p_val_adj) #
  # capstr <- paste0("FDR cutoff, ",FDRcutoff,"; logFC cutoff, ",logFCcutoff, "; nb genes signif, ",nbgenessig)
  # summary(as.factor(df$gene_type))
  
  # keyvals_colour <- factor(keyvals_colour, levels = unique(keyvals_colour))
  # unique(keyvals_colour)
  # names(keyvals_colour[1:3])
  df <- df[abs(df$logFC)>logFCcutoff & df$padj<pValuecutoff  & df$p_val<pValuecutoff,]
  # df$logFC <- sapply(df$logFC, function(x) replace(x, x > 3, 3))
  # df$logFC <- sapply(df$logFC, function(x) replace(x, x < (-3), -3))
  print(dim(df))
  df$gene_type <- ifelse(df$logFC>0,'Up-regulated','Down-regulated')
  st <- summary(as.factor(df$gene_type))
  # keyvals_shape <- ifelse(df$is_fitness_gene==T, 3, 16)
  # names(keyvals_shape) <- ifelse(df$is_fitness_gene==T,'Fitness genes','Others')
  
  # df$gene_type <- factor(df$gene_type, levels = unique(df$gene_type))
  if(capstr=='' & iscaption){
    # capstr <- paste0(capstr,'With abs(logFC)>0.25, FDR<0.01, pValue<0.05 \n')
    capstr <- paste0(capstr,names(st[1]),':',as.numeric(st[1]), ',  ')
    capstr <- paste0(capstr,names(st[2]),':',as.numeric(st[2]), ' \n')
    # for(i in rep(1:length(st),1)){
    #   # print(st[i])
    #   capstr <- paste0(capstr,names(st[i]),':',as.numeric(st[i]), ' ')
    # }
    
  }
  df$padj <- -log10(df$padj)
  df$padj <- sapply(df$padj, function(x) replace(x, is.infinite(x), 300))
  if(is.null(yl)){
    yl <- c(0, max(df$padj)+10)
  }
  # df$mlog10FDR <- -log10(df$FDR)
  # df$mlog10FDR <- sapply(df$mlog10FDR, function(x) replace(x, is.infinite(x), 300))
  # df$gt_alpha <- ifelse(df$gene_type=='In trans Up', 0.8,
  #                       ifelse(df$gene_type=='In cis Gain-Up',0.8,0.65))  
  # df$gt_alpha <- 0.8
  # df$gt_alpha <- ifelse(df$gene_type=='In trans Up' | df$gene_type=='In cis Gain-Up', 0.8, 0.9)  
  
  my_font <- "Helvetica"
  
  
  thesis_theme <- theme(  text=element_text(size = 8,family=my_font),
                          plot.title = element_text(color="black", size=11, hjust = 0, face = "bold", family=my_font),
                          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                          axis.text.y = element_text(size=7, hjust = 0.5, family=my_font),
                          axis.text.x = element_text(size=7, hjust = 0.5, family=my_font),
                          axis.title = element_text(size=8, hjust = 0.5, family=my_font),
                          plot.caption = element_text(size=8, hjust = 1, family=my_font),
                          legend.title = element_text(size=7, hjust = 0.5, family=my_font),
                          legend.text = element_text(size=7, hjust = 0, family=my_font),
                          strip.text.x = element_text(color="black",size=9, family=my_font),
                          strip.text.y = element_text(color="black",size=9, family=my_font),
                          legend.spacing.x = unit(0.1, 'mm'),
                          legend.spacing.y = unit(0.1, 'mm'),
                          legend.key.height=unit(1,"line"),
                          # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          legend.position = legend_verbose)
  
  p <- ggplot(df, aes(x = logFC, y = padj)) +
    geom_point(color='red', size=2, alpha=0.6) +#aes(color = gene_type)
    geom_vline(xintercept=c(-0.25,0.25), linetype="dotted")+
    # scale_color_manual(values = keyvals_colour, name = "Gene Type") +
    thesis_theme + 
    geom_text_repel(family = my_font,
                    data = df[df$gene_symb %in% topGenes, ], 
                    aes(label = gene_symb), size = 3.5, 
                    box.padding = unit(0.35, "lines"), 
                    point.padding = unit(0.3, "lines"),
                    max.overlaps = Inf,
                    min.segment.length = 0,  # draw segment lines, not matter how short they are)
                    color='black', segment.alpha = 0.2)
  # p <- p + labs(x= bquote(~Log[2] ~ ' Fold Change'), y=bquote(~-Log[10]~italic(FDR)),title = plttitle,
  #               caption = capstr)
  p <- p + labs(x= bquote(~Log[2] ~ ' Fold Change'), y=bquote(~-Log[10]~italic(P_adjusted_value)),title = plttitle,
                caption = capstr)
  if(!is.null(xl)){
    p <- p + xlim(xl[1],xl[2])
  }
  if(!is.null(yl)){
    p <- p + ylim(yl[1],yl[2])
  }
  # if(legend_verbose!='none'){
  #   return(p)
  # }
  
  # xdens <- cowplot::axis_canvas(p, axis = "x") +
  #   stat_density(data = df, aes(x = logFC, group = gene_type_classify,color =gene_type_classify),
  #                alpha = 1, size = 1, position="identity",geom="line") +
  #   scale_color_manual(values = keyvals_colour_classify, name = "Gene Classify") 
  
  
  # xdens <- cowplot::axis_canvas(p, axis = "x") +
  #          # geom_jitter(data = df, aes(x = logFC, fill = gene_type_classify), 
  #          #             shape=16, position=position_jitter(0.2)) + 
  #          geom_boxplot(data = df, aes(x = logFC, color = gene_type_classify),
  #                alpha = 0.6, size = 0.3) +
  #         # scale_color_grey()
  #         scale_color_manual(values = keyvals_colour_classify, name = "Gene Classify")
  #         
  # 
  # p_main <- cowplot::insert_xaxis_grob(p, xdens, grid::unit(.2, "null"), position = "top")
  # 
  
  # Just to get legend
  # t <- ggplot(df, aes(x = logFC, color = gene_type_classify)) +
  #      geom_boxplot(alpha = 0.6, size = 0.3) +
  #      scale_color_manual(values = keyvals_colour_classify, name = "Gene Classify") + 
  #      theme(text=element_text(size = 8,family=my_font),
  #           legend.title = element_text(size=7, hjust = 0.5, family=my_font),
  #           legend.text = element_text(size=7, hjust = 0, family=my_font),
  #           strip.text.x = element_text(color="black",size=9, family=my_font),
  #           strip.text.y = element_text(color="black",size=9, family=my_font),
  #           legend.spacing.x = unit(0.1, 'mm'),
  #           legend.spacing.y = unit(0.1, 'mm'),
  #           legend.key.height=unit(1,"line"))
  # lg <- cowplot::get_legend(t)
  # plg <- cowplot::ggdraw() + cowplot::draw_plot(lg)
  
  # blank_theme <- theme_minimal()+
  #   theme(
  #     axis.title.x = element_blank(),
  #     axis.title.y = element_blank(),
  #     panel.border = element_blank(),
  #     panel.grid=element_blank(),
  #     axis.ticks = element_blank(),
  #     plot.title=element_text(size=14, face="bold")
  #   )
  
  # p1 <- ggplot(df, aes(x="", y=gene_type, fill=gene_type))+
  #   geom_bar(width = 1, stat = "identity", alpha=0.8) +
  #   coord_polar("y", start=0) + 
  #   scale_fill_manual(values=keyvals_colour) + 
  #   blank_theme + 
  #   theme(axis.text.x=element_blank(), legend.position = "none")
  # 
  # caption = capstr
  # plg <- ggplot() +
  #   annotate("text", x = 8, y = 20, color="black", size=4.5, family=my_font, label = capstr) +
  #   theme_void()
  # p2
  # p1 <- cowplot::plot_grid(p2, p1, ncol=2, rel_widths = c(2,1))
  # p_total <- cowplot::plot_grid(p, p1, rel_heights = c(1,0.2), ncol=1)
  # p_total
  # tag <- ifelse(legend_verbose=='none','','_with_legend')
  # png(paste0(save_dir,"DE_",gsub(' ','_',plttitle),tag,".png"), height = 2*650, width=2*500, res = 2*72)
  # print(p_total)
  # dev.off()
  
  
  if(save_plot){
    plttitle <- gsub(':','_',plttitle)
    plttitle <- gsub(' ','_',plttitle)
    tag <- ''
    saveRDS(p, file=paste0(save_dir,"DE_",plttitle,".rds"))
    ggsave(paste0(save_dir,"DE_",gsub(' ','_',plttitle),tag,".pdf"),
           plot = p,
           height = 5.5,
           width = 6.5,
           useDingbats=F)
    ggsave(paste0(save_dir,"DE_",gsub(' ','_',plttitle),tag,".png"),
           plot = p,
           height = 5.5,
           width = 6.5,
           # useDingbats=F,
           type = "cairo-png",
           dpi=200)
  }
  # return(list(p_vol=p, p_prop=p1))
  return(p)
}

