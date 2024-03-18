
suppressPackageStartupMessages({
  require(RColorBrewer)
  # require(glue)
  # require(ggrepel)
  # require(tidyr)
  # require(fgsea)
  # require(scales)
  # require(ComplexHeatmap)
  # require(DOSE)
  require(enrichplot)
  require(stringr)
  require(ggraph)
  require(ggplot2)
  require(igraph)
  require(tidyverse)
})

# input_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/dlp_cnv/mapping_CNA_geneID/'
# df <- data.table::fread(paste0(input_dir, 'mapped_wholedata_SA1035_v2.csv.gz'))
# dim(df)
# head(df)
source('/home/htran/Projects/hakwoo_project/metastasis_material/scripts/bulk_rna/bulk_utils.R')

plot_upSetR <- function(){
  # Specific library
  # BiocManager::install("UpSetR")
  library(UpSetR)
  
  # Dataset
  input <- c(
    M.acuminata = 759,
    P.dactylifera = 769,
    "P.dactylifera&M.acuminata" = 467
  )  
  input <- c(
    B_met_vs_B_pri_upReg = 481,
    B_met_vs_B_pri_downReg = 259,
    C_met_vs_B_pri_upReg = 60,
    C_met_vs_B_pri_downReg = 41,
    "B_met_vs_B_pri_upReg&C_met_vs_B_pri_upReg" = 50,
    "B_met_vs_B_pri_downReg&C_met_vs_B_pri_downReg" = 24,
    "B_met_vs_B_pri_upReg&C_met_vs_B_pri_downReg" = 1,
    "B_met_vs_B_pri_downReg&C_met_vs_B_pri_upReg" = 0
  )  
  
  # Plot
  p_upset <- upset(fromExpression(input), 
        # nintersects = 40, 
        # nsets = 6, 
        order.by = "freq", 
        decreasing = T, 
        # mb.ratio = c(0.6, 0.4),
        number.angles = 0, 
        text.scale = 1.1, 
        point.size = 2.8, 
        line.size = 1
  )
  library(grid)
  ggsave(paste0(save_dir,"BB_CB_upset.png"),  
         plot =  print(p_upset),  
         height = 4,  
         width = 5.5,  
         # useDingbats=F,  
         type = "cairo-png",  
         dpi=150  
  )
  png(paste0(save_dir,"BB_CB_upset.png"), units = "in", height = 3.5, width = 5.5, res = 150)
  print(p_upset)                               ## use print here
  dev.off()
}
debug_DE_analysis <- function(){
  input_dir <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/bulkRNAseq/'
  bb <- data.table::fread(paste0(input_dir, 'SA919_DE_analysis/Passage7_CloneBMeta_vs_CloneBPrim_7Passage_DESeq_Sig_GS.csv'))
  bc <- data.table::fread(paste0(input_dir, 'SA919_DE_analysis/Passage7_CloneCMeta_vs_CloneBPrim_DESeq_Sig_GS.csv'))
  save_dir <- paste0(input_dir, 'SA919_DE_analysis/pathways_gprofiler/')
  # dir.create(save_dir)
  dim(bb)
  dim(bc)
  
  
  bb_up <- bb %>%
    dplyr::filter(log2FoldChange>=1 & padj<0.05)
  bb_down <- bb %>%
    dplyr::filter(log2FoldChange<=-1 & padj<0.05)
  bc_up <- bc %>%
    dplyr::filter(log2FoldChange>=1 & padj<0.05)
  bc_down <- bc %>%
    dplyr::filter(log2FoldChange<=-1 & padj<0.05)
  dim(bb_up)
  dim(bb_down)
  dim(bc_up)
  dim(bc_down)
  length(intersect(bb_up$Gene.name, bc_up$Gene.name))
  length(intersect(bb_down$Gene.name, bc_down$Gene.name))
  length(intersect(bb_down$Gene.name, bc_up$Gene.name))
  length(intersect(bb_up$Gene.name, bc_down$Gene.name))
  bc <- bc %>%
    # dplyr::filter(!Gene.name %in% obs_genes_symb) %>%
    dplyr::rename(logFC=log2FoldChange, gene_symbol=Gene.name)
  
  
  ## Redo DE analysis C met vs. B pri
  raw_counts <- data.table::fread(paste0(input_dir, 'SA919/SA919_total_raw_counts.csv.gz')) %>% as.data.frame()
  dim(raw_counts)
  sids <- colnames(raw_counts)[colnames(raw_counts) != 'ens_gene_id']
  t <- raw_counts[, sids]
  colSums(t)
}
get_connected_genes_network <- function(){
  ## First, for clone B, C. 
  input_dir <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/bulkRNAseq/'
  bb <- data.table::fread(paste0(input_dir, 'SA919_DE_analysis/Passage7_CloneBMeta_vs_CloneBPrim_7Passage_DESeq_Sig_GS.csv'))
  bc <- data.table::fread(paste0(input_dir, 'SA919_DE_analysis/Passage7_CloneCMeta_vs_CloneBPrim_DESeq_Sig_GS.csv'))
  save_dir <- paste0(input_dir, 'SA919_DE_analysis/pathways_gprofiler/')
  # dir.create(save_dir)
  bc <- bc %>%
    # dplyr::filter(!Gene.name %in% obs_genes_symb) %>%
    dplyr::rename(logFC=log2FoldChange, gene_symbol=Gene.name)
  
  obs_genes_symb <- intersect(bb$Gene.name, bc$Gene.name)
  # obs_genes_symb <- intersect(bb$g_ID, bc$g_ID)
  length(obs_genes_symb)
  
  cna <- data.table::fread(paste0(input_dir, 'SA919_DE_analysis/mapping_gene_cnv_SA919.csv.gz'))
  dim(cna)
  head(cna)
  # cna <- cna %>%
  #   dplyr::filter(chr %in% c(5, 7, 10))
  dim(cna)
  
  var_genes <- cna$gene_symbol[matrixStats::rowVars(as.matrix(cna[, 5:6]))>0]
  var_genes[1:30]
  length((var_genes))
  
  inf_graph <- read.delim(paste0(input_dir, 'SA919_DE_analysis/influence_graph.txt'), header = TRUE, sep = "\t", dec = ".")
  dim(inf_graph)
  head(inf_graph)
  inf_graph <- inf_graph %>%
    dplyr::mutate(desc=paste0(gene_a, '_',gene_b,'_',weight))
  direct_df <- inf_graph %>%
    dplyr::filter(gene_a %in% var_genes & gene_b %in% bc$gene_symbol)
  dim(direct_df)
  direct_df1 <- direct_df %>%
    dplyr::group_by(gene_b) %>%
    dplyr::summarise(desc2=max(weight)) %>%
    dplyr::mutate(desc2=paste0(gene_b, '_',desc2))
  direct_df <- direct_df %>%
    dplyr::mutate(desc2=paste0(gene_b, '_',weight)) %>%
    dplyr::filter(desc2 %in% direct_df1$desc2)
  
  dim(direct_df)
  head(direct_df)
  direct_df2 <- direct_df %>%
    dplyr::select(gene_b, desc) %>%
    dplyr::mutate(desc=paste0(desc, ': direct_1st_layer')) %>%
    dplyr::rename(C_vs_B_DE_gene=gene_b, connect_to_cis_gene=desc)
  head(direct_df2)
  remain_genes <- bc$gene_symbol[!bc$gene_symbol %in% direct_df$gene_b]
  length(remain_genes)
  head(direct_df)
  length(unique(direct_df$gene_b))
  indirect_df <- inf_graph %>%
    dplyr::filter(gene_a %in% var_genes & !desc %in% direct_df$desc) %>%
    dplyr::select(desc, gene_b) %>%
    dplyr::rename(desc1=desc, gene_1st=gene_b)
  indirect_df2 <- inf_graph %>%
    dplyr::filter(gene_a %in% unique(indirect_df$gene_1st) 
                  & gene_b %in% remain_genes)
    # dplyr::full_join(indirect_df, by=c('gene_a'='gene_1st'), relationship = "many-to-many") %>%
  indirect_df3 <- indirect_df2 %>%
    dplyr::group_by(gene_b) %>%
    dplyr::summarise(desc2=max(weight)) %>%
    dplyr::mutate(desc2=paste0(gene_b, '_',desc2))
  indirect_df2 <- indirect_df2 %>%
    dplyr::mutate(desc2=paste0(gene_b, '_',weight)) %>%
    dplyr::filter(desc2 %in% indirect_df3$desc2)
  
  dim(indirect_df2)
  length(unique(indirect_df2$C_vs_B_DE_gene))
  head(indirect_df2)
  indirect_df2 <- indirect_df2 %>%
    dplyr::select(gene_b, desc) %>%
    dplyr::mutate(desc=paste0(desc, ': indirect_2nd_layer')) %>%
    dplyr::rename(C_vs_B_DE_gene=gene_b, connect_to_cis_gene=desc)
  
  indirect_df2 <- indirect_df2[!duplicated(indirect_df2$C_vs_B_DE_gene),]
  direct_df2 <- direct_df2[!duplicated(direct_df2$C_vs_B_DE_gene),]
  total_df <- dplyr::bind_rows(direct_df2, indirect_df2)
  dim(total_df)
  head(total_df)
  dim(direct_df2)
  dim(indirect_df2)
  data.table::fwrite(total_df, paste0(save_dir, 'cloneC_vs_cloneB_DE_genes_cis_connected.csv'))
  
  pw_df <- data.table::fread(paste0(save_dir, 'pathways_SA919_common_genes_BC_2024Feb17_003335.csv.gz'))
  dim(pw_df)
  colnames(pw_df)
  pw_df <- as.data.frame(pw_df)
  pw_df$reference_set <- ifelse(pw_df$reference_set=="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_EMT",pw_df$reference_set)
  rownames(pw_df) <- pw_df$reference_set
  pattern_use <- '^hallmark_'
  # geneSets <- list()
  ref_genes <- tibble::tibble()
  for(p in rownames(pw_df)){
    
    gene_ls <- strsplit(as.character(pw_df[p,'signif_genes']),',')
    gene_ls <- unlist(gene_ls, use.names=FALSE)
    p <- gsub(pattern_use,'',tolower(p))
    # p <- gsub(pattern_use,'',toupper(p))
    tmp <- tibble::tibble(gene=gene_ls, significant_pathway=p)
    ref_genes <- dplyr::bind_rows(ref_genes, tmp)
    # geneSets[[p]] <- gene_ls
    # print(length(gene_ls))
  }
  dim(ref_genes)
  head(ref_genes)
  ref_genes <- ref_genes %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(significant_pathway=paste(significant_pathway, collapse = '_'))
  dim(ref_genes)
  head(ref_genes)
  total_df1 <- total_df
  head(total_df1)
  dim(total_df1)
  total_df1 <- total_df1 %>%
    left_join(ref_genes, by=c('C_vs_B_DE_gene'='gene'))
  View(total_df1)
  sum(!is.na(total_df1$significant_pathway))
  data.table::fwrite(total_df1, paste0(save_dir, 'cloneC_vs_cloneB_DE_genes_cis_connected.csv'))
}

get_pathway_viz <- function(){
  ## First, need gprofiler enrichment pathways
  pw_df <- get_gprofiler_pathways_obsgenes(obs_genes_symb, save_dir, datatag, 
                                           custom_id=NULL, pathway_fn=NULL, save_data=T)
  dim(pw_df)
  pw_df <- data.table::fread(paste0(save_dir, 'pathways_SA919_common_genes_BC_2024Feb17_003335.csv.gz'))
  dim(pw_df)
  pw_df <- as.data.frame(pw_df)
  pw_df$signif_genes[2]
  length(ref_set$HALLMARK_ANGIOGENESIS)
  
  # logFC_fn file should contain logFC, gene_symbol columns
  # pw_fn file should contain Term, ledge_genes, es, fdr. ledge_genes contains the list of genes, separated by ;
  colnames(bb)
  bb_common <- bb %>%
    dplyr::filter(Gene.name %in% obs_genes_symb) %>%
    dplyr::rename(logFC=log2FoldChange, gene_symbol=Gene.name)
  bc_common <- bc %>%
    dplyr::filter(Gene.name %in% obs_genes_symb) %>%
    dplyr::rename(logFC=log2FoldChange, gene_symbol=Gene.name)
  
  bb_common1 <- bb_common %>%
    dplyr::select(logFC, gene_symbol) %>%
    dplyr::rename(logFC_BB=logFC)
  bc_common1 <- bc_common %>%
    dplyr::select(logFC, gene_symbol) %>%
    dplyr::rename(logFC_BC=logFC)
  df <- bb_common1 %>%
    left_join(bc_common1, by='gene_symbol')
  sum(df$logFC_BB>0 & df$logFC_BC>0)
  sum(df$logFC_BB<0 & df$logFC_BC<0)
  t <- df %>%
    filter(logFC_BB>0 & logFC_BC<0)
  pw_df1 <- pw_df
  colnames(pw_df1)
  pw_df1 <- pw_df1 %>%
    dplyr::rename(Term=reference_set, ledge_genes=signif_genes) %>%
    dplyr::mutate(es=1, fdr=0.01)
  logFC_fn <- paste0(save_dir, 'logFC_common_BB.csv.gz')
  pw_fn <- paste0(save_dir, 'pw_BB_BC.csv.gz')
  data.table::fwrite(bb_common,logFC_fn)
  data.table::fwrite(pw_df1, pw_fn)
  
  
  pw_fn <- paste0(save_dir, 'pw_BB_BC.csv.gz')
  logFC_fn <- paste0(save_dir, 'logFC_common_BC.csv.gz')
  data.table::fwrite(bc_common,logFC_fn)
  datatag <- 'SA919_common_genes_BB_BC'
  tag <- 'logFC_BC'
  
  viz_genes_network(save_dir,
                    pw_fn,
                    logFC_fn,
                    datatag,
                    tag
  )
  
  summary(bb_common$logFC)
  
  for(i in seq(1:dim(pw_df1)[1])){
    gene_ls <- strsplit(as.character(pw_df1[i,'ledge_genes']),',')
    tmp <- bc_common %>%
      dplyr::filter(gene_symbol %in% unlist(gene_ls))
    print(dim(tmp))
    print(paste0(pw_df1$Term[i],': up-reg:',sum(tmp$logFC>0), ' ,down-reg:',sum(tmp$logFC<0)))
  }
  bb_common$gene_symbol[1:10]
  
}

get_norm_heatmap <- function(){
  norm_df <- data.table::fread(paste0(input_dir, 'tmp_Hoa/SA919_total_normalized_sizefactor_gene_symbols.csv.gz'))
  dim(norm_df)
  meta_samples <- data.table::fread(paste0(input_dir, 'SA919/library_groupings_bulk_SA919_cloneIds.csv'))
  meta_samples <- meta_samples %>%
    dplyr::select(mainsite, clone_id, bulk_sid) %>%
    dplyr::filter(bulk_sid %in% colnames(norm_df) & clone_id %in% c('B','C'))%>%
    dplyr::mutate(bulk_sid1 =stringr::str_sub(bulk_sid, nchar(bulk_sid)-8, nchar(bulk_sid)),
                  label=paste0(mainsite, '_',clone_id, '_',bulk_sid1))
  # colnames(norm_df) %in% meta_samples$bulk_sid
  meta_samples$label <- gsub('XB0','_',meta_samples$label)
  # meta_samples$label
  meta_samples <- as.data.frame(meta_samples)
  rownames(meta_samples) <- meta_samples$label
  order_sids <- c("Primary_B_X7_5402","Primary_B_X7_5378", 
                  "Metastasis_B_X7_5691", "Metastasis_B_X7_5692", "Metastasis_C_X7_5604")
  
  meta_samples <- meta_samples[order_sids,]
  meta_samples
  i <- 3
  pw_df1$Term[i]
  gene_ls <- unlist(strsplit(as.character(pw_df1[i,'ledge_genes']),','))
  
  sids <- meta_samples$label
  names(sids) <- meta_samples$bulk_sid
  
  norm_df_emt <- norm_df %>%
    dplyr::filter(symbol %in% gene_ls) %>%
    tibble::column_to_rownames('symbol') %>%
    dplyr::select(all_of(meta_samples$bulk_sid))
  dim(norm_df_emt)
  colnames(norm_df_emt) <- sids[colnames(norm_df_emt)]
  norm_df_emt <- log2(norm_df_emt + 1)
  
  col_mt <- c("Primary" = "darkgreen", "Metastasis" = "purple")
  mts <- meta_samples[colnames(norm_df_emt),'mainsite']
  names(mts) <- colnames(norm_df_emt)
  anno_col <- ComplexHeatmap::columnAnnotation(
    MainSite = mts, gap = unit(1, "points"),
    col = list(MainSite=col_mt), 
    simple_anno_size = unit(0.45, "cm"), show_legend = T#,
    # annotation_name_rot = 90,
    # annotation_name_gp= gpar(fontsize = 8)
  )
  # library(circlize)
  # col_fun = colorRamp2(c(0, 10, 20), c("blue", "white", "red"))
  p <- ComplexHeatmap::Heatmap(as.matrix(norm_df_emt), na_col = "white",
                               show_column_names=T,
                               show_row_names = T,
                               # col = col_fun,
                               # cluster_rows=cluster_within_group(t(exp_mtx), genes_clusters),
                               cluster_rows=T, 
                               cluster_columns=F,
                               name = "Normalized Exp", 
                               # row_order = sort(rownames(test)),
                               # row_split= genes_clusters,
                               # row_title_rot = 0,
                               row_gap = unit(1.2, "mm"),
                               # column_split = obs_cells_df, 
                               # column_title = "ddd",
                               column_gap = unit(1.2, "mm"),
                               column_names_gp = grid::gpar(fontsize = 13),
                               column_title_gp = grid::gpar(fontsize = 11),
                               row_names_gp = grid::gpar(fontsize = 13, hjust=0.5),
                               row_title_gp = grid::gpar(fontsize = 11, hjust=0.5),
                               row_title_rot = 0,
                               row_names_rot = 0,
                               top_annotation=anno_col,
                               # left_annotation = left_anno,
                               # right_annotation = anno_row,
                               # bottom_annotation = yaxis_label,
                               # cell_fun = cell_func,
                               # row_dend_reorder=F,
                               show_column_dend = T,
                               show_row_dend = T,
                               show_heatmap_legend = T,
                               # row_title = plttitle, 
                               heatmap_legend_param = list(direction = "horizontal")
                               # column_title = "Pseudotime"
                               # column_title_side = "bottom"
  )
  p
  length(gene_ls)
  tag <- 'angiogenesis'
  png(paste0(save_dir,tag,"_genes_hm.png"), height = 2*80*length(gene_ls), 
      width=2*600, res = 2*72)
  print(p)
  dev.off()
  
}
