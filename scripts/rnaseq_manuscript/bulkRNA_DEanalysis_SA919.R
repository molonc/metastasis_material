
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
source('/home/htran/Projects/hakwoo_project/metastasis_material/scripts/bulk_rna/bulk_utils.R')
## First, for clone B, C. 
input_dir <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/bulkRNAseq/'
bb <- data.table::fread(paste0(input_dir, 'SA919_DE_analysis/Passage7_CloneBMeta_vs_CloneBPrim_7Passage_DESeq_Sig_GS.csv'))
bc <- data.table::fread(paste0(input_dir, 'SA919_DE_analysis/Passage7_CloneCMeta_vs_CloneBPrim_DESeq_Sig_GS.csv'))
save_dir <- paste0(input_dir, 'SA919_DE_analysis/pathways_gprofiler/')
# dir.create(save_dir)

obs_genes_symb <- intersect(bb$Gene.name, bc$Gene.name)


## First, need gprofiler enrichment pathways
pw_df <- get_gprofiler_pathways_obsgenes(obs_genes_symb, save_dir, datatag, 
                                custom_id=NULL, pathway_fn=NULL, save_data=T)
dim(pw_df)
pw_df <- data.table::fread(paste0(save_dir, 'pathways_SA919_common_genes_BC_2024Feb17_003335.csv.gz'))
dim(pw_df)
pw_df <- as.data.frame(pw_df)
pw_df$signif_genes[3]
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
