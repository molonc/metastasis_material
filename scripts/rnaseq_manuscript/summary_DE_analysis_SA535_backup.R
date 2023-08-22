

# DE analysis summary for SA535
gc()
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  options(dplyr.summarise.inform = FALSE)
  options(tidyverse.quiet = TRUE)
  library(data.table)
})  

input_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/DE_analysis_introns/'
fds <- list.files(input_dir)
fds <- fds[fds != 'temp']
tag <- 'SA535'
metadata <- data.frame(de_comparison=fds, datatag=tag)
metadata$mouse_id <- c('M1_M3','M1_M3','M2','M2','M2','M2')
rownames(metadata) <- metadata$de_comparison
metadata$nb_DE_genes <- 0
de_genes <- tibble::tibble()
for(de in metadata$de_comparison){
  fn <- paste0(input_dir, de, '/de_significant_genes.csv.gz')
  print(de)
  tmp <- data.table::fread(fn)
  print(dim(tmp))
  tmp$de_comparison <- de
  metadata[de,'nb_DE_genes'] <- dim(tmp)[1]
  de_genes <- dplyr::bind_rows(de_genes, tmp)
}
data.table::fwrite(metadata, paste0(input_dir, 'summary_comparisons.csv'))

dim(de_genes)
head(de_genes)
de_genes <- de_genes %>%
  dplyr::left_join(metadata, by='de_comparison')

colnames(de_genes)


de_genes <- de_genes %>%
  dplyr::mutate(direction=
                  case_when(
                    avg_log2FC>=0 ~ 'pos',
                    TRUE ~ 'neg'
                  )) %>%
  dplyr::mutate(desc=paste0(gene_symb,'_',mouse_id,'_',direction))
t <- de_genes[!duplicated(de_genes$desc),]
dim(t)
t$desc[1:20]
head(t)
t1 <- t %>%
  dplyr::group_by(gene_symb) %>%
  dplyr::summarise(nb_occurrence=n()) %>%
  dplyr::filter(nb_occurrence>=2)

dim(t1)
t1$gene_symb
t <- de_genes %>%
  dplyr::group_by(gene_symb, mouse_id) %>%
  dplyr::summarise(nb_occurrence=n()) #
head(t)
sum(t$nb_occurrence>=3)
t1 <- t %>%
  dplyr::filter(nb_occurrence>2)
dim(t1)
head(t1[,1:3])
View(t1)
t2 <- t1 %>%
  dplyr::group_by(gene_symb) %>%
  dplyr::summarise(nb_occur=n()) %>%
  dplyr::filter(nb_occur>=2)
summary(as.factor(t1$direction))
head(t2)
dim(t2)
t2
table(t$gene_symb[1:50], t$mouse_id[1:50])
stat <- de_genes %>%
  dplyr::select(avg_log2FC, gene_symb, de_comparison) %>%
  tidyr::pivot_wider(names_from = 'gene_symb', values_from = 'avg_log2FC') %>%
  as.data.frame()
dim(stat)
length(unique(de_genes$gene_symb))
stat[1:3,1:2]
rownames(stat) <- stat$de_comparison
stat$de_comparison <- NULL
stat[1,]
p <- ComplexHeatmap::Heatmap(as.matrix(stat), na_col = "white",
                             show_column_names=F,
                             show_row_names = T,
                             cluster_rows=T,cluster_columns=F,
                             name = "Avg Exp", 
                             # row_order = sort(rownames(test)),
                             # row_split= obs_genes_df,
                             row_title_rot = 0,
                             row_gap = unit(2, "mm"),
                             # column_split = obs_cells_df, 
                             # column_title = "ddd",
                             column_gap = unit(2, "mm"),
                             # column_names_gp = grid::gpar(fontsize = 6),
                             row_names_gp = grid::gpar(fontsize = 15),
                             show_heatmap_legend = T,
                             # top_annotation=anno_col,
                             # left_annotation = left_anno,
                             # cell_fun = cell_func,
                             row_dend_reorder=F,
                             column_title = "Gene", 
                             column_title_side = "bottom")
p

sum(is.na(stat[6,]))
dim(stat)



pathway_fn = '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/pathway_set/h.all.v7.0.symbols.gmt'  
dim(t)
obs_genes_symb <- t$gene_symb
save_data <- F
tolower(stat$reference_set)
'PCLO' %in% de_genes$gene_symb
# [1] "hallmark_androgen_response"                
# [2] "hallmark_myc_targets_v1"                   
# [3] "hallmark_epithelial_mesenchymal_transition"
# [4] "hallmark_tnfa_signaling_via_nfkb"

ref_met_genes <- data.table::fread('/home/htran/Projects/hakwoo_project/metastasis_material/materials/ref_metastasis_related_genes/Metastasis_related_genes_from_Hakwoo.csv')
dim(ref_met_genes)
ref_met_genes <- ref_met_genes %>%
  dplyr::filter(gene_symbol!="")
sum(unique(de_genes$gene_symb) %in% unique(ref_met_genes$gene_symbol))
sum(t1$gene_symb %in% unique(ref_met_genes$gene_symbol))
length(t1$gene_symb)
