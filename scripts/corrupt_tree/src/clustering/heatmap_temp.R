dim(copynumber)
# clones$clone_id <- factor(clones$clone_id, levels = unique(clones$clone_id))

clones <- clones[gtools::mixedorder(as.character(clones$clone_id)),]
clones$clone_id[1:50]
copynumber <- copynumber[,clones$cell_id]
p <- ComplexHeatmap::Heatmap(as.matrix(copynumber), na_col = "black",
                                 show_column_names=F,
                                 show_row_names = F,
                                 cluster_rows=F,cluster_columns=F,
                                 name = "CNA", 
                                 # row_order = sort(rownames(test)),
                                 # row_split= obs_genes_df,
                                 row_title_rot = 0,
                                 # row_gap = unit(2, "mm"),
                                 # column_split = obs_cells_df, 
                                 # column_title = "ddd",
                                 # column_gap = unit(2, "mm"),
                                 column_names_gp = grid::gpar(fontsize = 6),
                                 row_names_gp = grid::gpar(fontsize = 7),
                                 show_heatmap_legend = T,
                                 # top_annotation=anno_col,
                                 # left_annotation = left_anno,
                                 # cell_fun = cell_func,
                                 row_dend_reorder=F,
                                 column_title = "Pseudotime", 
                                 column_title_side = "bottom")
# p
png(paste0(save_dir,"hm.png"), height = 2*1500, width=2*1200, res = 2*72)
print(p)
dev.off()
