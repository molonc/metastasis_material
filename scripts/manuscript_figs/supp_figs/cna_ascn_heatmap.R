suppressPackageStartupMessages({
  require("data.table")
  require("dplyr")
})


source('/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/scripts/corrupt_tree/src/tree_viz/make_cell_copynumber_tree_heatmap.R')

viz_cna_heatmap <- function(){
  datatag <- 'SA919'
  results_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/dlp_trees/SA919/'
  copynumber_fn <- paste0(results_dir,'total_merged_filtered_states.csv.gz')
  clones_fn <- paste0(results_dir,'cell_clones.csv.gz')
  
  grouping_file <- paste0(results_dir,'library_groupings.csv.gz')
  # tree <- read.tree(tree_fn)
  copynumber <- data.table::fread(copynumber_fn) %>% as.data.frame()
  rownames(copynumber) <- copynumber$V1
  copynumber$V1 <- NULL
  head(copynumber[1:3,1:3]) #a matrix with chr desc: chr_start_end in row names, and cell id in colnames: SA919X7XB05691-A96204B-R29-C07: sampleid-libraryid-rowid-colid
  
  
  results_pseudobk_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/results_CNA/ascn_SA919/'
  ascn_fn = paste0(results_pseudobk_dir,'ascn_inference.rds')
  ascn_raw <- readRDS(ascn_fn)
  ascn <- ascn_raw$data
  print(dim(ascn))
  print(colnames(ascn))
  
  clones <- data.table::fread(clones_fn)
  dim(copynumber)
  dim(clones)
  clones <- clones %>%
    dplyr::filter(!clone_id %in% c('Un','unassigned','None'))
  dim(clones)
  
  cells_used <- intersect(clones$cell_id, unique(ascn$cell_id))
  
  ascn <- ascn %>%
    dplyr::filter(cell_id %in% cells_used) %>%
    dplyr::mutate(chr_desc=paste0(chr,'_',start,'_',end))%>%
    dplyr::select(-chr,-start,-end)
  
  
  filtered_bins_used <- intersect(rownames(copynumber),unique(ascn$chr_desc))
  length(filtered_bins_used)
  copynumber <- copynumber[filtered_bins_used, cells_used]
  
  
  
  tree_fn <- paste0(results_dir, 'tree.newick')
  tree <- read.tree(tree_fn)
  tree$tip.label[1]
  # all_cells <- grep('cell', tree$tip.label, value = T)
  # all_cells <- grep('SA', tree$tip.label, value = T)
  all_cells <- tree$tip.label
  length(all_cells)
  none_cells <- clones$cell_id[!clones$cell_id %in% all_cells]
  # none_cells <- all_cells[!all_cells %in% paste0('cell_',clones$cell_id)]
  
  none_cells <- cells_used[!cells_used %in% all_cells]
  length(none_cells)
  # none_cells <- all_cells[!all_cells %in% paste0('cell_',clones$cell_id)]
  
  print(length(none_cells))
  if(length(none_cells)>0){
    tree <- ape::drop.tip(tree, none_cells, trim.internal =T, collapse.singles = T)  
  }
  
  # tree
  save_dir <- paste0(results_dir,'tree_viz/')
  dir.create(save_dir)
  
  png(paste0(save_dir,'cell_cn_tree_heatmap_',datatag,'.png'), height = 2*500, width=2*700, res = 2*72)
  make_cell_copynumber_tree_heatmap_demo(
    tree, copynumber, clones, NULL, grouping_file, save_dir, paste0(datatag, '_cna')
  )
  dev.off()
  
}    

viz_ascn_heatmap <- function(){
  
  datatag <- 'SA919'
    
  results_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/dlp_trees/SA919/'
  copynumber_fn <- paste0(results_dir,'total_merged_filtered_states.csv.gz')
  clones_fn <- paste0(results_dir,'cell_clones.csv.gz')
  
  grouping_file <- paste0(results_dir,'library_groupings.csv.gz')
  # tree <- read.tree(tree_fn)
  copynumber <- data.table::fread(copynumber_fn) %>% as.data.frame()
  rownames(copynumber) <- copynumber$V1
  copynumber$V1 <- NULL
  head(copynumber[1:3,1:3]) #a matrix with chr desc: chr_start_end in row names, and cell id in colnames: SA919X7XB05691-A96204B-R29-C07: sampleid-libraryid-rowid-colid
  
  clones <- data.table::fread(clones_fn)
  dim(copynumber)
  dim(clones)
  clones <- clones %>%
    dplyr::filter(!clone_id %in% c('Un','unassigned','None'))
  dim(clones)
  cells_used <- intersect(clones$cell_id, unique(ascn$cell_id))
  length(cells_used)
  
  
  
  # unique(ascn$state_phase)
  # [1] "A-LOH"    "Balanced" "A-Gained" "B-LOH"  
  ascn_state_phase <- ascn %>%
    # dplyr::filter(cell_id %in% cells_used) %>%
    # dplyr::select(chr, start, end, cell_id, state_phase) %>%
    dplyr::select(chr_desc, cell_id, state_phase) %>%
    # dplyr::mutate(chr_desc=paste0(chr,'_',start,'_',end))%>%
    # dplyr::select(-chr,-start,-end)%>%
    dplyr::filter(chr_desc %in% filtered_bins_used) %>%
    dplyr::mutate(state_phase=   #new term
                    case_when(
                      state_phase=='A-LOH' ~ 'A-Hom',
                      state_phase=='B-LOH' ~ 'B-Hom',
                      TRUE ~ state_phase
                    ))%>%
    tidyr::pivot_wider(names_from='cell_id', values_from = 'state_phase')%>%
    tibble::column_to_rownames('chr_desc')
  
  # ascn_state_phase <- ascn_state_phase %>%
  dim(ascn_state_phase)
  head(ascn_state_phase[1:3,1:3]) 
  
  
  tree_fn <- paste0(results_dir, 'tree.newick')
  tree <- read.tree(tree_fn)
  tree$tip.label[1]
  # all_cells <- grep('cell', tree$tip.label, value = T)
  all_cells <- tree$tip.label
  # all_cells <- grep('SA', tree$tip.label, value = T)
  length(all_cells)
  none_cells <- cells_used[!cells_used %in% all_cells]
  length(none_cells)
  # none_cells <- all_cells[!all_cells %in% paste0('cell_',clones$cell_id)]
  dim(ascn_state_phase)
  dim(clones)
  print(length(none_cells))
  if(length(none_cells)>0){
    tree <- ape::drop.tip(tree, none_cells, trim.internal =T, collapse.singles = T)  
  }
  
  save_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/figures/fig_cnv/'
  # dir.create(save_dir)
  png(paste0(save_dir,'cell_ascn_tree_heatmap_state_phase_',datatag,'.png'), height = 2*500, width=2*700, res = 2*72)
  make_allele_copynumber_tree_heatmap_demo(
    tree, ascn_state_phase, clones, plotcol='state_phase', NULL, grouping_file, save_dir, paste0(datatag, '_ascn_state_phase')
  )
  dev.off()
  
  
  
  # unique(ascn$state_BAF) # ratio of A, B
  # [1] 0.0 0.5 0.1 0.4 0.3 0.2 1.0
  ascn_state_BAF <- ascn %>%
    # dplyr::filter(cell_id %in% cells_used) %>%
    # dplyr::select(chr, start, end, cell_id, state_BAF) %>%
    dplyr::select(chr_desc, cell_id, state_BAF) %>%
    # dplyr::mutate(chr_desc=paste0(chr,'_',start,'_',end))%>%
    # dplyr::select(-chr,-start,-end)%>%
    dplyr::filter(chr_desc %in% filtered_bins_used) %>%
    tidyr::pivot_wider(names_from='cell_id', values_from = 'state_BAF')%>%
    tibble::column_to_rownames('chr_desc')
  
  # ascn_state_phase <- ascn_state_phase %>%
  dim(ascn_state_BAF)
  head(ascn_state_BAF[1:3,1:3]) 
  
  save_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/figures/fig_cnv/'
  # dir.create(save_dir)
  png(paste0(save_dir,'cell_ascn_state_BAF_tree_heatmap_',datatag,'.png'), height = 2*500, width=2*700, res = 2*72)
  make_allele_copynumber_tree_heatmap_demo(
    tree, ascn_state_BAF, clones, plotcol='state_BAF', NULL, grouping_file, save_dir, paste0(datatag, '_ascn_state_BAF')
  )
  dev.off()
  
  p1 <- grid.grabExpr(ComplexHeatmap::draw(p, padding = unit(c(1, 1, 1, 1), "mm")))
  
  p1
  
}

viz_supp_fig3 <- function(){
  save_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/figures/fig_cnv/'
  datatag <- 'SA919'
  pcna <- readRDS(paste0(save_dir, paste0(datatag, '_cna','_hm.rds')))
  pascn_state_phase <- readRDS(paste0(save_dir, paste0(datatag, '_ascn_state_phase','_hm.rds')))
  pascn_state_baf <- readRDS(paste0(save_dir, paste0(datatag, '_ascn_state_BAF','_hm.rds')))
  
  library(grid)
  pcna1 <- grid.grabExpr(ComplexHeatmap::draw(pcna, padding = unit(c(2, 2, 2, 2), "mm"), 
                         annotation_legend_side="right",
                         heatmap_legend_side="top"))
  pascn_state_phase1 <- grid.grabExpr(ComplexHeatmap::draw(pascn_state_phase, padding = unit(c(2, 2, 2, 2), "mm"),
                                      annotation_legend_side="right",
                                      heatmap_legend_side="top")) 
  pascn_state_baf1 <- grid.grabExpr(ComplexHeatmap::draw(pascn_state_baf, padding = unit(c(2, 2, 2, 2), "mm"), 
                                    annotation_legend_side="right",
                                    heatmap_legend_side="top"))

  p_total <- cowplot::plot_grid(pcna1, pascn_state_phase1, pascn_state_baf1, ncol=1)
  # p_total
  
  ggsave(paste0(save_dir,"SUPP_Fig3_ascn.svg"),
         plot = p_total,
         height = 15,
         width = 8,
         # useDingbats=F,
         dpi = 150
  )
  
  
  # ht_list = ht1 %v% ht2 %v% ht3
  # draw(ht_list)
  
}
    