library(dplyr)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(viridis)
library(igraph)
# library(plyr)
library(RColorBrewer)
library(ComplexHeatmap)
viz_marker <- function(cell_profiles_fn, datatag, save_dir,
                       xmax=NULL, ymax=NULL, invert_coord=TRUE){
  # library(dplyr)
  # library(ggplot2)
  
  # input_dir <- '/Users/hoatran/Documents/python_workspace/TissueJ4Proteins/testing_dataset/testing_macros/'
  # save_dir <- '/Users/hoatran/Documents/python_workspace/TissueJ4Proteins/testing_dataset/testing_macros/results/'
  # dir.create(save_dir)
  # datatag <- 'small_tissue'
  # cell_profiles_fn <- paste0(input_dir, 'C6-NUC_cell_profiles.csv')
  
  
  
  nodes_df <- data.table::fread(cell_profiles_fn) %>% as.data.frame()
  # nodes_df$ts <- NULL
  # data.table::fwrite(nodes_df, cell_profiles_fn)
  print(dim(nodes_df))
  # colnames(nodes_df)
  
  background_margin <- 10
  if(is.null(xmax)){
    xmax <- max(nodes_df$x) + background_margin
  }
  if(is.null(ymax)){
    ymax <- max(nodes_df$y) + background_margin
  }
  # invert_coord <- T
  if(invert_coord){
    nodes_df$y <- ymax - nodes_df$y  
  }
  pct_coverage_thres <- 15 # at least 15% of cell zone should be covered by a marker 
  
  ls_markers <- c('C1_BFP','C2_tSapphire','C3_venus','C4_tomato','C5_katushka')
  cns <- colnames(nodes_df)
  cns <- gsub('-','_', cns)
  colnames(nodes_df) <- cns
  
  plt_ls <- list()
  plt_ls[['xmax']] <- xmax
  plt_ls[['ymax']] <- ymax
  plt_ls[['datatag']] <- datatag
  
  dot_size <- 0.0001
  for(m in ls_markers){
    plt_color <- paste0(m,'_mean_intensity_nuc')##_cellzone
    plt_color <- cns[grepl(plt_color, cns)]
    
    plt_var <- paste0(m,'_pct_coverage')
    plt_var <- cns[grepl(plt_var, cns)]
    nodes_filtered_df <- nodes_df %>%
      dplyr::filter(!!sym(plt_var) > pct_coverage_thres)
    print(dim(nodes_filtered_df))
    
    dim(nodes_df)
    # sum(nodes_df1$C1_BFP_mean_corerage>0)
    
    p <- ggplot(nodes_df, aes(x = x, y = y)) + 
      geom_point(color='#C0C0C0', size=dot_size) + # drawing all cells in the tissue area with gray colors
      geom_point(data=nodes_filtered_df, aes_string(color=plt_color), size=dot_size) + 
      # geom_text(aes(label=celltype_desc))+
      # annotate('text', x = genes$xcoord, y = -log10(genes$ycoord), label = genes$celltype_desc)+
      scale_color_viridis(option = "D")+
      theme_bw() + 
      theme(legend.position = "bottom",
            panel.grid = element_blank(),
            axis.ticks = element_blank(),
            axis.text  = element_blank(),
            axis.title = element_blank(),
            plot.margin = unit(c(0.001, 0.001, 0.001, 0.001), "null"),
            panel.spacing = unit(c(0, 0, 0, 0), "null"),
            legend.text = element_text(size=9),
            legend.title = element_text(size=11))+
      guides(color = guide_colourbar(barwidth = 5, barheight = 0.5, title=m, title.position = "left"))
    # p
    
    plg <- cowplot::get_legend(p)
    p <- p + theme(legend.position = "none")
    plt_ls[[m]] <- p
    plt_ls[[paste0(m,'_legend')]] <- plg
    
  }  
  
  ## Plotting nuc
  nodes_df$cell_id <- as.factor(nodes_df$cell_id)
  p_nuc <- ggplot(nodes_df, aes(x = x, y = y)) + 
    geom_point(aes(color=cell_id), size=dot_size, shape=1)+
    theme_bw() + 
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text  = element_blank(),
          axis.title = element_blank(),
          plot.margin = unit(c(0.001, 0.001, 0.001, 0.001), "null"),
          panel.spacing = unit(c(0, 0, 0, 0), "null"))
  
  wd <- plt_ls$xmax/1000
  ht <- plt_ls$ymax/1000
  ## Small datasets
  # plt_labels <- cowplot::plot_grid(NULL, NULL, NULL, NULL, NULL, labels = ls_markers, 
  #                                  nrow = 1, rel_widths = rep(wd, length(ls_markers)),
  #                                  label_size = 11)
  # 
  # plt_main <- cowplot::plot_grid(plt_ls$C1_BFP, plt_ls$C2_tSapphire, plt_ls$C3_venus, 
  #                                plt_ls$C4_tomato, plt_ls$C5_katushka,NULL,  
  #                                nrow = 2, rel_widths = rep(wd, length(ls_markers)))
  # plt_total_lgs <- cowplot::plot_grid(plt_ls$C1_BFP_legend, plt_ls$C2_tSapphire_legend, 
  #                                     plt_ls$C3_venus_legend, plt_ls$C4_tomato_legend, plt_ls$C5_katushka_legend, 
  #                                     nrow = 1, rel_widths = rep(wd, length(ls_markers)))
  plt_labels <- cowplot::plot_grid(NULL, NULL, NULL, NULL, NULL, labels = ls_markers, 
                                   nrow = 2, #rel_widths = rep(wd, length(ls_markers)),
                                   label_size = 11)
  
  plt_main <- cowplot::plot_grid(plt_ls$C1_BFP, plt_ls$C2_tSapphire, plt_ls$C3_venus, 
                                 plt_ls$C4_tomato, plt_ls$C5_katushka, p_nuc,  
                                 nrow = 2) #rel_widths = rep(wd, length(ls_markers)))
  plt_total_lgs <- cowplot::plot_grid(plt_ls$C1_BFP_legend, plt_ls$C2_tSapphire_legend, 
                                      plt_ls$C3_venus_legend, plt_ls$C4_tomato_legend, plt_ls$C5_katushka_legend, 
                                      nrow = 2)#, rel_widths = rep(wd, length(ls_markers)))
  
  p_markers_total <- cowplot::plot_grid(plt_labels, plt_main, NULL, plt_total_lgs, ncol=1,
                                        rel_heights = c(0.2, ht, 0.05, 0.2))
  # p_markers_total
  png(paste0(save_dir,'marker_protein_exp.png'),
      height = 2*ht*10*2+50, width=2*wd*10*3+10, res = 2*200)
  print(p_markers_total)
  dev.off()
  dim(nodes_df)
  png(paste0(save_dir,'marker_protein_exp.png'),
      height = 2*ht*15*2, width=2*wd*15*3, res = 2*200)
  print(plt_main)
  dev.off()
  
  
  return(plt_ls)
}

viz_wholetissue <- function(nodes_df, save_dir, datatag, xmax=NULL, ymax=NULL, 
                            col=NULL, invert_coord=F){
  if(is.null(col)){
    # col <- colorRampPalette(brewer.pal(8, "Set2"))(dim(clone_df)[1])
    # clone_df <- data.table::fread(paste0(save_dir,'predefined_clones.csv')) %>% as.data.frame()
    # names(col) <- clone_df$clone_id
    script_dir <- '/Users/hoatran/Documents/jean_project/data/script/'
    metaclone_df <- data.table::fread(paste0(script_dir,'predefined_clones_v2.csv')) %>% as.data.frame()
    # head(metaclone_df)
    col <- metaclone_df$cluster_color
    names(col) <- metaclone_df$clone_id
  }  
  
  nodes_df <- nodes_df %>%
    dplyr::filter(cluster_label!='Clone_0')
  dim(nodes_df)
  
  
  # summary(nodes_df$y)
  # unique(nodes_df$cluster_label)
  if(is.null(xmax)){
    xmax <- max(nodes_df$x) + 10
  }
  if(is.null(ymax)){
    ymax <- max(nodes_df$y) + 10
  }
  if(invert_coord){
    nodes_df$y <- ymax - nodes_df$y  
  }
  nodes_df <- nodes_df[gtools::mixedorder(nodes_df$cluster_label),]
  if(dim(nodes_df)[1]<10000){
    cell_size <- 2
    ht <- round(ymax/3,0)
    wd <- round(xmax/3,0)
  }else{
    cell_size <- 0.001
    ht <- round(ymax/10,0)
    wd <- round(xmax/10,0)
  }
  print(unique(nodes_df$cluster_label))
  p <- ggplot(nodes_df, aes(x = x, y = y)) + 
    geom_point(aes(color = cluster_label), size=cell_size) + 
    # geom_text(aes(label=celltype_desc))+
    # annotate('text', x = genes$xcoord, y = -log10(genes$ycoord), label = genes$celltype_desc)+
    scale_color_manual(values = col) +
    theme_bw(base_size = 12) + 
    ylim(0, ymax) + 
    xlim(0, xmax) + 
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text  = element_blank(),
          axis.title = element_blank(),
          panel.background = element_rect(fill = "black",
                                          colour = "black"))
  # p
  png(paste0(save_dir,datatag,'_nodes.png'), 
      height = 2*ht, width=2*wd,res = 2*72)
  print(p)
  dev.off()
  
  nodes_stat <- nodes_df %>% 
    dplyr::group_by(cluster_label) %>% 
    dplyr::summarise(nb_cells=n())
  names(col) <- gsub('_','',names(col))
  nodes_stat$cluster_label <- gsub('_','',nodes_stat$cluster_label)
  nodes_stat <- nodes_stat[gtools::mixedorder(nodes_stat$cluster_label),]
  nodes_stat$cluster_label <- factor(nodes_stat$cluster_label, levels = nodes_stat$cluster_label)
  nodes_stat$clone_info <- paste0(nodes_stat$cluster_label,' (',nodes_stat$nb_cells,')')
  # nodes_stat$clone_info
  # nodes_stat$cluster_label
  map_clones <- nodes_stat$clone_info
  names(map_clones) <- nodes_stat$cluster_label
  cls <- nodes_stat$cluster_label
  nodes_stat$number_of_cells <- ifelse(nodes_stat$nb_cells>0, log10(nodes_stat$nb_cells), 0.01)
  p1 <- ggplot(nodes_stat, aes(x=factor(cluster_label, levels = cls), 
                               y=number_of_cells, fill=cluster_label)) +
    geom_bar(stat="identity", width = 0.5)+
    theme_bw(base_size = 12) + 
    scale_fill_manual(values = col, labels=map_clones[names(col)]) + 
    # scale_y_continuous(breaks=c(0,100,500, 5000, 10000, 15000)) + 
    scale_y_continuous(breaks=scales::breaks_pretty(n = 8)) + 
    theme(axis.text.x = element_text(size=9, angle=90),
          legend.position = 'bottom',
          legend.title = element_blank(),
          legend.text = element_text(size=7)) + 
    labs(x='Clone Id', y= 'log10(# cells)', title=paste0('Summary cell type in - ',datatag))
  # p1
  png(paste0(save_dir,datatag, 'cell_types_summary.png'), 
      height = 2*850, width=2*850,res = 2*72)
  print(p1)
  dev.off()
  metaclone_df$clone_id <- gsub('_','',metaclone_df$clone_id)
  nodes_stat$clone_info <- NULL
  nodes_stat <- nodes_stat %>% left_join(metaclone_df, by=c('cluster_label'='clone_id'))
  data.table::fwrite(nodes_stat, paste0(save_dir,datatag, '_cell_types_summary.csv'))
  
  
  nodes_stat$clone_desc <- paste0(nodes_stat$cluster_label, ' (',nodes_stat$nb_cells,'): ',nodes_stat$clone_desc)
  map_clones2 <- nodes_stat$clone_desc
  names(map_clones2) <- nodes_stat$cluster_label
  cls <- nodes_stat$cluster_label
  plg <- ggplot(nodes_stat, aes(x=factor(cluster_label, levels = cls), 
                               y=log10(nb_cells), fill=cluster_label)) +
    geom_bar(stat="identity", width = 0.5)+
    theme_bw(base_size = 12) + 
    scale_fill_manual(values = col, labels=map_clones2[names(col)])+ 
    scale_y_continuous(breaks=scales::breaks_pretty(n = 8)) + 
    theme(axis.text.x = element_text(size=9, angle=90),
          legend.position = 'bottom',
          legend.title = element_blank(),
          legend.text = element_text(size=9)) + 
    labs(x='Clone Id', y= 'log10(# cells)', title=paste0('# cells - ',datatag))
  png(paste0(save_dir,datatag, 'cell_types_summary_clone_labels.png'), 
      height = 2*850, width=2*1400,res = 2*72)
  print(plg)
  dev.off()
  # plg11 <- cowplot::ggdraw() + cowplot::draw_plot(cowplot::get_legend(plg))
  
  ls_markers <- c('C1_BFP','C2_tSapphire','C3_venus','C4_tomato','C5_katushka')
  for(m in ls_markers){
    plt_color <- paste0(m,'_mean_intensity_cellzone')
    plt_var <- paste0(m,'_pct_coverage_cellzone')
    nodes_df1 <- nodes_df
    nodes_df1[,plt_color] <- ifelse(nodes_df1[,plt_var]>0,nodes_df1[,plt_color],0)
    
    # sum(nodes_df1$C1_BFP_mean_corerage>0)
    p <- ggplot(nodes_df1, aes(x = x, y = y)) + 
      geom_point(aes_string(color = plt_color), size=cell_size) + 
      # geom_text(aes(label=celltype_desc))+
      # annotate('text', x = genes$xcoord, y = -log10(genes$ycoord), label = genes$celltype_desc)+
      scale_color_viridis(option = "D")+
      theme_bw(base_size = 12) + 
      theme(legend.position = "bottom",
            panel.grid = element_blank(),
            axis.ticks = element_blank(),
            axis.text  = element_blank(),
            axis.title = element_blank())
    # p
    png(paste0(save_dir,str_sub(plt_color, 1, 6),'_exp.png'), 
        height = 2*ht+70, width=2*wd,res = 2*72)
    print(p)
    dev.off()
  }
  
  
}
get_cell_type_stat <- function(nodes_df1, tag, thres_vol=0.3, min_cells_counts=5){
  nodes_df1 <- nodes_df1 %>%
    dplyr::filter(pixvol>=50)
  dim(nodes_df1)
  nodes_df1 <- get_celltype_v2(nodes_df1, save_dir, datatag, 
                               save_data=T, thres_vol_marker=thres_vol)
  # dim(nodes_df1)
  # summary(as.factor(nodes_df1$cluster_label))
  clones_use <- nodes_df1 %>%
    group_by(cluster_label) %>%
    dplyr::summarise(nb_cells=n()) %>%
    dplyr::filter(nb_cells>=min_cells_counts) %>%
    dplyr::pull(cluster_label)
  nodes_df1 <- nodes_df1 %>%
    dplyr::filter(cluster_label %in% clones_use)
  nodes_df1$pcell_id <- paste0(tag,'_', nodes_df1$cell_id)
  return(nodes_df1)
}

# Get Spearman correlation between proportion of cells clone labels in image and in facs method
get_correlation <- function(nodes_df, facs_df){
  facs <- data.table::fread(facs_fn) %>% as.data.frame()
  facs <- facs %>%
    dplyr::mutate(clone_id=paste0('Clone_',clone_id))
  cr <- cor(facs$pct, ct$pct, method = 'spearman')
  return(cr)
}

# pct_ct <- runif(10)
# set.seed(123)
# mat = matrix(rnorm(100), 10)
# rownames(mat) = paste0("R", 1:10)
# colnames(mat) = paste0("R", 1:10)
# hm <- viz_heatmap_cellphenotype(mat, pct_ct)
# mat: matrix of cell-cell interaction
# pct_ct: value for barplot percentage of each cell type
viz_heatmap_cellphenotype <- function(mat, pct_ct, save_dir, col=NULL){
  # library(ComplexHeatmap)
  # library(circlize)
  # library(RColorBrewer)
  nbcelltypes=dim(mat)[1]
  if(is.null(col)){
    col <- colorRampPalette(brewer.pal(8, "Set2"))(nbcelltypes)
    # names(col) <- paste0("R", 1:10)  
    names(col) <- rownames(mat)
  }
  # pct_ct <- runif(10)
  col_top = columnAnnotation(cell_phenotype = rownames(mat), col = list(cell_phenotype = col),
                             show_legend = F, show_annotation_name=F)
  row_right = rowAnnotation(bar2 = anno_barplot(pct_ct, width = unit(3, "cm"), 
                                                gp = gpar(fill = col)),annotation_label='Log2(cell counts)')
  row_left = rowAnnotation(cell_phenotype = rownames(mat), 
                           col = list(cell_phenotype = col),show_annotation_name=F)
  # draw(row_left)
  hm <- Heatmap(mat, name = "Log2(cell-cell contact + 1)", 
                column_title = "Cell Phenotype",
                row_title = "Cell Phenotype",
                left_annotation = row_left,
                right_annotation = row_right, 
                top_annotation = col_top,
                cluster_rows = F, cluster_columns = F)
  png(paste0(save_dir,"hm_cell2cell_connection.png"), height = 2*400, width=2*700,res = 2*72)
  print(hm)
  dev.off()
  return(hm)
}
# g <- readRDS(paste0(save_dir,'graph_ct.rds'))
# res <- get_edges_nodes_from_igraph(g)
# edges <- res$edges
# nodes <- res$nodes
# Based on cell type
compute_cell_cell_interaction <- function(edges, nodes, save_dir){
  cts <- unique(nodes$cellphenotype)
  cts <- gtools::mixedsort(cts)
  connect_mtx <- matrix(0, nrow = length(cts), ncol = length(cts))  
  dim(connect_mtx)
  colnames(connect_mtx) <- cts
  rownames(connect_mtx) <- cts
  
  for(ct in cts){
    obs_nodes <- nodes %>%
      dplyr::filter(cellphenotype==ct)%>%
      dplyr::pull(name)
    edges_tmp <- edges %>%
      dplyr::filter(from %in% obs_nodes)
    cts2 <- cts[cts!=ct]
    for(ct2 in cts2){
      obs_nodes2 <- nodes %>%
        dplyr::filter(cellphenotype==ct2)%>%
        dplyr::pull(name)
      edges_tmp2 <- edges_tmp %>%
        dplyr::filter(to %in% obs_nodes2)
      connect_mtx[ct,ct2] <- dim(edges_tmp2)[1]
    }  
  }
  print(connect_mtx)
  # View(connect_mtx)
  data.table::fwrite(connect_mtx, paste0(save_dir,'connect_mtx.csv'), quote=F)
  return(connect_mtx)
}


# get_clone_label(1, 0, 0, 1, 0)

get_clone_label <- function(e, ts, v, td, k){
  markers <- c('eBFP2','tSapphire','Venus','tdTomato','Katushka')
  # vals <- c(0,1)
  # ls_clones <- list()  
  lb <- ''
  s <- 0
  if(e==1){
    lb <- paste0(lb, markers[1])
    s <- s + 2^0
  }
  if(ts==1){
    lb <- paste0(lb,', ',markers[2])
    s <- s + 2^1
  }
  if(v==1){
    lb <- paste0(lb,', ',markers[3])
    s <- s + 2^2
  }
  if(td==1){
    lb <- paste0(lb,', ',markers[4])
    s <- s + 2^3
  }
  if(k==1){
    lb <- paste0(lb,', ',markers[5])
    s <- s + 2^4
  }
  if(grepl('^, ', lb)){
    lb <- str_sub(lb, 3, str_length(lb))
  }
  if(s==0 & lb==''){
    lb <- 'unlabeled'
  }
  # res <- list(clone_label=paste0('Clone_',s),clone_desc=lb)
  return(paste0('Clone_',s))
}




viz_wholetissue <- function(nodes_df, save_dir, xmax=NULL, ymax=NULL, col=NULL, invert_coord=F){
  if(is.null(col)){
    # col <- colorRampPalette(brewer.pal(8, "Set2"))(dim(clone_df)[1])
    clone_df <- data.table::fread(paste0(save_dir,'predefined_clones_v2.csv')) %>% as.data.frame()
    col <- clone_df$cluster_color
    names(col) <- clone_df$clone_id
  }  
  
  nodes_df <- nodes_df %>%
    dplyr::filter(cluster_label!='Clone_0')
  dim(nodes_df)
  
  
  # summary(nodes_df$y)
  # unique(nodes_df$cluster_label)
  if(is.null(xmax)){
    xmax <- max(nodes_df$x) + 10
  }
  if(is.null(ymax)){
    ymax <- max(nodes_df$y) + 10
  }
  if(invert_coord){
    nodes_df$y <- ymax - nodes_df$y  
  }
  nodes_df <- nodes_df[gtools::mixedorder(nodes_df$cluster_label),]
  p <- ggplot(nodes_df, aes(x = x, y = y)) + 
    geom_point(aes(color = cluster_label)) + 
    # geom_text(aes(label=celltype_desc))+
    # annotate('text', x = genes$xcoord, y = -log10(genes$ycoord), label = genes$celltype_desc)+
    scale_color_manual(values = col) +
    theme_bw(base_size = 12) + 
    ylim(0, ymax) + 
    xlim(0, xmax) + 
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text  = element_blank(),
          axis.title = element_blank())
  # p
  png(paste0(save_dir,datatag,'_celltypes.png'), 
      height = 2*round(ymax/3,0), width=2*round(xmax/3,0),res = 2*72)
  print(p)
  dev.off()
  
  cluster_ls <- gtools::mixedsort(unique(nodes_df$cluster_label))
  nodes_df$cluster_label <- factor(nodes_df$cluster_label, levels = cluster_ls)
  
  p1 <- ggplot(nodes_df, aes(x = x, y = y)) + 
    geom_point(aes(color = cluster_label)) + 
    # geom_text(aes(label=celltype_desc))+
    # annotate('text', x = genes$xcoord, y = -log10(genes$ycoord), label = genes$celltype_desc)+
    scale_color_manual(values = col) +
    theme_bw(base_size = 12) + 
    theme(legend.position = "bottom",
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text  = element_blank(),
          axis.title = element_blank()) +
    guides(color = guide_legend(title="",override.aes = list(shape = 15, size=3), nrow = 3))
  # p1
  plg <- cowplot::get_legend(p1)
  png(paste0(save_dir,datatag,'_celltypes_legend.png'), height = 150, width=900,res = 2*72)
  print(cowplot::ggdraw(plg))
  dev.off()
  
  
  
  
  
  
  colnames(nodes_df1)
  ls_markers <- c('C1_BFP','C2_tSapphire','C3_Venus','C4_Tomato','C5_Katushka')
  for(m in ls_markers){
    plt_color <- paste0(m,'_mean_intensity')
    plt_var <- gsub('intensity','corerage',plt_color)
    nodes_df1 <- nodes_df
    nodes_df1[,plt_color] <- ifelse(nodes_df1[,plt_var]>0,nodes_df1[,plt_color],0)
    
    # sum(nodes_df1$C1_BFP_mean_corerage>0)
    p <- ggplot(nodes_df1, aes(x = x, y = y)) + 
      geom_point(aes_string(color = plt_color)) + 
      # geom_text(aes(label=celltype_desc))+
      # annotate('text', x = genes$xcoord, y = -log10(genes$ycoord), label = genes$celltype_desc)+
      scale_color_viridis(option = "D")+
      theme_bw(base_size = 12) + 
      theme(legend.position = "bottom",
            panel.grid = element_blank(),
            axis.ticks = element_blank(),
            axis.text  = element_blank(),
            axis.title = element_blank())
    # p
    png(paste0(input_dir,'testing_cell_type/CT/',str_sub(plt_color, 1, 6),'_exp.png'), 
        height = 2*round(ymax/3,0)+70, width=2*round(xmax/3,0),res = 2*72)
    print(p)
    dev.off()
  }
  
  
}


# g_backup <- g
# TO DO: add average coordinates here 
collapse_cluster <- function(g, save_dir, col=NULL){
  if(is.null(V(g)$name)){
    g <- g %>%
      set_vertex_attr("name", value = paste0('V',rep(1:length(V(g)),1)))
  }
  edges <- igraph::as_data_frame(g, what="edges")
  # edges = get.edgelist(g) %>% as.data.frame()
  # colnames(edges) <- c('from','to')
  print(head(edges))
  nodes_df <- igraph::as_data_frame(g, what = c("vertices"))
  print(head(nodes_df))
  print(colnames(nodes_df))
  dim(nodes_df)
  
  # nodes_df <- data.frame(node_name=V(g)$name, cluster=V(g)$cluster_label,
  #                        celltype=V(g)$cellphenotype, stringsAsFactors=F)
  clusters <- unique(nodes_df$subcluster_label)
  clusters <- clusters[!is.na(clusters)]
  print(summary(as.factor(nodes_df$subcluster_label)))
  nodes_ls <- tibble::tibble()
  edges_ls <- tibble::tibble()
  
  
  # comp <- gtools::combinations(length(clusters),2,clusters)
  # class(comp)
  # dim(comp)
  for(obs_cls in clusters){
    nodes_source <- nodes_df %>%
      dplyr::filter(subcluster_label==obs_cls) 
    # %>%  dplyr::pull(node_name)
    s <- nodes_source %>%   # collapse a cluster to a node with mean coordinate values
      dplyr::summarise(x=mean(x),
                       y=mean(y))
    ntmp <- tibble(node_name=paste0('C_',obs_cls),node_size=length(nodes_source$name), 
                   celltype=nodes_source$cellphenotype[1], cluster=obs_cls, x=s$x, y = s$y)
    nodes_ls <- dplyr::bind_rows(nodes_ls, ntmp)
    # edges_obs <- edges %>%
    #   dplyr::filter(from %in% nodes_source)
    # for(cls in clusters[clusters!=obs_cls]){
    #   nodes_target <- nodes_df %>%
    #     dplyr::filter(subcluster_label==cls) %>%
    #     dplyr::pull(name)
    #   edges_tmp <- edges_obs %>%
    #     dplyr::filter(to %in% nodes_target)
    #   if(dim(edges_tmp)[1]>0){
    #     etmp <- tibble(from=paste0('C_',obs_cls),to=paste0('C_',cls),edge_width=dim(edges_tmp)[1])
    #     edges_ls <- dplyr::bind_rows(edges_ls, etmp)
    #   }
    # }
    
  }
  nodes_ls <- as.data.frame(nodes_ls)
  # edges_ls <- as.data.frame(edges_ls)
  dim(nodes_ls)
  head(nodes_ls)
  rownames(nodes_ls) <- nodes_ls$node_name
  summary(nodes_ls$node_size)
  dim(nodes_ls)
  data.table::fwrite(nodes_ls, paste0(save_dir,'collapse_',datatag,'_nodes.csv'))
  
  # collapsed_g <- igraph::graph_from_edgelist(as.matrix(edges_ls[,c('from','to')]), directed = F)#
  # V(collapsed_g)$id <- seq(vcount(collapsed_g))
  # V(collapsed_g)$nb_cells <- nodes_ls[V(collapsed_g)$name,'node_size']
  # V(collapsed_g)$celltype <- nodes_ls[V(collapsed_g)$name,'celltype']
  # V(collapsed_g)$cluster <- nodes_ls[V(collapsed_g)$name,'cluster']
  # celltypes <- unique(V(collapsed_g)$celltype)
  # if(is.null(col)){
  #   col <- colorRampPalette(brewer.pal(8, "Set2"))(length(celltypes))
  #   names(col) <- celltypes
  # }  
  # edges_list <- get.edgelist(collapsed_g) %>% as.data.frame()
  # edges_list$desc <- paste0(edges_list$V1, edges_list$V2)
  # edges_ls$desc <- paste0(edges_ls$from, edges_ls$to)
  # rownames(edges_ls) <- edges_ls$desc
  # edges_list$edge_width <- edges_ls[edges_list$desc,'edge_width'] # a bit manual but ok
  # # E(g) - better way to include info is based on this vector
  # png(paste0(save_dir,"test_graph.png"), height = 2*400, width=2*400,res = 2*72)
  # plot(collapsed_g, vertex.color=col[V(collapsed_g)$celltype], vertex.size=3*log2(V(collapsed_g)$nb_cells), 
  #      vertex.label=NA, edge.width=edges_list$edge_width, edge.arrow.size=0, vertex.label.dist=3)
  # dev.off()
  # 
  # 
  # return(g)
  ## TO DO: return nodes and centroid coordinates here
}  


# cells with same phenotype and connected together form a sub population
# nodes_df: cluster and node_name

get_clusters <- function(edges_df, nodes_df, save_dir, col=NULL){
  # rownames(nodes_df) <- nodes_df$node_name
  # V(g)$cellphenotype <- nodes_df[V(g)$name,'cluster']
  
  t <- summary(as.factor(nodes_df$cluster_label))
  cluster_ls <- names(t)[t>1]  # get only communities with more than 1 cell
  print(cluster_ls)
  
  nodes_df <- nodes_df %>%
    dplyr::filter(cluster_label %in% cluster_ls)
  dim(nodes_df)
  rownames(nodes_df) <- nodes_df$cell_id
  edges_df <- edges_df %>%
    dplyr::filter(from %in% nodes_df$cell_id & to %in% nodes_df$cell_id)  # TO DO: need to think about connection here
  dim(edges_df)
  g <- read_image_tree(edges_df)
  nodes_df <- nodes_df %>%
    dplyr::filter(cell_id %in% V(g)$name)
  dim(nodes_df)
  # summary(as.factor(nodes_df$cluster_label))
  idx_cls <- 0
  nodes_ls <- tibble::tibble()
  # obs_cls <- 'Clone_2'
  for(obs_cls in cluster_ls){
    nodes_tmp <- nodes_df %>%
      dplyr::filter(cluster_label==obs_cls)%>%
      dplyr::pull(cell_id)
    print(length(nodes_tmp))
    print(obs_cls)
    # g1 <- induced_subgraph(g, paste0('V',nodes_tmp))
    g1 <- induced_subgraph(g, nodes_tmp)  # get cells with same phenotype
    cm <- components(g1) # get all connected components
    
    if(idx_cls==0){
      ntmp <- tibble(cell_id=V(g1)$name,subcluster_label=cm$membership)
    }else{
      ntmp <- tibble(cell_id=V(g1)$name,subcluster_label=idx_cls+cm$membership)
    }
    print(dim(ntmp))
    # head(ntmp)
    nodes_ls <- dplyr::bind_rows(nodes_ls, ntmp)
    idx_cls <- idx_cls + cm$no
  }
  nodes_ls <- nodes_ls %>% as.data.frame()
  print(dim(nodes_ls))
  rownames(nodes_ls) <- nodes_ls$cell_id
  print(summary(as.factor(nodes_ls$subcluster_label)))
  V(g)$subcluster_label <- nodes_ls[V(g)$name,'subcluster_label']
  V(g)$cellphenotype <- nodes_df[V(g)$name,'cluster_label']
  V(g)$cellphenotype_desc <- nodes_df[V(g)$name,'clone_desc']
  V(g)$x <- nodes_df[V(g)$name,'x']
  V(g)$y <- nodes_df[V(g)$name,'y']
  
  # if(is.null(col)){
  #   col <- colorRampPalette(brewer.pal(8, "Set2"))(length(unique(nodes_ls$subcluster_label)))
  #   names(col) <- unique(nodes_ls$subcluster_label)
  # }
  # # col_use <- col[nodes_ls$cluster_label]
  # 
  # png(paste0(save_dir,"test_graph.png"), height = 2*400, width=2*400,res = 2*72)
  # # plot(wc, g)
  # # plot(nodes_df$cluster, g)
  # plot.igraph(g, vertex.color=col[V(g)$cluster_label])
  # dev.off()
  saveRDS(g, paste0(save_dir,'graph_ct.rds'))
  # g <- readRDS(paste0(save_dir,'graph_ct.rds'))
  return(g)
}  

read_image_tree <- function(edge_list) {
  # Find the root
  g <- igraph::graph_from_edgelist(as.matrix(edge_list))
  V(g)$id <- seq(vcount(g))
  print(g)
  return(g)
}

# get_celltype_v2 <- function(nodes_df, save_dir, datatag, save_data=T, thres_vol_marker=0.2){
#   if(!file.exists(save_dir)){
#     file.create(save_dir)
#   }
#   # thrsBFP <- 20
#   # thrstSapphire <- 50
#   # thrsVenus <- 50
#   # thrsTomato <- 70
#   # thrsKatushka <- 20
#   
#   # summary(nodes_df$C5_Katushka_mean_intensity)
#   # No intensity info here
#   nodes_df$C1_BFP_mean_corerage <- ifelse(nodes_df$C1_BFP_mean_corerage>=thres_vol_marker,1,0)
#   nodes_df$C2_tSapphire_mean_corerage <- ifelse(nodes_df$C2_tSapphire_mean_corerage>=thres_vol_marker,1,0)
#   nodes_df$C3_Venus_mean_corerage <- ifelse(nodes_df$C3_Venus_mean_corerage>=thres_vol_marker,1,0)
#   nodes_df$C4_Tomato_mean_corerage <- ifelse(nodes_df$C4_Tomato_mean_corerage>=thres_vol_marker,1,0)
#   nodes_df$C5_Katushka_mean_corerage <- ifelse(nodes_df$C5_Katushka_mean_corerage>=thres_vol_marker,1,0)
#   
#   # nodes_backup <- nodes_df  
#   # nodes_df <- nodes_backup
#   nodes_df$cluster_label <- 'unlabeled'
#   nodes_df$clone_desc <- 'unlabeled'
#   nbcores <- 4
#   labels <- parallel::mclapply(1:dim(nodes_df)[1], function(i) {
#     clone_label <- get_clone_label(nodes_df[i,'C1_BFP_mean_corerage'], nodes_df[i,'C2_tSapphire_mean_corerage'], 
#                                    nodes_df[i,'C3_Venus_mean_corerage'], nodes_df[i,'C4_Tomato_mean_corerage'], 
#                                    nodes_df[i,'C5_Katushka_mean_corerage'])
#     return(clone_label)
#   }, mc.cores = nbcores)
#   nodes_df$cluster_label <- unlist(labels)
#   # for(i in 1:dim(nodes_df)[1]){
#   #   # get_clone_label(e, ts, v, td, k) ##c('eBFP2','tSapphire','Venus','tdTomato','Katushka')
#   #   res <- get_clone_label(nodes_df[i,'C1_BFP_mean_corerage'], nodes_df[i,'C2_tSapphire_mean_corerage'], 
#   #                          nodes_df[i,'C3_Venus_mean_corerage'], nodes_df[i,'C4_Tomato_mean_corerage'], 
#   #                          nodes_df[i,'C5_Katushka_mean_corerage'])
#   #   nodes_df[i,'cluster_label'] <- res$clone_label
#   #   nodes_df[i,'clone_desc'] <- res$clone_desc
#   # }
#   print(summary(as.factor(nodes_df$cluster_label)))
#   
#   # data.table::fwrite(nodes_df, paste0(save_dir, datatag,'_nodes_allcells.csv'))
#   # nodes_df <- data.table::fread(paste0(save_dir, datatag,'_nodes_allcells.csv')) %>% as.data.frame()
#   nodes_df$cell_id <- paste0('V',nodes_df$cell_id)
#   nodes_df <- nodes_df %>%
#     dplyr::filter(cluster_label!='Clone_0')  # 0 is unlabelled cells
#   print(dim(nodes_df))
#   
#   
#   if(save_data){
#     data.table::fwrite(nodes_df, paste0(save_dir, datatag,'_nodes_ct.csv.gz'))
#   }
#   return(nodes_df)
#   
# }

# save_dir <- paste0(input_dir,'testing_cell_type/CT/')
# nodes_df <- data.table::fread(paste0(save_dir,datatag,'_nodes.csv')) %>% as.data.frame()
# edges_df <- data.table::fread(paste0(save_dir,datatag,'_edges.csv')) %>% as.data.frame()
# head(edges_df)

# get_celltype <- function(nodes_df, edges_df, 
#                          save_dir, datatag, 
#                          save_data=T, thres_vol_marker=0.2){
#   if(!file.exists(save_dir)){
#     file.create(save_dir)
#   }
#   thrsBFP <- 20
#   thrstSapphire <- 50
#   thrsVenus <- 50
#   thrsTomato <- 70
#   thrsKatushka <- 20
#   # summary(nodes_df$C5_Katushka_mean_intensity)
#   nodes_df$C1_BFP_mean_corerage <- ifelse(nodes_df$C1_BFP_mean_corerage>=thres_vol_marker &
#                                           nodes_df$C1_BFP_mean_intensity>thrsBFP,1,0)
#   nodes_df$C2_tSapphire_mean_corerage <- ifelse(nodes_df$C2_tSapphire_mean_corerage>=thres_vol_marker &
#                                                   nodes_df$C2_tSapphire_mean_intensity>thrstSapphire
#                                                 ,1,0)
#   nodes_df$C3_Venus_mean_corerage <- ifelse(nodes_df$C3_Venus_mean_corerage>=thres_vol_marker &
#                                               nodes_df$C3_Venus_mean_intensity>thrsVenus,1,0)
#   nodes_df$C4_Tomato_mean_corerage <- ifelse(nodes_df$C4_Tomato_mean_corerage>=thres_vol_marker &
#                                                nodes_df$C4_Tomato_mean_intensity>thrsTomato,1,0)
#   nodes_df$C5_Katushka_mean_corerage <- ifelse(nodes_df$C5_Katushka_mean_corerage>=thres_vol_marker &
#                                                  nodes_df$C5_Katushka_mean_intensity>thrsKatushka,1,0)
#   
#   # nodes_backup <- nodes_df  
#   # nodes_df <- nodes_backup
#   nodes_df$cluster_label <- 'unlabeled'
#   nodes_df$clone_desc <- 'unlabeled'
#   nbcores <- 8
#   labels <- parallel::mclapply(1:dim(nodes_df)[1], function(i) {
#     clone_label <- get_clone_label(nodes_df[i,'C1_BFP_mean_corerage'], nodes_df[i,'C2_tSapphire_mean_corerage'], 
#                     nodes_df[i,'C3_Venus_mean_corerage'], nodes_df[i,'C4_Tomato_mean_corerage'], 
#                     nodes_df[i,'C5_Katushka_mean_corerage'])
#     return(clone_label)
#   }, mc.cores = nbcores)
#   nodes_df$cluster_label <- unlist(labels)
#   # for(i in 1:dim(nodes_df)[1]){
#   #   # get_clone_label(e, ts, v, td, k) ##c('eBFP2','tSapphire','Venus','tdTomato','Katushka')
#   #   res <- get_clone_label(nodes_df[i,'C1_BFP_mean_corerage'], nodes_df[i,'C2_tSapphire_mean_corerage'], 
#   #                          nodes_df[i,'C3_Venus_mean_corerage'], nodes_df[i,'C4_Tomato_mean_corerage'], 
#   #                          nodes_df[i,'C5_Katushka_mean_corerage'])
#   #   nodes_df[i,'cluster_label'] <- res$clone_label
#   #   nodes_df[i,'clone_desc'] <- res$clone_desc
#   # }
#   print(summary(as.factor(nodes_df$cluster_label)))
#   
#   # data.table::fwrite(nodes_df, paste0(save_dir, datatag,'_nodes_allcells.csv'))
#   # nodes_df <- data.table::fread(paste0(save_dir, datatag,'_nodes_allcells.csv')) %>% as.data.frame()
#   nodes_df$cell_id <- paste0('V',nodes_df$cell_id)
#   nodes_df <- nodes_df %>%
#     dplyr::filter(cluster_label!='Clone_0')  # 0 is unlabelled cells
#   print(dim(nodes_df))
#   
#   
#   edges_df$from <- paste0('V',edges_df$from)
#   edges_df$to <- paste0('V',edges_df$to)
#   
#   edges_df <- edges_df %>%
#     dplyr::filter(from %in% nodes_df$cell_id & to %in% nodes_df$cell_id)
#   dim(edges_df)
#   res <- list(nodes=nodes_df, edges=edges_df)
#   if(save_data){
#     data.table::fwrite(nodes_df, paste0(save_dir, datatag,'_nodes_ct.csv.gz'))
#     data.table::fwrite(edges_df, paste0(save_dir, datatag,'_edges_ct.csv.gz'))
#   }
#   return(res)
#   
# }
get_celltype_v4 <- function(nodes_df, edges_df, 
                         save_dir, datatag, 
                         celltype_by_pct=T, save_data=T, thres_vol_marker=0.2, 
                         thrsBFP=20, thrstSapphire=50, thrsVenus=50,
                         thrsTomato=70, thrsKatushka=20){
  if(!file.exists(save_dir)){
    file.create(save_dir)
  }
  thres_vol_marker=15 # at least 15% of the marker that cover the a given cell zone area
  # sum(nodes_df$C1_BFP_pct_coverage_cellzone>15)
  # thrsBFP <- 20
  # thrstSapphire <- 50
  # thrsVenus <- 50
  # thrsTomato <- 70
  # thrsKatushka <- 20
  # summary(nodes_df$C5_Katushka_mean_intensity)
  # nodes_df$C1_BFP_mean_corerage <- ifelse(nodes_df$C1_BFP_mean_corerage>=thres_vol_marker &
  #                                           nodes_df$C1_BFP_mean_intensity>thrsBFP,1,0)
  # nodes_df$C2_tSapphire_mean_corerage <- ifelse(nodes_df$C2_tSapphire_mean_corerage>=thres_vol_marker &
  #                                                 nodes_df$C2_tSapphire_mean_intensity>thrstSapphire
  #                                               ,1,0)
  # nodes_df$C3_Venus_mean_corerage <- ifelse(nodes_df$C3_Venus_mean_corerage>=thres_vol_marker &
  #                                             nodes_df$C3_Venus_mean_intensity>thrsVenus,1,0)
  # nodes_df$C4_Tomato_mean_corerage <- ifelse(nodes_df$C4_Tomato_mean_corerage>=thres_vol_marker &
  #                                              nodes_df$C4_Tomato_mean_intensity>thrsTomato,1,0)
  # nodes_df$C5_Katushka_mean_corerage <- ifelse(nodes_df$C5_Katushka_mean_corerage>=thres_vol_marker &
  #                                                nodes_df$C5_Katushka_mean_intensity>thrsKatushka,1,0)
  colnames(nodes_df) <- gsub('-','_',colnames(nodes_df))
  print(paste0("Total number of cells in tissue: ", dim(nodes_df)[1]))
  
  if(celltype_by_pct==TRUE){
    print("Detecting cell type using percentage coverage of each marker within cell zone")
    nodes_df <- nodes_df %>%
      mutate(
        C1_BFP_status = case_when(
          C1_BFP_pct_coverage_cellzone >= thres_vol_marker ~ 1,
          TRUE ~ 0
        ),
        C2_tSapphire_status = case_when(
          C2_tSapphire_pct_coverage_cellzone >= thres_vol_marker ~ 1,
          TRUE ~ 0
        ),
        C3_venus_status = case_when(
          C3_venus_pct_coverage_cellzone >= thres_vol_marker ~ 1,
          TRUE ~ 0
        ),
        C4_tomato_status = case_when(
          C4_tomato_pct_coverage_cellzone >= thres_vol_marker ~ 1,
          TRUE ~ 0
        ),
        C5_katushka_status = case_when(
          C5_katushka_pct_coverage_cellzone >= thres_vol_marker ~ 1,
          TRUE ~ 0
        )
      )  
  }else{
    print("Detecting cell type using percentage coverage of each marker within cell zone, 
          and the mean intensity of each marker within cell zone")
    if(thrsBFP<=0 || thrstSapphire<=0 || thrsVenus<=0 || thrsTomato<=0 || thrsKatushka<=0){
      stop("The threshold values should greater than 0, please double check input data")
    }
    nodes_df %>%
      mutate(
        C1_BFP_status = case_when(
          C1_BFP_pct_coverage_cellzone >= thres_vol_marker & C1_BFP_mean_intensity_cellzone >=thrsBFP ~ 1,
          TRUE ~ 0
        ),
        C2_tSapphire_status = case_when(
          C2_tSapphire_pct_coverage_cellzone >= thres_vol_marker & C2_tSapphire_mean_intensity_cellzone >=thrstSapphire ~ 1,
          TRUE ~ 0
        ),
        C3_venus_status = case_when(
          C3_venus_pct_coverage_cellzone >= thres_vol_marker & C3_venus_mean_intensity_cellzone >=thrsVenus  ~ 1,
          TRUE ~ 0
        ),
        C4_tomato_status = case_when(
          C4_tomato_pct_coverage_cellzone >= thres_vol_marker & C4_tomato_mean_intensity_cellzone >=thrsTomato ~ 1,
          TRUE ~ 0
        ),
        C5_katushka_status = case_when(
          C5_katushka_pct_coverage_cellzone >= thres_vol_marker & C5_katushka_mean_intensity_cellzone >=thrsKatushka ~ 1,
          TRUE ~ 0
        )
      )  
  }
  
  
  # nodes_backup <- nodes_df  
  # nodes_df <- nodes_backup
  nodes_df$cluster_label <- 'unlabeled'
  nodes_df$clone_desc <- 'unlabeled'
  nbcores <- 8
  labels <- parallel::mclapply(1:dim(nodes_df)[1], function(i) {
    clone_label <- get_clone_label(nodes_df[i,'C1_BFP_status'], nodes_df[i,'C2_tSapphire_status'], 
                                   nodes_df[i,'C3_venus_status'], nodes_df[i,'C4_tomato_status'], 
                                   nodes_df[i,'C5_katushka_status'])
    return(clone_label)
  }, mc.cores = nbcores)
  nodes_df$cluster_label <- unlist(labels)
  # for(i in 1:dim(nodes_df)[1]){
  #   # get_clone_label(e, ts, v, td, k) ##c('eBFP2','tSapphire','Venus','tdTomato','Katushka')
  #   res <- get_clone_label(nodes_df[i,'C1_BFP_mean_corerage'], nodes_df[i,'C2_tSapphire_mean_corerage'], 
  #                          nodes_df[i,'C3_Venus_mean_corerage'], nodes_df[i,'C4_Tomato_mean_corerage'], 
  #                          nodes_df[i,'C5_Katushka_mean_corerage'])
  #   nodes_df[i,'cluster_label'] <- res$clone_label
  #   nodes_df[i,'clone_desc'] <- res$clone_desc
  # }
  print(summary(as.factor(nodes_df$cluster_label)))
  
  # data.table::fwrite(nodes_df, paste0(save_dir, datatag,'_nodes_allcells.csv'))
  # nodes_df <- data.table::fread(paste0(save_dir, datatag,'_nodes_allcells.csv')) %>% as.data.frame()
  nodes_df$cell_id[1]
  nodes_df$cell_id <- paste0('C',nodes_df$cell_id)
  excluded_cols <- c('C1_BFP_status','C2_tSapphire_status','C3_venus_status',
                     'C4_tomato_status', 'C5_katushka_status')
  nodes_df <- nodes_df %>%
    dplyr::select(-all_of(excluded_cols))
  
  nodes_df <- nodes_df %>%
    dplyr::filter(cluster_label!='Clone_0') %>%  # 0 is unlabelled cells
    dplyr::select(cell_id, x, y, cluster_label, clone_desc, everything())
  print(dim(nodes_df))
  
  
  # edges_df$from <- paste0('C',edges_df$from)
  # edges_df$to <- paste0('C',edges_df$to)
  # 
  # edges_df <- edges_df %>%
  #   dplyr::filter(from %in% nodes_df$cell_id & to %in% nodes_df$cell_id)
  # dim(edges_df)
  # res <- list(nodes=nodes_df, edges=edges_df)
  if(save_data){
    data.table::fwrite(nodes_df, paste0(save_dir, datatag,'_nodes_celltype.csv.gz'))
    # data.table::fwrite(edges_df, paste0(save_dir, datatag,'_edges_celltype.csv.gz'))
  }
  return(res)
  
}
get_cell_type_v3 <- function(nodes_df, edges_df, 
                         save_dir, datatag, 
                         save_data=T, thres_vol_marker=15){
  if(!file.exists(save_dir)){
    file.create(save_dir)
  }
  # BFP: 8/5
  # Sapphire: 20/10
  # Venus : 13/10
  # Tomato: 13/10
  # Katushka: 13/10
  thrsBFP <- 5
  thrstSapphire <- 10
  thrsVenus <- 10
  thrsTomato <- 10
  thrsKatushka <- 10
  # summary(nodes_df$C5_Katushka_mean_intensity)
  colnames(nodes_df) <- gsub('-','_',colnames(nodes_df))
  # C1_BFP_mean_coverage
  # summary(nodes_df$`C1-BFP_pct_coverage`)
  
  nodes_df <- nodes_df %>%
    mutate(
      C1_BFP_status = case_when(
        (C1_BFP_pct_coverage >= thres_vol_marker) & (C1_BFP_mean_intensity_nuc >= thrsBFP) ~ 1,
        TRUE ~ 0
      ),
      C2_tSapphire_status = case_when(
        (C2_tSapphire_pct_coverage >= thres_vol_marker) & (C2_tSapphire_mean_intensity_nuc >= thrstSapphire) ~ 1,
        TRUE ~ 0
      ),
      C3_venus_status = case_when(
        (C3_venus_pct_coverage >= thres_vol_marker) & (C3_venus_mean_intensity_nuc >= thrsVenus) ~ 1,
        TRUE ~ 0
      ),
      C4_tomato_status = case_when(
        (C4_tomato_pct_coverage >= thres_vol_marker) & (C4_tomato_mean_intensity_nuc >= thrsTomato) ~ 1,
        TRUE ~ 0
      ),
      C5_katushka_status = case_when(
        (C5_katushka_pct_coverage >= thres_vol_marker) & (C5_katushka_mean_intensity_nuc >= thrsKatushka) ~ 1,
        TRUE ~ 0
      )
    )
  # summary(as.factor(nodes_df$C1_BFP_status))  
  # summary(as.factor(nodes_df$C2_tSapphire_status))  
  # summary(as.factor(nodes_df$C3_venus_status))  
  # summary(as.factor(nodes_df$C4_tomato_status))  
  # summary(as.factor(nodes_df$C5_katushka_status))  
  nodes_df$cluster_label <- 'unlabeled'
  nodes_df$clone_desc <- 'unlabeled'
  
  nbcores <- 12
  labels <- parallel::mclapply(1:dim(nodes_df)[1], function(i) {
    clone_label <- get_clone_label(nodes_df[i,'C1_BFP_status'], 
                                   nodes_df[i,'C2_tSapphire_status'], 
                                   nodes_df[i,'C3_venus_status'], 
                                   nodes_df[i,'C4_tomato_status'], 
                                   nodes_df[i,'C5_katushka_status'])
    return(clone_label)
  }, mc.cores = nbcores)
  nodes_df$cluster_label <- unlist(labels)
  print(summary(as.factor(nodes_df$cluster_label)))
  dim(nodes_df)
  excluded_cols <- c('C1_BFP_status','C2_tSapphire_status','C3_venus_status',
                     'C4_tomato_status', 'C5_katushka_status')
  nodes_df <- nodes_df %>%
    dplyr::select(-all_of(excluded_cols))
  data.table::fwrite(nodes_df, paste0(save_dir,'filtered_cell_profiles.csv.gz'))
  
}
# annotate_prob=0.85: if a sub population with more cells than others in 0.85 quantile --> shown
# xmax: width of image
# ymax: height of image
viz_downsample_phenotype_v2 <- function(df, meta_clones, 
                                        xmax, ymax,col=NULL,
                                        annotate_prob=0.85 
){
  library(ggplot2)
  library(ggrepel)
  
  df <- df %>% inner_join(meta_clones, by=c('celltype'='clone_id'))
  if(is.null(col)){
    # cts <- gtools::mixedsort(unique(df$celltype))
    # meta_clones <- data.table::fread(paste0(script_dir,'predefined_clones_v2.csv')) %>% as.data.frame()
    meta_clones <- meta_clones %>%
      dplyr::filter(clone_id %in% unique(df$celltype))
    col <- meta_clones$cluster_color
    names(col) <- meta_clones$clone_id
    # nbcelltypes <- length(cts)
    # col <- colorRampPalette(brewer.pal(8, "Set2"))(nbcelltypes)
    # names(col) <- cts
  }
  df <- df %>%
    dplyr::filter(node_size>1)
  dim(df)
  df$node_size
  # df$nb_cells <- round(log2(sample(10:1000, dim(df)[1])),2)
  df_annotated <- df
  summary(df_annotated$node_size)
  df_annotated <- df_annotated %>% 
    dplyr::filter(node_size>quantile(df_annotated$node_size,probs=annotate_prob))
  dim(df_annotated)
  df$ct <- gsub('Clone_','',df$celltype)
  df$y <- ymax - df$y
  p <- ggplot(df, aes(x = x, y = y)) + 
    geom_point(aes(color = celltype, size=log2(node_size))) + 
    # geom_text(aes(label=celltype_desc))+
    annotate('text', x = df$x, y = df$y,
             label = df$ct, size=2)+
    scale_color_manual(values = col) +
    theme_bw(base_size = 12) + 
    theme(legend.position = "bottom",
          legend.box = "vertical",
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text  = element_blank(),
          axis.title = element_blank()) + 
    geom_text_repel( data = df_annotated, aes(label = clone_desc), max.overlaps = Inf,
                     size = 3.5, box.padding = unit(0.35, "lines"), 
                     point.padding = unit(0.3, "lines"),
                     min.segment.length = 0) + 
    guides(color = FALSE,size=guide_legend(title="log2(cell counts)"))
  # p
  png(paste0(save_dir,"summary_population.png"), height = ymax+100, width=xmax,res = 2*72)
  print(p)
  dev.off()
  
  return(p)
  
}
get_edges_nodes_from_igraph <- function(g){
  library(igraph)
  if(is.null(V(g)$name)){
    g <- g %>%
      set_vertex_attr("name", value = paste0('V',rep(1:length(V(g)),1)))
  }
  edges <- igraph::as_data_frame(g, what="edges")
  # edges = get.edgelist(g) %>% as.data.frame()
  # colnames(edges) <- c('from','to')
  print(head(edges))
  nodes_df <- igraph::as_data_frame(g, what = c("vertices"))
  print(head(nodes_df))
  print(colnames(nodes_df))
  print(dim(nodes_df))
  return(list(nodes=nodes_df, edges=edges))
}


# script_dir <- '/Users/hoatran/Documents/jean_project/data/script/'
# meta_clones <- get_reference_clones(script_dir)
# meta_clones$cluster_color <- get_color_clone(meta_clones$clone_desc)
# data.table::fwrite(meta_clones, paste0(script_dir,'predefined_clones_v2.csv'))
# meta_clones <- data.table::fread(paste0(script_dir,'predefined_clones_v2.csv')) %>% as.data.frame()

# get_reference_clones <- function(save_dir){
#   markers <- c('eBFP2','tSapphire','Venus','tdTomato','Katushka') # Protein channel
#   vals <- c(0,1)
#   ls_clones <- list()
#   for(e in vals){
#     for(ts in vals){
#       for(v in vals){
#         for(td in vals){
#           for(k in vals){
#             lb <- ''
#             s <- 0
#             if(e==1){
#               lb <- paste0(lb, markers[1])
#               s <- s + 2^0
#             }
#             if(ts==1){
#               lb <- paste0(lb,', ',markers[2])
#               s <- s + 2^1
#             }
#             if(v==1){
#               lb <- paste0(lb,', ',markers[3])
#               s <- s + 2^2
#             }
#             if(td==1){
#               lb <- paste0(lb,', ',markers[4])
#               s <- s + 2^3
#             }
#             if(k==1){
#               lb <- paste0(lb,', ',markers[5])
#               s <- s + 2^4
#             }
#             # if(grepl('^, ', lb)){
#             #   lb <- str_sub(lb, 3, str_length(lb))
#             # }
#             ls_clones[[paste0('Clone_',s)]] <- lb
#           } 
#         } 
#       }  
#     }
#   }
#   clone_df <- data.frame(clone_desc=unlist(ls_clones), clone_id=names(ls_clones))
#   clone_df <- clone_df %>%
#     dplyr::filter(clone_id!='Clone_0')
#   clone_df$clone_desc <- gsub('^, ','',clone_df$clone_desc)
#   # View(clone_df)
#   clone_df <- clone_df[gtools::mixedsort(clone_df$clone_id),]
#   data.table::fwrite(clone_df, paste0(save_dir,'predefined_clones.csv'))
#   return(clone_df)
# }  

get_color_clone <- function(clones){
  col <- c("eBFP2"="#fbeb14",                                      
           "tSapphire"="#e32636",                                  
           "eBFP2, tSapphire"="#7fbf7f",                           
           "Venus"="#008040",                                      
           "eBFP2, Venus"="#c06000",                               
           "tSapphire, Venus"="#cc99a2",                           
           "eBFP2, tSapphire, Venus"="#ffc0cb",                    
           "tdTomato"="#408000",                                   
           "eBFP2, tdTomato"="#b266b2",                            
           "tSapphire, tdTomato" ="#ffc04c",                       
           "eBFP2, tSapphire, tdTomato"="#ff6666",                 
           "Venus, tdTomato"="#198c19",                            
           "eBFP2, Venus, tdTomato"="#7f7fff",                     
           "tSapphire, Venus, tdTomato"="#ffe4b2",                 
           "eBFP2, tSapphire, Venus, tdTomato"="#540a87",          
           "Katushka"="#e31f50",                                   
           "eBFP2, Katushka"="#af9210",                            
           "tSapphire, Katushka"="#7a660b",                        
           "eBFP2, tSapphire, Katushka"="#463a06",                 
           "Venus, Katushka"="#d7c887",                            
           "eBFP2, Venus, Katushka"="#8fe73e",                     
           "tSapphire, Venus, Katushka"="#57617",                 
           "eBFP2, tSapphire, Venus, Katushka"="#db9eeb",          
           "tdTomato, Katushka"="#400f4d",                         
           "eBFP2, tdTomato, Katushka"="#73c8da",                  
           "tSapphire, tdTomato, Katushka"="#9e1209",              
           "eBFP2, tSapphire, tdTomato, Katushka"="#84ed7c",       
           "Venus, tdTomato, Katushka"="#eb9710",                  
           "eBFP2, Venus, tdTomato, Katushka"="#9e6317",           
           "tSapphire, Venus, tdTomato, Katushka"="#07a70b",       
           "eBFP2, tSapphire, Venus, tdTomato, Katushka"="#1109d8")
  return(col[clones])
  
}

# sce <- readRDS('~/Downloads/SA535_clonealign/SA535X10XB03693.rds')
# dim(sce)
# sce$ml_params



get_reference_clones <- function(meta_clones_fn=''){
  # meta_clones_fn <- paste0(save_dir,'predefined_clones_v2.csv')
  
  if(file.exists(meta_clones_fn)){  #predefined clones labels
    clone_df <- data.table::fread(meta_clones_fn) %>% as.data.frame()
    print(dim(clone_df))
    return(clone_df)
  }else{# if do not exist predefined clones labels, define it here
    markers <- c('eBFP2','tSapphire','Venus','tdTomato','Katushka') # Protein channel
    vals <- c(0,1)
    ls_clones <- list()
    for(e in vals){
      for(ts in vals){
        for(v in vals){
          for(td in vals){
            for(k in vals){
              lb <- ''
              s <- 0
              if(e==1){
                lb <- paste0(lb, markers[1])
                s <- s + 2^0
              }
              if(ts==1){
                lb <- paste0(lb,', ',markers[2])
                s <- s + 2^1
              }
              if(v==1){
                lb <- paste0(lb,', ',markers[3])
                s <- s + 2^2
              }
              if(td==1){
                lb <- paste0(lb,', ',markers[4])
                s <- s + 2^3
              }
              if(k==1){
                lb <- paste0(lb,', ',markers[5])
                s <- s + 2^4
              }
              # if(grepl('^, ', lb)){
              #   lb <- str_sub(lb, 3, str_length(lb))
              # }
              ls_clones[[paste0('Clone_',s)]] <- lb
            } 
          } 
        }  
      }
    }
    clone_df <- data.frame(clone_desc=unlist(ls_clones), clone_id=names(ls_clones))
    clone_df <- clone_df %>%
      dplyr::filter(clone_id!='Clone_0')
    clone_df$clone_desc <- gsub('^, ','',clone_df$clone_desc)
    # View(clone_df)
    clone_df <- clone_df[gtools::mixedsort(clone_df$clone_id),]
    # data.table::fwrite(clone_df, meta_clones_fn)
    # dim(meta_clones)
    # View(head(meta_clones))
    # data.table::fwrite(meta_clones, paste0(script_dir,'predefined_clones.csv'), quote=F)
    # meta_clones <- data.table::fread(paste0(script_dir,'predefined_clones.csv')) %>% as.data.frame()
    # dim(meta_clones)
    
    # col <- colorRampPalette(brewer.pal(8, "Set2"))(dim(meta_clones)[1])
    # names(col) <- meta_clones$clone_id
    # meta_clones$cluster_color <- col  # just temporary colors, will change it later
    clone_df$cluster_color <- get_color_clone(clone_df$clone_desc)
    data.table::fwrite(clone_df, meta_clones_fn)
    
    return(clone_df)
  }
} 


get_clone_label_v1 <- function(e, ts, v, td, k){
  markers <- c('eBFP2','tSapphire','Venus','tdTomato','Katushka')
  # vals <- c(0,1)
  # ls_clones <- list()  
  lb <- ''
  s <- 0
  if(e==1){
    lb <- paste0(lb, markers[1])
    s <- s + 2^0
  }
  if(ts==1){
    lb <- paste0(lb,', ',markers[2])
    s <- s + 2^1
  }
  if(v==1){
    lb <- paste0(lb,', ',markers[3])
    s <- s + 2^2
  }
  if(td==1){
    lb <- paste0(lb,', ',markers[4])
    s <- s + 2^3
  }
  if(k==1){
    lb <- paste0(lb,', ',markers[5])
    s <- s + 2^4
  }
  if(grepl('^, ', lb)){
    lb <- str_sub(lb, 3, str_length(lb))
  }
  if(s==0 & lb==''){
    lb <- 'unlabeled'
  }
  res <- list(clone_label=paste0('Clone_',s),clone_desc=lb)
  return(res)
}

viz_facs_imaging_summary_results <- function(){
  script_dir <- '/Users/hoatran/Documents/jean_project/data/script/'
  metaclone_df <- data.table::fread(paste0(script_dir,'predefined_clones_v2.csv')) %>% as.data.frame()
  head(metaclone_df)
  metaclone_df <- metaclone_df %>%
    dplyr::select(-cluster_color)%>%
    dplyr::rename(cluster_label=clone_id)
  df <- df %>% inner_join(metaclone_df, by=c("cluster_label"))
  data.table::fwrite(df, paste0(save_dir,'total_nodes_celltype_v3.csv'))
  
  df1 <- df
  df1 <- df %>%
    dplyr::filter(!grepl('tdTomato',clone_desc))
  cr <- cor(df1$FACS_counts, df1$imaging_counts, method = 'pearson')
  cr
  df$cluster_label
  dim(df)
  sum_images <- sum(df$imaging_counts)
  metaclone_df <- metaclone_df %>%
    dplyr::filter(clone_id %in% unique(df$cluster_label))
  col <- metaclone_df$cluster_color
  names(col) <- paste0(gsub('Clone_','',metaclone_df$clone_id),':',metaclone_df$clone_desc)
  
  df$clone <- gsub('Clone_','',df$cluster_label)
  
  # df <- df %>%
  #   dplyr::mutate(pct_imaging=imaging_counts/sum_images)
  df$clone_desc <- paste0(df$clone,':',df$clone_desc)
  p <- ggplot(df,aes(x = FACS_counts, y = imaging_counts, color=clone_desc)) + 
    geom_point(size=5, alpha=0.6)+
    scale_color_manual(values = col) +
    annotate('text', x = df$FACS_counts, y = df$imaging_counts,
             label = df$clone, size=2.5)+
    theme_bw(base_size = 8) + 
    theme(legend.position = "bottom",
          legend.text = element_text(size=5),
          legend.box = "vertical",
          panel.grid = element_blank()) +
    guides(color = guide_legend(title="",
                                override.aes = list(size=1), ncol = 3))
  p
  png(paste0(save_dir,'vitro_celltypes_correlation_v1.png'), 
      height = 1100, width=750,res = 2*72)
  print(p)
  dev.off()
  
  p <- ggplot(df,aes(x = FACS_Percent, y = pct_imaging, color=cluster_label)) + 
    geom_point()
  
  
  data.table::fwrite(df, paste0(save_dir,'summary_facs_imaging_thrs_vol_200.csv'))
  data.table::fwrite(df, paste0(save_dir,'summary_facs_imaging_thrs_vol_50.csv'))
  
}