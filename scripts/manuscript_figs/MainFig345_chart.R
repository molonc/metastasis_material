


# install.packages("ggplot2")
# install.packages("ggrepel")
# install.packages("tidyverse")
library(ggplot2)
library(dplyr)
library(ggrepel)
library(tidyverse)

viz_Fig3_clone_prevalence_samples <- function(){
  
  input_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/dlp_trees/clonal_prevalence_per_sample/'
  datatag <- 'SA535'
  save_fig_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/figures/barchart_Fig345/'
  dir.create(save_fig_dir)
  colors_df <- data.table::fread(paste0(input_dir, 'color_code_',datatag,'.csv'))
  colors_code <- colors_df$color
  names(colors_code) <- colors_df$clone_id
  
  prevalence_df <- data.table::fread(paste0(input_dir,datatag, '_cells_proportions.csv'))
  View(head(prevalence_df))
  unique(prevalence_df$sample_id)
  ## Primary or Metastasis
  mts <- sapply(strsplit(prevalence_df$sample_id, '_'), function(x){
    return(x[3])
  })
  mouse_ids <- sapply(strsplit(prevalence_df$sample_id, '_'), function(x){
    # return(gsub('mouse','',x[2]))
    return(gsub('mouse','',paste0(x[1],'_', x[2])))
  })
  mouse_ids1 <- sapply(strsplit(prevalence_df$sample_id, '_'), function(x){
    # return(gsub('mouse','',x[2]))
    return(gsub('mouse','',x[2]))
  })
  prevalence_df$mainsite <- as.character(mts)
  prevalence_df$mouse_id <- as.character(mouse_ids)
  prevalence_df$mouse_id <- as.character(mouse_ids1)
  unique(prevalence_df$mainsite)
  unique(prevalence_df$mouse_id)
  
  # t <- prevalence_df %>%
  #   dplyr::group_by(mouse_id, clone_id, mainsite) %>%
  #   dplyr::summarise(nb_cells_per_clone=sum(nb_cells)) %>%
  #   dplyr::filter(mouse_id=='M2')
  # 
  # st <- t %>%
  #   dplyr::group_by(mainsite) %>%
  #   dplyr::summarise(nb_cells_total=sum(nb_cells_per_clone))
  # t <- t %>%
  #   dplyr::left_join(st, by=c('mainsite'))
  # t$pct_cells <- round(100*t$nb_cells_per_clone/t$nb_cells_total, 2)
  # View(t)
  # Primary samples
  prev_primary_df <- prevalence_df %>%
    dplyr::filter(mainsite=='Primary')
  primary_plt_list <- list()
  for(s in unique(prev_primary_df$sample_id)){
    df <- prev_primary_df %>%
      dplyr::filter(sample_id==s)
    print(sum(df$pct_cells))
    mouseId <- unique(df$mouse_id)
    print(mouseId)
    cols_use <- colors_code[unique(df$clone_id)]
    print(cols_use)
    p <- clone_prevalence_piechart(df, mouseId, save_fig_dir, cols_use)
    primary_plt_list[[mouseId]] <- p
  }
  primary_plt_list[['M1']]
  primary_plt_list[['M3']]
  primary_plt_list[['M4']]
  p_piechart_lg <- clone_prevalence_piechart_legend(save_fig_dir, cols_use=NULL)
  
  # Met samples
  prev_met_df <- prevalence_df %>%
    dplyr::filter(mainsite=='Metastasis')
  met_plt_list <- list()
  for(s in unique(prev_met_df$sample_id)){
    df <- prev_met_df %>%
      dplyr::filter(sample_id==s)
    print(sum(df$pct_cells))
    mouseId <- unique(df$mouse_id)
    print(mouseId)
    cols_use <- colors_code[unique(df$clone_id)]
    print(cols_use)
    p <- clone_prevalence_donutchart(df, paste0(mouseId), save_fig_dir, cols_use)
    met_plt_list[[mouseId]] <- p
  }
  names(met_plt_list)
  
  met_plt_list[["SA535X4XB05667_M3"]]
  p_total_pri <- cowplot::plot_grid(primary_plt_list[['M1']], primary_plt_list[['M2']], primary_plt_list[['M3']], primary_plt_list[['M4']], nrow=1)
  
  p_total_met <- cowplot::plot_grid(met_plt_list[['SA535X4XB05723_M1']], 
                                    met_plt_list[['SA535X4XB09113_M2']],met_plt_list[['SA535X4XB09109_M2']],met_plt_list[['SA535X4XB09106_M2']],
                                    met_plt_list[['SA535X4XB05667_M3']],
                                    met_plt_list[['SA535X4XB05989_M4']], nrow=1)
  p_total <- cowplot::plot_grid(p_total_pri, p_total_met, ncol=1, rel_heights = c(1,0.3))
  
}

clone_prevalence_piechart <- function(df, datatag, save_dir, cols_use=NULL){
  
  
  # df <- data.frame(pct_cells = c(15, 25, 32, 28),
  #                  clone_id = paste0("G", 1:4))
  # cols_use <- c('red','blue','cyan','pink')
  # names(cols_use) <- df$clone_id
  # Get the positions
  df2 <- df %>% 
    mutate(csum = rev(cumsum(rev(pct_cells))), 
           pos = pct_cells/2 + lead(csum, 1),
           pos = if_else(is.na(pos), pct_cells/2, pos))
  
  p <- ggplot(df, aes(x = "" , y = pct_cells, fill = fct_inorder(clone_id))) +
    geom_col(width = 1, color = 1) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = cols_use) + 
    # scale_fill_brewer(palette = "Pastel1") +
    geom_text_repel(data = df2,
                    aes(y = pos, label = paste0(clone_id,' ',pct_cells, "%")),
                    size = 4.5, nudge_x = 1, show.legend = FALSE, segment.color = NA) +
    # guides(fill = guide_legend(title = "clone_id")) +
    theme_void() + 
    theme(legend.position = 'none') + 
    labs(title = datatag)
  
  
  saveRDS(p, paste0(save_dir, datatag, '_piechart.rds'))
  return(p)
}

clone_prevalence_piechart_legend <- function(save_dir, cols_use=NULL){
  
  df <- data.frame(pct_cells = c(50,50),
                   clone_id = paste0("G", 1:2))
  # df <- data.frame(pct_cells = c(15, 25, 32, 28),
  #                  clone_id = paste0("G", 1:4))
  # cols_use <- c('red','blue','cyan','pink')
  # names(cols_use) <- df$clone_id
  # Get the positions
  df2 <- df %>% 
    mutate(csum = rev(cumsum(rev(pct_cells))), 
           pos = pct_cells/2 + lead(csum, 1),
           pos = if_else(is.na(pos), pct_cells/2, pos))
  
  
  cols_use <- rep('white', length(unique(df$clone_id)))
  names(cols_use) <- df$clone_id
  p_legend <- ggplot(df, aes(x = "" , y = pct_cells, fill= fct_inorder(clone_id))) + #, fill = fct_inorder(clone_id))
    geom_col(width = 1, color = 1) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = cols_use)+
    theme_void() + 
    theme(legend.position = 'none')
  # p_legend
  saveRDS(p_legend, paste0(save_dir, 'piechart_legend_shape.rds'))  
  
  return(p_legend)
}
clone_prevalence_donutchart_legend <- function(save_dir, cols_use=NULL){
  
  df <- data.frame(pct_cells = c(50,50),
                   clone_id = paste0("G", 1:2))
  # df <- data.frame(pct_cells = c(15, 25, 32, 28),
  #                  clone_id = paste0("G", 1:4))
  # cols_use <- c('red','blue','cyan','pink')
  # names(cols_use) <- df$clone_id
  # Get the positions
  df2 <- df %>% 
    mutate(csum = rev(cumsum(rev(pct_cells))), 
           pos = pct_cells/2 + lead(csum, 1),
           pos = if_else(is.na(pos), pct_cells/2, pos))
  
  
  cols_use <- rep('white', length(unique(df$clone_id)))
  names(cols_use) <- df$clone_id
  p_legend <- ggplot(df, aes(x = "" , y = pct_cells, fill= fct_inorder(clone_id))) + #, fill = fct_inorder(clone_id))
    geom_col(width = 1, color = 1) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = cols_use)+
    theme_void() + 
    theme(legend.position = 'none')
  # p_legend
  saveRDS(p_legend, paste0(save_dir, 'piechart_legend_shape.rds'))  
  
  return(p_legend)
}

clone_prevalence_donutchart <- function(df, datatag, save_dir, cols_use){
  # df <- data.frame(value = c(10, 30, 32, 28),
  #                  clone_id = paste0("G", 1:4))
  
  # Hole size
  hsize <- 1.5
  
  df <- df %>% 
    mutate(x = hsize)
  df2 <- df %>% 
    mutate(csum = rev(cumsum(rev(pct_cells))), 
           pos = pct_cells/2 + lead(csum, 1),
           pos = if_else(is.na(pos), pct_cells/2, pos))
  
  p <- ggplot(df, aes(x = hsize, y = pct_cells, fill = clone_id)) +
    geom_col(color = "black") +
    # geom_text(aes(label = pct_cells),
    #           position = position_stack(vjust = 0.5)) +
    coord_polar(theta = "y") +
    # scale_fill_brewer(palette = "GnBu") +
    scale_fill_manual(values = cols_use) + 
    xlim(c(0.2, hsize + 0.5)) +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.position = 'none')
  # p
  
  saveRDS(p, paste0(save_dir, datatag, '_donutchart.rds'))
  return(p)
}

clone_prevalence_piechart <- function(df, datatag, save_dir, cols_use){
  
  
  df <- data.frame(pct_cells = c(15, 25, 32, 28),
                   clone_id = paste0("G", 1:4))
  cols_use <- c('red','blue','cyan','pink')
  names(cols_use) <- df$clone_id
  # Get the positions
  df2 <- df %>% 
    mutate(csum = rev(cumsum(rev(pct_cells))), 
           pos = pct_cells/2 + lead(csum, 1),
           pos = if_else(is.na(pos), pct_cells/2, pos))
  
  p <- ggplot(df, aes(x = "" , y = pct_cells, fill = fct_inorder(clone_id))) +
    geom_col(width = 1, color = 1) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = cols_use) + 
    # scale_fill_brewer(palette = "Pastel1") +
    geom_text_repel(data = df2,
                    aes(y = pos, label = paste0(clone_id,' ',pct_cells, "%")),
                    size = 4.5, nudge_x = 1, show.legend = FALSE, segment.color = NA) +
    # guides(fill = guide_legend(title = "clone_id")) +
    theme_void() + 
    theme(legend.position = 'none')
  p
  cols_use <- rep('white', 4)
  names(cols_use) <- df$clone_id
  p_legend <- ggplot(df, aes(x = "" , y = pct_cells, fill= fct_inorder(clone_id))) + #, fill = fct_inorder(clone_id))
    geom_col(width = 1, color = 1) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = cols_use)+
    theme_void() + 
    theme(legend.position = 'none')
  # p_legend
  saveRDS(p, paste0(save_dir, datatag, '_piechart_legend_shape.rds'))
  saveRDS(p, paste0(save_dir, datatag, '_piechart.rds'))
  
}


donutchart_with_label <- function(){
  mdat <- data.frame(
    category = c("C", "E", "I", "L", "Mi", "Mo", 
                 "O", "Q", "S", "V"), 
    ct = c(147, 275, 431, 967, 121, 105, 17, 186, 620, 42))
  
  mdat$category <- factor(mdat$category, levels = mdat$category)
  
  # Fractions
  mdat$fraction <- mdat$ct / sum(mdat$ct)
  
  # Cumulative fractions; this forms the top of each rectangle
  mdat$ymax <- cumsum(mdat$fraction)
  
  # This will the the bottom of the rectangle
  mdat$ymin <- c(0, head(mdat$ymax, n = -1))
  
  # Label position - this isn't right
  mdat$labelPosition <- ((mdat$ymax + mdat$ymin) / 2) 
  
  # Labels
  mdat$label <- paste0(mdat$category, " Fraction: \n", 
                       round(mdat$ct/sum(mdat$ct), 4) * 100, "%")
  
  # The donut chart is just a stacked bar chart in polar coordinates. 
  # So if you choose your x within geom_label() below the lower boundary within xlim(), 
  # the labels wander off to the opposite side of the plot.
  g <- ggplot(mdat, aes(ymax = ymax, ymin = ymin, xmax = 11, xmin = 10, 
                        fill = category))
  g <- g + geom_rect()
  g <- g + geom_label(x = 12.1, aes(y = labelPosition, label = label), size = 3)
  g <- g + scale_fill_brewer(palette = "Set3")
  g <- g + scale_color_brewer(palette = "Set3")
  g <- g + coord_polar(theta = "y")
  g <- g + xlim(c(7, 12))
  g <- g + theme_void()
  g <- g + theme(legend.position = "none")
  g
}
