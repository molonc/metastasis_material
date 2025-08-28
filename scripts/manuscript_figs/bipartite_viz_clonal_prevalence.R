
# Figure 6
suppressPackageStartupMessages({
  library("dplyr")
  library("ggplot2")
  library("ggbump")
  library("cowplot")
})

## Loading utility functions
script_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/scripts/'
source(paste0(script_dir, 'corrupt_tree/src/stat/prevalence_utils.R'))


viz_bipartite_graph_SA535 <- function(){
  input_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/dlp_trees/'
  datatag <- 'SA535'
  save_fig_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/figures/fig_cnv/'
  cellclones_fn <- paste0(input_dir, datatag, '/', 'cell_clones.csv.gz')
  cell_clones <- data.table::fread(cellclones_fn)
  excluded_clones <- c('Un','unassigned','None')
  cell_clones <- cell_clones[!cell_clones$clone_id %in% excluded_clones,]
  print(dim(cell_clones))
  print(summary(as.factor(cell_clones$clone_id)))
  
  library_grouping_fn <- paste0(input_dir, datatag, '/' ,'library_groupings.csv.gz')
  meta_data <- get_meta_data(cell_clones, results_dir, library_grouping_fn)
  print(dim(meta_data))
  colnames(meta_data)
  View(head(meta_data))
  
  prevalence_df <- meta_data %>%
    dplyr::filter(origin!='Tumor_Recur')%>%
    dplyr::group_by(clone_id, pdxid, mainsite) %>%
    dplyr::summarise(nb_cells=n())
  # View(prevalence_df)
  
  color_codes <- get_color_clone_v2(unique(cell_clones$clone_id), datatag)
  # clone, mainsite, clonal_prevalence
  plt_ls <- list()
  for(pdx in unique(prevalence_df$pdxid)){
    print(pdx)
    df_pri <- prevalence_df %>%
      dplyr::filter(pdxid==pdx & mainsite=='Primary') # & nb_cells>=5
    total_cells <- sum(df_pri$nb_cells)
    print(total_cells)
    df_pri <- df_pri %>%
      dplyr::mutate(clonal_prevalence=round(100*nb_cells/total_cells,1))
    
    df_met <- prevalence_df %>%
      dplyr::filter(pdxid==pdx & mainsite=='Metastasis') # & nb_cells>=5
    total_cells <- sum(df_met$nb_cells)
    print(total_cells)
    df_met <- df_met %>%
      dplyr::mutate(clonal_prevalence=round(100*nb_cells/total_cells,1))
    df <- dplyr::bind_rows(df_pri, df_met)
    cols_use <- color_codes[unique(df$clone_id)]
    p <- viz_bipartite_graph(df, save_fig_dir, datatag=paste0('M',stringr::str_sub(pdx, nchar(pdx), nchar(pdx))), cols_use)
    # p
    plt_ls[[pdx]] <- p
  }
  plt_ls1 <- list(plt_ls[["X0011_2361"]],plt_ls[["X0011_2362"]],plt_ls[["X0011_2363"]],plt_ls[["X0011_2364"]])
  p_total <- cowplot::plot_grid(plotlist = plt_ls1, nrow=1) +
    theme(plot.background = element_rect(fill = "white", colour = "white"))
  ggsave(  
    filename = paste0(save_fig_dir, datatag,"_total_graph.svg"),  
    plot = p_total,  
    height = 7, width = 16, dpi = 150)
  
  
}


viz_bipartite_graph <- function(df, save_fig_dir, datatag='',cols_use=NULL){
  # df <- tibble(clone = c("A", "B", "C", "A", "B", "C"),
  #              mainsite = c("Primary", "Primary", "Primary", "Metastasis", "Metastasis", "Metastasis"),
  #              clonal_prevalence = c(0.6, 0.3, 0.1, 0.2, 0.25, 0.55))
  # if(is.null(cols_use) & datatag=='SA919'){
  #   
  # }else if(is.null(cols_use) & datatag=='SA535'){
  #   
  # }else{
  #   print('')
  # }
  
  ## just for order in x axis
  df <- df %>%
    dplyr::mutate(xval=
                    case_when(mainsite=='Primary' ~ 1,
                              TRUE ~ 2)
    )
  df_pri <- df %>% 
    dplyr::filter(mainsite=="Primary") %>% 
    dplyr::arrange(-clonal_prevalence)
  df_pri$rank <- rep(1:dim(df_pri)[1],1)
  
  df_met <- df %>% 
    dplyr::filter(mainsite=="Metastasis") %>% 
    dplyr::arrange(-clonal_prevalence)
  df_met$rank <- rep(1:dim(df_met)[1],1)
  
  df <- dplyr::bind_rows(df_pri, df_met)
  df$desc <- paste0(df$clone_id,' (',df$clonal_prevalence,')')
  
  p <- ggplot(df, aes(xval, rank, color = clone_id)) +
    geom_point(size = 7) +
    geom_text(data = df %>% filter(xval == min(xval)),
              aes(x = xval - .1, label = desc), size = 5, hjust = 1) +
    geom_text(data = df %>% filter(xval == max(xval)),
              aes(x = xval + .1, label = desc), size = 5, hjust = 0) +
    ggbump::geom_bump(size = 2, smooth = 8) +
    scale_x_continuous(limits = c(min(df$xval)-0.5, max(df$xval)+0.5),
                       breaks = c(min(df$xval), max(df$xval)),
                       labels = c('Primary',"Metastasis")) + #seq(2011, 2013, 1)
    cowplot::theme_minimal_grid(font_size = 14, line_size = 0) +
    theme(legend.position = "none",
          axis.text.x = element_text(size=14),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank()) +
    labs(y=NULL,#y = "Clonal Prevalence",
         x = NULL,
         title = datatag) +
    scale_y_reverse() +
    scale_color_manual(values = cols_use)
  
  # png(paste0(save_fig_dir,datatag,"_bipartite_graph.png"), 
  #     height = 550, width=2*450, res = 2*50)
  # print(p)
  # dev.off()
  return(p)
  
}

get_color_clone_v2 <- function(clone_levels, datatag){
  clone_levels <- gtools::mixedsort(clone_levels[!grepl("None", clone_levels)])
  
  ## Meta colors based on package inlmisc::GetColors(n) -- problem at installation, so keep values of color codes in csv file
  # colorcode_fn <- "/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_wholedata/config/colorcode.csv"
  # colorcode_fn <- "/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_wholedata/config/colorcode_v2.csv"
  # colorcode_fn <- paste0('/home/htran/Projects/hakwoo_project/metastasis_material/materials/dlp_trees/colorCode_clones/color_code_',datatag,'.csv')
  input_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/dlp_trees/'
  datatag <- 'SA535'
  colorcode_fn <- paste0(input_dir, 'colorCode_clones/color_code_',datatag, '.csv')
  if(file.exists(colorcode_fn)){
    print("Loading color codes from file for plotting...")    
    color_df <- data.table::fread(colorcode_fn)
    # color_df <- color_df %>%
    #   dplyr::filter(series==datatag)
    clone_pal <- color_df$color
    names(clone_pal) <- color_df$clone_id
    if(sum(clone_levels %in% color_df$clone)==length(clone_levels)){
      clone_pal <- clone_pal[clone_levels]
    }else{
      print("Color code file do not contain full list of clones and its colors")  
      print("Generating the color codes for plotting...")    
      clone_pal <- make_clone_palette(clone_levels)
    }
    
  }else{
    print("Generating the color codes for plotting...")    
    clone_pal <- make_clone_palette(clone_levels)
  }
  
  for(i in names(clone_pal)){
    if(clone_pal[i]=="#E7EBFA"){
      print(i)
      clone_pal[i] <- "#C0C0C0"
    }
  }
  
  print(clone_pal)
  return(clone_pal)
}  


