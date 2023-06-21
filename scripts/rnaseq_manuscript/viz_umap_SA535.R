suppressPackageStartupMessages({
  require("SingleCellExperiment")
  require("stringr")
  require("tidyverse")
  require("Seurat")
  # require("sctransform")
  require("dplyr")
  require("inlmisc")
})

# library(extrafont)
# font_import(prompt=F, paths ='/usr/share/fonts/truetype/myfonts/') # import Helvetica font
# fonts()
my_font <- "Helvetica"
# Just modify function for predefined clone colors. 
main(){
  meta_data <- data.table::fread('/home/htran/Projects/hakwoo_project/metastasis_material/materials/10x/SA535_10x_metadata_passage_X4_full_infos.csv')
  
  # meta$nb_cells_introns
  # meta <- meta %>% 
  #   dplyr::filter(mouse_id=='SA535X4XB05462')
  meta_data <- meta_data %>% 
    dplyr::filter(pdxid!='M2364')
  colnames(meta_data)
  dim(meta_data)
  # View(meta_data)
  base_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/'
  output_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/normalized_introns/'
  datatag <- 'SA535' 
  
  
}
get_clone_labels <- function(output_dir, meta_data){
  # Get clone labels first
  
  clonealign_dir <- paste0(base_dir,'clonealign_introns/TreeAlign_result_08/')

  fns <- list.files(clonealign_dir)    
  fns <- fns[grepl('*_clone_assign.csv.gz',fns)]
  fns <- gsub('_clone_assign.csv.gz','',fns)
  fns <- fns[fns %in% meta_data$mouse_id]
  fns <- fns[fns != 'SA535X4XB05989']
  clonealign_stat <- tibble::tibble()
  for(f in fns){
    # First clonealign output
    sdf <- data.table::fread(paste0(clonealign_dir,f,'_clone_assign.csv.gz')) %>% as.data.frame()
    print(dim(sdf))
    if(dim(sdf)[1]>0){
      sdf$sample_id <- f
      clonealign_stat <- dplyr::bind_rows(clonealign_stat, sdf)  
    }
    
  }  
  dim(clonealign_stat)
  head(clonealign_stat)
    
  cells_df <- data.table::fread(paste0(output_dir,'meta_cells.csv.gz'))
  dim(cells_df)
  colnames(cells_df)
  # sids[!sids %in% unique(clonealign_stat$sample_id)]
  # sum(clonealign_stat$cell_id %in% cells_df$Barcode)
  clonealign_stat <- clonealign_stat %>%
    dplyr::inner_join(cells_df, by = c('sample_id'='mouse_id','cell_id'='Barcode'))
  dim(clonealign_stat)
  
  summary(as.factor(clonealign_stat$clone_id))
  sum(clonealign_stat$clone_id=='')
  clonealign_stat <- clonealign_stat %>%
    dplyr::mutate(clone_id=
                    case_when(
                      clone_id=='' ~ 'unassigned',
                      TRUE ~ clone_id
                    ))
  
  
  colnames(meta_data)
  colnames(clonealign_stat)
  meta_data1 <- meta_data %>%
    dplyr::select(mouse_id, pdxid, Site_origin, Grouping,
                  library_id)
  meta_data1$library_id[1]
  clonealign_stat$library_id[1]
  clonealign_stat <- clonealign_stat %>%
    dplyr::left_join(meta_data1, by = c('sample_id'='mouse_id',
                                        'library_id'))
  data.table::fwrite(clonealign_stat, paste0(output_dir,'total_clonal_assignment_08.csv.gz'))
  
  umap_df <- data.table::fread(paste0(output_dir, datatag, "_scTransform_umap_filtered_cells.csv.gz"))
  dim(umap_df)
  sum(clonealign_stat$cell_id %in% umap_df$Barcode)
  dim(clonealign_stat)
  clonealign_stat <- clonealign_stat %>%
    dplyr::select(cell_id, library_id, clone_id) #%>%
    # dplyr::filter()
  dim(umap_df)
  umap_df <- umap_df %>% 
    left_join(clonealign_stat, by=c('Barcode'='cell_id',
                                    'library_id'))
  
  umap_df <- umap_df %>% 
    left_join(meta_data1, by=c('library_id'))
  colnames(meta_data)
  meta_data1$library_id
  data.table::fwrite(umap_df, paste0(output_dir, datatag, "_scTransform_umap_filtered_cells_clones.csv.gz"))
  umap_df <- data.table::fread(paste0(output_dir, datatag, "_scTransform_umap_filtered_cells_clones.csv.gz"))
  
  # clonealign_stat$library_id <- gsub('.cache/','',clonealign_stat$Sample[1])
  # clonealign_stat$library_id <- gsub('/filtered_feature_bc_matrix','',clonealign_stat$library_id)
  # clonealign_stat$Sample <- NULL
  # clonealign_stat$datatag <- stringr::str_sub(clonealign_stat$id,1,6)
  # clonealign_stat$datatag <- gsub('X$','',clonealign_stat$datatag)
  # unique(clonealign_stat$datatag)
  # clonealign_stat$lcell_id <- paste0(clonealign_stat$library_id,'_',clonealign_stat$Barcode)
  # clonealign_stat$unique_clone <- get_unique_clone_id(clonealign_stat$clone)
  # output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/umap_figs/'
  # data.table::fwrite(clonealign_stat, paste0(output_dir,'clonealign_labels.csv'))
  # View(head(clonealign_stat))
  cls <- unique(umap_df$clone_id)
  unassign_clones <- c('unassigned','Unassigned','Un','un','None')
  cls <- cls[!cls %in% unassign_clones]
  print(cls)
  
  obs_clones <- sort(unique(umap_df$clone_id))
  unassign_clones <- c('unassigned','Unassigned','Un','un','None')
  obs_clones <- obs_clones[!obs_clones %in% unassign_clones]
  print(obs_clones)
  cols_use <- make_clone_palette(obs_clones)
  cols_use <- cols_use[obs_clones]
  other_clones <- '#F0F0F0'
  names(other_clones) <- 'Others'
  cols_use <- c(cols_use,other_clones)
  rna_SA535 <- list()
  for(obs_given_clone in cls){
    tmp1 <- umap_df %>%
      dplyr::filter(clone_id==obs_given_clone)
    for(st in unique(tmp1$Grouping)){
      res <- viz_umap_given_clones(umap_df, datatag, output_dir,
                                   obs_given_clone=obs_given_clone, 
                                   obs_site=st,
                                   cols_use=cols_use,
                                   plottitle='', plotlegend=F)
      rna_SA535[[paste0(obs_given_clone,'_',st)]] <- res  
    }
    
    
  }
  names(rna_SA535)
  p535_total <- cowplot::plot_grid(rna_SA535[['A_Primary']]$p,rna_SA535[['C_Primary']]$p,
                                   rna_SA535[['D_Primary']]$p,rna_SA535[['E_Tumor_Recur']]$p,
                                   rna_SA535[['E_Metastasis']]$p,rna_SA535[['F_Tumor_Recur']]$p,
                                   rna_SA535[['F_Metastasis']]$p,rna_SA535[['F_Metastasis']]$p,
                                   rna_SA535[['G_Primary']]$p, rna_SA535[['G_Tumor_Recur']]$p,
                                   ncol=5, align = 'hv')
  
  p535_total <- cowplot::plot_grid(rna_SA535$A$p,rna_SA535$C$p,rna_SA535$D$p,rna_SA535$E$p,rna_SA535[['F']]$p,
                                   rna_SA535$G$p, rna_SA535$M$p, rna_SA535$N$p, rna_SA535$O$p, rna_SA535$P$p,
                                   rna_SA535$Q$p, rna_SA535$R$p, rna_SA535$S$p, rna_SA535[['T']]$p, rna_SA535$J$p,
                                   ncol=5, align = 'hv')
  # p535_total1 <- cowplot::plot_grid(res_SA535_UnRx$plg, p535_total, ncol=1, rel_heights = c(1,10), labels = c('SA535',''))
  png(paste0(output_dir,datatag,".png"), height = 2*500, width=2*1250,res = 2*72)
  print(p535_total)
  dev.off()
}  
viz_umap_given_clones <- function(umap_df, datatag, output_dir, 
                                obs_given_clone=NULL, obs_site=NULL,cols_use=NULL, 
                                plottitle='', plotlegend=F){
  
  xl <- c(min(umap_df$UMAP_1),max(umap_df$UMAP_1))
  yl <- c(min(umap_df$UMAP_2),max(umap_df$UMAP_2))
  
  my_font <- "Helvetica"
  
  umap_df <- umap_df %>%
    dplyr::mutate(clone_id=replace(clone_id, clone_id=='','None'))
  if(!is.null(obs_given_clone)){
    tmp <- umap_df %>%
      dplyr::filter(clone_id==obs_given_clone)
    plottitle=paste0(plottitle, ' Clone-',obs_given_clone)
  }else{
    tmp <- umap_df 
  }  
  if(!is.null(obs_site)){
    tmp <- tmp %>%
      dplyr::filter(Grouping==obs_site)
    plottitle=paste0(plottitle, ', Site-', obs_site)
  }else{
    print('Observe all sites')
  }  
  
  if(dim(tmp)[1]==0){
    print('No cell, check input data')
    return(NULL)
  }
  
  print(dim(tmp))
  obs_clones <- sort(unique(umap_df$clone_id))
  unassign_clones <- c('unassigned','Unassigned','Un','un','None')
  obs_clones <- obs_clones[!obs_clones %in% unassign_clones]
  print(obs_clones)
  if(is.null(cols_use)){
    cols_use <- make_clone_palette(obs_clones)
    cols_use <- cols_use[obs_clones]
    other_clones <- '#F0F0F0'
    names(other_clones) <- 'Others'
    cols_use <- c(cols_use,other_clones)
  }
  
  # print(cols_use)
  # umap_df <- umap_df %>%
  #   dplyr::mutate(clone=replace(clone, !clone %in% obs_clones,'Others'))
  # umap_df$cell_id[1]
  rownames(umap_df) <- umap_df$cell_id
  # cells_excluded <- umap_df$cell_id[!umap_df$cell_id %in% tmp$cell_id]
  # umap_df[cells_excluded,'clone'] <- 'Others'
  # umap_df$alpha_val <- 1
  # umap_df[cells_excluded,'alpha_val'] <- 0.03
  # umap_df <- umap_df %>%
  #   dplyr::filter(cell_id %in% tmp$cell_id)
  # print(summary(as.factor(umap_df$clone)))
  
  thesis_theme <- get_theme()
  umap_df$clone_id <- factor(umap_df$clone_id, levels = sort(unique(umap_df$clone_id)))
  summary(as.factor(umap_df$clone_id))
  p <- ggplot(umap_df, aes(UMAP_1, UMAP_2)) + 
    geom_point(size=0.5, shape=1, color='#e0e0e0') +
    geom_point(data=tmp, size=0.5, shape=1, aes(color=clone_id)) +  # color='#e0e0e0', grey color for all cells landscape displaying in background
    # geom_point(data=tmp, aes(color=clone), size=0.5, shape=1) +
    # facet_wrap(~ clone_id, nrow=3) +
    scale_color_manual(values=cols_use, name='') +
    labs(title = plottitle, x=' ') +
    xlim(xl[1], xl[2]) + 
    ylim(yl[1], yl[2])
  
  p <- p + thesis_theme
  # p
  # lg <- cowplot::get_legend(p + guides(color = guide_legend(nrow = 1, 
  #                                                           title.position = "left", 
  #                                                           override.aes = list(size=2))))
  # plg <- cowplot::ggdraw() + cowplot::draw_plot(lg)
  if(!plotlegend){
    lg_pos <- "none"
    p <- p + ggplot2::theme(legend.position = lg_pos)  
  }else{
    lg_pos <- "right"
    p <- p + ggplot2::theme(legend.position = lg_pos)  
  }
  
  results <- list(p=p, cols_use=cols_use, obs_clone=obs_given_clone)
  # basename <- paste0(datatag,"_", gsub(' ','',plottitle))
  # saveRDS(results, paste0(output_dir, basename, ".rds"))
  # 
  # png(paste0(output_dir,basename,".png"), height = 2*350, width=2*500,res = 2*72)
  # print(p)
  # dev.off()
  return(results)
  
}

viz_umap_obs_clones <- function(umap_df, cols_use, datatag, output_dir, 
                                obs_treatment, obs_passage=NULL, plottitle='', plotlegend=F){
  xl <- c(min(umap_df$UMAP_1),max(umap_df$UMAP_1))
  yl <- c(min(umap_df$UMAP_2),max(umap_df$UMAP_2))
  
  my_font <- "Helvetica"
  if(!'treatmentSt' %in% colnames(umap_df)){
    if('treat' %in% colnames(umap_df)){
      umap_df <- umap_df %>%
        dplyr::rename(treatmentSt=treat)
    }else{
      stop('Please check input treatmentSt')
    }
  }
  
  umap_df <- umap_df %>%
    dplyr::mutate(clone=replace(clone, is.na(clone),'None'))
  # umap_df$treatmentSt <- umap_df$treat
  if(!is.null(obs_passage)){
    tmp <- umap_df %>%
      dplyr::filter(timepoint==obs_passage)
  }else{
    tmp <- umap_df
  }
  
  if(dim(tmp)[1]==0){
    print('Double check input passage info, exit')
    print('No cell at this passage')
    return(NULL)
  }
  
  print(dim(tmp))
  obs_clones <- sort(unique(tmp$clone))
  unassign_clones <- c('unassigned','Unassigned','Un','un','None')
  obs_clones <- obs_clones[!obs_clones %in% unassign_clones]
  # print(obs_clones)
  # if(is.null(obs_clones)){
  #   
  # }
  cols_use <- cols_use[obs_clones]
  other_clones <- '#F0F0F0'
  names(other_clones) <- 'Others'
  cols_use <- c(cols_use,other_clones)
  # print(cols_use)
  # umap_df <- umap_df %>%
  #   dplyr::mutate(clone=replace(clone, !clone %in% obs_clones,'Others'))
  
  rownames(umap_df) <- umap_df$cell_id
  cells_excluded <- umap_df$cell_id[!umap_df$cell_id %in% tmp$cell_id]
  # umap_df[cells_excluded,'clone'] <- 'Others'
  # umap_df$alpha_val <- 1
  # umap_df[cells_excluded,'alpha_val'] <- 0.03
  # umap_df <- umap_df %>%
  #   dplyr::filter(cell_id %in% tmp$cell_id)
  # print(summary(as.factor(umap_df$clone)))
  
  thesis_theme <- get_theme()
  umap_df$clone <- factor(umap_df$clone, levels = sort(unique(umap_df$clone)))
  p <- ggplot(umap_df, aes(UMAP_1, UMAP_2)) + 
    geom_point(color='#e0e0e0', size=0.5, shape=1) +  # grey color for all cells landscape displaying in background
    geom_point(data=tmp, aes(color=clone), size=0.5, shape=1) + 
    scale_color_manual(values=cols_use, name='') + 
    labs(title = plottitle, x=' ') +
    xlim(xl[1], xl[2]) + 
    ylim(yl[1], yl[2])
  
  p <- p + thesis_theme
  # lg <- cowplot::get_legend(p + guides(color = guide_legend(nrow = 1, 
  #                                                           title.position = "left", 
  #                                                           override.aes = list(size=2))))
  # plg <- cowplot::ggdraw() + cowplot::draw_plot(lg)
  if(!plotlegend){
    lg_pos <- "none"
    p <- p + ggplot2::theme(legend.position = lg_pos)  
  }else{
    lg_pos <- "right"
    p <- p + ggplot2::theme(legend.position = lg_pos)  
  }
  
  results <- list(p=p, cols_use=cols_use, df=tmp)
  # basename <- paste0(datatag,"_", gsub(' ','',plottitle))
  # saveRDS(results, paste0(output_dir, basename, ".rds"))
  # 
  # png(paste0(output_dir,basename,".png"), height = 2*350, width=2*500,res = 2*72)
  # print(p)
  # dev.off()
  return(results)
  
}

make_clone_palette <- function(levels) {
  # install.packages("inlmisc", dependencies = TRUE)  # TO DO: check this package
  clone_names <- sort(levels)
  pal <- as.character(inlmisc::GetColors(length(clone_names)))  
  # if (length(levels) <= 12 & length(levels)>8) {
  #   pal <- brewer.pal(max(length(levels), 3), "Set3")
  # } else if (length(levels) <= 20 & length(levels) > 12) {
  #   pal <- clone_palette_20
  # } else {
  #   pal <- colorRampPalette(brewer.pal(8, "Set2"))(length(levels))
  #   print("WARNING: more clones than palette can accomodate!")
  # }
  pal[pal=='#E7EBFA']<- '#ABB9ED' # too bright, can not see well
  names(pal) <- clone_names
  pal <- pal[levels]
  return(pal)
}

get_unique_clone_id <- function(clone_labels){
  set.seed(42) 
  cls <- sapply(strsplit(clone_labels,'_'), function(x){
    if(length(x)==1){
      return(x[1])
    }else if(length(x)==2){
      idx <- sample(c(1,2),1)
      return(x[idx])
    }else{
      idx <- sample(c(1,2,3),1)
      return(x[idx])
    }
    
  })
  return(as.character(cls))
}
get_color_clones <- function(tag, color_fn){
  # if(datatag=='SA609'){
  #   datatag <-'SA1000'
  # }
  # input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/umap_figs/'
  # color_df <- data.table::fread(paste0(input_dir,'all_clones_new_colors.csv')) %>% as.data.frame()
  # # color_df$cache_path <- NULL
  # color_df <- color_df %>%
  #   dplyr::filter(family == grep(paste0("*",datatag,"*"),unique(color_df$ref_datatag), value=T))#%>%
  # 
  # color_df <- color_df %>%
  #   dplyr::filter(family == grep(paste0("*",datatag,"*"),unique(color_df$ref_datatag), value=T))
  #   # dplyr::select(clone_id, datatag, family, colour)
  # unique(color_df$ref_datatag)
  
  color_df <- data.table::fread(color_fn)
  # color_df <- data.table::fread(paste0(output_dir,'colorcode_total.csv'))
  color_df <- color_df %>%
    dplyr::filter(datatag==tag)
  
  if(dim(color_df)[1] > 0){
    cols_use <- color_df$colour
    names(cols_use) <- color_df$clone_id
    cols_use['None'] <- '#D3D3D3' # unassigned
    # if('None' %in% names(cols_use)){
    #     
    # }
    
  }else{
    # cols_use <- make_clone_palette(obs_clones)
    stop('Error, check color mapping')
  }
  return(cols_use)
}

get_theme <- function(){
  thesis_theme <- ggplot2::theme(
    text = element_text(color="black",size = 8, hjust = 0.5, family=my_font),
    # axis.title.x = element_text(color="black",size=8, hjust = 0.5, family=my_font),
    # axis.title.y = element_text(color="black",size=8, hjust = 0.5, family=my_font),
    # axis.text.x = element_text(color="black",size=7, hjust = 0.5, family=my_font, angle = 90),
    # axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    strip.placement = "outside",
    # axis.line = element_line(colour = "black"),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(color="black",size=11, face="bold",family=my_font, hjust = 0.5),
    # legend.title=element_text(color="black",size=7, hjust = 0.5, family=my_font),
    # legend.text=element_text(color="black",size=7, hjust = 0.5, family=my_font),
    strip.text.x = element_text(color="black",size=11, family=my_font),
    strip.text.y = element_text(color="black",size=11, family=my_font),
    # legend.position = lg_pos,
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-2,-2,-2,-2),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "grey50", fill=NA, size=0.5),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  )
  return(thesis_theme)  
}
plot_all_clones_legends <- function(){
  input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/umap_figs/figs/'
  
  colorcode <- data.table::fread(paste0(input_dir,'colorcode_total.csv')) %>% as.data.frame()
  
  unassign_clones <- c('unassigned','Unassigned','Un','un','None')
  
  for(pdx in unique(colorcode$datatag)){
    df <- colorcode %>%
      dplyr::filter(datatag==pdx)
    obs_clones <- unique(df$clone_id)
    obs_clones <- obs_clones[!obs_clones %in% unassign_clones]
    df <- df %>%
      dplyr::filter(clone_id %in% obs_clones)
    plot_legends(df, input_dir, pdx)
  }
}
