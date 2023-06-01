suppressPackageStartupMessages({
  require("dplyr")
  require("RColorBrewer")
  require("ggplot2")
  require("ComplexHeatmap")
})
clone_palette_20 <- c(
  "#be5f72", "#d74058", "#dc4229", "#a6552c", "#df956f", "#e47a33",
  "#d49f34", "#836e2c", "#b2ad5a", "#92b539", "#4c7d38", "#4dc041",
  "#5dba7f", "#47b8c3", "#6280ca", "#7b57db", "#ce8bd1", "#934f94",
  "#cb48cb", "#d74391"
)
data_frame_to_list <- function(clustering) {
  clusts <- unique(clustering$clone_id)
  res <- list()
  for (cl in clusts) {
    res[[cl]] <- clustering$cell_id[clustering$clone_id == cl]
  }
  return(res)
}
get_color_clone <- function(cell_clones){
  aug_cut <- data_frame_to_list(cell_clones)
  myColors <- get_cluster_colours(length(aug_cut))
  if(length(myColors)>length(aug_cut)){  # in case there are only 2 levels but brewer.pal(nClusters, "Set2") generate at least 3 levels
    myColors <- myColors[1:length(aug_cut)]
  }
  names(myColors) <- names(aug_cut)
  print('DEBUG: myColors')
  print(myColors)
  # TODO: Since 0 is shared with the heatmap CNV state, use something else for the Un-assigned ones
  # myColors <- c(myColors, "0" = "#000000")
  clone_meta <- data.frame(clone_id=names(myColors),color=as.character(myColors), stringsAsFactors = FALSE)
  return(clone_meta)
}

make_clone_palette <- function(levels) {
  # install.packages("inlmisc", dependencies = TRUE)  # TO DO: check this package
  # clone_names <- sort(levels)
  # pal <- as.character(inlmisc::GetColors(length(clone_names)))
  pal <- colorRampPalette(brewer.pal(8, "Set2"))(length(clone_levels))
  names(pal) <- clone_levels
  if (length(levels) <= 12 & length(levels)>8) {
    pal <- brewer.pal(max(length(levels), 3), "Set3")
  } else if (length(levels) <= 20 & length(levels) > 12) {
    pal <- clone_palette_20
  } else {
    # pal <- colorRampPalette(brewer.pal(8, "Set2"))(length(levels))
    pal <- brewer.pal(length(levels), "Set2")
    print("WARNING: more clones than palette can accomodate!")
  }
  # names(pal) <- levels
  names(pal) <- clone_names
  
  
  
  pal <- pal[levels]
  return(pal)
}
# clone_levels <- unique(cell_clones$clone_id)
get_color_clone_v2 <- function(clone_levels, datatag){
  clone_levels <- gtools::mixedsort(clone_levels[!grepl("None", clone_levels)])
  
  ## Meta colors based on package inlmisc::GetColors(n) -- problem at installation, so keep values of color codes in csv file
  # colorcode_fn <- "/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_wholedata/config/colorcode.csv"
  # colorcode_fn <- "/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_wholedata/config/colorcode_v2.csv"
  colorcode_fn <- paste0('/home/htran/Projects/hakwoo_project/metastasis_material/materials/dlp_trees/colorCode_clones/color_code_',datatag,'.csv')
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
  print(clone_pal)
  return(clone_pal)
}  

get_proportion_mathmodel <- function(meta_data, datatag, save_dir){
    input_dir <- '/home/htran/Projects/hakwoo_project/rscript/math_model/data_met_proj_v3/'
    dir.create(input_dir)
  
    meta_data <- meta_data %>%
      dplyr::filter(origin!="Tumor_Recur")
    
    dim(meta_data)
    
    meta_data$passage <- gsub(datatag,'',meta_data$sample_id)
    meta_data$passage <- stringr::str_sub(meta_data$passage, 1, 2)
    unique(meta_data$passage)
    unique(meta_data$mainsite)
    # colnames(meta_data)+
    # summary(as.factor(meta_data$clone_id))
    # View(head(meta_data))
    meta_data$origin <- paste0(meta_data$origin,'_',meta_data$passage)  # add passage to tumour origin site, for SA535: only one passage
    meta_data$origin <- paste0(meta_data$origin,'_',meta_data$pdxid)
    unique(meta_data$origin)
    # meta_data$origin <- ifelse(grepl('Primary',meta_data$origin),'Primary', meta_data$origin)    
    for(m in unique(meta_data$pdxid)){
      # df <- meta_data %>% 
      #   tidyr::pivot_wider(id_cols=c('origin'), 
      #                    names_from='clone_id', values_from = 'origin',
      #                    values_fn = list(origin = length))
      # t <- with(meta_data, table(origin, clone_id))
      # df <- with(meta_data, prop.table(origin, clone_id))
      # df <- prop.table(meta_data$origin, meta_data$clone_id)
      df <- meta_data %>%
        dplyr::filter(pdxid==m) %>%
        dplyr::group_by(origin,clone_id) %>%
        dplyr::summarize(sum=n()) %>%
        tidyr::pivot_wider(id_cols="origin",names_from="clone_id",values_from="sum")%>%
        tibble::column_to_rownames('origin')
      # df[is.na(df)] <- 0
      df[is.na(df)] <- 1
      
      df <- df/rowSums(df)
      # df1[df1==0] <- 0 # just formatting
      print(dim(df))
      # 
      df$clone_id <- rownames(df)
      df <- df %>%
        dplyr::select(clone_id, everything())
      # View(df)
      # m <- 'total'
      # data.table::fwrite(df, paste0(input_dir, datatag, '_', m, '_raw.csv'))
      data.table::fwrite(df, paste0(input_dir, datatag, '_', m, '.csv'))
      
      # df <- data.table::fread(paste0(input_dir, datatag, '_', m, '.csv')) %>% as.data.frame()
      # View(df)
      # df <- df %>% tibble::column_to_rownames('clone_id')
    }
    
}
viz_heatmap_frequencies <- function(df, datatag, save_dir){
  datatag
  save_dir <- input_dir
  cell_func = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.3f", df[i, j]), x, y, gp = gpar(fontsize = 8))
  }
  # library(ComplexHeatmap)
  p <- ComplexHeatmap::Heatmap(as.matrix(df), na_col = "white",
                               show_column_names=T,
                               show_row_names = T,
                               cluster_rows=F,cluster_columns=F,
                               name = "% origin-clones", 
                               # row_order = sort(rownames(df)),
                               # row_split= annotation_row,
                               row_title_rot = 0,
                               row_gap = unit(2, "mm"),
                               # column_split = annotation_col, 
                               column_title = paste0(datatag, ": frequencies"),
                               column_gap = unit(2, "mm"),
                               column_names_gp = grid::gpar(fontsize = 15),
                               row_names_gp = grid::gpar(fontsize = 12),
                               show_heatmap_legend = T,
                               # top_annotation=top_anno,
                               # left_annotation = left_anno,
                               cell_fun = cell_func,
                               row_dend_reorder=F)
  png(paste0(save_dir,datatag,"_clone_prevalence_hm.png"), height = 2*(30*dim(df)[1]+50), width=2*(50*dim(df)[2]+50), res = 2*72)
  print(p)
  dev.off()
  
}
plot_legend <- function(clone_color){
  my_font <- "Helvetica"
  df <- data.frame(color_code=clone_color, clone_id=names(clone_color))
  p <- ggplot(df, aes(fill=clone_id, y=1, x=seq(length(clone_color)))) + 
    geom_bar(position="fill", stat="identity",width=0.55) +
    scale_fill_manual(values = clone_color)
  p + theme_bw() + 
    theme(legend.title=element_blank(),
          # legend.title=element_text(color="black",size=9, hjust = 0.5, family=my_font),
          legend.text=element_text(color="black",face="bold",size=16, hjust = 0.5, family=my_font),
          # legend.key.height= unit(0.2, 'cm'),
          # legend.key.width= unit(0.5, 'cm'),
          legend.key = element_blank()) +
    guides(fill = guide_legend(nrow=3, title.position = "left", 
                                 override.aes = list(shape = 0, size=0.4)))
  plg <- cowplot::ggdraw() + cowplot::draw_plot(cowplot::get_legend(p))
  return(plg)
}
assign_color_clone <- function(cell_clones){
  # cell_clones <- read.delim(cellclones, stringsAsFactors = FALSE, sep = ",")
  aug_cut <- data_frame_to_list(cell_clones)
  myColors <- get_cluster_colours(length(aug_cut))
  names(myColors) <- names(aug_cut)
  # TODO: Since 0 is shared with the heatmap CNV state, use something else for the Un-assigned ones
  # myColors <- c(myColors, "0" = "#000000")
  return(myColors)
}

get_cluster_colours <- function(nClusters) {
  if (nClusters > 8) {
    clust_colours <- colorRampPalette(brewer.pal(8, "Set2"))(nClusters)
  } else {
    clust_colours <- brewer.pal(nClusters, "Set2")
  }
  return(clust_colours)
}
# get_cluster_colours <- function(nClusters) {
#   if (nClusters > 8) {
#     #clust_colours <- colorRampPalette(brewer.pal(12, "Set3"))(nClusters)
#     clust_colours <- colorRampPalette(brewer.pal(8, "Set2"))(nClusters)
#   } else {
#     clust_colours <- brewer.pal(nClusters, "Set2")
#   }
#   # Todo: the 5th and 8th (E and H) have similar colours, change them
#   if (nClusters == 9) {
#     clust_colours[[5]] <- '#564527'
#     
#     # change I to pink too
#     clust_colours[[9]] <- '#FFC0CB'
#   }
#   
#   clust_colours
#   
# }


plot_prevalences_origin <- function(cell_clones, results_dir, save_dir, 
                             datatag, meta_data=NULL, clone_color=NULL){
  
  tag <- ''
  if (!file.exists(save_dir)){
    dir.create(save_dir)
  }
  if(is.null(meta_data)){
    meta_data <- get_meta_data(cell_clones, results_dir)
  }
  table(meta_data$mainsite, meta_data$clone_id)
  # if(is.null(clone_color)){
  #   # use same color as in dream tree 
    # input_fn <- paste0(results_dir,'tree_viz_dream/',datatag,'_meta_data_ls.rds')
  #   if(file.exists(input_fn)){
  #     meta_data_ls <- readRDS(input_fn)
  #     clone_color = meta_data_ls$clone_meta
  #     clone_color <- clone_color[clone_color$clone_id!='None',]
  #     clone_color <- clone_color[order(clone_color$clone_id),]
  #     dim(clone_color)
  #     rownames(clone_color) <- clone_color$clone_id
  #   } else{
  #     clone_color <- get_color_clone(cell_clones)
  #     rownames(clone_color) <- clone_color$clone_id
  #   }
  # }
  clone_color <- get_color_clone_v2(unique(cell_clones$clone_id), datatag)
  # View(clone_color)
  colorcode <- clone_color
  get_pct_origin_v2(meta_data, colorcode, save_dir,
                 save_fn='primary_origin',
                 datatag, xlabel=NULL,
                 ylabel="Clone Fraction",
                 plottitle=NULL)
  
  
  # meta_data_tmp <- meta_data[meta_data$mainsite=='Primary',]
  # print('Debug')
  # print(dim(meta_data_tmp))
  # if(dim(meta_data_tmp)[1]>0){
  #   clone_color_pri <- clone_color[gtools::mixedsort(unique(meta_data_tmp$clone_id))]
  #   colorcode <- clone_color_pri
  #   # clone_color_pri <- clone_color[unique(meta_data_tmp$clone_id),]
  #   # clone_color_pri <- clone_color_pri[order(clone_color_pri$clone_id),]
  #   # View(clone_color_pri)
  #   # colorcode <- clone_color_pri$color
  #   # get_pct_origin(meta_data_tmp, colorcode, save_dir,
  #   #                save_fn='primary_origin',
  #   #                datatag, xlabel="Origin",
  #   #                ylabel="Relative Clone Prevalence",
  #   #                plottitle=paste0("Primary - Origin ",tag))
  #   get_pct_origin(meta_data_tmp, colorcode, save_dir,
  #                  save_fn='primary_origin',
  #                  datatag, xlabel=NULL,
  #                  ylabel="Clone Fraction",
  #                  plottitle=paste0("Primary cells",tag))
  #   
  #   # get_pct_pdx(meta_data_tmp, colorcode, save_dir,
  #   #             save_fn='primary_pdx',
  #   #             datatag, xlabel="Pdx",
  #   #             ylabel="Relative Clone Prevalence",
  #   #             plottitle=paste0("Primary - Pdx ",tag))
  # }
  # 
  # 
  # meta_data_tmp <- meta_data[meta_data$mainsite=='Metastasis',]
  # print('Debug')
  # print(dim(meta_data_tmp))
  # colnames(meta_data_tmp)
  # summary(as.factor(meta_data_tmp$pdxid))
  # if(dim(meta_data_tmp)[1]>0){
  #   # clone_color_meta <- clone_color[unique(meta_data_tmp$clone_id),]
  #   # clone_color_meta <- clone_color_meta[order(clone_color_meta$clone_id),]
  #   # colorcode <- clone_color_meta$color
  #   clone_color_meta <- clone_color[gtools::mixedsort(unique(meta_data_tmp$clone_id))]
  #   colorcode <- clone_color_pri
  #   get_pct_origin(meta_data_tmp, colorcode, save_dir,
  #                  save_fn='metastasis_origin',
  #                  datatag, xlabel=NULL,
  #                  ylabel="Clone Fraction",
  #                  plottitle=paste0("Metastasis cells",tag))
  #   # get_pct_origin(meta_data_tmp, colorcode, save_dir,
  #   #                save_fn='metastasis_origin',
  #   #                datatag, xlabel="Origin",
  #   #                ylabel="Relative Clone Prevalence",
  #   #                plottitle=paste0("Metastasis - Origin ",tag))
  #   
  #   # get_pct_pdx(meta_data_tmp, colorcode, save_dir,
  #   #             save_fn='metastasis_pdx',
  #   #             datatag, xlabel="Pdx",
  #   #             ylabel="Relative Clone Prevalence",
  #   #             plottitle=paste0("Metastasis - Pdx ",tag))
  #   
  #   for(pd in unique(meta_data_tmp$pdxid)){
  #     meta_pdx <- meta_data_tmp[meta_data_tmp$pdxid==pd,]
  #     clone_color_meta_pd <- clone_color[unique(meta_pdx$clone_id),]
  #     clone_color_meta_pd <- clone_color_meta_pd[order(clone_color_meta_pd$clone_id),]
  #     colorcode_pd <- clone_color_meta_pd$color
  #     
  #     
  #     get_pct_origin(meta_pdx, colorcode_pd, save_dir,
  #                    save_fn=paste0('metastasis_origin_',pd),
  #                    datatag, xlabel="Metastasis Origin",
  #                    ylabel="Relative Clone Prevalence",
  #                    plottitle=paste0(" ",pd))
  #   }
  # }  
  
  
}

get_library_id <- function(cell_ids) {
  labels <- sapply(strsplit(cell_ids, "-"), function(x) {
    return(x[2])
  })
  return(as.character(labels))
}
get_sample_id <- function(cell_ids) {
  labels <- sapply(strsplit(cell_ids, "-"), function(x) {
    return(x[1])
  })
  return(as.character(labels))
}

get_meta_data <- function(cell_clones, results_dir, library_grouping_fn=NULL){
  
  # cell_clones <- read.delim(cellclones, stringsAsFactors = FALSE, sep = ",")
  # dim(cell_clones)
  if(is.null(library_grouping_fn)){
    library_grouping_fn <- paste0(results_dir,'library_groupings.csv')
  }
  
  grouping_df <- read.csv(library_grouping_fn, header=T,check.names=F, 
                          stringsAsFactors = F)
  # colnames(grouping_df)[which(names(grouping_df) == "passage")] <- "timepoint"
  # colnames(grouping_df)[which(names(grouping_df) == "mainsite")] <- "mainsite"
  colnames(grouping_df)[which(colnames(grouping_df) == "grouping")] <- "library_id"
  # dim(grouping_df)
  # length(unique(grouping_df$library_id))
  print(colnames(grouping_df))
  # lib_ids <- strsplit(as.character(cell_clones$cell_id),"-")
  # length(lib_ids)
  # nbcores <- 5
  # ps <- c()
  # lids <- parallel::mclapply(lib_ids, function(f) {
  #   ps <- c(ps, as.character(f[2]))
  # }, mc.cores = nbcores)
  # length(lids)
  # length(unique(cell_clones$library_id))
  # sum(unique(cell_clones$library_id) %in% unique(grouping_df$library_id))
  cell_clones$library_id <- get_library_id(cell_clones$cell_id)
  cell_clones$sample_id <- get_sample_id(cell_clones$cell_id)
  cell_clones <- cell_clones %>% left_join(grouping_df, by = c("library_id","sample_id"))
  
  # tm_df <- data.frame(treatmentSt=c('UT','UTT','UTTT','UTTTT',
  #                                   'U','UU','UUU','UUUU','UUUUU',
  #                                   'UTU','UTTU','UTTTU'),
  #                     passage=c(rep('absolute_treatment',4),
  #                               rep('untreated',5),
  #                               rep('drug_vacant',3)), stringsAsFactors = F)
  
  
  
  # cell_clones <- cell_clones %>% left_join(tm_df, by = "treatmentSt")
  return(cell_clones)
}

plot_function <- function(meta_data, xstring, ystring, plottype, 
                          plottitle="Cluster 0", 
                          colorcode="blue", xlabel='time point', ylabel='pct',
                          yl=NULL) {
  
  p <- ggplot(meta_data, aes_string(x=xstring, y=ystring, group = 1)) +  #, 
    # geom_line(aes_string(color=plottype)) +   #,color=colorcode linetype=plottype, 
    geom_line(color=colorcode) +   #,
    geom_point(aes_string(shape=xstring))  #+aes_string(shape=xstring)
  # scale_color_manual(values = colorcode) #+ labs(colour=legendtitle)
  
  if(!is.null(yl)){
    p <- p + ylim(yl[1],yl[2])
  }
  p <- p + labs(x=xlabel,y=ylabel,title=plottitle)
  p <- p + theme(legend.title = element_blank(), 
                 plot.title = element_text(color="black", size=13,hjust = 0.5),
                 legend.position = "none", 
                 # axis.line = element_blank(), 
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                 panel.border = element_blank(),
                 
                 axis.text.x = element_text(color="black", size=9, hjust = 0.5),
                 axis.text.y = element_text(color="black", size=9, hjust = 0.5)
                 #axis.title = element_blank(),
                 # axis.ticks = element_blank()
  )
  # p
  return(p)
  
}  

legend_function <- function(meta_data, xstring, ystring, plottype, 
                            plottitle="Cluster 0", colorcode="blue", 
                            legendtitle="Cluster_Prevalence", legendlabels=NULL,
                            xlabel='Treatment Status', ylabel='Percentage') {
  p <- ggplot(meta_data, aes_string(x=xstring, y=ystring, group=plottype)) +
    geom_line(aes_string(color=plottype)) +  #,color=colorcode
    # geom_line(aes(color=colorcode)) 
    scale_color_manual(values = colorcode) + labs(colour=legendtitle)
  
  p <- p + labs(x=xlabel,y=ylabel,title=plottitle)
  p <- p + theme(legend.title = element_blank(), #element_text(color="black", size=9,hjust = 0.5),
                 legend.text = element_blank(),#element_text(color="black", size=9,hjust = 0.5),
                 plot.title = element_text(color="black", size=13,hjust = 0.5),
                 legend.position = "none",
                 # axis.line = element_blank(), 
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                 panel.border = element_blank(),
                 axis.text.x = element_text(color="black", size=9, hjust = 0.5),
                 axis.text.y = element_text(color="black", size=9, hjust = 0.5)
                 #axis.title = element_blank(),
                 # axis.ticks = element_blank()
  )
  # p <- p + theme(
  #   legend.title = element_text(size=17), 
  #   legend.text = element_text(size=16), 
  #   plot.title = element_text(color="black", size=20, hjust = 0.5)
  # )
  return(p)
}

# meta_data <- cls_tp_df

plot_passage_clusters <- function(output_fn, meta_data, color_ls, datatag, tag, save_dir, 
                                  xstring="passage",ystring="clone_freq",
                                  plttitle="Treatment Prevalence",
                                  lgtitle="Treatment Status",xlb='Treatment Status'){
  
  
  
  
  plottype <- "clone_id"
  plots <- list()
  idx = 0
  lscls <- sort(unique(meta_data$clone_id))
  
  # Plot treament status percentage for each cluster
  for(c in lscls){
    meta_data_tmp <- meta_data[meta_data$clone_id==c,]
    idx <- idx+1
    p <- plot_function(meta_data_tmp, xstring, ystring, plottype,
                       plottitle=paste0("Clone ",c), colorcode=as.character(color_ls[c]), 
                       xlabel=' ', ylabel='pct', yl=c(0,0.36))
    
    plots[[idx]] <- p
  }
  # length(plots)
  # plots[[1]]
  # Adding legend
  plg <- legend_function(meta_data, xstring, ystring, plottype, 
                         plottitle=plttitle, colorcode=color_ls, 
                         legendtitle=lgtitle, xlabel=xlb, ylabel='Percentage')
  
  
  # legend1 <- cowplot::get_legend(plg)
  plots[[idx+1]] <- plg
  
  # library(gridExtra)
  # # install.packages("gridExtra")
  # layout <- rbind(c(1,2,3,4),
  #                 c(5,6,7,NA))
  # select_grobs <- function(lay) {
  #   id <- unique(c(t(lay)))
  #   id[!is.na(id)]
  # }
  p <- cowplot::plot_grid(plotlist = plots, align = "hv", ncol=4)
  
  png(output_fn, height = 2*700, width=2*1200,res = 2*72)
  # pdf(paste(save_dir,"cluster_prevalence.pdf",sep=""), height=6, width=13)
  # print(grid.arrange(grobs = plots[select_grobs(layout)], layout_matrix = layout,
  #                    bottom=" ",right=" "))
  print(p)
  dev.off()
  
  # pc <- plot_grid(plotlist = plots, ncol = 4, align = 'hv')
  # png(paste(save_dir,datatag,'_',tag,"cluster_prevalence.png",sep=""),height = 2*800, width=2*1000,res = 2*72)
  # # pdf(paste(save_dir,"cluster_prevalence.pdf",sep=""), height=6, width=13)
  # print(pc)
  # 
  # dev.off()
  
  # pdf(paste(save_dir,"cluster_prevalence.pdf",sep=""), height=5, width=11)
  # print(grid.arrange(grobs = plots[select_grobs(layout)], layout_matrix = layout,
  #                    bottom=" ",right=" "))
  # 
  # dev.off()
  
} 

get_pct_clusters <- function(meta_data, save_dir=NULL,
                             datatag='SA1035',xlabel="Treatment Status",
                             ylabel="Relative Clone Prevalence",
                             plottitle=" "){
  
  if (!file.exists(save_dir)){
    dir.create(save_dir)
  }
  # colnames(meta_data)
  tab <- with(meta_data, table(passage, clone_id))
  rs <- prop.table(tab, margin = 1)
  cls_df <- as.data.frame(rs)
  rownames(cls_df) <- paste0(cls_df$clone_id, "_", cls_df$passage)
  cls_df$desc <- paste0(cls_df$clone_id, "_", cls_df$passage)
  print(head(cls_df))
  
  # legend_title <- "Cluster"
  
  # pdf(paste0(save_dir,"cluster_percent.pdf"),width=8, height=5)
  p <- ggplot(cls_df, aes(fill=clone_id, y=Freq, x=passage)) + 
    geom_bar(position="fill", stat="identity", width=0.7)
  p <- p + labs(x=xlabel, y=ylabel, title=plottitle)
  p <- p + theme(legend.title = element_blank(), 
                 legend.text = element_text(size=20),
                 plot.title = element_text(color="black", size=20,hjust = 0.5),
                 # legend.position = "none", 
                 # axis.line = element_blank(), 
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                 panel.border = element_blank(),
                 # axis.text = element_blank(),
                 #axis.title = element_blank(),
                 axis.ticks = element_blank()
  )
  png(paste(save_dir,"clone_prevalence.png",sep="/"), height = 2*500, width=2*800, res = 2*72)
  print(p)
  dev.off()
  return(cls_df) 
}

plot_fill_barplot <- function(cls_df, colorcode, xlabel, ylabel, plottitle, fa='clone_id', xa='origin'){
  p <- ggplot(cls_df, aes_string(fill=fa, y='Freq', x=xa)) + 
    geom_bar(position="fill", stat="identity") + 
    scale_fill_manual(values = colorcode)
  p <- p + labs(x=xlabel, y=ylabel, title=plottitle)
  p <- p + theme(legend.title = element_text(size=12), 
                 legend.text = element_text(size=12),
                 plot.title = element_text(color="black", size=18, hjust = 0.5),
                 legend.position = "top",
                 # axis.line = element_blank(), 
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                 panel.border = element_blank(),
                 axis.text.x = element_text(color="black", size=11, hjust=0.5, angle = 90),
                 axis.text.y = element_text(color="black", size=11, hjust = 0.5),
                 axis.title = element_text(color="black", size=15, hjust = 0.5),
                 axis.ticks = element_blank()
  )
  return(p)
}  

# include mainsite as a facet_grid
plot_fill_barplot_v2 <- function(cls_df, colorcode, xlabel, ylabel, plottitle, fa='clone_id', xa='origin'){
  p <- ggplot(cls_df, aes_string(fill=fa, y='Freq', x=xa)) + 
    geom_bar(position="fill", stat="identity",width =0.65) + 
    facet_grid(rows = vars(mainsite)) + 
    scale_fill_manual(values = colorcode)
  p <- p + labs(x=xlabel, y=ylabel, title=plottitle)
  p <- p + theme(legend.title = element_text(size=11), 
                 legend.text = element_text(size=12),
                 plot.title = element_text(color="black", size=18, hjust = 0.5),
                 legend.position = "top",
                 # axis.line = element_blank(), 
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                 panel.border = element_blank(),
                 axis.text.x = element_text(color="black", size=13, hjust=0.5, angle = 90),
                 axis.text.y = element_text(color="black", size=11, hjust = 0.5),
                 axis.title = element_text(color="black", size=15, hjust = 0.5),
                 axis.ticks = element_blank(),
                 strip.text.y = element_text(size=13, color="black")
  )
  return(p)
}  

plot_stack_barplot <- function(cls_df, colorcode, xlabel, ylabel, plottitle, fa='clone_id', xa='origin'){
  # df_cumsum <- ddply(cls_df, "timepoint",
  #                    transform, 
  #                    label_ypos=cumsum(Freq) - 0.5*Freq + 1)
  # 
  p <- ggplot(cls_df, aes_string(fill=fa, y='Freq', x=xa)) + 
    geom_bar(position="dodge", stat="identity") + 
    # geom_text(aes(y=label_ypos, label=Freq), vjust=1.6, 
    #           color="white", size=3.5) +
    scale_fill_manual(values = colorcode)
  p <- p + labs(x=xlabel, y=ylabel, title=plottitle)
  # p <- p + facet_wrap(~timepoint)
  p <- p + theme(legend.title = element_blank(), 
                 legend.text = element_text(size=12),
                 plot.title = element_text(color="black", size=15,hjust = 0.5),
                 # legend.position = "none", 
                 # axis.line = element_blank(), 
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                 panel.border = element_blank(),
                 axis.text.x = element_text(color="black", size=11, hjust=0.5, angle = 90),
                 axis.text.y = element_text(color="black", size=11, hjust = 0.5),
                 axis.title = element_text(color="black", size=13, hjust = 0.5),
                 axis.ticks = element_blank()
  )
  return(p)
}
get_pct_pdx <- function(meta_data, colorcode, save_dir=NULL,
                        save_fn='drug_holiday',
                        datatag='SA1035', xlabel="Time Point",
                        ylabel="Relative Clone Prevalence",
                        plottitle=" "){
  # colnames(meta_data)
  tab <- with(meta_data, table(pdxid, clone_id))
  rs <- prop.table(tab, margin = 1)
  cls_df <- as.data.frame(rs)
  cls_df_count <- as.data.frame(tab)
  rownames(cls_df) <- paste0(cls_df$clone_id, "_", cls_df$pdxid)
  cls_df$desc <- paste0(cls_df$clone_id, "_", cls_df$pdxid)
  print(head(cls_df))
  nb_pdx <- length(unique(cls_df$pdxid))
  
  p <- plot_stack_barplot(cls_df_count, colorcode, xlabel, ylabel, plottitle, fa='clone_id', xa='pdxid')
  
  p1 <- plot_fill_barplot(cls_df, colorcode, xlabel, ylabel, plottitle, fa='clone_id', xa='pdxid')
  
  write.csv(cls_df_count, paste0(save_dir,datatag,'_',save_fn,"_clone_prevalence_counts.csv"), quote=F, row.names=T)
  write.csv(cls_df, paste0(save_dir,datatag,'_',save_fn,"_clone_prevalence.csv"), quote=F, row.names=T)
  
  # plottitle <- "All Drug Holiday Samples"
  # legend_title <- "Cluster"
  
  # pdf(paste0(save_dir,"cluster_percent.pdf"),width=8, height=5)
  
  
  wunit1 <- 160
  wunit2 <- 210
  png(paste0(save_dir,datatag,'_',save_fn,"_clone_prevalence.png"), height = 2*520, width=2*(wunit1*nb_pdx+50), res = 2*72)
  print(p1)
  dev.off()
  png(paste0(save_dir,datatag,'_',save_fn,"_clone_prevalence_counts.png"), height = 2*430, width=2*(wunit2*nb_pdx+50), res = 2*72)
  print(p)
  dev.off()
  return(cls_df) 
}

get_pct_origin_v2 <- function(meta_data, colorcode, save_dir=NULL,
                           save_fn='drug_holiday',
                           datatag='SA1035', xlabel="Time Point",
                           ylabel="Relative Clone Prevalence",
                           plottitle=" "){
  # colnames(meta_data)
  tab <- with(meta_data, table(mainsite, origin, clone_id))
  rs <- prop.table(tab, margin = 1)
  cls_df <- as.data.frame(rs)
  cls_df$mainsite <- factor(cls_df$mainsite, levels=c('Primary','Metastasis'))
  sites <- unique(cls_df$origin)
  sites <- sites[sites!='Primary']
  cls_df$origin <- factor(cls_df$origin, levels=c('Primary',as.character(sites)))
  cls_df_count <- as.data.frame(tab)
  rownames(cls_df) <- paste0(cls_df$clone_id, "_", cls_df$origin,'_',cls_df$mainsite)
  cls_df$desc <- paste0(cls_df$clone_id, "_", cls_df$origin)
  # print(head(cls_df))
  nb_origin <- length(unique(cls_df$origin))
  
  # p <- plot_stack_barplot(cls_df_count, colorcode, xlabel, ylabel, plottitle, fa='clone_id', xa='origin')
  
  p1 <- plot_fill_barplot_v2(cls_df, colorcode, xlabel, ylabel, plottitle, fa='clone_id', xa='origin')
  
  # write.csv(cls_df_count, paste0(save_dir,datatag,'_',save_fn,"_clone_prevalence_counts.csv"), quote=F, row.names=T)
  # write.csv(cls_df, paste0(save_dir,datatag,'_',save_fn,"_clone_prevalence.csv"), quote=F, row.names=T)
  
  # plottitle <- "All Drug Holiday Samples"
  # legend_title <- "Cluster"
  
  # pdf(paste0(save_dir,"cluster_percent.pdf"),width=8, height=5)
  
  
  wunit1 <- 50
  wunit2 <- 210
  png(paste0(save_dir,datatag,'_',save_fn,"_clone_prevalence.png"), height = 2*500, width=2*(wunit1*nb_origin+50), res = 2*72)
  print(p1)
  dev.off()
  # png(paste0(save_dir,datatag,'_',save_fn,"_clone_prevalence_counts.png"), height = 2*430, width=2*(wunit2*nb_origin+50), res = 2*72)
  # print(p)
  # dev.off()
  # saveRDS()
  # return(cls_df) 
}

get_pct_origin <- function(meta_data, colorcode, save_dir=NULL,
                              save_fn='drug_holiday',
                              datatag='SA1035', xlabel="Time Point",
                              ylabel="Relative Clone Prevalence",
                              plottitle=" "){
  # colnames(meta_data)
  tab <- with(meta_data, table(origin, clone_id))
  rs <- prop.table(tab, margin = 1)
  cls_df <- as.data.frame(rs)
  cls_df_count <- as.data.frame(tab)
  rownames(cls_df) <- paste0(cls_df$clone_id, "_", cls_df$origin)
  cls_df$desc <- paste0(cls_df$clone_id, "_", cls_df$origin)
  print(head(cls_df))
  nb_origin <- length(unique(cls_df$origin))
  
  # p <- plot_stack_barplot(cls_df_count, colorcode, xlabel, ylabel, plottitle, fa='clone_id', xa='origin')
  
  p1 <- plot_fill_barplot(cls_df, colorcode, xlabel, ylabel, plottitle, fa='clone_id', xa='origin')
  
  write.csv(cls_df_count, paste0(save_dir,datatag,'_',save_fn,"_clone_prevalence_counts.csv"), quote=F, row.names=T)
  write.csv(cls_df, paste0(save_dir,datatag,'_',save_fn,"_clone_prevalence.csv"), quote=F, row.names=T)
  
  # plottitle <- "All Drug Holiday Samples"
  # legend_title <- "Cluster"
  
  # pdf(paste0(save_dir,"cluster_percent.pdf"),width=8, height=5)
  
  
  wunit1 <- 160
  wunit2 <- 210
  png(paste0(save_dir,datatag,'_',save_fn,"_clone_prevalence.png"), height = 2*520, width=2*(wunit1*nb_origin+50), res = 2*72)
  print(p1)
  dev.off()
  # png(paste0(save_dir,datatag,'_',save_fn,"_clone_prevalence_counts.png"), height = 2*430, width=2*(wunit2*nb_origin+50), res = 2*72)
  # print(p)
  # dev.off()
  # saveRDS()
  return(cls_df) 
}

get_proportion_plt <- function(data, plottitle='', cols=NULL){
  
  # Remove small cells clusters
  # data <- data %>%
  #   dplyr::filter(count>=10)
  if(is.null(cols)){
    cols <- get_color_clone_v2(as.character(unique(data$label)))
  }
  
  # Compute percentages
  # data$label <- factor(data$label, levels = names(cols))
  rownames(data) <- data$label
  data <- data[names(cols)[names(cols) %in% as.character(data$label)],]
  data$fraction <- data$count / sum(data$count)
  
  # Compute the cumulative percentages (top of each rectangle)
  data$ymax <- cumsum(data$fraction)
  
  # Compute the bottom of each rectangle
  data$ymin <- c(0, head(data$ymax, n=-1))
  
  # Compute label position
  data$labelPosition <- (data$ymax + data$ymin) / 2
  
  # Compute a good label
  # data$label <- data$category
  # data$label <- c('1Rx','2Rx','3Rx','1RxH')
  cols_use <- cols[as.character(unique(data$label))]
  
  p <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=label)) +
    geom_rect() +
    # geom_label(x=4.5, aes(y=labelPosition, label=label), size=4, label.size=0) +
    scale_fill_manual(values = cols_use) +
    coord_polar(theta="y") +
    xlim(c(2, 6)) +
    theme_void() +
    theme(legend.position = "none",
          plot.title = element_text(color="black", size=13, hjust = 0.5, face = "bold"),
          panel.background = element_rect(fill = "transparent",
                                          colour = NA), # necessary to avoid drawing panel outline
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          plot.background = element_rect(fill = "transparent",
                                         colour = NA)) # necessary to avoid drawing plot outline)
  # p
  p <- p + labs(title=plottitle)   
  
  return(p)
  
}
make_discrete_palette <- function(pal_name, levels) {
  if (length(levels) <= 8) {
    pal <- brewer.pal(max(length(levels), 3), pal_name)
  } else if (length(levels) <= 12) {
    pal <- brewer.pal(max(length(levels), 3), "Set3")
  } else if (length(levels) <= 20 & length(levels)>8) {
    pal <- clone_palette_20
  } else {
    pal <- clone_palette_20
    print("WARNING: more clones than palette can accomodate!")
  }
  names(pal) <- levels
  pal <- pal[levels]
  return(pal)
}

# Donut plots
plot_clones_proportion <- function(results_dir, outputfile=NULL, 
                                   cellclones=NULL, datatag = ''){
  print("Generating donut plots per mouse id")
  if(!endsWith(results_dir,'/')){
    results_dir <- paste0(results_dir,'/')  
  }
  if(is.null(cellclones)){
    cellclones <- paste0(results_dir, 'cell_clones.csv')
  }
  save_dir <- paste0(dirname(outputfile),'/')
  if (!file.exists(save_dir)){
    dir.create(save_dir)
  }
  cell_clones <- read.delim(cellclones, stringsAsFactors = FALSE, sep = ",")
  excluded_clones <- c('Un','unassigned','Unassigned','None')
  cell_clones <- cell_clones[!cell_clones$clone_id %in% excluded_clones,]
  print(dim(cell_clones))
  print(summary(as.factor(cell_clones$clone_id)))
  # cell_clones_backup <- cell_clones
  library_grouping_fn <- paste0(results_dir,'library_groupings.csv')
  meta_data <- get_meta_data(cell_clones, results_dir, library_grouping_fn)
  print(dim(meta_data))
  # unique(meta_data$pdxid), originx
  colnames(meta_data)
  unique(meta_data$origin)
  
  # if(datatag=='SA919'){
  #   # clone_color <- make_discrete_palette("Set2", gtools::mixedsort(unique(cell_clones$clone_id)))
  #   clone_color <- get_color_clone_v2(unique(cell_clones$clone_id), datatag)  
  # }
  clone_color <- get_color_clone_v2(unique(cell_clones$clone_id), datatag)  
  
  plg <- plot_legend(clone_color)
  ggsave(
    plot = plg,
    height = 10,
    width = 1,
    # To Do: add height, width here
    filename = paste0(save_dir,"clones_legend.png"),
    bg = "transparent"
  )  
  
  mainsite_levels <- gtools::mixedsort(unique(meta_data$mainsite))
  mainsites_cols <- make_discrete_palette("Set1", mainsite_levels)
  plttype <- 'origin'
  save_dir_fg <- paste0(save_dir,plttype,'/')
  if (!file.exists(save_dir_fg)){
    dir.create(save_dir_fg)
  }
  unique(meta_data$pdxid)
  for(pd in unique(meta_data$pdxid)){
    for(or in unique(meta_data$origin)){
      tmp <- meta_data %>%
        dplyr::filter(origin==or & pdxid==pd) %>%
        dplyr::select(origin,clone_id)
      
      # tmp2 <- meta_data %>%
      #   dplyr::filter(origin==or)%>%
      #   dplyr::select(origin,mainsite)
      
      if(dim(tmp)[1]>0){
        print(paste0(pd,': ',or))
        print(dim(tmp))
        print(summary(as.factor(tmp$clone_id)))
        out <- table(tmp$clone_id, tmp$origin)  %>% as.data.frame()
        colnames(out) <- c('label','cluster', 'count')
        # print(out)
        # View(clone_color)
        colorcode <- clone_color[as.character(out$label)]
        plttitle <- gsub('_',' ',or)
        # if(plttitle=='VentralSpinal'){
        #   plttitle <- 'Ventral Spinal'
        # }else if(plttitle=='SupraSpinal'){
        #   plttitle <- 'Supra Spinal'
        # }
        # if(out$count[1]>0){
          p <- get_proportion_plt(out, NULL, colorcode)
          # png(paste0(save_dir_fg,pd,"_",or,"_",plttype, "_clone.png"),height = 2*200, width=2*200,res = 2*72)
          # print(p)
          # dev.off()
          ggsave(
            plot = p,
            height = 0.5,
            width =  0.5,
            # To Do: add height, width here
            filename = paste0(save_dir_fg,pd,"_",or,"_",plttype, "_clone.png"),
            bg = "transparent"
          )  
        # }
        
        # out2 <- table(tmp2$mainsite, tmp2$origin)  %>% as.data.frame()
        # colnames(out2) <- c('label','cluster', 'count')
        # print(out2)
        # plttitle <- gsub('_',' ',or)
        # p2 <- get_proportion_plt(out2, plttitle, mainsites_cols)
        # png(paste0(save_dir_fg,or,"_",plttype, "_mainsite.png"),height = 2*250, width=2*250,res = 2*72)
        # print(p2)
        # dev.off()
      }
    }  
    
  }
  
  
  
  plttype <- 'pdxid'
  save_dir_fg <- paste0(save_dir,plttype,'/')
  if (!file.exists(save_dir_fg)){
    dir.create(save_dir_fg)
  }
  for(pdx in unique(meta_data$pdxid)){
    tmp <- meta_data %>%
      dplyr::filter(pdxid==pdx)%>%
      dplyr::select(pdxid, clone_id)
    
    tmp2 <- meta_data %>%
      dplyr::filter(pdxid==pdx)%>%
      dplyr::select(pdxid, mainsite)
    print(pdx)
    print(dim(tmp))
    if(dim(tmp)[1]>0){
      out <- table(tmp$clone_id, tmp$pdxid)  %>% as.data.frame()
      colnames(out) <- c('label','cluster', 'count')
      # print(out)
      # View(clone_color)
      colorcode <- clone_color[as.character(out$label)]
      plttitle <- gsub('_',' ',pdx)
      p <- get_proportion_plt(out, plttitle, colorcode)
      png(paste0(save_dir_fg, pdx,"_",plttype, "_clone.png"),height = 2*250, width=2*250,res = 2*72)
      print(p)
      dev.off()
      
      out2 <- table(tmp2$mainsite, tmp2$pdxid)  %>% as.data.frame()
      colnames(out2) <- c('label','cluster', 'count')
      print(out2)
      p2 <- get_proportion_plt(out2, plttitle, mainsites_cols)
      png(paste0(save_dir_fg, pdx,"_",plttype, "_mainsite.png"),height = 2*250, width=2*250,res = 2*72)
      print(p2)
      dev.off()
    }
    
  }
  
  plttype <- 'pdxid_origin_SA919'
  save_dir_fg <- paste0(save_dir,plttype,'/')
  if (!file.exists(save_dir_fg)){
    dir.create(save_dir_fg)
  }
  for(pdx in unique(meta_data$pdxid)){
    tmp <- meta_data %>%
      dplyr::filter(pdxid==pdx)#%>%
      # dplyr::select(pdxid, clone_id)
    
    # tmp2 <- tmp %>%
    #   dplyr::filter(pdxid==pdx)%>%
    #   dplyr::select(pdxid, mainsite)
    print(pdx)
    print(dim(tmp))
    if(dim(tmp)[1]>0){
      for(or in unique(tmp$origin)){
        
        tmp2 <- tmp %>%
          dplyr::filter(origin==or)
        print(or)
        out <- table(tmp2$clone_id, tmp2$pdxid)  %>% as.data.frame()
        colnames(out) <- c('label','cluster', 'count')
        # print(out)
        # View(clone_color)
        colorcode <- clone_color[as.character(out$label)]
        plttitle <- paste0(gsub('X0011_','',pdx),'_',or)
        p <- get_proportion_plt(out, plttitle, colorcode)
        png(paste0(save_dir_fg, plttitle,"_",plttype, ".png"),height = 2*250, width=2*250,res = 2*72)
        print(p)
        dev.off()
        
      }
     
    }
    
  }
  
  
  plttype <- 'pdxid_mainsites'
  save_dir_fg <- paste0(save_dir,plttype,'/')
  if (!file.exists(save_dir_fg)){
    dir.create(save_dir_fg)
  }
  for(pdx in unique(meta_data$pdxid)){
    meta_data1 <- meta_data %>%
      dplyr::filter(pdxid==pdx)
    
    tmp <- meta_data1 %>%
      dplyr::filter(mainsite=='Primary')%>%
      dplyr::select(pdxid,clone_id)
    
    tmp2 <- meta_data1 %>%
      dplyr::filter(mainsite=='Metastasis')%>%
      dplyr::select(pdxid,clone_id)
    print(pdx)
    print(dim(tmp))
    if(dim(tmp)[1]>0){
      out <- table(tmp$clone_id, tmp$pdxid)  %>% as.data.frame()
      if(dim(out)[1]>0){
        colnames(out) <- c('label','cluster', 'count')
        # print(out)
        # View(clone_color)
        colorcode <- clone_color[as.character(out$label)]
        plttitle <- gsub('_',' ',pdx)
        plttitle <- paste0(plttitle, ' ', 'Primary')
        p <- get_proportion_plt(out, plttitle, colorcode)
        png(paste0(save_dir_fg, gsub(' ','_',plttitle),"_",plttype, "_clone.png"),height = 2*250, width=2*250,res = 2*72)
        print(p)
        dev.off()
        
      }
      
      out2 <- table(tmp2$clone_id, tmp2$pdxid)  %>% as.data.frame()
      if(dim(out2)[1]>0){
        colnames(out2) <- c('label','cluster', 'count')
        # print(out)
        # View(clone_color)
        colorcode <- clone_color[as.character(out2$label)]
        plttitle <- gsub('_',' ',pdx)
        plttitle <- paste0(plttitle, ' ', 'Metastasis')
        p <- get_proportion_plt(out2, plttitle, colorcode)
        png(paste0(save_dir_fg, gsub(' ','_',plttitle),"_",plttype, "_clone.png"),height = 2*250, width=2*250,res = 2*72)
        print(p)
        dev.off()
      }  
    }
    
  }
  # in case of SA919, more than 1 passage
  meta_data$passage <- stringr::str_sub(meta_data$sample_id, 6, 7)
  unique(meta_data$passage)
  plttype <- 'passage'
  save_dir_fg <- paste0(save_dir,plttype,'/')
  if (!file.exists(save_dir_fg)){
    dir.create(save_dir_fg)
  }
  for(obs in unique(meta_data$passage)){
    tmp <- meta_data %>%
      dplyr::filter(passage==obs)%>%
      dplyr::select(passage,clone_id)
    
    tmp2 <- meta_data %>%
      dplyr::filter(passage==obs)%>%
      dplyr::select(passage,mainsite)
    print(obs)
    print(dim(tmp))
    if(dim(tmp)[1]>0){
      out <- table(tmp$clone_id, tmp$passage)  %>% as.data.frame()
      colnames(out) <- c('label','cluster', 'count')
      # print(out)
      # View(clone_color)
      colorcode <- clone_color[as.character(out$label)]
      plttitle <- gsub('_',' ',obs)
      p <- get_proportion_plt(out, plttitle, colorcode)
      png(paste0(save_dir_fg, obs,"_",plttype, "_clone.png"),height = 2*250, width=2*250,res = 2*72)
      print(p)
      dev.off()
      
      out2 <- table(tmp2$mainsite, tmp2$passage)  %>% as.data.frame()
      colnames(out2) <- c('label','cluster', 'count')
      print(out2)
      p2 <- get_proportion_plt(out2, plttitle, mainsites_cols)
      png(paste0(save_dir_fg, obs,"_",plttype, "_mainsite.png"),height = 2*250, width=2*250,res = 2*72)
      print(p2)
      dev.off()
    }
    
  }
}