# source_dir <- '/home/htran/Projects/farhia_project/rscript/dlp/visualize_tree/'
# source(paste0(source_dir, "make_cell_copynumber_tree_heatmap.R"))
# source(paste0(source_dir, "tree_viz.R"))
# source(paste0(source_dir, "utils.R"))
# source(paste0(source_dir,'fitclone_edge_list_utils.R'))
# library(network)
# library(ggnet)
suppressMessages(require("network"))
suppressMessages(require("ggnet"))
suppressMessages(require("ape"))
suppressMessages(require("tools"))
suppressMessages(require("phytools"))
suppressMessages(require("igraph"))
suppressMessages(require("RColorBrewer"))
suppressMessages(require("outliers"))
suppressMessages(require("data.table"))
suppressMessages(require("dplyr"))
suppressMessages(require("optparse"))
suppressMessages(require("tidyr"))
suppressMessages(require("tibble"))
suppressMessages(require("ggplot2"))
suppressMessages(require("cowplot"))
suppressMessages(require("ggtree"))
suppressMessages(require("optparse"))


convert_newick2graphml <- function(results_dir){
  save_dir <- paste0(results_dir,'graph_cut/')
  if (!file.exists(save_dir)){
    dir.create(save_dir)
  }
  newick <- paste0(results_dir, 'tree.newick')
  tree <- read.tree(newick)
  edge_list <- tree_2_edge_list(tree)
  g <- read_ltm_tree(edge_list)
  igraph::write.graph(g, paste0(save_dir,'tree.graphml'),format = 'graphml')
  
}
# source_dir <- '/home/htran/Projects/farhia_project/rscript/dlp/fitness/'
# source(paste0(source_dir,'fitclone_edge_list_utils.R'))

plot_trees <- function(outputfile, cellclones, tree_info, results_dir, save_dir, 
                       rm_threshold=-1, datatag='SA1035', save_data=T){
  h <- tree_info$height
  if(rm_threshold==-1){
    qt <- as.integer(quantile(h$dist, probs = .25))
    print(qt)
    sub <- 55
    if((qt-sub) > 0){
      rm_threshold <- qt - sub
    } else{
      rm_threshold <- qt
    }
  }
  print(rm_threshold)
  
  # Using rm_threshold to chop down the loci at the beginning of tree, near root
  # If layout is break-down, please change this threshold to a smaller or larger number 
  rem <- h$id[h$dist > rm_threshold]
  bb <- tree_info$edge_list
  bb <- bb[bb$source %in% rem & bb$target %in% rem, ]
  
  # Find a first locus to re-attach the root
  first_locus <- grep('locus_', rev(rem), value=T)[1]
  bb <- rbind(bb, data.frame(source = 'root', target = first_locus))
  net <- network(as.matrix(bb), directed = T)
  #ggnet2(net, directed = T, size=2)
  v_names <- network.vertex.names(net)
  print(length(v_names))
  
  
  # Get metasample: time point, treatment status,.. for each cell
  # Please change this part based on your metadata in "library_grouping.csv" file
  node_names <- unique(c(bb$source, bb$target))
  length(grep('locus_', node_names, value=T))
  meta_sample <- data.frame(cell_id=node_names,node_number=rep(1:length(node_names),1), 
                            row.names = node_names, stringsAsFactors = F)
  meta_sample <- meta_sample[v_names,]
  print(dim(meta_sample))
  print('Debug')
  # meta_sample <- get_clones(meta_sample, cellclones, results_dir)
  meta_data_ls <- get_meta_data(meta_sample, cellclones, results_dir)
  meta_sample <- meta_data_ls$meta_sample
  
  # sum(meta_sample$cell_id==v_names)
  rownames(meta_sample) <- meta_sample$cell_id
  meta_sample$node_size <- rep(1,length(meta_sample$cell_id))
  locus_ids <- grep('locus_', meta_sample$cell_id, value=T)
  
  # Get non-assigned cells - black color cells in heatmap plot
  none_cells <- meta_sample$cell_id[meta_sample$clone_id=='None']
  # grep('locus_', none_cells, value=T)
  none_cells <- none_cells[!none_cells %in% c(locus_ids, 'root')]
  length(none_cells)
  meta_sample[none_cells,'node_size'] <- 0
  meta_sample[locus_ids,'node_size'] <- 0
  
  
  # Process root
  rid <- grep('root',meta_sample$cell_id)
  meta_sample[rid,'node_size'] <- 2
  # meta_sample[rid,'clone_id'] <- 'Root'
  meta_sample[rid,'clone_color'] <- '#FF0000'  
  meta_sample[rid,'passage_color'] <- '#FF0000'
  meta_sample[rid,'treatment_color'] <- '#FF0000'
  meta_sample[rid,'color_tp'] <- '#FF0000'
  meta_sample[none_cells,'treatment_color'] <- '#000000'
  meta_sample[none_cells,'passage_color'] <- '#000000'
  meta_sample[none_cells,'color_tp'] <- '#000000'
  
  
  # Generate layout for dream tree
  # If layout is break-down, please change rv_threshold to a smaller or larger number 
  filtered_g <- graph_from_edgelist(as.matrix(bb))
  layoutkk <- fitclone_set_dream_tree_layout(net, filtered_g, paste0(datatag,'_dreamtree_'), save_dir)
  
  
  # Plot clone color
  pcl <- plot_tree(net, meta_sample$clone_color, meta_sample$node_size, layoutkk,  
                   datatag, plottype='clone', save_dir)
  
  # Plot timepoint color
  ptp <- plot_tree(net, meta_sample$color_tp, meta_sample$node_size, layoutkk,  
            datatag, plottype='time_point', save_dir)
  
  # Plot passage : abs treatment, untreated, drug vacant
  pps <- plot_tree(net, meta_sample$passage_color, meta_sample$node_size, layoutkk,  
                   datatag, plottype='passage', save_dir)
  
  # color by treatment 
  pts <- plot_tree(net, meta_sample$treatment_color, meta_sample$node_size, layoutkk,  
                   datatag, plottype='treatment', save_dir)
  
  plot_ls <- list(pcl, ptp, pps, pts)
  
  if(save_data){
    saveRDS(meta_data_ls,file = paste0(save_dir,datatag,'_meta_data_ls.rds'))
    saveRDS(meta_sample,file = paste0(save_dir,datatag,'_meta_sample.rds'))
    saveRDS(net,file = paste0(save_dir,datatag,'_net.rds'))
    saveRDS(bb,file = paste0(save_dir,datatag,'_bb_edge_list.rds'))
  }
  plot_metadata(meta_data_ls, save_dir, datatag)
  
  saveRDS(plot_ls, outputfile)
}


# datatag <- 'SA535_Tyler'
# save_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535/SA535_Tyler_v2/dream_tree/'
# meta_sample <- readRDS(paste0(save_dir,datatag,'_meta_sample.rds'))
# net <- readRDS(paste0(save_dir,datatag,'_net.rds'))
# layoutkk <- readRDS(paste0(save_dir,"main_layoutkk",datatag,'_dreamtree_','.rds'))
# outputfile <- paste0(save_dir,'plots_dream_tree.rds')
# meta_data_ls <- readRDS(paste0(save_dir,datatag,'_meta_data_ls.rds'))


read_ltm_tree_from_ls <- function(edge_list) {
  # Find the root
  g <- igraph::graph_from_edgelist(as.matrix(edge_list))
  V(g)$id <- seq(vcount(g))
  return(g)
}
get_tree_info <- function(results_dir, datatag){
  # Load tree
  newick <- paste0(results_dir, 'tree.newick')
  tree <- read.tree(newick)
  edge_list <- tree_2_edge_list(tree)
  g <- read_ltm_tree_from_ls(edge_list)
  # Get height dat
  h <- get_height_dat(g)
  h <- h[order(h$dist, decreasing = T), ]
  print(summary(h$dist))
  # Remove all < rm_threshold to shorten the root
  #rm_threshold <- 70
  # rm_threshold <- 280
  # if (datatag == 'SA532') {
  #   rm_threshold <- 86
  # }
  # rm_threshold <- 20 # nb cells > 3855
  return(list(edge_list=edge_list,height=h))
}

plot_metadata <- function(meta_data_ls=NULL, save_dir='', datatag='SA1035'){
  
  # meta_sample <- readRDS(paste0(save_dir,datatag,'_meta_sample.rds')) 
  if(is.null(meta_data_ls)){
    meta_data_ls <- readRDS(paste0(save_dir,datatag,'_meta_data_ls.rds'))
  }
  
  # treament status
  meta_data = meta_data_ls$treatment_meta
  meta_data <- meta_data[meta_data$treatmentSt!='Locus',]
  meta_data <- meta_data[order(meta_data$treatmentSt),]
  plottype <- 'treatmentSt'
  xstring <- 'treatmentSt'
  ystring <- 'treatment_color'
  plg_ts <- legend_function(meta_data, xstring, ystring, plottype,
                            plottitle="Treament Status", 
                            colorcode=meta_data$treatment_color, 
                            legendtitle="Treament Status", legendlabels=NULL)
  
  p_ts <- ggdraw() + draw_plot(cowplot::get_legend(plg_ts))
  # lg <- cowplot::get_legend(plg_ts)
  png(paste(save_dir,"lg_",plottype,'_',datatag,".png",sep=""),height = 2*300, width=2*200,res = 2*72)
  # print(grid.arrange(grobs = plots[select_grobs(layout)], layout_matrix = layout,
  #                    bottom=" ",right=" "))
  print(p_ts)
  dev.off()
  
  
  ## Passage status
  # rm(meta_data)
  meta_data = meta_data_ls$treatment_meta
  meta_data <- meta_data[meta_data$passage!='Locus',]
  meta_data <- meta_data[order(meta_data$passage),]
  meta_data$passage <- ifelse(meta_data$passage=="abs_treatment",'absolute_treatment',
                              ifelse(meta_data$passage=="drug_vacant",'drug_holiday',meta_data$passage))
  xstring <- 'passage'
  ystring <- 'passage_color'
  plottype <- 'passage'
  
  plg_ps <- legend_function(meta_data, xstring, ystring, plottype,
                            plottitle="Passage Status", 
                            colorcode=unique(meta_data$passage_color), 
                            legendtitle="Passage Status", legendlabels=NULL)
  
  p_ps <- ggdraw() + draw_plot(cowplot::get_legend(plg_ps))
  png(paste(save_dir,"lg_",plottype,'_',datatag,".png",sep=""),height = 2*200, width=2*220,res = 2*72)
  print(p_ps)
  dev.off()
  
  # clone id
  # rm(meta_data)
  meta_data = meta_data_ls$clone_meta
  meta_data <- meta_data[meta_data$clone_id !='None',]
  meta_data <- meta_data[order(meta_data$clone_id),]
  plottype <- 'clone_id'
  xstring <- 'clone_id'
  ystring <- 'color'
  plg_cl <- legend_function(meta_data, xstring, ystring, plottype,
                            plottitle="Clones", 
                            colorcode=meta_data$color, 
                            legendtitle="Clone Id", legendlabels=NULL)
  
  
  p_cl <- ggdraw() + draw_plot(cowplot::get_legend(plg_cl))
  png(paste(save_dir,"lg_",plottype,'_',datatag,".png",sep=""),height = 2*500, width=2*160,res = 2*72)
  # print(grid.arrange(grobs = plots[select_grobs(layout)], layout_matrix = layout,
  #                    bottom=" ",right=" "))
  print(p_cl)
  dev.off()
  
  # timepoint
  # rm(meta_data)
  meta_data = meta_data_ls$timepoint_meta
  meta_data <- meta_data[meta_data$timepoint != 'Locus',]
  meta_data <- meta_data[order(meta_data$timepoint),]
  plottype <- 'timepoint'
  xstring <- 'timepoint'
  ystring <- 'color_tp'
  plg_tp <- legend_function(meta_data, xstring, ystring, plottype,
                            plottitle="Time Point", 
                            colorcode=meta_data$color_tp, 
                            legendtitle="Time Point", legendlabels=NULL)
  
  
  p_tp <- ggdraw() + draw_plot(cowplot::get_legend(plg_tp))
  png(paste(save_dir,"lg_",plottype,'_',datatag,".png",sep=""),height = 2*200, width=2*160,res = 2*72)
  # print(grid.arrange(grobs = plots[select_grobs(layout)], layout_matrix = layout,
  #                    bottom=" ",right=" "))
  print(p_tp)
  dev.off()
  plots <- list(plg_ts, plg_ps, plg_cl, plg_tp)
  saveRDS(plots, file = paste0(save_dir,datatag,'_plots_list.rds'))
  
}
plot_tree <- function(net, color_nodes, node_sizes, layoutkk=NULL,  
                      datatag ='SA919_Sohrab', plottype='clone', 
                      save_dir=NULL, ht=400, wd=500){
  
  # need to get same layout for 2 plots
  # color_nodes <- meta_sample$clone_color
  # node_sizes <- meta_sample$node_size
  if(is.null(layoutkk)){
    p1 <- ggnet2(net, size = node_sizes, mode = "kamadakawai", color = color_nodes) 
  } else{
    # p1 <- ggnet2(net, node.size = node_sizes * 10, mode = layoutkk, color = color_nodes)
    mode <- layoutkk
    net %v% "xx" <- mode[, 1]
    net %v% "yy" <- -mode[, 2]
    mode <- c("xx", "yy")
    p1 <- ggnet2(net,  node.shape = 16, size=1.5,
                 color = color_nodes, mode = mode)
  }
  p1 <- p1 + geom_point(color = 'black', shape = 19, alpha = 1/1000) #, size=3
  
  
  p1 <- p1 + ggtitle(paste0(datatag, ' ',plottype)) + 
    guides(size = FALSE) +
    theme(legend.text=element_blank(), 
          plot.title = element_text(size = 17,hjust = 0.5))
  # if(plottype=='clone'){
  #   cls_df <- get_color_clones(meta_sample)
  # } else if(plottype=='time_point'){
  #   cls_df <- get_color_timepoint(meta_sample)
  # } else if(plottype=='passage'){
  #   cls_df <- get_color_passage(meta_sample)
  # }
  
  # leg_cl <- legend_plot(clone_dic = cls_df, vertical = T,font.size = 3, point_size = 3.5)
  p_1 <- ggdraw() + draw_plot(p1)
  png(paste0(save_dir,"trajectory_",plottype,'_',datatag,".png"),height = 2*ht, width=2*wd,res = 2*72)
  # pdf(paste(save_dir,"cluster_prevalence.pdf",sep=""), height=6, width=13)
  print(p_1)
  dev.off()
  return(p_1)
  # png(paste(save_dir,"lg_",plottype,'_',datatag,".png",sep=""),height = 2*200, width=2*200,res = 2*72)
  # # pdf(paste(save_dir,"cluster_prevalence.pdf",sep=""), height=6, width=13)
  # print(leg_cl)
  # dev.off()
}  

# plot_dream_tree <- function(net, meta_sample, layoutkk=NULL, datatag ='SA919_Sohrab'){
#   
#   # need to get same layout for 2 plots
#   if(is.null(layoutkk)){
#     p1 <- ggnet2(net, size = meta_sample$size, mode = "kamadakawai", color = meta_sample$color) 
#   } else{
#     p1 <- ggnet2(net, node.size = meta_sample$size * 20, mode = layoutkk, color = meta_sample$color)
#   }
#   
#   p12 <- p1
#   if(is.null(layoutkk)){
#     p2 <- ggnet2(net, size = meta_sample$size, mode = "kamadakawai", color = meta_sample$color_group) 
#   } else{
#     p2 <- ggnet2(net, node.size = meta_sample$size * 20, mode = layoutkk, color = meta_sample$color_group) 
#   }
#   
#   p22 <- p2
#   
#   p12 <- p12 + ggtitle(paste0(datatag)) + 
#     guides(size = FALSE) +
#     theme(legend.text=element_blank(), 
#           plot.title = element_text(size = 17,hjust = 0.5))
#   
#   p22 <- p22 + ggtitle(paste0(datatag)) + 
#     guides(size = FALSE) +
#     theme(legend.text=element_blank(), 
#           plot.title = element_text(size = 17,hjust = 0.5))
#   # p2
#   cls_df_cl <- get_color_clones(meta_sample)
#   cls_df_grp <- get_color_timepoint(meta_sample)
#   
#   
#   leg_cl <- legend_plot(clone_dic = cls_df_cl, vertical = T,font.size = 5, point_size = 3.5)
#   leg_grp <- legend_plot(clone_dic = cls_df_grp, vertical = T,font.size = 5, point_size = 3.5)
#   
#   # leg
#   
#   p_1 <- ggdraw() +
#     draw_plot(p12) #+
#   # draw_plot(leg, x = -0.35, y = 0.12, scale = .3)
#   # draw_plot(leg, x = 0.4, y = -0.35, scale = .3)
#   # draw_plot(leg_cl, x = 0.32, y = 0.2, scale = .3)
#   # p_1
#   
#   p_2 <- ggdraw() +
#     draw_plot(p22) #+
#   # draw_plot(leg, x = -0.35, y = 0.12, scale = .3)
#   # draw_plot(leg, x = 0.4, y = -0.35, scale = .3)
#   # draw_plot(leg_grp, x = 0.32, y = 0.2, scale = .3)
#   
#   
#   png(paste(save_dir,"tree_trajectory_clone_",datatag,".png",sep=""),height = 2*700, width=2*450,res = 2*72)
#   # pdf(paste(save_dir,"cluster_prevalence.pdf",sep=""), height=6, width=13)
#   print(p_1)
#   dev.off()
#   
#   png(paste(save_dir,"legend_clone_",datatag,".png",sep=""),height = 2*150, width=2*200,res = 2*72)
#   # pdf(paste(save_dir,"cluster_prevalence.pdf",sep=""), height=6, width=13)
#   print(leg_cl)
#   dev.off()
#   
#   png(paste(save_dir,"tree_trajectory_timepoint_",datatag,".png",sep=""),height = 2*700, width=2*450,res = 2*72)
#   # pdf(paste(save_dir,"cluster_prevalence.pdf",sep=""), height=6, width=13)
#   print(p_2)
#   dev.off()
#   
#   png(paste(save_dir,"legend_timepoint_",datatag,".png",sep=""),height = 2*150, width=2*200,res = 2*72)
#   # pdf(paste(save_dir,"cluster_prevalence.pdf",sep=""), height=6, width=13)
#   print(leg_grp)
#   dev.off()
# }
# 

fitclone_set_dream_tree_layout <- function(net, the_g, datatag, save_dir) {
  # the_g <- g_sans_chains
  out_path <- paste0(save_dir,"main_layoutkk",datatag,'.rds')
  if (file.exists(out_path)) {
    return(readRDS(out_path))
  }
  x <- layout.kamada.kawai(the_g)
  ord <- seq(nrow(x))
  names(ord) <- V(the_g)$name
  x <- x[unname(ord[network.vertex.names(net)]), ]
  saveRDS(x, file = out_path)
  return(x)
}

# fitclone_set_dream_tree_layout <- function(net, the_g) {
#   # the_g <- g_sans_chains
#   x <- layout.kamada.kawai(the_g)
#   ord <- seq(nrow(x))
#   names(ord) <- V(the_g)$name
#   x <- x[unname(ord[network.vertex.names(net)]), ]
#   
#   return(x)
# }

get_clones <- function(meta_sample, cellclones, results_dir){
  # node_names <- res$node_names
  cell_clones <- read.delim(cellclones, stringsAsFactors = FALSE, sep = ",")
  aug_cut <- data_frame_to_list(cell_clones)
  myColors <- get_cluster_colours(length(aug_cut))
  names(myColors) <- names(aug_cut)
  
  # TODO: Since 0 is shared with the heatmap CNV state, use something else for the Un-assigned ones
  myColors <- c(myColors, "None" = "#000000")
  
  meta_sample <- meta_sample %>% left_join(cell_clones, by = "cell_id")
  rownames(meta_sample) <- meta_sample$cell_id
  meta_sample$clone_id <- ifelse(is.na(meta_sample$clone_id),'None',meta_sample$clone_id)
  color_df <- data.frame(clone_id=names(myColors),color=as.character(myColors), stringsAsFactors = FALSE)
  meta_sample <- meta_sample %>% left_join(color_df, by = "clone_id")
  colnames(meta_sample)[which(names(meta_sample) == "color")] <- "clone_color"
  return(meta_sample)
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


get_meta_data <- function(meta_sample, cellclones, results_dir){
  
  cell_clones <- read.delim(cellclones, stringsAsFactors = FALSE, sep = ",")
  aug_cut <- data_frame_to_list(cell_clones)
  myColors <- get_cluster_colours(length(aug_cut))
  names(myColors) <- names(aug_cut)
  
  # TODO: Since 0 is shared with the heatmap CNV state, use something else for the Un-assigned ones
  myColors <- c(myColors, "None" = "#000000")
  # meta_sample$cell_id <- as.character(meta_sample$cell_id)
  meta_sample <- meta_sample %>% left_join(cell_clones, by = "cell_id")
  rownames(meta_sample) <- meta_sample$cell_id
  meta_sample$clone_id <- ifelse(is.na(meta_sample$clone_id),'None',meta_sample$clone_id)
  # meta_sample$clone_id <- as.character(meta_sample$clone_id)
  clone_color_df <- data.frame(clone_id=names(myColors),color=as.character(myColors), stringsAsFactors = FALSE)
  meta_sample <- meta_sample %>% left_join(clone_color_df, by = "clone_id")
  colnames(meta_sample)[which(names(meta_sample) == "color")] <- "clone_color"
  
  grouping_df <- read.csv(paste0(results_dir,'library_groupings.csv'),header=T, check.names=F, row.names = 1, stringsAsFactors=F)
  colnames(grouping_df)[which(names(grouping_df) == "passage")] <- "timepoint"
  colnames(grouping_df)[which(names(grouping_df) == "treatment_st")] <- "treatmentSt"
  colnames(grouping_df)[which(names(grouping_df) == "grouping")] <- "library_id"
  
  meta_sample$library_id <- get_library_id(as.character(meta_sample$cell_id))
  meta_sample$sample_id <- get_sample_id(as.character(meta_sample$cell_id))
  meta_sample <- meta_sample %>% left_join(grouping_df, by = c("library_id","sample_id"))
  # colnames(meta_data)[which(names(meta_data) == "mainsite")] <- "treatmentSt"                                     
  meta_sample$treatmentSt <- ifelse(is.na(meta_sample$treatmentSt),'Locus',meta_sample$treatmentSt)
  meta_sample$timepoint <- ifelse(is.na(meta_sample$timepoint),'Locus',meta_sample$timepoint)
  
  tm_df <- data.frame(treatmentSt=c('UT','UTT','UTTT','UTTTT',
                                    'U','UU','UUU','UUUU','UUUUU',
                                    'UTU','UTTU','UTTTU','Locus'),
                      passage=c(rep('absolute_treatment',4),
                                rep('untreated',5),
                                rep('drug_vacant',3),'Locus'),
                      passage_color=c(rep('#458B00',4),
                                 rep('#CD3333',5),
                                 rep('#008B8B',3),'#0F0F0F'),
                      stringsAsFactors=F)
  
  library(RColorBrewer)
  tm_df$treatment_color <- colorRampPalette(brewer.pal(8, "Set2"))(length(tm_df$treatmentSt))
  # Define the number of colors you want
  lstp <- unique(meta_sample$timepoint)
  # lstp <- lstp[!lstp %in% 'Locus'] 
  timepoint_color <- colorRampPalette(brewer.pal(8, "Set2"))(length(lstp))
  
  color_tp_df <- data.frame(timepoint=lstp, color_tp=as.character(timepoint_color), stringsAsFactors = FALSE)
  
  meta_sample <- meta_sample %>% left_join(color_tp_df, by = "timepoint")
  # meta_sample$treatmentSt <- as.character(meta_sample$treatmentSt)
  meta_sample <- meta_sample %>% left_join(tm_df, by = "treatmentSt")
  return(list(meta_sample=meta_sample,treatment_meta=tm_df,timepoint_meta=color_tp_df, clone_meta=clone_color_df))
}

get_timepoint <- function(meta_sample, results_dir){
  
  grouping_df <- read.csv(paste0(results_dir,'library_groupings.csv'),header=T,check.names=F, row.names = 1)
  grouping_df <- grouping_df[,c('grouping','mainsite')]
  lib_ids <- strsplit(as.character(meta_sample$cell_id),"-")
  length(lib_ids)
  nbcores <- 5
  ps <- c()
  passage_id <- parallel::mclapply(lib_ids, function(f) {
    ps <- c(ps, as.character(f[2]))
  }, mc.cores = nbcores)
  length(passage_id)
  # passage_id[1:5]
  
  meta_sample$library_labels <- as.character(passage_id)
  meta_sample <- meta_sample %>% left_join(grouping_df, by = c("library_labels"="grouping"))
  # colnames(meta_data)[which(names(meta_data) == "mainsite")] <- "treatmentSt"                                     
  meta_sample$mainsite <- ifelse(is.na(meta_sample$mainsite),'Locus',meta_sample$mainsite)
  
  mt_colors <- structure(
    c(
      "#FF3232", "#6495ED", "#A0A0A0"
    ),
    names=c("Metastasis", "Primary", "Locus")
  )
  color_group <- data.frame(mainsite=names(mt_colors),color_group=as.character(mt_colors), stringsAsFactors = FALSE)
  meta_sample <- meta_sample %>% left_join(color_group, by = "mainsite")
  return(meta_sample)
}


# TODO: set the sorting for clones...
legend_plot <- function(clone_dic, vertical = TRUE, font.size = 3, point_size = 3, main = NULL) {
  if (vertical) {
    ll.dat <- data.frame(x = c(rep(1, length(clone_dic$letters))), y = rev(seq_along(clone_dic$letters)), 
                         label = clone_dic$letters, color=clone_dic$color, stringsAsFactors = F)
  } else {
    ll.dat <- data.frame(x = seq_along(clone_dic$letters), y = c(rep(1, length(clone_dic$letters))), label = clone_dic$letters, stringsAsFactors = F)
  }
  
  # Better results with spaces
  ll.dat$label <- paste0('  ', ll.dat$label)
  clone_dic$letters <- paste0('  ', clone_dic$letters)
  
  nudge_x = 0; nudge_y = 0
  if (vertical) {
    #nudge_y = .5
    nudge_x = 0
  } else {
    nudge_x = .5
  }
  
  myColors <- clone_dic$color
  names(myColors) <- clone_dic$letters
  pg <- ggplot(data = ll.dat, mapping = aes(x = x, y = y, color = label)) + 
    geom_point(aes(color=label)) + 
    geom_text(aes(label = label), colour = 'black', nudge_x = nudge_x, nudge_y = nudge_y, fontface = "bold", size = font.size, hjust = 0) + 
    fitclone_get_theme_no_grid() + 
    fitclone_get_theme_no_axis() + 
    fitclone_get_theme_no_legend() + 
    scale_color_manual(values = as.character(myColors)) + 
    theme(panel.background = element_blank(),
          panel.border=element_blank(),
          plot.margin = margin(0, 0, 0, 0, "cm"))
  if (point_size > 0) {
    # ll.dat.point <- ll.dat[-c(1), ]
    ll.dat.point <- ll.dat
    pg <- pg + geom_point(data = ll.dat.point, size = point_size, shape = 15)
    
  }
  
  if (!is.null(main)) {
    pg <- pg + ggtitle(label = '', subtitle = main)
  }
  
  if (vertical) {
    #pg <- pg + xlim(c(1-.1, 1+.1))
    pg <- pg + xlim(c(1-nudge_x, 1+nudge_x*2))
  } else {
    pg <- pg + ylim(c(1-.1, 1+.1))
    # pg <- pg + ylim(c(0.3-.1, 0.3+.1))
  }
  
  return(pg)
}

get_unique_values <- function(meta_sample){
  tab <- with(meta_sample, table(letters, color))
  rs <- prop.table(tab, margin = NULL)
  cls_df <- as.data.frame(rs)
  cls_df <- cls_df[cls_df$Freq>0,]
  cls_df <- cls_df[order(cls_df$letters),]
  rownames(cls_df) <- rep(1:dim(cls_df)[1])
  return(cls_df)
}
get_color_clones <- function(meta_sample){
  tab <- with(meta_sample, table(clone_id, color))
  rs <- prop.table(tab, margin = NULL)
  cls_df <- as.data.frame(rs)
  cls_df <- cls_df[cls_df$Freq>0,]
  colnames(cls_df)[which(colnames(cls_df) == "clone_id")] <- "letters"
  cls_df <- cls_df[order(cls_df$letters),]
  rownames(cls_df) <- rep(1:dim(cls_df)[1])
  return(cls_df)
}

get_color_timepoint <- function(meta_sample){
  tab <- with(meta_sample, table(mainsite, color_group))
  rs <- prop.table(tab, margin = NULL)
  cls_df <- as.data.frame(rs)
  cls_df <- cls_df[cls_df$Freq>0,]
  colnames(cls_df)[which(colnames(cls_df) == "mainsite")] <- "letters"
  cls_df <- cls_df[order(cls_df$letters),]
  rownames(cls_df) <- rep(1:dim(cls_df)[1])
  return(cls_df)
}


tree_2_edge_list_path <- function(newick, results_dir, save_dir){
  # Load tree
  tree <- read.tree(newick)
  edge_list <- tree_2_edge_list(tree)
  edge_list_path <- paste0(save_dir,'edge_list.csv')
  write.csv(edge_list, file = edge_list_path, quote=F, row.names = F)
  if(file.exists(edge_list_path)){
    return(edge_list_path)
  } else{
    return(NULL)
  }
  
}
get_clones_subgraph <- function(results_dir, cellclones){
  # Load tree
  newick <- paste0(results_dir, 'tree.newick')
  tree <- read.tree(newick)
  edge_list <- tree_2_edge_list(tree)
  g <- read_ltm_tree(edge_list)
  
  all_cells <- gsub('cell_', '', grep('cell', tree$tip.label, value = T))
  length(all_cells)
  cell_clones <- read.delim(cellclones, stringsAsFactors = FALSE, sep = ",")
  none_cells <- all_cells[!all_cells %in% cell_clones$cell_id]
  length(none_cells)
  # none_cells <- none_cells[!none_cells %in% grep('locus_', none_cells, value=T)]
  length(none_cells)
  cg_dir <- paste0(results_dir,'corrupt_grow/')
  if (!file.exists(cg_dir)){
    dir.create(cg_dir)
  }
  
  edge_list <- edge_list[!edge_list$target %in% none_cells, ]
  # edge_list_filtered_path <- paste0(cg_dir,'edges_filtered_SA919.csv')
  # write.csv(edge_list, file = edge_list_filtered_path, quote=F, row.names = F)
  
  # test corrupt grow
  res_cg <- convert_edge_list_to_ape(as.matrix(edge_list))
  tree_cg <- res_cg$tree
  tree_cg$tip.label <- paste0('cell_',tree_cg$tip.label)
  write.tree(tree_cg, file = paste0(cg_dir,'tree.newick.clones'), append = FALSE, digits = 10, 
             tree.names = FALSE)
  
  # newick_clones <- paste0(cg_dir, 'tree.newick.clones')
  # tree_clones <- read.tree(newick_clones)
  # p <- ggtree(tree_clones)
  # p
}


legend_function <- function(meta_data, xstring, ystring, plottype,
                            plottitle="Treament label", 
                            colorcode="blue", legendtitle="Treament label", legendlabels=NULL) {
  p <- ggplot(meta_data, aes_string(x=xstring, y=ystring, group=plottype)) +
    # geom_line(aes_string(color=plottype)) +  #,color=colorcode
  geom_point(aes_string(color=plottype),size=8) +
  scale_color_manual(values = colorcode) + labs(colour=legendtitle)
  
  p <- p + theme(
    legend.title = element_text(size=17), 
    legend.text = element_text(size=16), 
    plot.title = element_text(color="black", size=20, hjust = 0.5)
  )
  return(p)
}


plot_tree_timepoint <- function(meta_data, group_use,
                                datatag, plttype='clone', save_dir){
  meta_data$timepoint
  group_use <- c('UT','UTT','UTTT','UTTTT','U','UU','UUU','UUUU','UUUUU')
  for(g in group_use){
    meta_data_tmp <- meta_data
    meta_data_tmp[meta_data_tmp$treatmentSt==g,'node_size'] <- 2
    meta_data_tmp[meta_data_tmp$treatmentSt!=g,'node_size'] <- 0
    meta_data_tmp[meta_data_tmp$treatmentSt!=g,'clone_color'] <- '#000000'
    rid <- grep('root',meta_data_tmp$cell_id)
    
    meta_sample[rid,'node_size'] <- 2
    # meta_sample[rid,'clone_id'] <- 'Root'
    meta_sample[rid,'clone_color'] <- '#FF0000'  
    rv_cells <- meta_data_tmp[meta_data_tmp$treatmentSt!=g,'cell_id']
    locus_ids <- grep('locus_', rv_cells, value=T)
    rv_cells <- rv_cells[!rv_cells %in% c(locus_ids, 'root')]
    v_names <- network.vertex.names(net)
    length(v_names)
    vertices <- data.frame(cell_id=v_names,cell_idx=rep(1:length(v_names),1), row.names = rep(1:length(v_names),1))
    rv_idx <- vertices[vertices$cell_id %in% rv_cells,'cell_idx']
    
    
    bb1 <- bb
    bb1 <- bb1[!bb1$target %in% rv_cells,]
    net <- network(as.matrix(bb1), directed = T)
    filtered_g <- graph_from_edgelist(as.matrix(bb1))
    layoutkk <- fitclone_set_dream_tree_layout(net, filtered_g)
    meta_data_tmp <- meta_data_tmp[network.vertex.names(net),]
    
    pcl <- plot_tree(net, meta_data_tmp$clone_color, meta_data_tmp$node_size, layoutkk,  
                     paste0(datatag,'_',g), plottype=paste0(datatag,'_',g), save_dir)
    
    p1 <- ggnet2(net, node.size = meta_data_tmp$node_size * 10, mode = layoutkk, color = meta_data_tmp$clone_color)
    p1
    
  }
  
}
