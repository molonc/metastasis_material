# source_dir <- '/home/htran/Projects/farhia_project/rscript/dlp/fitness/'
# source(paste0(source_dir,'fitclone_corrupt_utils.R'))
# source(paste0(source_dir,'fitclone_tree_utils.R'))
# source(paste0(source_dir,'fitclone_autocut.R'))
# source(paste0(source_dir,'fitclone_tree.R'))
# source(paste0(source_dir,'fitclone_heatmaps.R'))
# source(paste0(source_dir,'fitclone_heatmap_utils.R'))
# source(paste0(source_dir,'fitclone_s.R'))
# source(paste0(source_dir,'fitclone_trajectories.R'))
# source(paste0(source_dir,'fitclone_private_cnv.R'))
# source(paste0(source_dir,'fitclone_paths.R'))
# source(paste0(source_dir,'fitclone_cell_assignment_utils.R'))
# source(paste0(source_dir,'fitclone_genotype_cut.R'))
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
suppressMessages(require("tidyverse"))
suppressMessages(require("stringr"))



get_loci <- function(str_arry) grep('locus_', str_arry, value = T)

g.get.loci <- function(g=NULL, edge_list_path=NULL) {
  if (is.null(g)) g <- read_ltm_tree(edge_list_path = edge_list_path)
  get_loci(V(g)$name)
}



get_cell_count_per_edge <- function(g) {
  edges <- c(g.get.loci(g), 'root')
  cell_counts <- c()
  for (edge in edges) {
    #print(edge)
    # edge <- edges[[50]]
    immediate_children <- igraph::neighborhood(graph = g, nodes=edge, mode='out')[[1]]$name[-c(1)]
    tmp_count <- sum(!grepl('locus', immediate_children))
    cell_counts <- c(cell_counts, tmp_count)
  }
  
  data.frame(edges = edges, cell_counts = cell_counts, stringsAsFactors = F)
}

# From original newick tree, get loci, cells and 
# nb of cells at each loci (nb neighbors cells related to this loci)
get_fixed_ebola_tree <- function(exp_path, edge_list_path, datatag="SA1035") {
  out_path <- file.path(exp_path, sprintf('ebola_tree_fixed_%s.rds',datatag))
  print(out_path)
  if (file.exists(out_path)) {
    return(readRDS(out_path))
  }
  
  g <- read_ltm_tree(edge_list_path)
  # all_cells <- g.get.cells(g)
  # keep_cells <- c('root')
  # # Don't remove root
  # g <- igraph::delete.vertices(g, setdiff(all_cells, keep_cells))
  dat <- get_cell_count_per_edge(g)
  res <- list(dat = dat, g = g)
  saveRDS(res, out_path)
  return(res)
}


# Get main tree: remove cells leaves, get key loci, count nb cells at each loci
get_main_tree <- function(save_dir, save_fn, edge_list_path, datatag){
  
  # Get tree and cell counts per loci
  shared_res <- get_fixed_ebola_tree(save_dir, edge_list_path, datatag)
  g <- shared_res$g
  # print(summary(shared_res$dat$cell_counts))
  dat <- shared_res$dat
  # ggtree(g)
  
  
  # Remove cells leaves from tree, keep only loci
  cells_sans_root <- g.get.cells(g)[g.get.cells(g) != 'root']
  g_sans_cells <- igraph::delete.vertices(g, cells_sans_root)
  vcount(g_sans_cells)  # update indices after removing cells
  
  # Remove the chains
  # Remove loci, keep only important loci, shorten tree trajectory
  g_sans_chains <- collapse_chains(the_g = g_sans_cells, cell_count_dat = dat, 
                                   cell_count_treshold = 4)
  dat <- dat[dat$edges %in% V(g_sans_chains)$name, ]
  print(summary(dat$cell_counts))
  
  
  res <- list(fixed_tree=shared_res, g_sans_chains = g_sans_chains, dat = dat)
  saveRDS(object = res, file = paste0(save_dir, save_fn))
  
  return(res)
}
get_cells_by_g_mainsite <- function(results_dir, g, datatag, mainsite=NULL) {

  all_cells <- g.get.cells(g)
  all_cells <- all_cells[!all_cells %in% 'root']
  length(all_cells)
  all_libids <- get_library_id(all_cells)
  sids <- get_sample_id(all_cells)
  # grouping_df <- get_timepoint(results_dir)
  
  grouping_df <- read.csv(paste0(results_dir,'library_groupings.csv'),header=T, check.names=F, 
                          row.names = 1, stringsAsFactors=F)
  colnames(grouping_df)[which(names(grouping_df) == "grouping")] <- "library_id"
  
  
  cell.dat <- data.frame(sample_id=sids,library_id = all_libids, cell_id = all_cells, stringsAsFactors = F)
 
  cell.dat <- dplyr::left_join(cell.dat, grouping_df, by=c('library_id','sample_id'))
  # rownames(cell.dat) <- cell.dat$cell_id
  # t <- cell.dat[which(cell.dat$timepoint==timepoint & cell.dat$treatmentSt=='UUU'),]
  # sum(cell.dat$timepoint==timepoint & cell.dat$treatmentSt=='UUU')
  # View(t)
  if (is.null(mainsite)) {
    # cells_use <- all_cells[all_libids %in% cell.dat$library_id]
    cells_use <- all_cells
  } else {
    # cells_use <- all_cells[all_libids %in% cell.dat$library_id[cell.dat$timepoint %in% timepoint]]
    cells_use <- cell.dat[cell.dat$mainsite==mainsite,'cell_id']
    print(length(cells_use))
  }
  
  # View(get_library_id(cells_use))
  print(paste0('Main site: ',mainsite, ' library id is: '))
  print(unique(get_library_id(cells_use)))
  return(cells_use)
}


get_cells_by_g_pdx <- function(results_dir, g, datatag, pdxid=NULL) {
  
  all_cells <- g.get.cells(g)
  all_cells <- all_cells[!all_cells %in% 'root']
  length(all_cells)
  all_libids <- get_library_id(all_cells)
  # grouping_df <- get_timepoint(results_dir)
  
  grouping_df <- read.csv(paste0(results_dir,'library_groupings.csv'),header=T, check.names=F, 
                          row.names = 1, stringsAsFactors=F)
  colnames(grouping_df)[which(names(grouping_df) == "grouping")] <- "library_id"
  
  
  cell.dat <- data.frame(library_id = all_libids, cell_id = all_cells, stringsAsFactors = F)
  
  cell.dat <- dplyr::left_join(cell.dat, grouping_df, by=c('library_id'))
  # rownames(cell.dat) <- cell.dat$cell_id
  # t <- cell.dat[which(cell.dat$timepoint==timepoint & cell.dat$treatmentSt=='UUU'),]
  # sum(cell.dat$timepoint==timepoint & cell.dat$treatmentSt=='UUU')
  # View(t)
  if (is.null(pdxid)) {
    # cells_use <- all_cells[all_libids %in% cell.dat$library_id]
    cells_use <- all_cells
  } else {
    # cells_use <- all_cells[all_libids %in% cell.dat$library_id[cell.dat$timepoint %in% timepoint]]
    cells_use <- cell.dat[cell.dat$pdxid==pdxid,'cell_id']
    print(length(cells_use))
  }
  
  # View(get_library_id(cells_use))
  print(paste0('Pdx id: ',pdxid, ' library id is: '))
  print(unique(get_library_id(cells_use)))
  return(cells_use)
}

get_cells_by_g_origin <- function(results_dir, g, datatag, origin=NULL) {
  
  all_cells <- g.get.cells(g)
  all_cells <- all_cells[!all_cells %in% 'root']
  length(all_cells)
  all_libids <- get_library_id(all_cells)
  # grouping_df <- get_timepoint(results_dir)
  
  grouping_df <- read.csv(paste0(results_dir,'library_groupings.csv'),header=T, check.names=F, 
                          row.names = 1, stringsAsFactors=F)
  colnames(grouping_df)[which(names(grouping_df) == "grouping")] <- "library_id"
  
  
  cell.dat <- data.frame(library_id = all_libids, cell_id = all_cells, stringsAsFactors = F)
  
  cell.dat <- dplyr::left_join(cell.dat, grouping_df, by=c('library_id'))
  # rownames(cell.dat) <- cell.dat$cell_id
  # t <- cell.dat[which(cell.dat$timepoint==timepoint & cell.dat$treatmentSt=='UUU'),]
  # sum(cell.dat$timepoint==timepoint & cell.dat$treatmentSt=='UUU')
  # View(t)
  if (is.null(origin)) {
    # cells_use <- all_cells[all_libids %in% cell.dat$library_id]
    cells_use <- all_cells
  } else {
    # cells_use <- all_cells[all_libids %in% cell.dat$library_id[cell.dat$timepoint %in% timepoint]]
    cells_use <- cell.dat[cell.dat$origin==origin,'cell_id']
    print(length(cells_use))
  }
  
  # View(get_library_id(cells_use))
  print(paste0('Origin id: ',origin, ' library id is: '))
  print(unique(get_library_id(cells_use)))
  return(cells_use)
}

fitclone_get_ebola_tree_conditional <- function(res, results_dir, 
                                                save_dir, datatag, 
                                                obs=NULL, labeltype=1){
  # Compute a secondary cell_coount dat for the filtered down graph
  g.temp <- res$fixed_tree$g
  g <- res$fixed_tree$g
  
  if(labeltype==1 & !is.null(obs)){ # mainsite
    keep_cells <- c(get_cells_by_g_mainsite(results_dir, g, datatag, obs), 'root')
  } else if(labeltype==2 & !is.null(obs)){  # pdx
    keep_cells <- c(get_cells_by_g_pdx(results_dir, g, datatag, obs), 'root')
  } else if(labeltype==3 & !is.null(obs)){  # origin
    keep_cells <- c(get_cells_by_g_origin(results_dir, g, datatag, obs), 'root')
  } else{
    keep_cells <- g.get.cells(g.temp)
  }
  
  g.temp <- igraph::delete.vertices(g.temp, setdiff(g.get.cells(g.temp), keep_cells))
  dat_cond <- get_cell_count_per_edge(g.temp)
  # head(dat_cond)
  
  gtree_cond <- list(g_sans_chains = res$g_sans_chains, dat = dat_cond)
  # saveRDS(object = gtree_cond, file = paste0(save_dir,'gtree_cond_',timepoint,'_',datatag,'.rds'))
  return(gtree_cond)
}

fitclone_get_ebola_tree_conditional_pdx <- function(res, results_dir, 
                                                save_dir, datatag, 
                                                pdxid='2164'){
  # Compute a secondary cell_coount dat for the filtered down graph
  g.temp <- res$fixed_tree$g
  g <- res$fixed_tree$g
  keep_cells <- c(get_cells_by_g_pdx(results_dir, g, datatag, pdxid), 'root')
  g.temp <- igraph::delete.vertices(g.temp, setdiff(g.get.cells(g.temp), keep_cells))
  dat_cond <- get_cell_count_per_edge(g.temp)
  # head(dat_cond)
  
  gtree_cond <- list(g_sans_chains = res$g_sans_chains, dat = dat_cond)
  # saveRDS(object = gtree_cond, file = paste0(save_dir,'gtree_cond_',timepoint,'_',datatag,'.rds'))
  return(gtree_cond)
}


# fitclone_get_ebola_tree_conditional <- function(res, results_dir, 
#                                                     save_dir, datatag, 
#                                                     obs=NULL, labeltype=1){
#   # Compute a secondary cell_coount dat for the filtered down graph
#   g.temp <- res$fixed_tree$g
#   g <- res$fixed_tree$g
#   if(labeltype)
#   keep_cells <- c(get_cells_by_g_origin(results_dir, g, datatag, origin), 'root')
#   g.temp <- igraph::delete.vertices(g.temp, setdiff(g.get.cells(g.temp), keep_cells))
#   dat_cond <- get_cell_count_per_edge(g.temp)
#   # head(dat_cond)
#   
#   gtree_cond <- list(g_sans_chains = res$g_sans_chains, dat = dat_cond)
#   # saveRDS(object = gtree_cond, file = paste0(save_dir,'gtree_cond_',timepoint,'_',datatag,'.rds'))
#   return(gtree_cond)
# }


# results_dir <- '/home/htran/storage/datasets/metastasis_results/SA535_new_encoding/SA535_wholedata_v2/'
# outputfile <- paste0(results_dir,'tree_viz_dream/','tree_viz_dream.png')
# cellclones <- paste0(results_dir,'cell_clones.csv')
# datatag <- 'SA535'

# results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_wholedata/'
# outputfile <- paste0(results_dir,'tree_viz_dream/','tree_viz_dream.png')
# cellclones <- paste0(results_dir,'cell_clones.csv')
# datatag <- 'SA919'
plot_ebola_tree_condition_Shaocheng <- function(datatag, results_dir, 
                                                outputfile, cellclones=NULL){
  
  results_dir <- paste0(results_dir,'/')
  newick <- paste0(results_dir, 'tree.newick')
  save_dir <- paste0(results_dir,'tree_viz_dream/')
  # save_dir <- paste0(results_dir,'tree_viz_dream_testing/')
  
  if (!file.exists(save_dir)){
    dir.create(save_dir)
  }
  if (!file.exists(dirname(outputfile))){
    dir.create(dirname(outputfile))
  }
  if(is.null(cellclones)){
    cellclones <- paste0(results_dir, 'cell_clones.csv')
  }
  print("Get clone metadata")
  cell_clones <- read.delim(cellclones, stringsAsFactors = FALSE, sep = ",", check.names = F)
  # cell_clones <- cell_clones[cell_clones$clone_id!='None',]
  clone_meta <- get_color_clone(cell_clones)
  plot_clone_colour(clone_meta, save_dir, datatag)
  print(clone_meta)
  
  # Get summary tree, collapse cells in 1 clone to 1 node 
  cloneColours <- clone_meta$color
  names(cloneColours) <- clone_meta$clone_id
  p1 <- fast_get_summary_tree(cell_clones = cell_clones, 
                              newick = newick, 
                              cloneColours = cloneColours, 
                              edge_list = NULL, 
                              point_size_1 = 16, point_size_2 = 14)
  saveRDS(p1, paste0(save_dir,'summary_tree.rds'))
  png(paste0(save_dir,'summary_tree.png'), height = 2*400, width=2*800,res = 2*72)
  # pdf(paste(save_dir,"cluster_prevalence.pdf",sep=""), height=6, width=13)
  print(p1)
  dev.off()
  
  
}  


plot_ebola_tree_condition <- function(datatag, results_dir, outputfile, cellclones=NULL){
  
  results_dir <- paste0(results_dir,'/')
  newick <- paste0(results_dir, 'tree.newick')
  save_dir <- paste0(results_dir,'tree_viz_dream/')
  if (!file.exists(save_dir)){
    dir.create(save_dir)
  }
  if (!file.exists(dirname(outputfile))){
    dir.create(dirname(outputfile))
  }
  if(is.null(cellclones)){
    cellclones <- paste0(results_dir, 'cell_clones.csv')
  }
  print("Get clone metadata")
  cell_clones <- read.delim(cellclones, stringsAsFactors = FALSE, sep = ",", check.names = F)
  # excluded_clones <- c('Un','unassigned','None')
  # cell_clones <- cell_clones[!cell_clones$clone_id %in% excluded_clones,]
  print(dim(cell_clones))
  
  # cell_clones <- cell_clones[cell_clones$clone_id!='None',]
  # clone_meta <- get_color_clone(cell_clones)
  clone_meta <- get_color_clone_v2(unique(cell_clones$clone_id))
  # plot_clone_colour(clone_meta, save_dir, datatag)
  print(clone_meta)
  
  # Get summary tree, collapse cells in 1 clone to 1 node 
  cloneColours <- clone_meta$color
  names(cloneColours) <- clone_meta$clone_id
  # p1 <- fast_get_summary_tree(cell_clones = cell_clones, 
  #                             newick = newick, 
  #                             cloneColours = cloneColours, 
  #                             edge_list = NULL, 
  #                             point_size_1 = 16, point_size_2 = 14)
  # saveRDS(p1, paste0(save_dir,'summary_tree.rds'))
  # png(paste0(save_dir,'summary_tree.png'), height = 2*400, width=2*800,res = 2*72)
  # # pdf(paste(save_dir,"cluster_prevalence.pdf",sep=""), height=6, width=13)
  # print(p1)
  # dev.off()
  
  
  print("Get main tree trajectory")
  
  edge_list_path <- tree_2_edge_list_path(newick, results_dir, save_dir)
  save_fn <- 'g_san_chain_ls.rds'
  res <- get_main_tree(save_dir, save_fn, edge_list_path, datatag)
  # dat <- res$dat
  
  
  print("Constructing network ")
  # load layout
  edge_list <- as.matrix(igraph::as_edgelist(res$g_sans_chains))
  net <- network::network(edge_list, directed = T)
  net %e% 'names' <- 1:nrow(edge_list)
  
  print("Generating main layout...")
  layoutkk <- fitclone_set_dream_tree_layout(net, res$g_sans_chains, datatag, save_dir)
  v_names <- network.vertex.names(net)
  
  print("Get color cluster assignment for each loci based on clustering")
  locus_group <- get_color_cluster_assignment(results_dir, save_dir, datatag, cell_clones, clone_meta)
  # summary(locus_group$group)
  dim(locus_group)
  locus_group <- locus_group[v_names,]
  
  # summary(as.factor(locus_group$color))
  # unassigned_loci <- locus_group[locus_group$group=='0' && locus_group$label!='root','label']
  # 
  # unassigned_loci <- unassigned_loci[unassigned_loci %in% v_names]
  # length(unassigned_loci)
  # if(length(unassigned_loci)>0){
  #   locus_group[unassigned_loci,'color'] <- "#000000"  
  # }
  grouping_df <- read.csv(paste0(results_dir,'library_groupings.csv'),header=T, check.names=F, stringsAsFactors=F)
  colnames(grouping_df)[which(names(grouping_df) == "grouping")] <- "library_id"
  for(mt in unique(grouping_df$mainsite)){
    # net_tmp = net
    print(mt)
    plot_tree_condition(clone_meta, net, layoutkk, v_names, locus_group, 
                        res, results_dir, save_dir, datatag, mt, 1,
                        ht=350, wd=470)
  }
  
  # for(pdx in unique(grouping_df$pdxid)){
  #   # net_tmp = net
  #   print(pdx)
  #   plot_tree_condition(clone_meta, net, layoutkk, v_names, locus_group, 
  #                       res, results_dir, save_dir, datatag, pdx, 2, 
  #                       ht=350, wd=470)
  # }
  # for(origin in unique(grouping_df$origin)){
  #   # net_tmp = net
  #   print(origin)
  #   plot_tree_condition(clone_meta, net, layoutkk, v_names, locus_group, 
  #                           res, results_dir, save_dir, datatag, origin, 3, 
  #                           ht=350, wd=470)
  # }
  
  # For each tree, keep same main layout and list of loci
  # Update nb cells counted at each loci based on time points
  # Each loci represented by a cirle and size based on nb cells count
  # print("Generate ebola tree at each time point")
  # grouping_df <- get_timepoint(results_dir)
  # for(timepoint in unique(grouping_df$timepoint)){
  #   # timepoint = "X6"
  #   # "X4" "X5" "X6" "X7" "X8"
  #   print(timepoint)
  #   ts_ls <- unique(grouping_df[grouping_df$timepoint==timepoint,'treatmentSt'])
  #   for(treatmentSt in ts_ls){
  #     print(treatmentSt)
  #     net_tmp = net
  #     plot_tree_condition(clone_meta, net_tmp, layoutkk, v_names, locus_group, 
  #                        res, results_dir, save_dir, datatag, 
  #                        timepoint, treatmentSt,
  #                        ht=350, wd=470)
  #   }
  #   
  # }
  
  # # Plot whole tree 
  # p <- plot_whole_tree_clones(clone_meta, net, layoutkk, v_names, locus_group, 
  #                     res, results_dir, save_dir, datatag, 
  #                     ht=350, wd=470)
  # 
  # png(outputfile, height = 2*400, width=2*470,res = 2*72)
  # # pdf(paste(save_dir,"cluster_prevalence.pdf",sep=""), height=6, width=13)
  # print(p)
  # dev.off()
}

# plot_whole_tree <- function(layoutkk, v_names, locus_group, 
#                                 dat, results_dir, save_dir, datatag){
#   rownames(dat) <- dat$edges
#   "root" %in% rownames(dat)
#   dat["root",'cell_counts'] <- 200
#   dat <- dat[v_names,]
#   locus_group_tmp <- locus_group[v_names,]
#   size_node <- dat$cell_counts
#   length(size_node)
#   length(locus_group_tmp$color)
#   net_tmp <- net
#   pcl <- plot_tree(net_tmp, locus_group_tmp$color, size_node, layoutkk,
#                    datatag, plottype=paste0('clone_whole'), save_dir)
#   
#   
# }  

# plot_tree_condition(clone_meta, net_tmp, layoutkk, v_names, locus_group, 
#                     res, results_dir, save_dir, datatag, mt, 1,
#                     ht=350, wd=470)
plot_tree_condition <- function(clone_meta, net_tmp, layoutkk, v_names, locus_group, 
                                res, results_dir, save_dir, datatag, 
                                obs=NULL, labeltype=1,
                                ht=350, wd=470){
  plottype <- 'clone'
  if(labeltype==1 & !is.null(obs)){  # mainsite label
    plottype=paste0('main site ', obs)
    print(paste0('Plot tree trajectory at the condition: ',plottype))
    gtree_cond <- fitclone_get_ebola_tree_conditional(res, results_dir, save_dir, datatag, obs, 1)
    dat <- gtree_cond$dat
  } else if(labeltype==2 & !is.null(obs)){   # pdx label
    plottype=paste0('pdx ', obs)
    print(paste0('Plot tree trajectory at the condition: ',plottype))
    gtree_cond <- fitclone_get_ebola_tree_conditional(res, results_dir, save_dir, datatag, obs, 2)
    dat <- gtree_cond$dat
  } else if(labeltype==3 & !is.null(obs)){   # pdx label
    plottype=paste0('origin ', obs)
    print(paste0('Plot tree trajectory at the condition: ',plottype))
    gtree_cond <- fitclone_get_ebola_tree_conditional(res, results_dir, save_dir, datatag, obs, 3)
    dat <- gtree_cond$dat
  } else{
    print('Plot whole tree trajectory')
    plottype <- 'clone_wholedata'
    dat <- res$dat
  }
  
  rownames(dat) <- dat$edges
  # print("Contain root?")
  # print("root" %in% rownames(dat))
  dat <- dat[v_names,]
  # if(length(unassigned_loci)>0){
  #   dat[unassigned_loci,'cell_counts'] <- 0
  # }
  dat["root",'cell_counts'] <- 20
  # 
  sum(v_names %in% dat$edges) == length(v_names)
  
  # size_node <- dat$cell_counts
  # size_node <- ceiling(log(size_node))
  # size_node <- ifelse(!is.finite(size_node),0,size_node)
  
  locus_group_tmp <- locus_group[v_names,]
  
  # If at this time point, no cells present in this branch, change color of loci to black
  rv_loci <- dat[dat$cell_counts==0,'edges']
  locus_group_tmp[as.character(rv_loci),'color'] <- "#000000"
  print(summary(as.factor(locus_group_tmp$color)))
  print("DEBUG 1")
  # Node size by dat$cell_count
  ord <- ceiling(log(dat$cell_counts))
  ord[!is.finite(ord)] <- 0.0
  names(ord) <- dat$edges
  # ord[which(names(ord) == 'root')] <- as.character(max(ord))
  ord <- unname(ord[v_names])
  net_tmp %v% 'direct_cells_log' <- ord
  point_size_scale <- 0.8
  size_vals <- as.numeric(unique(net_tmp %v% 'direct_cells_log'))/point_size_scale
  # Exponentiate the names
  names(size_vals) <- sprintf('%.0f', as.numeric(unique(net_tmp %v% 'direct_cells_log'))**exp(1))
  net_tmp %v% 'direct_cells_log' <- sprintf('%.0f', as.numeric(net_tmp %v% 'direct_cells_log')**exp(1))
  
  size_vals <- size_vals[order(as.numeric(size_vals))]
  # Layout
  mode <- layoutkk
  net_tmp %v% "xx" <- -mode[, 1]
  net_tmp %v% "yy" <- -mode[, 2]
  mode <- c("xx", "yy")
  
  cols <- locus_group_tmp$color
  names(cols) <- locus_group_tmp$label
  net_tmp %v% 'color' <- cols
  
  # cloneColors <- clone_meta$color
  # names(cloneColors) <- clone_meta$clone_id
  #cloneColors <- c(cloneColors,  'root' = 'red', 'Un' = 'black')
  # cloneColors <- c(cloneColors,  'root' = 'red')
  
  # pcl <- plot_tree(net, locus_group_tmp$color, size_node, layoutkk,  
  #                  datatag, plottype=paste0('clone_tp_',timepoint,'_ts_',treatmentSt), save_dir)
  
  # Without label nb cells
  print("DEBUG 2")
  g1 <- ggnet2(net = net_tmp, label.size = 2, edge.size = 0.8, edge.color = 'black', 
               node.shape = 16, size = 'direct_cells_log', mode = mode, color = 'color') + 
    geom_point(aes(size = size), color = 'black', shape = 19) + 
    geom_point(aes(color = factor(color), size = size, fill = factor(color)), shape = 16) + 
    scale_size_manual(breaks = names(size_vals), labels = names(size_vals), limits = names(size_vals), values = unname(size_vals)) 
  
  g1 <- g1 + ggtitle(paste0(datatag,' ',plottype)) +
    guides(size = FALSE) +
    theme(legend.text=element_blank(),
          plot.title = element_text(size = 17,hjust = 0.5),
          legend.position = "none",
          panel.border = element_blank()) 
  
  # Label nb cells
  g2 <-
    ggnet2(net = net_tmp, label.size = 2, edge.size = 0.8,  edge.color = 'black', node.shape = 16, size = 'direct_cells_log', mode = mode, color = 'color') + 
    geom_point(aes(size = size), color = 'black', shape = 19) + 
    geom_point(aes(color = factor(color), size = size, fill = factor(color)), shape = 16) + 
    scale_size_manual(breaks = names(size_vals), labels = names(size_vals), limits = names(size_vals), values = unname(size_vals)) +
    theme(legend.position = 'top') + 
    # transition_time(factor(size)) + 
    guides(size = guide_legend(title = 'Nb. of cells'), fill = FALSE)
  
  # if (!show_colour_guide) {
  #   g1 <- g1 + scale_color_manual(values = cloneColors, guide = FALSE) 
  # }  else {
  #   g1 <- g1 + scale_color_manual(values = cloneColors) + labs(color = 'Clones')
  # }
  
  
  
  p_1 <- ggdraw() + draw_plot(g1)
  plottype <- str_replace_all(plottype, "([ ])", "_")
  png(paste0(save_dir,"trajectory_",plottype,'_',datatag,".png"),height = 2*ht, width=2*wd,res = 2*72)
  # pdf(paste(save_dir,"cluster_prevalence.pdf",sep=""), height=6, width=13)
  print(p_1)
  dev.off()
  
  p_2 <- ggdraw() + draw_plot(g2)
  
  png(paste0(save_dir,"trajectory_",plottype,'_',datatag,"_with_label.png"),height = 2*(ht+50), width=2*wd,res = 2*72)
  # pdf(paste(save_dir,"cluster_prevalence.pdf",sep=""), height=6, width=13)
  print(p_2)
  dev.off()
  
  return(p_2)
  
}

# plot_tree_condition_pdx <- function(clone_meta, net_tmp, layoutkk, v_names, locus_group, 
#                                 res, results_dir, save_dir, datatag, 
#                                 pdxid,
#                                 ht=350, wd=470){
#   plottype <- 'clone'
#   if(!is.null(mainsite)){
#     plottype=paste0('clone_pdx_',pdxid)
#     print(paste0('Plot tree trajectory at the condition: ',plottype))
#     gtree_cond <- fitclone_get_ebola_tree_conditional_pdx(res, results_dir, save_dir, datatag, 
#                                                       pdxid=pdxid)
#     dat <- gtree_cond$dat
#   } else{
#     print('Plot whole tree trajectory')
#     plottype <- 'clone_wholedata'
#     dat <- res$dat
#   }
#   
#   rownames(dat) <- dat$edges
#   # print("Contain root?")
#   # print("root" %in% rownames(dat))
#   dat <- dat[v_names,]
#   # if(length(unassigned_loci)>0){
#   #   dat[unassigned_loci,'cell_counts'] <- 0
#   # }
#   dat["root",'cell_counts'] <- 20
#   # 
#   sum(v_names %in% dat$edges) == length(v_names)
#   
#   # size_node <- dat$cell_counts
#   # size_node <- ceiling(log(size_node))
#   # size_node <- ifelse(!is.finite(size_node),0,size_node)
#   
#   locus_group_tmp <- locus_group[v_names,]
#   
#   # If at this time point, no cells present in this branch, change color of loci to black
#   rv_loci <- dat[dat$cell_counts==0,'edges']
#   locus_group_tmp[as.character(rv_loci),'color'] <- "#000000"
#   
#   # Node size by dat$cell_count
#   ord <- ceiling(log(dat$cell_counts))
#   ord[!is.finite(ord)] <- 0.0
#   names(ord) <- dat$edges
#   # ord[which(names(ord) == 'root')] <- as.character(max(ord))
#   ord <- unname(ord[v_names])
#   net_tmp %v% 'direct_cells_log' <- ord
#   point_size_scale <- 1.0
#   size_vals <- as.numeric(unique(net_tmp %v% 'direct_cells_log'))/point_size_scale
#   # Exponentiate the names
#   names(size_vals) <- sprintf('%.0f', as.numeric(unique(net_tmp %v% 'direct_cells_log'))**exp(1))
#   net_tmp %v% 'direct_cells_log' <- sprintf('%.0f', as.numeric(net_tmp %v% 'direct_cells_log')**exp(1))
#   
#   size_vals <- size_vals[order(as.numeric(size_vals))]
#   # Layout
#   mode <- layoutkk
#   net_tmp %v% "xx" <- -mode[, 1]
#   net_tmp %v% "yy" <- -mode[, 2]
#   mode <- c("xx", "yy")
#   
#   cols <- locus_group_tmp$color
#   names(cols) <- locus_group_tmp$label
#   net_tmp %v% 'color' <- cols
#   
#   # cloneColors <- clone_meta$color
#   # names(cloneColors) <- clone_meta$clone_id
#   #cloneColors <- c(cloneColors,  'root' = 'red', 'Un' = 'black')
#   # cloneColors <- c(cloneColors,  'root' = 'red')
#   
#   # pcl <- plot_tree(net, locus_group_tmp$color, size_node, layoutkk,  
#   #                  datatag, plottype=paste0('clone_tp_',timepoint,'_ts_',treatmentSt), save_dir)
#   
#   # Without label nb cells
#   g1 <- ggnet2(net = net_tmp, label.size = 2, edge.size = .3,  edge.color = 'black', node.shape = 16, size = 'direct_cells_log', mode = mode, color = 'color') + 
#     geom_point(aes(size = size), color = 'black', shape = 19) + 
#     geom_point(aes(color = factor(color), size = size, fill = factor(color)), shape = 16) + 
#     scale_size_manual(breaks = names(size_vals), labels = names(size_vals), limits = names(size_vals), values = unname(size_vals)) 
#   g1 <- g1 + ggtitle(paste0(datatag,' ',plottype)) +
#     guides(size = FALSE) +
#     theme(legend.text=element_blank(),
#           plot.title = element_text(size = 17,hjust = 0.5),
#           legend.position = "none",
#           panel.border = element_blank()) 
#   
#   # Label nb cells
#   g2 <-
#     ggnet2(net = net_tmp, label.size = 2, edge.size = .3,  edge.color = 'black', node.shape = 16, size = 'direct_cells_log', mode = mode, color = 'color') + 
#     geom_point(aes(size = size), color = 'black', shape = 19) + 
#     geom_point(aes(color = factor(color), size = size, fill = factor(color)), shape = 16) + 
#     scale_size_manual(breaks = names(size_vals), labels = names(size_vals), limits = names(size_vals), values = unname(size_vals)) +
#     theme(legend.position = 'top') + 
#     # transition_time(factor(size)) + 
#     guides(size = guide_legend(title = 'Number of cells'), fill = FALSE)
#   
#   # if (!show_colour_guide) {
#   #   g1 <- g1 + scale_color_manual(values = cloneColors, guide = FALSE) 
#   # }  else {
#   #   g1 <- g1 + scale_color_manual(values = cloneColors) + labs(color = 'Clones')
#   # }
#   
#   
#   
#   p_1 <- ggdraw() + draw_plot(g1)
#   png(paste0(save_dir,"trajectory_",plottype,'_',datatag,".png"),height = 2*ht, width=2*wd,res = 2*72)
#   # pdf(paste(save_dir,"cluster_prevalence.pdf",sep=""), height=6, width=13)
#   print(p_1)
#   dev.off()
#   
#   p_2 <- ggdraw() + draw_plot(g2)
#   png(paste0(save_dir,"trajectory_",plottype,'_',datatag,"_with_label.png"),height = 2*(ht+50), width=2*wd,res = 2*72)
#   # pdf(paste(save_dir,"cluster_prevalence.pdf",sep=""), height=6, width=13)
#   print(p_2)
#   dev.off()
#   
#   return(p_2)
#   
# }


plot_whole_tree_clones <- function(clone_meta, net_tmp, layoutkk, v_names, locus_group, 
                                res, results_dir, save_dir, datatag, 
                                ht=350, wd=470){
  plottype <- 'clone'
  print('Plot whole tree trajectory')
  plottype <- 'clone_wholedata'
  dat <- res$dat
  
  rownames(dat) <- dat$edges
  # print("Contain root?")
  # print("root" %in% rownames(dat))
  dat <- dat[v_names,]
  # if(length(unassigned_loci)>0){
  #   dat[unassigned_loci,'cell_counts'] <- 0
  # }
  dat["root",'cell_counts'] <- 20
  # 
  sum(v_names %in% dat$edges) == length(v_names)
  
  # size_node <- dat$cell_counts
  # size_node <- ceiling(log(size_node))
  # size_node <- ifelse(!is.finite(size_node),0,size_node)
  
  locus_group_tmp <- locus_group[v_names,]
  
  # If at this time point, no cells present in this branch, change color of loci to black
  rv_loci <- dat[dat$cell_counts==0,'edges']
  locus_group_tmp[as.character(rv_loci),'color'] <- "#000000"
  
  # Node size by dat$cell_count
  ord <- ceiling(log(dat$cell_counts))
  ord[!is.finite(ord)] <- 0.0
  names(ord) <- dat$edges
  # ord[which(names(ord) == 'root')] <- as.character(max(ord))
  ord <- unname(ord[v_names])
  net_tmp %v% 'direct_cells_log' <- ord
  point_size_scale <- 0.8
  size_vals <- as.numeric(unique(net_tmp %v% 'direct_cells_log'))/point_size_scale
  # Exponentiate the names
  names(size_vals) <- sprintf('%.0f', as.numeric(unique(net_tmp %v% 'direct_cells_log'))**exp(1))
  net_tmp %v% 'direct_cells_log' <- sprintf('%.0f', as.numeric(net_tmp %v% 'direct_cells_log')**exp(1))
  
  size_vals <- size_vals[order(as.numeric(size_vals))]
  # Layout
  mode <- layoutkk
  net_tmp %v% "xx" <- -mode[, 1]
  net_tmp %v% "yy" <- -mode[, 2]
  mode <- c("xx", "yy")
  
  cols <- locus_group_tmp$color
  names(cols) <- locus_group_tmp$label
  net_tmp %v% 'color' <- cols
  
  # cloneColors <- clone_meta$color
  # names(cloneColors) <- clone_meta$clone_id
  #cloneColors <- c(cloneColors,  'root' = 'red', 'Un' = 'black')
  # cloneColors <- c(cloneColors,  'root' = 'red')
  
  # pcl <- plot_tree(net, locus_group_tmp$color, size_node, layoutkk,  
  #                  datatag, plottype=paste0('clone_tp_',timepoint,'_ts_',treatmentSt), save_dir)
  
  # Without label nb cells
  g1 <- ggnet2(net = net_tmp, label.size = 2, edge.size = 0.8,  edge.color = 'black', node.shape = 16, size = 'direct_cells_log', mode = mode, color = 'color') + 
        geom_point(aes(size = size), color = 'black', shape = 19) + 
        geom_point(aes(color = factor(color), size = size, fill = factor(color)), shape = 16) + 
        scale_size_manual(breaks = names(size_vals), labels = names(size_vals), limits = names(size_vals), values = unname(size_vals)) 
        # geom_text(aes(label = label), color = "black", size = 5)
  g1 <- g1 + ggtitle(paste0(datatag,' ',plottype)) +
        guides(size = FALSE) +
        theme(legend.text=element_blank(),
              plot.title = element_text(size = 17,hjust = 0.5),
              legend.position = "none",
              panel.border = element_blank()) 
  
  # Label nb cells
  g2 <-
    ggnet2(net = net_tmp, label.size = 2, edge.size = 0.8,  edge.color = 'black', node.shape = 16, size = 'direct_cells_log', mode = mode, color = 'color') + 
    geom_point(aes(size = size), color = 'black', shape = 19) + 
    geom_point(aes(color = factor(color), size = size, fill = factor(color)), shape = 16) + 
    scale_size_manual(breaks = names(size_vals), labels = names(size_vals), limits = names(size_vals), values = unname(size_vals)) +
    theme(legend.position = 'top') + 
    # transition_time(factor(size)) + 
    guides(size = guide_legend(title = 'Nb. of cells'), fill = FALSE)
  
  # if (!show_colour_guide) {
  #   g1 <- g1 + scale_color_manual(values = cloneColors, guide = FALSE) 
  # }  else {
  #   g1 <- g1 + scale_color_manual(values = cloneColors) + labs(color = 'Clones')
  # }
  
  
  
  p_1 <- ggdraw() + draw_plot(g1)
  png(paste0(save_dir,"trajectory_",plottype,'_',datatag,".png"),height = 2*ht, width=2*wd,res = 2*72)
  # pdf(paste(save_dir,"cluster_prevalence.pdf",sep=""), height=6, width=13)
  print(p_1)
  dev.off()
  
  p_2 <- ggdraw() + draw_plot(g2)
  png(paste0(save_dir,"trajectory_",plottype,'_',datatag,"_with_label.png"),height = 2*(ht+50), width=2*wd,res = 2*72)
  # pdf(paste(save_dir,"cluster_prevalence.pdf",sep=""), height=6, width=13)
  print(p_2)
  dev.off()
  
  return(p_2)
  
}



plot_main_tree <- function(results_dir, save_dir, datatag, cellclones=NULL){
  if(is.null(cellclones)){
    cellclones <- paste0(results_dir, 'cell_clones.csv')
  }
  edge_list_path <- tree_2_edge_list_path(results_dir, save_dir)
  save_fn <- paste0('g_san_chain_ls_',datatag,'.rds')
  res <- get_main_tree(save_dir, save_fn, edge_list_path, datatag)
  
  
  edge_list <- as.matrix(igraph::as_edgelist(res$g_sans_chains))
  net <- network::network(edge_list, directed = T)
  net %e% 'names' <- 1:nrow(edge_list)
  layoutkk <- fitclone_set_dream_tree_layout(net, g_sans_chains, datatag, save_dir)
  # saveRDS(layoutkk, file = paste0(save_dir,"main_layoutkk",datatag,'.rds'))
  
  v_names <- network.vertex.names(net)
  cols <- rep('red',length(v_names))
  siz <- rep(2,length(v_names))
  dat <- res$dat
  rownames(dat) <- dat$edges
  sum(v_names %in% dat$edges) == length(v_names)
  # size_node <- ceiling(log(dat[v_names,'cell_counts']))
  size_node <- dat[v_names,'cell_counts']
  size_node <- ifelse(!is.finite(size_node),0,size_node)
  locus_group <- get_color_cluster_assignment(results_dir, save_dir, datatag, cellclones)
  # View(colnames(locus_group))
  
  locus_group <- locus_group %>% add_row(label = "root", group = "0",color="#FF0000")
  # locus_group[nrow(locus_group) + 1,] = c("root","0","#FF0000")
  "root" %in% locus_group$label
  rownames(locus_group) <- locus_group$label
  locus_group1 <- locus_group[v_names,]
  pcl <- plot_tree(net, locus_group1$color, size_node, layoutkk,  
                   datatag, plottype='clone_testing', save_dir)
}
get_clone_members <- function(clones) {
  clone_members <- list()
  for(c in unique(clones$clone_id)) {
    # if(c != "None") {
      clone_members[[c]] <- clones[clones$clone_id == c, "cell_id"]
    # }
  }
  return(clone_members)
}
data_frame_to_list <- function(clustering) {
  clusts <- unique(clustering$clone_id)
  res <- list()
  for (cl in clusts) {
    res[[cl]] <- clustering$cell_id[clustering$clone_id == cl]
  }
  return(res)
}

format_tree <- function(tree, brlen) {
  locus_tips <- grep('locus', tree$tip.label, value=TRUE)
  tree <- drop.tip(tree, locus_tips)
  
  if (!is.null(brlen)) {
    tree <- compute.brlen(tree, brlen)
  }
  
  tree$tip.label <- gsub('cell_', '', tree$tip.label)
  
  return(tree)
}
get_color_clone_v2 <- function(clone_levels){
  clone_levels <- gtools::mixedsort(clone_levels[!grepl("None", clone_levels)])
  
  ## Meta colors based on package inlmisc::GetColors(n) -- problem at installation, so keep values of color codes in csv file
  colorcode_fn <- "/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_wholedata/config/colorcode.csv"
  if(file.exists(colorcode_fn)){
    color_df <- data.table::fread(colorcode_fn)
    color_df <- color_df %>%
      dplyr::filter(nb_clones==length(clone_levels))
    clone_pal <- color_df$color
    names(clone_pal) <- clone_levels
  }else{
    clone_pal <- make_clone_palette(clone_levels)
  }
  # clone_pal <- make_clone_palette(clone_levels)
  print(clone_pal)
  # return(clone_pal)
  clone_meta <- data.frame(clone_id=names(clone_pal),color=as.character(clone_pal), stringsAsFactors = FALSE)
  clone_meta[nrow(clone_meta)+1,] <- c("None","#000000")
  return(clone_meta)
  
}  
get_color_clone <- function(cell_clones){
  aug_cut <- data_frame_to_list(cell_clones)
  print(length(aug_cut))
  myColors <- get_cluster_colours(length(aug_cut))
  if(length(myColors)>length(aug_cut)){  # in case there are only 2 levels but brewer.pal(nClusters, "Set2") generate at least 3 levels
    myColors <- myColors[1:length(aug_cut)]
  }
  names(myColors) <- names(aug_cut)
  # TODO: Since 0 is shared with the heatmap CNV state, use something else for the Un-assigned ones
  # myColors <- c(myColors, "0" = "#000000")  #unassigned cell = clone 0
  # myColors <- c(myColors, "None" = "#000000")  #unassigned cell = clone 0
  for(c in names(myColors)){
    if(c=="None"){
      print('Exist none clone')
      myColors["None"] = "#000000"
    }
  }
  
  clone_meta <- data.frame(clone_id=names(myColors),color=as.character(myColors), stringsAsFactors = FALSE)
  return(clone_meta)
}

get_cluster_colours <- function(nClusters) {
  if (nClusters > 8) {
    clust_colours <- colorRampPalette(brewer.pal(8, "Set2"))(nClusters)
  } else {
    clust_colours <- brewer.pal(nClusters, "Set2")
  }
  return(clust_colours)
}

# get_cluster_colours_v2 <- function(nClusters) {
#   if (nClusters > 8) {
#     #clust_colours <- colorRampPalette(brewer.pal(12, "Set3"))(nClusters)
#     clust_colours <- colorRampPalette(brewer.pal(8, "Set2"))(nClusters)
#   } else {
#     clust_colours <- brewer.pal(nClusters, "Set2")
#   }
#   
#   
#   # Todo: the 5th and 8th (E and H) have similar colours, change them
#   if (nClusters == 9) {
#     clust_colours[[5]] <- '#564527'
#     
#     # change I to pink too
#     clust_colours[[9]] <- '#FFC0CB'
#   }
#   #viz_colour_spectrum(myColors)
#   clust_colours
# }

get_color_cluster_assignment <- function(results_dir, save_dir, datatag, cell_clones=NULL, clone_meta=NULL){
  # if (file.exists(paste0(save_dir, "locus_group_",datatag,".rds"))) {
  #   return(readRDS(paste0(save_dir, "locus_group_",datatag,".rds")))
  # }
  newick <- paste0(results_dir, 'tree.newick')
  tree <- read.tree(newick)
  if(is.null(cell_clones)){
    cellclones <- paste0(results_dir, 'cell_clones.csv')
    cell_clones <- read.delim(cellclones, stringsAsFactors = FALSE, sep = ",", check.names = F)
    print("DEBUG*******")
    print(summary(as.factor(cell_clones$clone_id)))
  }
  
  tree <- format_tree(tree, brlen=NULL)
  clone_members <- get_clone_members(cell_clones)
  tree <- ggtree::groupOTU(tree, clone_members) 
  x <- as_tibble(tree)
  locus_ids <- grep('locus_', c(tree$tip.label,tree$node.label), value=T)
  print(length(tree$tip.label))
  print(length(tree$node.label))
  locus_group <- x[x$label %in% locus_ids, c("label","group")]
  locus_group <- as.data.frame(locus_group)
  summary(locus_group$group)
  locus_group$group <- as.character(locus_group$group)
  if(is.null(clone_meta)){
    # clone_meta <- get_color_clone(cell_clones)
    clone_meta <- get_color_clone_v2(unique(cell_clones$clone_id))
    # clone_meta[nrow(clone_meta)+1,] <- c("None","#000000")
    
    
    print('DEBUG: plot_clone_colour....')
    plot_clone_colour(clone_meta, save_dir, datatag)
  }
  print('DEBUG 11')
  print(clone_meta)
  # clone_meta$clone_id <- ifelse(clone_meta$clone_id=="None",'0',clone_meta$clone_id)
  
  locus_group <- locus_group %>% left_join(clone_meta, by = c("group"="clone_id"))
  print(summary(as.factor(locus_group$color)))
  # View(head(locus_group))
  library(tidyverse)
  locus_group <- locus_group %>% add_row(label = "root", group = "0",color="#FF0000")
  # locus_group[nrow(locus_group) + 1,] = c("root","0","#FF0000")
  print("root" %in% locus_group$label)
  rownames(locus_group) <- locus_group$label
  summary(locus_group$color)
  dim(locus_group)
  saveRDS(object = locus_group, file = paste0(save_dir, "locus_group_",datatag,".rds"))
  return(locus_group)
}

collapse_chains <- function(the_g, cell_count_dat, cell_count_treshold = 3) {
  # the_g <- g_sans_cells ; cell_count_dat <- dat
  node_names <- V(the_g)$name
  node_degrees <- igraph::degree(graph = the_g, v = cell_count_dat$edges, mode = 'out')
  df <- data.frame(names = cell_count_dat$edges, count = cell_count_dat$cell_counts, out_degree = node_degrees, stringsAsFactors = F)
  
  # Keep cells if (i) cell_count > 1 (ii) out_degree > 1
  
  keep_nodes <- df$names[(df$count > cell_count_treshold | df$out_degree > 1 | df$names == 'root')]
  del_nodes <- setdiff(node_names, keep_nodes)
  
  # Draw a gaph using the keep_nodes
  ggnet2(net = network::network(as.matrix(igraph::as_edgelist(the_g)), directed = T), label = NULL, size = 1, label.size = 2, edge.size = 0.8, edge.color = 'black', node.shape = 16)
  #plot(the_g)
  tmp_g <- the_g
  for (node in del_nodes) {
    # node = del_nodes[[2]]
    out_degree <- node_degrees[[node]]
    #print(df[df$names == node, ])
    if (out_degree == 0) { # terminal
      tmp_g <- igraph::delete.vertices(graph = tmp_g, v = node)
    } else if (out_degree == 1) {
      par <- igraph::neighborhood(graph = tmp_g, nodes=node, mode='in')[[1]]$name[[2]]
      child <- igraph::neighborhood(graph = tmp_g, nodes=node, mode='out')[[1]]$name[[2]]
      #print(sprintf('par=%s -> node=%s -> child=%s', par, node, child))
      tmp_g <- igraph::delete.vertices(graph = tmp_g, v = node)
      tmp_g <- igraph::add.edges(graph = tmp_g, edges = c(par, child))
    } else {
      stop(sprintf('An out_degree %d is not supported. It is not 0 or 1.', out_degree))
    }
  }
  
  ggnet2(net = network::network(as.matrix(igraph::as_edgelist(tmp_g)), directed = T), label = NULL, size = 1, label.size = 2, edge.size = 0.8, edge.color = 'black', node.shape = 16)
  
  tmp_g
}




get_libid_from_cell_names <- function(cell_names) {
  # TODO: add the condition for SA906 and SA666
  gsub('SA([0-9]|[A-Z]|[a-z])+-(A([0-9]|[A-Z])+)-.*', '\\2', cell_names)
}

get_library_id <- function(cell_ids) {
  labels <- sapply(strsplit(cell_ids, "-"), function(x) {
    return(x[2])
  })
  return(labels)
}

get_metadata <- function(results_dir){
  grouping_df <- read.csv(paste0(results_dir,'library_groupings.csv'),header=T, check.names=F, 
                          row.names = 1, stringsAsFactors=F)
  # colnames(grouping_df)[which(names(grouping_df) == "passage")] <- "timepoint"
  # colnames(grouping_df)[which(names(grouping_df) == "treatment_st")] <- "treatmentSt"
  colnames(grouping_df)[which(names(grouping_df) == "grouping")] <- "library_id"
  return(grouping_df)
}
get_timepoint <- function(results_dir){
  grouping_df <- read.csv(paste0(results_dir,'library_groupings.csv'),header=T, check.names=F, 
                          row.names = 1, stringsAsFactors=F)
  colnames(grouping_df)[which(names(grouping_df) == "passage")] <- "timepoint"
  colnames(grouping_df)[which(names(grouping_df) == "treatment_st")] <- "treatmentSt"
  colnames(grouping_df)[which(names(grouping_df) == "grouping")] <- "library_id"
  return(grouping_df)
}

parse_cell_names <- function(cell_names, results_dir) {
  # Find libids
  #libids <- libid_from_cell_id(cell_names)
  libids <- get_libid_from_cell_names(cell_names)
  cell.dat <- data.frame(library_id = libids, cellids = cell_names, stringsAsFactors = F)

  # cell.dat <- dplyr::left_join(cell.dat, core_dat[, c('timepoint', 'library_id')], by=c('library_id'))
  # stopifnot(all(cell_names == cell.dat$cellids, na.rm = T))

  grouping_df <- get_timepoint(results_dir)
  cell.dat <- dplyr::left_join(cell.dat, grouping_df[, c('timepoint', 'library_id')], by=c('library_id'))

  # Set SA928 to itself
  # cell.dat$timepoint[grepl('SA928', cell.dat$cellids)] <- cell.dat$cellids[grepl('SA928', cell.dat$cellids)]

  cell.dat$timepoint
}


plot_clone_colour <- function(meta_data, save_dir, datatag){
  # meta_data <- data.frame(clone_id=clone_id, color=color)
  
  meta_data <- meta_data[!meta_data$clone_id %in% c('0','None','Unassigned','Un'),]
  meta_data <- meta_data[!is.na(meta_data$clone_id),]
  meta_data <- meta_data[order(meta_data$clone_id),]
  write.csv(meta_data, paste0(save_dir,'clone_colors.csv'), quote=F, row.names = F)
  plottype <- 'clone_id'
  xstring <- 'clone_id'
  ystring <- 'color'
  plg_cl <- legend_function(meta_data, xstring, ystring, plottype,
                            plottitle="Clones", 
                            colorcode=meta_data$color, 
                            legendtitle="Clone Id", legendlabels=NULL)
  
  
  lg <- cowplot::get_legend(plg_cl)
  saveRDS(lg, file=paste0(save_dir,"lg_",plottype,".rds"))
  # lg <- readRDS(paste0(save_dir,"lg_",plottype,".rds"))
  
  pc <- ggdraw() + draw_plot(lg)
  png(paste(save_dir,"lg_",plottype,'_',datatag,".png",sep=""),height = 2*500, width=2*150,res = 2*72)
  # print(grid.arrange(grobs = plots[select_grobs(layout)], layout_matrix = layout,
  #                    bottom=" ",right=" "))
  print(pc)
  dev.off()
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

# retrieve_time_point_per_clade <- function(results_dir, the_graph, decendents=NULL, data_tag=NULL, 
#                                           trim_cell_name=TRUE, include_loci=FALSE) {
#   if (is.numeric(unlist(decendents)))
#     decendents <- node_number_to_label(the_graph, decendents)
#   res <- list()
#   for (clade_root in names(decendents)) {
#     #print(clade_root)
#     # clade_root = names(decendents)[[1]]
#     cells <- decendents[[clade_root]]
#     # cells <- grep('SA609', cells, value = T)
#     if (trim_cell_name) {
#       # cells <- parse_cell_names(cell_names = cells)
#       cells <- parse_cell_names(cell_names = cells, results_dir)
#     }
#     res[[clade_root]] <- cells
#   }
#   res
# }

# get_edge_data <- function(the_graph) {
#   # Number of edges=Number of nodes - number of leaves
#   # TODO: generalise this...
#   
#   # TODO: @Sohrab: remove this
#   # out_ts <- sprintf('%s_%s_ts_dat.rds', './outputs/', edge_density(g))
#   # if (file.exists(out_ts)) return(readRDS(out_ts))
#   
#   n_nodes <- length(V(the_graph)$name)
#   n_leaves <- length(V(the_graph)$name[!grepl("root", V(the_graph)$name) & !grepl("locus", V(the_graph)$name)])
#   n_edges <- n_nodes - n_leaves
#   
#   # internal_nodes
#   i_nodes <- V(the_graph)$name[grepl("root", V(the_graph)$name) | grepl("locus", V(the_graph)$name)]
#   
#   # Account for GM cells
#   # TODO REMOVE
#   # GM_nodes <- grep("SA928", V(the_graph)$name, value = TRUE)
#   
#   # Sanity check
#   # stopifnot(sum(tcc$Freq) == (n_leaves - length(GM_nodes)))
#   
#   # Find all descendents
#   desc <- get_decendents(clone_roots = i_nodes, the_graph = the_graph, min_cell_per_clone = 1)
#   # Find out which one is a cells and which is a locus (will be NA)
#   res <- retrieve_time_point_per_clade(the_graph = the_graph, 
#                                        decendents = desc)
#   
#   # Remove all NAs entries in res
#   for (cn in names(res)) {
#     if (all(is.na(res[[cn]]))) {
#       res[[cn]] <- NULL
#     } else {
#       res[[cn]] <- res[[cn]][!is.na(res[[cn]])]
#     }
#   }
#   
#   # Count the number of cells in each subgroup
#   ts <- lapply(names(res), function(x) {length(res[[x]])})
#   ts <- data.frame(clone_id = names(res), N = unlist(ts), stringsAsFactors = FALSE)
#   stopifnot(ts$N[ts$clone_id == "root"] == n_leaves)
#   
#   ts$frac <- ts$N / n_leaves
#   stopifnot(max(ts$frac) <= 1 & min(ts$frac) >= 0)
#   
#   # Add height that
#   height_dat <- get_height_dat(the_graph)
#   height_dat <- height_dat[height_dat$id %in% ts$clone_id, ]
#   colnames(height_dat) <- c("clone_id", "height")
#   
#   ts <- dplyr::right_join(ts, height_dat, by = c("clone_id"))
#   
#   ts <- ts[order(ts$height), ]
#   ts$id <- paste0("(height=", ts$height, ", frac=", format(ts$frac, digits=2), ")")
#   
#   # TODO:@Sohrab remove this
#   # saveRDS(ts, out_ts)
#   
#   return(ts)
# }

# get_the_cut <- function(the_graph, minimum_fraction, maximum_fraction) {
#   ts <- get_edge_data(the_graph)
#   
#   candiate_edges <- ts[ts$frac > minimum_fraction,]
#   
#   the_cut <- h1_auto_cut(the_graph = the_graph,
#                          candiate_edges = candiate_edges,
#                          minimum_fraction = minimum_fraction,
#                          maximum_fraction = maximum_fraction
#   )
#   
#   return(the_cut)
# }





get_tcc <- function(g, results_dir) {
  tcc <- as.data.frame(table(parse_cell_names(cell_names = V(g)$name[!grepl('root', V(g)$name) & !grepl('locus', V(g)$name)],results_dir)))
  # Remove GM cells
  tcc <- tcc[grep('X[0-9][0-9]*', tcc$Var1), ]
  tcc[order(as.numeric(gsub('X', '', tcc$Var1))), ]
}




read_ltm_tree <- function(edge_list_path) {
  if (edge_list_path == '') {
    return(NULL)
  }
  ss <- read.csv(edge_list_path, stringsAsFactors = F, header=T)
  # Find the root
  dd <- (as.matrix(ss))
  g <- graph_from_edgelist(dd)
  V(g)$id <- seq(vcount(g))
  g
}

# create_cell_assignment_dat <- function(results_dir, datatag, edge_list_path, clone_dic, g = NULL, the_cut = NULL, aug_cut = NULL) {
#   # All assigned cells
#   if (is.null(aug_cut)) {
#     desc <- get_decendents(the_cut, the_graph = g, min_cell_per_clone = 1)
#     # res <- retrieve_time_point_per_clade(g, desc, datatag, FALSE)
#     res <- retrieve_time_point_per_clade(g, desc,FALSE)
#   } else {
#     # res <- retrieve_time_point_per_clade(g, aug_cut, datatag, FALSE)
#     res <- retrieve_time_point_per_clade(g, aug_cut,FALSE)
#   }
#   
#   #names(res) <- gsub('locus_', '', names(res))
#   
#   dat <- ltm_clust_to_df(res, names(res))
#   dat$time <- parse_cell_names(dat$single_cell_id, results_dir)
#   dat$genotype <- as.character(dat$genotype)
#   clone_dic$old_K <- as.character(clone_dic$old_K)
#   dat <- dplyr::right_join(dat, clone_dic, by=c('genotype'='old_K'))
#   
#   # For unassigned cells
#   all_cells <- V(g)$name[!grepl('root', V(g)$name) & !grepl('locus_', V(g)$name)]
#   # Remove GM cells
#   # TODO: centralised GM cells removal
#   if (datatag %in% c('SA906a', 'SA906b')) {
#     all_cells <- all_cells[grepl(datatag, all_cells)]
#   } else if (datatag == 'SA922') {
#     # Just do nothing
#   } else {
#     all_cells <- all_cells[grep('X[0-9]+', all_cells)]
#   }
#   
#   un_assigned_cell <- setdiff(all_cells, dat$single_cell_id)
#   # any(un_assigned_cell %in% dat$single_cell_id)
#   #dat1 = data.frame(single_cell_id=un_assigned_cell, genotype='unassigned', time=parse_cell_names(un_assigned_cell), K=-1, is_ref=NA, letters = 'Un', pretty_names='')
#   if (length(un_assigned_cell) > 0) {
#     dat1 <- data.frame(single_cell_id=un_assigned_cell, genotype='unassigned', 
#                        time=parse_cell_names(un_assigned_cell, results_dir), 
#                        is_ref=NA, letters = 'Un', pretty_names='', K = NA, 
#                        stringsAsFactors = F)
#     
#     dat <- dplyr::bind_rows(dat, dat1)
#   }
#   
#   dat
# }

addrow <- function(orig, tmp) {
  if (is.null(orig)) {
    orig <- tmp
  } else {
    orig <- rbind(orig, tmp)
  }
  orig
}


find_parent <- function(the_g, node) {
  igraph::neighborhood(the_g, nodes=node, mode='in')[[1]]$name[[2]]
}

