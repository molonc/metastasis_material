project_dir <- '/home/htran/Projects/hakwoo_project/corrupt_tree/'
source(paste0(project_dir, "src/tree_viz/make_cell_copynumber_tree_heatmap.R"))
source(paste0(project_dir, "src/tree_viz/tree_viz.R"))
source(paste0(project_dir, "src/tree_viz/utils.R"))



# Test Shaocheng
# test_dir <- '/home/htran/Projects/hakwoo_project/rscript/'
# newick <- paste0(test_dir, 'outtree_TFRIPAIR8_merged')
# tree <- read.tree(newick)
# p <- ggtree(tree)
# p
# edge_list <- tree_2_edge_list(tree)

results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_2252_global_97/'

# copynumber <- paste0(results_dir, 'total_merged_filtered_states.csv')
cellclones <- paste0(results_dir, 'tree_cut_out/cell_clones_if_0.1_af_0.55_p0.75_e0.04.csv')
# copy_number <- read.csv(copynumber, header=T, row.names=1, check.names = F,stringsAsFactors = FALSE)
# copynumber <- read.delim(copy_number, check.names = FALSE, stringsAsFactors = FALSE, sep = ",")

# copy_number <- format_copynumber_matrix(copy_number)
cell_clones <- read.delim(cellclones, stringsAsFactors = FALSE, sep = ",")
ex_clones <- c('I','J')

rownames(cell_clones) <- cell_clones$cell_id

results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_new_encode/SA1035_whole_local/'
plot_title <- 'ct_v1_8clones'
cellclones <- paste0(results_dir,'tree_cut_out/cell_clones_if_0.1_af_0.65_p0.75_e0.04.csv')
plot_tree(cellclones, results_dir, plot_title)


cellclones <- paste0(results_dir,'tree_cut_out/cell_clones_if_0.02_af_0.45_p0.75_e0.04.csv')
plot_title <- 'ct_v2_13clones'
plot_tree(cellclones, results_dir, plot_title)

plot_tree <- function(cellclones, results_dir, plot_title='tree_trajectory'){

    # The tree
  newick <- paste0(results_dir, 'tree.newick')
  tree <- read.tree(newick)
  the_graph <- read_ltm_tree(tree_2_edge_list(tree))
  res <- convert_edge_list_to_ape(edge_list_mat = igraph::as_edgelist(the_graph))
  
  
  # tree$edge.length=rep(5, dim(tree$edge)[1]) nicer tree
  # ctree <- res$tree
  node_names <- res$node_names
  cell_clones <- read.delim(cellclones, stringsAsFactors = FALSE, sep = ",")
  aug_cut <- data_frame_to_list(cell_clones)
  myColors <- get_cluster_colours(length(aug_cut))
  names(myColors) <- names(aug_cut)
  
  # TODO: Since 0 is shared with the heatmap CNV state, use something else for the Un-assigned ones
  myColors <- c(myColors, "None" = "#000000")
  
  node_names <- node_names %>% left_join(cell_clones, by = c("node_name" = "cell_id"))
  sum(is.na(node_names$clone_id))
  # unique(node_names$clone_id)
  rownames(node_names) <- node_names$node_number
  node_names$clone_id <- ifelse(is.na(node_names$clone_id),'None',node_names$clone_id)
  color_df <- data.frame(clone_id=names(myColors),color=as.character(myColors), stringsAsFactors = FALSE)
  node_names <- node_names %>% inner_join(color_df, by = "clone_id")
 
  sum(is.na(node_names$clone_id))
  # the_graph <- read_ltm_tree(ctree$edge)
  # node_names <- node_names[ctree$tip.label,]
  V(the_graph)$node_name <- node_names[V(the_graph)$id, 'node_name']
  V(the_graph)$clone_id <- node_names[V(the_graph)$id, 'clone_id']
  V(the_graph)$color <- node_names[V(the_graph)$id, 'color']
  V(the_graph)$label.color <- "black"
  V(the_graph)$label <- NA
  V(the_graph)$name <- node_names[V(the_graph)$id, 'node_name']
  # V(the_graph)$size <- 0.5
  E(the_graph)$weight <- 10
  
  # Set edge width based on weight:
  # E(net)$width <- E(net)$weight/6
  # E(net)$edge.color <- "gray80"
  the_root <- find_root(the_graph)
  root_name <- node_names[node_names$node_number==the_root,'node_name']
  ht <- get_height_dat(the_graph)
  rownames(ht) <- ht$id
  ht1 <- ht
  
  none_cells <- node_names[node_names$clone_id=='None','node_name']
  idx <- grep('locus_', ht1$id, value = T)
  ht1[idx,'dist'] <- 0
  V(the_graph)[the_root]$color <- '#FF0000'
  
  ht1[root_name,'dist'] <- 60
  
  ht1 <- ht1 %>% left_join(node_names, by = c("id" = "node_name"))
  rownames(ht1) <- ht1$id
  ht1$dist <- ifelse(ht1$clone_id=='None',0,ht1$dist)
  
  length(unique(ht1$id))
  dim(ht1)
  # the_graph1 <- the_graph
  the_graph1 <- igraph::delete.vertices(the_graph, c(idx, none_cells))
  ht2 <- ht1[V(the_graph1)$name,]
  
  
  # edge_list_mat = igraph::as_edgelist(the_graph1)
  # edge_list <- as.data.frame(edge_list_mat)
  # colnames(edge_list) <- c("source", "target")
  # View(head(edge_list))
  # edge_list$none_clone <- ifelse(edge_list$target %in% none_cells,TRUE,FALSE)
  # sum(edge_list$none_clone==T)
  # 
  # sum(edge_list$source %in% none_cells)
  # lfr <- layout_with_gem(the_graph)
  # l <- layout_with_lgl(the_graph)
  # l <- layout_nicely(the_graph, dim = 2)
  # l <- layout_as_tree(graph = the_graph, root = c(the_root))
  # l <- layout_with_lgl(the_graph, root = the_root)
  vcount(the_graph1)
  # l <- layout_with_lgl(the_graph)
  # l1 <- layout_with_lgl(the_graph1)
  # l <- layout_with_fr(the_graph)
  l1 <- layout_with_kk(the_graph1)
  
  # pdf(paste0(results_dir,'trajectory_plot.pdf'),width = 10, height = 9)
  # png(paste0(results_dir,plot_title,'.png'),width = 1050, height = 800)
  # plot(the_graph, rescale=T, layout=l,  #layout_with_lgl, l*4, 
  #      vertex.size=5  #ht1$dist*0.07,  #change color to metastasis label, vertex.size=ht$dist*0.4
  # ) #edge.curved=.1, edge.color='black'
  # legend(x=-1.2, y=-0.2, names(myColors), pch=21,
  #        
  #        col="#777777", pt.bg=myColors, pt.cex=2.5, cex=1, bty="n", ncol=1)
  # dev.off()
  
  png(paste0(results_dir,plot_title,'.png'),width = 1050, height = 800)
  plot(the_graph1, rescale=T, layout=l1,  #layout_with_lgl, l*4, 
       vertex.size=ht2$dist*0.03,edge.arrow.size=.01  #change color to metastasis label, vertex.size=ht$dist*0.4
  ) #edge.curved=.1, edge.color='black'
  legend(x=-1.2, y=-0.2, names(myColors), pch=21,
         
         col="#777777", pt.bg=myColors, pt.cex=2.5, cex=1, bty="n", ncol=1)
  dev.off()
  
}


# # Visualize tree
# ctree$edge.length <- rep(4,length(ctree$edge.length))
# p <- ggtree(ctree, right = TRUE)
# plot(p)
# p <- p %<+% node_names + geom_tippoint(aes(color=color), size=1)
# plot(p)




# png(paste0(results_dir,'trajectory_plot.png'),width = 1000, height = 800)
# plot(the_graph, rescale=T, layout=lfr * 3,  #layout_with_lgl, l*4
#      vertex.color=node_names$color,vertex.size=ht$dist*0.4,
#      edge.color='black') #edge.curved=.1,
# legend(x=-1.2, y=-0.2, c("B","A","C","Non-assign"), pch=21,
#        
#        col="#777777", pt.bg=myColors, pt.cex=3, cex=1.1, bty="n", ncol=1)
# dev.off()


# p <- ggtree(tree, right = TRUE)
# plot(p)
# 
# 
# 
# tip.loci <- grep("locus", tree$tip.label, value = TRUE)
# tree <- drop.tip(tree, tip.loci)
# 
# 
# # Plot
# aug_cut <- data_frame_to_list(cell_clones)
# 
# #browser()
# 
# xtree <- ggtree::groupOTU(tree, aug_cut)
# xtree$edge.length=rep(5, dim(tree$edge)[1])
# p <- ggtree(xtree, aes(color = group))
# p
# # p <- sos_heat(p = p,
# #               data = genotype,
# #               offset = offset,
# #               width = 8,
# #               colnames = TRUE,
# #               color = NULL,
# #               colnames_level = colnames_level,
# #               colnames_filter = colnames_filter,
# #               font.size = chr_font_size,
# #               colnames_offset_y = -40
# # ) +
# #   scale_fill_manual(breaks = breaks,
# #                     values = cn_colours,
# #                     guide = guide_legend(nrow = 1, direction = "horizontal", label.position = "bottom")
# #   )
# 
# # Add the cluster annotations
# p <- NULL
# p <- annoate_tree_aug(p = p,
#                       the_graph = the_graph,
#                       candiate_edges = NULL,
#                       aug_cut = aug_cut,
#                       the_height = the_height,
#                       show_guide = FALSE,
#                       show_label = FALSE
# )
# 
# 
# 
# candiate_edges <- get_dummy_candid_edge(the_graph, aug_cut)
# 
# # A hack to make colouring work
# orig_p <- NULL
# 
# # TODO: Sohrab, sort this out...
# res <- graph_to_tree_dic(the_graph)
# 
# tree <- res$tree
# node_names <- res$node_names
# 
# 
# # Cut the tail
# if (!is.null(the_height)) {
#   tree_res <- trim_tree_before_height(tree, node_names, the_graph, the_height)
#   tree <- tree_res$tree
# }
# 
# myColors <- get_cluster_colours(length(aug_cut))
# names(myColors) <- names(aug_cut)
# 
# # TODO: Since 0 is shared with the heatmap CNV state, use something else for the Un-assigned ones
# myColors <- c(myColors, "0" = "black")
# 
# cdat <- data.frame(clone_id = names(myColors), colours = myColors, stringsAsFactors = FALSE)
# 
# # Group 0, the default of the tree should be black
# if (is.null(p)) {
#   xtree <- ggtree::groupOTU(tree, aug_cut)
#   p <- ggtree(xtree, aes(color = group), size = 1)
# }
# 
# p <- p + scale_color_manual(values = myColors, guide = FALSE)
# p
# 
# V(the_graph)$color <- 
#   
# V(the_graph)$label <- "" 
# l <- layout_with_fr(the_graph)
# 
# l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
# plot(the_graph, rescale=F, layout=l*2, )
# ggtree(tree, aes(color = group), size = 1)
# 
# E(the_graph)
# 
# xtree$tip.label[1:3]
# length(xtree$edge)
# xtree$edge.leng
