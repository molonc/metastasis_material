get_pretty_names_for_loci <- function(str_array) {
  # str_array = candiate_edges$clone_id
  str_array <- gsub('locus_', '', str_array)
  chr = gsub('([0-9]+|X|Y)_.*', '\\1', str_array)
  str_array <- gsub('.*_([0-9]+)_.*', '\\1', str_array)
  str_array <- substr(str_array, 1, 4)
  str_array <- paste0('chr', chr, '_', str_array, '')
  gsub('chrroot_root', 'root', str_array)
}

collapse_tree_by_edge <- function(the_g, nodes, do_plot=FALSE, annotate_nodes=NULL) {
  # the_g = g
  # nodes = candiate_edges$clone_id
  h <- get_height_dat(the_g)
  h <- h[h$id %in% nodes, ]
  h <- h[order(h$dist, decreasing = T), ]
  
  # The neighbour with the smallest height is the parent
  edge_list <- matrix(NA, nrow=nrow(h), ncol=2, dimnames = list(c(), c('source', 'target')))
  for (i in seq(nrow(h))) {
    # i = 1
    p <- find_all_parents(the_g, h$id[i], h$dist[i])
    parent <- 'root'
    if (nrow(h[h$id %in% p, ]) > 0)
      parent <- h$id[h$id %in% p][1]
    
    edge_list[i, ] <- c(parent, h$id[i])
  }
  
  orig_edge_list = edge_list
  edge_list[, 1] <- get_pretty_names_for_loci(edge_list[, 1])
  edge_list[, 2] <- get_pretty_names_for_loci(edge_list[, 2])
  
  # Convert cell name to clone number
  clone_g <- igraph::graph_from_edgelist(el = edge_list)
  #plot(clone_g)
  #saveRDS(list(el=el, g=clone_g), paste0('/Users/sohrabsalehi/projects/fitclone/presentations/July_31/', datatag, '/clone_graph5_clones.rds'))
  # if (do_plot) {
  #   out_path = file.path('~/projects/fitclone/presentations/August28/edge_cut/', paste0(datatag, '_clone_tree.png'))
  #   plot_graph(ft = as.matrix(edge_list), out_path = out_path, annotate_nodes = get_pretty_names_for_loci(annotate_nodes))
  #   print(sprintf('Results in %s', out_path))
  # }
  
  gg <- igraph::graph_from_edgelist(el = orig_edge_list)
  V(gg)$id <- seq(vcount(gg))
  list(edge_list=orig_edge_list, graph=gg)
}

# Most recent common ancestor
get_mrca_for_clone_roots <- function(the_g, aug_cut, h = NULL, check_graph_integrity = TRUE) {
  if (is.null(h)) {
    h <- get_height_dat(the_g)
  }
  h <- h[order(h$dist, decreasing = T), ]
  
  clone_roots <- list()
  for (cn in names(aug_cut)) {
    # cn = names(aug_cut)[[1]]
    # cn = 'B_A'
    print(cn)
    nodes <- aug_cut[[cn]]
    
    # Only look at nodes with minimum height
    sub.h <- h[h$id %in% nodes, ]
    if (length(nodes) > 10) {
      nodes <- sub.h$id[sub.h$dist == min(sub.h$dist)]
    } 
    
    #sub.h[sub.h$dist == min(sub.h$dist), ]
    all_parents <- c()
    for (node in nodes) {
      # node = nodes[[1]]
      # node = nodes[[2]]
      #temp <- find_all_parents(the_g = the_g, node = node, max_height = 1)
      temp <- find_all_parents(the_g = the_g, node = node, max_height = max(h$dist))
      if (length(all_parents) == 0) {
        all_parents <- temp
      } else {
        all_parents <- intersect(all_parents, temp)
      }
    }
    # Pick the MRCA
    #clone_roots[[cn]] <- rev(h$id[h$id %in% all_parents])[1]
    clone_roots[[cn]] <- h$id[h$id %in% all_parents][1]
    
    # Sanity check - does this parent contain all the nodes?
    # annotate_nodes_on_tree(p = NULL, edge_list_path = NULL, node_names = c(aug_cut[[cn]], clone_roots[[cn]]), tree_node_dic = graph_to_tree_dic(the_g), show_label = F)
    
    desc <- find_all_descendents(the_g, clone_roots[[cn]], max(h$dist))
    # Set to false if you're using a modified graph, say after cutting finger or tail
    if (check_graph_integrity) {
      stopifnot(all(aug_cut[[cn]] %in% desc))
    }
    
    # max_height <- max(max(sub.h$dist) - min(sub.h$dist), 1)
    # all_parents <- c()
    # for (node in nodes) {
    #   temp <- find_all_parents(the_g = the_g, node = node, max_height = max_height)
    #   if (length(all_parents) == 0) {
    #     all_parents <- temp
    #   } else {
    #     all_parents <- intersect(all_parents, temp)
    #   }
    # }
    # # Pick the MRCA
    # clone_roots[[cn]] <- rev(h$id[h$id %in% all_parents])[1]
  }
  
  unlist(unname(clone_roots))
}


MRCA_based_collapse_tree <- function(the_g, aug_cut) {
  # Sort nodes by their height and merge pairs from the highest to lowest
  #the_g = g
  h <- get_height_dat(the_g)
  h <- h[order(h$dist, decreasing = T), ]
  
  # Get the mrca for all the nodes in the clone
  the_nodes <- get_mrca_for_clone_roots(the_g, aug_cut, h)
  
  # annotate_nodes_on_tree(p = NULL, edge_list_path = NULL, node_names = the_nodes, tree_node_dic = graph_to_tree_dic(the_g), show_label = F)
  edge_node_dic <- names(aug_cut)
  names(edge_node_dic) <- the_nodes
  
  node_heights <- h[h$id %in% the_nodes, ]
  
  # Find most recent common ancestor with all the current ones in the queue and pick the highest one (further from the root)
  additional_nodes <- list()
  
  node_queue_bag <- node_heights$id
  for (i in seq(nrow(node_heights)-1)) {
    # i = 1
    print(i)
    MRCA_list <- c()
    for (j in 2:length(node_queue_bag)) {
      current_nodes <- node_queue_bag[c(1,j)]
      print(c('Joining ', current_nodes))
      MRCA_list[[node_queue_bag[[j]]]] <- find_MRCA(the_g, current_nodes, h)
    }
    
    MRCA <- h$id[h$id %in%  unname(unlist(MRCA_list))][1]
    joining_companion_index <- unname(which(MRCA_list == MRCA))
    joining_companion <- names(MRCA_list)[joining_companion_index]
    
    node_queue_bag <- node_queue_bag[-c(1,joining_companion_index)]
    
    print(c('MRCA was ', MRCA))
    #annotate_nodes_on_tree(p = NULL, edge_list_path = NULL, node_names = current_nodes, tree_node_dic = graph_to_tree_dic(the_g), show_label = F)
    
    additional_nodes <- append(additional_nodes, MRCA)
    
    # Add MCRA to the queue
    node_queue_bag <- c(MRCA, node_queue_bag)
  }
  
  # node_queue <- node_heights$id
  # for (i in seq(nrow(node_heights)-1)) {
  #   current_nodes <- node_queue[c(1,2)]
  #   node_queue <- node_queue[-c(1,2)]
  #   print(c('Joining ', current_nodes))
  #   MRCA <- find_MRCA(the_g, current_nodes, h)
  #   print(c('MRCA was ', MRCA))
  #   #annotate_nodes_on_tree(p = NULL, edge_list_path = NULL, node_names = current_nodes, tree_node_dic = graph_to_tree_dic(the_g), show_label = F)
  #   
  #   additional_nodes <- append(additional_nodes, MRCA)
  #   
  #   # Add MCRA to the queue
  #   node_queue <- c(MRCA, node_queue)
  # }
  
  all_nodes <- c(unlist(additional_nodes), node_heights$id)
  # annotate_nodes_on_tree(p = NULL, edge_list_path = NULL, node_names = all_nodes, tree_node_dic = graph_to_tree_dic(the_g), show_label = F)
  
  res <- collapse_tree_by_edge(the_g = the_g, nodes = all_nodes, do_plot = FALSE, annotate_nodes = node_heights$id)
  
  # Updata the clone root names using aug_cut
  res$edge_list[, 1] <- translate(s_arry = res$edge_list[, 1], the_dic = edge_node_dic)
  res$edge_list[, 2] <- translate(s_arry = res$edge_list[, 2], the_dic = edge_node_dic)
  
  res$graph <- igraph::graph_from_edgelist(el = res$edge_list, directed = T)
  
  res
}

find_all_descendents <- function(the_g, node, max_height) {
  igraph::neighborhood(the_g, nodes=node, mode='out', order = max_height)[[1]]$name[-c(1)]
}
# |the_dic|: a named list, where names(the_dic) are terms in s_array
translate <- function(s_arry, the_dic, trans_func = NULL) {
  # Looks up the terms in s_array in the src column of the_dic, and returns the value
  r1 <- unname(the_dic[s_arry])
  
  # Handle null values
  if (!is.null(trans_func)) {
    r1[is.na(r1)] <- trans_func(s_arry[which(is.na(r1))])
  } else {
    r1[is.na(r1)] <- s_arry[which(is.na(r1))]
  }
  
  r1
}
setup_graph <- function(g=NULL) {
  V(g)$id <- seq(vcount(g))
  h <- get_height_dat(g)
  stopifnot(all(h$id == V(g)$name))
  V(g)$height <- as.integer(h$dist)
  g
}
fast_get_summary_edges <- function(cell_clones, newick) {
  #aug_cut <- get_aug_cut(NULL, datatag, edge_list_path); str(aug_cut)
  #g <- read_ltm_tree(edge_list_path)
  aug_cut <- data_frame_to_list(clustering = cell_clones[, c('cell_id', 'clone_id')])  # %>% dplyr::rename(cluster = clone_id, single_cell_id = cell_id)
  tree <- read.newick(newick)
  # all_cells <- grep('cell', tree$tip.label, value = T)
  # print(length(all_cells)) 
  # none_cells <- all_cells[!all_cells %in% cell_clones$cell_id]
  # tree <- ape::drop.tip(tree, none_cells, trim.internal =T, collapse.singles = T)
  g <- fast_read_ltm_tree(tree_2_edge_list(tree))
  
  edge_list <- MRCA_based_collapse_tree(the_g = g, aug_cut = aug_cut)
  edge_mat <- edge_list$edge_list
  
  edge_mat
}

fast_get_summary_tree <- function(cell_clones, newick, cloneColours = NULL, 
                                  edge_list = NULL, point_size_1 = 11, point_size_2 = 9) {
  # Add the summary phylogeny to the side
  #edge_list <- get_summary_edges(datatag, edge_list_path)
  if (is.null(edge_list)) edge_list <- fast_get_summary_edges(cell_clones, newick)
  print(dim(edge_list))
  # View(edge_list)
  # fwrite(edge_list, '~/projects/fitness_material/reviews/tables/SA609_summary_tree_draft.csv')
  net <- network::network(as.matrix(edge_list), directed = T)
  
  if (!is.matrix(edge_list)) edge_list <- as.matrix(edge_list %>% as.data.frame())
  
  # Generate a tree loyout
  summary_g <- graph_from_edgelist(el = edge_list, directed = T)
  summary_g <- setup_graph(summary_g)
  # plot(summary_g)
  x <- layout.reingold.tilford(summary_g)
  
  # Sort the names to match the network object
  ord <- seq(nrow(x))
  names(ord) <- V(summary_g)$name
  x <- x[unname(ord[network.vertex.names(net)]), ]
  
  net %v% "xx" <- -x[, 1]
  net %v% "yy" <- -x[, 2]
  
  # Add the colours 
  if (is.null(cloneColours)) stop()
  
  
  myColors <- cloneColours
  # Add empty clone for all loci (e.g., root) not in the clone_dic. These constitue the connecting edges
  all_nodes <- unique(c(edge_list[, 1], edge_list[, 2]))
  tmp_nodes <- all_nodes[!(all_nodes %in%  names(cloneColours))]
  
  # Set the temp nodes to black
  myColors <- as.list(myColors)
  for (tmp_node in tmp_nodes) { 
    #myColors[[tmp_node]] <- 'black'
    myColors[[tmp_node]] <- 'white'
  }
  myColors <- unlist(myColors)
  
  net %v% 'letter' <- myColors[network.vertex.names(net)]
  
  # Set the ones that are black to nothing...
  # Unless they are clone A...
  tmp.labs <- network.vertex.names(net)
  tmp.labs[(net %v% 'letter' == 'black') & (tmp.labs != 'A') ] <- ''
  tmp.labs[tmp.labs == 'None'] <- ''
  net %v% 'labels' <- tmp.labs
  
  # Geom_text actually reads from net %v% 'vertex.names', so updaet that too...
  tmp.labs <- network.vertex.names(net)
  tmp.labs[net %v% 'letter' == 'black'  & (tmp.labs != 'A')] <- ''
  tmp.labs[tmp.labs == 'None'] <- ''
  net %v% 'vertex.names' <- tmp.labs
  
  g1 <- ggnet2(net = net, label = 'labels', mode = c("xx", "yy"), color = 'letter', size = 0, label.size = 7, edge.size = 1, edge.color = 'black', node.shape = 16) +
    geom_point(aes(color = color), size = point_size_1, color = "black", shape = 16) +
    geom_point(aes(color = color), size = point_size_2, shape = 16) +
    #geom_text(aes(label = label), color = "black", size = 5)
    geom_text(aes(label = label), color = "white", size = 5) + 
    theme(
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
    )  
  # Add extra margins; It should be easier than this...
  coef <- .1
  plot_width <- max(net %v% "yy") - min(net %v% "yy") 
  lim_top_x <- abs(min(net %v% "yy")) + coef*plot_width
  lim_bottom_x <- max(net %v% "yy") + coef*plot_width
  
  plot_height <- max(net %v% "xx") - min(net %v% "xx") 
  lim_top_y <- abs(min(net %v% "xx")) + coef*plot_height
  lim_bottom_y <- max(net %v% "xx") + coef*plot_height
  
  g1 + coord_flip(ylim = c(-lim_top_x, lim_bottom_x), xlim = c(-lim_top_y, lim_bottom_y)) 
}

# viz_summary_tree <- function(datatag){
#   cc <- ftt('/Users/sohrabsalehi/projects/fitness_material/fitness_code_repo/fitclone/data/SA535_combined_CISPLATIN_CX5461_to_Sohrab/processed/lumberjack/JANKV6689P/cell_clones_updated.csv')
#   clone_names <- sort(u(cc %>% dplyr::pull('clone_id')))
#   # cloneColours <- as.character(inlmisc::GetColors(length(clone_names)))
#   names(cloneColours) <- clone_names
#   cloneColours[names(cloneColours) == 'A'] <- 'black'
#   cloneColours[names(cloneColours) == 'G'] <- '#E67F33'
#   cloneColours[names(cloneColours) == 'H'] <- '#D0B440'
#   
#   
#   
#   p2 <- fast_get_summary_tree(cell_clones = cc, 
#                               newick = '/Users/sohrabsalehi/projects/fitness_material/fitness_code_repo/fitclone/data/SA535_combined_CISPLATIN_CX5461_to_Sohrab/processed/sa535_cisplatin.newick', 
#                               cloneColours = cloneColours, edge_list = NULL, point_size_1 = 16, point_size_2 = 14)
#   
# }




