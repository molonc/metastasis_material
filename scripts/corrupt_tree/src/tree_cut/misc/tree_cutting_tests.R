require("RColorBrewer")
require('ggtree')
require('ggplot2')

sample_run <- function() {
  #--minimum-fraction  --maximum-fraction .3 /Users/sohrabsalehi/projects/tree_cutting/input_trees/tree_tyler.newick /Users/sohrabsalehi/projects/tree_cutting/ ./outputs
  # Test SA1091
  
  argv <- list()
  argv$newick <- './input_trees/tree_1091.newick'
  argv$copynumber <- './supplementary_input/tyler_intput/filtered_cell_cn_1091.tsv'
  argv$minimum_fraction <- .02
  argv$maximum_fraction <- .38
  argv$output <- './outputs/cell_cnv_1091.tsv'
  
  tree <- read.tree(argv$newick)
  g <- read_ltm_tree(tree_2_edge_list(tree))
  
  copynumber <- read.delim(
    argv$copynumber, check.names=FALSE, stringsAsFactors=FALSE)
  copynumber <- format_copynumber_matrix(copynumber)
  
  # TODO: export diagnostics, including the internal fist cut, and each of the indiv cuts.
  
  tree_split <- split_by_genotype_driver(
    g=g,
    minimum_fraction=argv$minimum_fraction,
    maximum_fraction=argv$maximum_fraction,
    mat=copynumber, 
    min_merge_back_fraction = .001
  )
  
  # Viz results
  fast_tree_aug(g = g, aug_cut = tree_split$new_aug, outdir = dirname(argv$output))
  

  cell_clones <- list_to_data_frame(tree_split$new_aug)
  cell_clones <- rename_clones(cell_clones)
  
  # Viz resulats
  aug_cut <- data_frame_to_list(cell_clones)  
  fast_tree_aug(g = g, aug_cut = aug_cut)
  
  
  write.table(
    cell_clones, argv$output, quote=FALSE, sep="\t", row.names=FALSE)
}

fast_tree_aug <- function(g, aug_cut, show_guide = TRUE, outdir = NULL, title = NULL) {
  label_offset = .8; angle=0; fontsize=8
  candiate_edges <- get_dummy_candid_edge(g, aug_cut)
  
  # Remove leaf loci
  degrees <- igraph::degree(graph = g, v = V(g)$name, mode = 'out')
  degree.dat <- data.frame(node_name = names(degrees), degree = unname(degrees), stringsAsFactors = F)
  leaf_loci <- degree.dat$node_name[degree.dat$degree == 0 & grepl('locus', degree.dat$node_name)]
  g <- igraph::delete.vertices(graph = g, v = leaf_loci)
  
  res <- graph_to_tree_dic(g)
  
  tree <- res$tree
  node_names <- res$node_names
  
  myColors <- get_cluster_colours(length(aug_cut))
  candiate_edges$colours <- myColors
  candiate_edges$letters <- candiate_edges$clone_id
  names(myColors) <- names(aug_cut)
  
  # TODO: Since 0 is shared with the heatmap CNV state, use something else for the Un-assigned ones
  myColors <- c(myColors, '0'='black')
  
  cdat <- data.frame(clone_id=names(myColors), colours=myColors, stringsAsFactors = F)
  
  # Group 0, the default of the tree should be black
  xtree <- groupOTU(tree, .node=aug_cut)
  p <- ggtree(xtree, aes(color=group), size = 1)
  

  
  # TODO: where are the f..ing clone names? Why aren't they shown 
  p <- p + scale_color_manual(values = myColors, guide=show_guide)
  
  qq <- p$data
  qq <- qq[qq$isTip == TRUE, ]
  
  compute_mean_clade_size <- function(candiate_edges, qq) {
    clade_size <- c()
    for (i in seq_along(candiate_edges$clone_id)) {
      the_key = candiate_edges$clone_id[[i]]
      sq = qq[qq$label %in% aug_cut[[the_key]], ]
      clade_size <- c(clade_size, max(sq$y) - min(sq$y))
    }
    min(clade_size)
  }
  
  ave_clade_size <- compute_mean_clade_size(candiate_edges, qq)
  
  for (i in seq_along(candiate_edges$clone_id)) {
    # Find the two ends of this pseudo-clade
    # TODO: handle the case where the two ends of the pseudo-clade run over another clade
    the_key = candiate_edges$clone_id[[i]]
    sq = qq[qq$label %in% aug_cut[[the_key]], ]
    max_node = sq$node[which.max(sq$y)]
    min_node = sq$node[which.min(sq$y)]
    
    # Find the edge of the left-most side of the tree-box 
    y_cent <- min(sq$y) + (max(sq$y) - min(sq$y))/2
    if (ave_clade_size > .05 * max(qq$y)) {
      ave_clade_size <- .05 * max(qq$y)
    }
    
    # where we want the label to be, i.e., close to its clade, but not overlapping the tree
    max_local_x <- max(qq$x[qq$y < (y_cent + ave_clade_size/2) & qq$y > (y_cent - ave_clade_size/2)]) 
    max_x <- max(qq$x) # Where the label will be shown by default
    
    #the_x_offset <- max_local_x - .9 * max_x
    the_x_offset <- max_local_x - max_x
    
    # Don't show the label with the heatmap
    clone_name_label <- paste0(candiate_edges$letters[[i]])
    clone_name_label <- paste0(clone_name_label, ' (', format(candiate_edges$frac[[i]], digits = 2), ')')
    
    
    p <- p + geom_strip(max_node, min_node, barsize=0, 
                        label = clone_name_label, 
                        color=candiate_edges$colours[[i]], 
                        offset=(label_offset + the_x_offset), fontsize = fontsize, angle = angle, hjust = 0, offset.text = 1)
  }
  
  if (!is.null(title)) {
    p <- p + ggtitle(label = title)
  }
  
  if (!is.null(outdir)) {
    dir.create(outdir, showWarnings = F, recursive = T)
    ggsave(filename = file.path(outdir, sprintf('sub_g_%s.png', generate_random_str())), plot = p)
  }

  p
}


annotate_nodes_on_tree <- function(p = NULL, tree_node_dic = NULL, edge_list_path, node_names, show_label=TRUE, leaf_shape=21, leaf_colour='steelblue', leaf_size=3, equal_branch=FALSE, use_pretty_names=F) {
  if (is.null(tree_node_dic)) {
    tree_node_dic <- convert_edge_list_to_ape(edge_list_path)
  }
  
  tree <- tree_node_dic$tree
  vertex_names <- tree_node_dic$node_names
  p_local <- ggtree(tree) + theme_tree2()
  # Get node labels, it is differently stored for tip names and internal nodes
  qq <- p_local$data
  tip_names <- node_names[node_names %in% qq$label[qq$isTip == TRUE]]
  tip_numbers <- qq$node[qq$label %in% tip_names]
  internal_nodes <- setdiff(node_names, tip_names)
  internal_numbers <- vertex_names$node_number[vertex_names$node_name %in% internal_nodes]
  
  # 1. annotate tips
  xtree <- groupOTU(tree, .node=list(tiplabs=tip_names))
  # Shape 16 (black circle), 23 (lozi), 21 (fillable circle)
  
  if (equal_branch) {
    p <- ggtree(xtree, branch.length="none")
  } else {
    p <- ggtree(xtree)
  }
  
  p <- p + geom_point2(aes(subset=(node %in% tip_numbers)), size=leaf_size, shape=leaf_shape, fill=leaf_colour, color=leaf_colour)
  
  if (show_label) p <- p + geom_tiplab(aes(subset=(group=='tiplabs')))
  
  # 2. annotate internanl nodes
  p <- p + geom_point2(aes(subset=(node %in% internal_numbers)), size=3, shape=23, fill="red") 
  p$data <- dplyr::left_join(p$data, vertex_names, by=c('node'='node_number'))
  #p <- p + geom_text2(aes(subset=(node %in% internal_numbers), label=node_name), hjust=-.3, color = 'red', position = position_jitter(height=150)) 
  if (use_pretty_names) {
    p$data$node_name <- get_pretty_names_for_loci(p$data$node_name)
  }
  p <- p + geom_text2(aes(subset=(node %in% internal_numbers), label=node_name), hjust=-.3, color = '#ED7953FF', angle = 45, face = "bold") 
  p  
  
  # How to annotate tips
  #p1 + geom_tiplab(aes(subset=(group=='baza')), align=TRUE, linesize = .01) + 
  # geom_point2(aes(subset=(node %in% nodes)), size=3, shape=21, fill="steelblue") + 
  #  theme_tree2() + 
  #  xlim(0, 400)
  p
}

generate_random_str <- function(n = 1) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}

generate_time_str <- function(n = 1) {
  #format(Sys.time(), "%Y%m%d %a %b %d %X %Y")
  format(Sys.time(), "%Y%m%d%H-%M-%S")
}


get_dummy_candid_edge <- function(g, aug_cut) {
  ntotal <- g.count.cells(g)
  frac <- unlist(lapply(names(aug_cut), function(x) length(aug_cut[[x]])/ntotal))
  data.frame(clone_id = names(aug_cut), frac = frac, stringsAsFactors = F)
}

graph_to_tree_dic <- function(the_graph) {
  edge_list <- as.data.frame(igraph::as_edgelist(the_graph), stringsAsFactors = F)
  colnames(edge_list) <- c('source', 'target')
  res <- convert_edge_list_to_ape(edge_list_mat = edge_list)
  res
}


convert_edge_list_to_ape <- function(edge_list_mat=NULL) {
  options(stringsAsFactors = F)
  edge_list <- as.data.frame(edge_list_mat)
  colnames(edge_list) <- c('source', 'target')
  
  leaf_names <- get_leaves_names(edge_list)
  root_name <- get_root_name(edge_list)
  all_nodes <- unique(c(edge_list$source, edge_list$target))
  
  internal_nodes <- setdiff(all_nodes, leaf_names)
  internal_nodes <- setdiff(internal_nodes, root_name)
  
  n_leaves <- length(leaf_names)
  n_internal <- length(internal_nodes) + 1
  
  d1 <- data.frame(node_name=leaf_names, node_number=seq(n_leaves))
  d2 <- data.frame(node_name=root_name, node_number=n_leaves+1)
  
  # Account for a start where the only loci is the root
  if (length(internal_nodes) == 0) {
    d <- rbind(d1, d2)
  } else {
    d3 <- data.frame(node_name=internal_nodes, node_number=seq(n_internal-1))
    d3$node_number <- d3$node_number + n_leaves + 1
    d <- rbind(d1, d2, d3)
  }
  
  edges <- edge_list
  
  for (node in all_nodes) {
    edges$source[edges$source == node] <- d$node_number[d$node_name == node]
    edges$target[edges$target == node] <- d$node_number[d$node_name == node]
  }
  
  edges$source <- as.numeric(edges$source)
  edges$target <- as.numeric(edges$target)
  
  edges = as.matrix(edges)
  colnames(edges) = NULL
  
  d <- d[order(d$node_number), ]
  tr <- list(edge = edges, tip.label = d$node_name[1:n_leaves], Nnode = n_internal)
  class(tr) <- "phylo"
  
  Nedge <- nrow(tr$edge)
  tr$edge.length <- rep(1, Nedge)
  list(tree=tr, node_names=d)
}


get_leaves_names <- function(edge_list) {
  leaves <- unique(edge_list$target)[!(unique(edge_list$target) %in%  unique(edge_list$source))]
  target_freq <- as.data.frame(table(edge_list$target))
  stopifnot(all (target_freq$Freq[target_freq$Var1 %in% leaves] == 1))
  leaves
}

get_root_name <- function(edge_list) {
  unique(edge_list$source)[!(unique(edge_list$source) %in% unique(edge_list$target))]
}

get_cluster_colours <- function(nClusters) {
  if (nClusters > 8) {
    clust_colours <- colorRampPalette(brewer.pal(8, "Set2"))(nClusters)
  } else {
    clust_colours <- brewer.pal(nClusters, "Set2")
  }
  clust_colours
}

setup_graph <- function(g=NULL) {
  V(g)$id <- seq(vcount(g))
  h <- get_height_dat(g)
  stopifnot(all(h$id == V(g)$name))
  V(g)$height <- as.integer(h$dist)
  g
}