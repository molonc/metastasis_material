suppressPackageStartupMessages({
  require("tools")
  require("phytools")
  require("igraph")
  require("RColorBrewer")
  require("outliers")
  require("data.table")
  require("dplyr")
  # require("optparse")
  require('tidyr')
  require('tibble')
  require('ggplot2')
  require('cowplot')
  require('ggtree')
  library("yarrr")
  library('ggrepel')
})

# option_list <- list(make_option(c("-i", "--input"), type="character", default=NULL, help="filtered_cnv_input", metavar="character"),
#                     make_option(c("-o", "--outfile"), type="character", default=NULL, help="tip_probability_output", metavar="character"),
#                     make_option(c("-f", "--filterthrs"), type="character", default='0.9', help="outlier_threshold", metavar="character")
# )
# 
# opt_parser <- OptionParser(option_list=option_list)
# opt <- parse_args(opt_parser)


## Constants
# 1. Check if chromosome linings are included X
# 2. Remove the bins spanning multiple bins
get_single_feature_bins <- function(dat, mat_delta) {
  #original_bins <- rownames(sort_mat_by_bins(load_new_cn_data(datatag))) 
  original_bins <- row.names(dat)
  delta_bin_names <- row.names(mat_delta)
  delta_bin_dat <- data.frame(delta_bin_names = delta_bin_names, left_bin = original_bins[-c(length(original_bins))], right_bin = original_bins[-c(1)], stringsAsFactors = F)
  stopifnot(all(delta_bin_dat$delta_bin_names == delta_bin_dat$left_bin))
  
  cnv_txt_left <- parse_bin_names(delta_bin_dat$left_bin)
  cnv_txt_right <- parse_bin_names(delta_bin_dat$right_bin)
  
  delta_bin_dat$dist <- cnv_txt_right$start -  cnv_txt_left$start
  # Set to inf for those that span multiple chromosomes
  delta_bin_dat$dist <- abs(ifelse(cnv_txt_left$chr == cnv_txt_right$chr, 1, Inf) * delta_bin_dat$dist)
  
  stopifnot(min(delta_bin_dat$dist, na.rm = T) == 500000)
  
  # Keep bins that only cover 500K bases
  single_feature_bins <- delta_bin_dat$delta_bin_names[delta_bin_dat$dist == 500000]
  length(single_feature_bins)
  # nrow(delta_bin_dat)
  single_feature_bins
}


get_pretty_names_for_loci <- function(str_array) {
  str_array <- gsub('locus_', '', str_array)
  chr = gsub('([0-9]+|X|Y)_.*', '\\1', str_array)
  str_array <- gsub('.*_([0-9]+)_.*', '\\1', str_array)
  str_array <- substr(str_array, 1, 4)
  str_array <- paste0('chr', chr, '_', str_array, '')
  gsub('chrroot_root', 'root', str_array)
}

addrow <- function(orig, tmp) {
  if (is.null(orig)) {
    orig <- tmp
  } else {
    orig <- rbind(orig, tmp)
  }
  orig
}


format_copynumber_matrix <- function(copynumber) {
  # bin_ids <- paste(
  #   copynumber$chr, copynumber$start, copynumber$end, sep="_")
  
  # rownames(copynumber) <- bin_ids
  # copynumber <- subset(copynumber, select=-c(chr, start, end, width))
  copynumber <- as.matrix(copynumber)
  
  return(copynumber)
}


tree_2_edge_list <- function(tree) {
  tree$node.label[1] <- "root"
  node_names <- c(tree$tip.label, tree$node.label)
  
  edges <- tree$edge
  edges <- data.frame(
    source=node_names[edges[, 1]],
    target=node_names[edges[, 2]],
    stringsAsFactors=FALSE
  )
  edges$source <- gsub("cell_", "", edges$source)
  edges$target <- gsub("cell_", "", edges$target)
  
  return(edges)
}


fitclone_plot_tree_heatmap <- function(g, cell_clones = NULL, mat = mat, tag='', output=NULL, just_tree_loci=TRUE, chr_filter = c(1,2), ...) {
  fitclone_plot_tree_heatmap_more(g = g, cell_clones = cell_clones, mat = mat, tag = tag, output = output, just_tree_loci=just_tree_loci, chr_filter = chr_filter, ...)
}

fitclone_plot_tree_heatmap_more <- function(g, cell_clones, mat = mat, tag='', output=NULL, aug_cut=NULL, chr_filter=NULL, just_tree_loci=NULL, chr_font_size = 3.8, ...) {
  all_cells <- g.count.cells(g)
  covered_cells <- nrow(cell_clones)
  title_str <- sprintf('%d/%d', covered_cells, all_cells)
  
  print('Generating heatmap object...')
  pp <- tree_and_heatmap(g = g, cell_clones = cell_clones, mat = mat, chr_filter=chr_filter, just_tree_loci = just_tree_loci, chr_font_size = chr_font_size, ...)
  # pp <- tree_and_heatmap(g, cell_clones, mat, chr_filter, just_tree_loci, chr_font_size = chr_font_size)
  
  #p_res <- pp
  
  pp <- pp + 
    ggtitle(title_str) + 
    theme(legend.text=element_text(size=15, angle = 45), 
          plot.title = element_text(size = 28), 
          plot.subtitle=element_text(size=20, face="italic", color="black"))
  
  print('Saving the heatmap object...')
  coef = 1
  height = 8*coef
  #width = 14 * coef
  width = 20 * coef
  ggsave(plot = pp, filename = output,  width = width, height = height, units = "in", limitsize = FALSE)
  print(output)
}

tree_and_heatmap <- function(g, cell_clones, mat, just_tree_loci = FALSE, chr_filter = NULL, chr_font_size = 3.8, drop_loci_names = FALSE) {
  if (!is.null(chr_filter)) {
    chr <- gsub('([0-9]+|X|Y)_[0-9]+_[0-9]+', '\\1', rownames(mat))
    mat <- mat[chr %in% chr_filter, ]
  }
  
  mat <- sort_mat_by_bins(the_mat = mat)  
  #mat <- remove_pad_bins(mat)
  
  mat <- mat[, colnames(mat) %in% cell_clones$cell_id]
  mat <- mat[, match(cell_clones$cell_id, colnames(mat))]  
  
  
  # The tree
  tree_res <- convert_edge_list_to_ape(edge_list_mat = igraph::as_edgelist(g))
  
  tree <- tree_res$tree
  node_names <- tree_res$node_names
  the_height <- find_tail_height(g = g)$h_cut
  
  if (!is.null(the_height)) {
    tree_res <- trim_tree_before_height(tree, node_names, g, the_height, 20)
    tree <- tree_res$tree
  }
  
  # Add the heatmap
  sub_mat <- heatmap_mat_from_bin_cellID(mat, update_colnames = FALSE)
  
  if (just_tree_loci) {
    the_loci <- grep('locus', node_names$node_name, value = T)
    the_loci <- gsub('locus_', '', the_loci)
    if (!any(colnames(sub_mat) %in% the_loci)) {
      stop('Non of the selected chromosoms are in the tree_loci.')
    }
    sub_mat <- sub_mat[, colnames(sub_mat) %in% the_loci]
  }
  
  
  colnames_filter <- strip_chr_names(colnames(sub_mat))
  ss <- which(colnames_filter != '')
  colnames_filter <- data.frame(from = colnames(sub_mat)[ss], new_label=colnames_filter[ss])
  
  # Add white borders between chromosomes
  genotype <- as_tibble(sub_mat)
  white_border_size <- 10
  if (!is.null(just_tree_loci)) white_border_size <- 10
  
  if (white_border_size > 0) {
    for (cname in colnames_filter$from[-c(1)]) {
      for (i in 1:white_border_size) {
        new_name <- gsub('([0-9]+|X|Y)_(.*)', paste0('\\1_', i-1, '_\\2'), cname)
        genotype <- add_column(genotype, !!new_name := NA, .before = cname)
      }
    }
  }
  
  # Convert to character format
  d <- dim(genotype)
  cn <- colnames(genotype)
  rn <- rownames(sub_mat)
  genotype <- as.matrix(genotype)
  genotype <- as.character(genotype)
  genotype <- matrix(data=genotype, nrow = d[1], ncol = d[2])
  colnames(genotype) <- cn
  rownames(genotype) <- rn
  
  colnames_level <- unique(colnames(genotype))
  
  # Add colours for CNV value
  cn_colours <- legacy_colours_for_CN_heatmap()
  names(cn_colours) <- paste0(seq_along(cn_colours)-1)
  
  # Add extra colours for the input to corrupt (filtered.csv): shades of green, black, darkblue
  cn_colours <- c(cn_colours, "40" = 'black', "41" = 'darkgreen')
  
  # Add extra colours for the posterior of corrupt (average.csv): shades of blue
  pos_valuse <- get_corrupt_posterior_values()
  pos_colours  <- colorRampPalette(brewer.pal(8, "Greys"))(len(pos_valuse))
  names(pos_colours) <- pos_valuse
  cn_colours <- c(cn_colours, pos_colours)
  
  
  # Show library_ids for annotation for now
  # colAnnotMat = NULL, colAnnotWidthMult = 1) 
  lib_ids <- get_libid_from_cell_names(rownames(genotype))
  colAnnotMat <- matrix(data = lib_ids, nrow = length(rownames(genotype)), ncol = 1)
  rownames(colAnnotMat) <- rownames(genotype)
  
  # Augment the colours for library ids
  lib_colours <- colorRampPalette(colors = piratepal(palette = "basel"))(lenu(lib_ids))
  names(lib_colours) <- unique(lib_ids)
  cn_colours <- c(cn_colours, lib_colours)
  
  tip.loci <- grep('locus', tree$tip.label, value= T)
  tree <- drop.tip(tree, tip.loci)
  
  ## Plot
  aug_cut <- data_frame_to_list(cell_clones)
  
  
  xtree <- ggtree::groupOTU(tree, aug_cut)
  p <- ggtree(xtree, aes(color=group))
  
  
  # Annotate loci names 
  the_loci <- NULL
  if (!drop_loci_names) {
    the_loci <- grep('locus', node_names$node_name, value = T)
    # Just keep the loci shown in the matrix
    if (!is.null(chr_filter)) {
      lchr <- gsub('([0-9]+|X|Y)_[0-9]+_[0-9]+', '\\1', gsub('locus_', '', the_loci))
      the_loci <- the_loci[lchr %in% chr_filter]
    }
    #the_loci <- the_loci[gsub('locus_', '', the_loci) %in% colnames(sub_mat)]
    p <- annotate_nodes_on_tree(p = p, 
                                edge_list_path = NULL, 
                                node_names = the_loci, 
                                tree_node_dic = tree_res, 
                                show_label = F)
    # For the heatmap
    # For enhanced loci, use their end coordinates
    # Find them and replace their start with the one that exists in the mat
    lbins <- parse_bin_names(the_loci)
    lbins$start <- lbins$start + 2
    tmp_loci <- paste0('locus_', collapse_bin_names(lbins))
    if (any(tmp_loci %in% colnames(sub_mat))) {
      # We're showing the enhanced bins
      the_loci <- tmp_loci
    }
    the_loci <- data.frame(from = gsub('locus_', '', the_loci), new_label = the_loci)
  }
  
  
  
  
  
  # which(levels(mapping_l$from) == '1_101166670_101166668')
  offset <- 10
  breaks <- factor(x = names(cn_colours), levels = names(cn_colours), ordered = T)
  # The names in the legend
  labels <- c(as.character(1:11-1), '11+', 'F0', 'F1', 
              '0.0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0', 
              tail(levels(breaks), lenu(lib_ids)))
  
  p <- sos_heat(p = p, 
                data = genotype, 
                offset = offset, 
                width = 8, 
                colnames = TRUE, 
                color = NULL, 
                colnames_level = colnames_level, 
                colnames_filter = colnames_filter, 
                loci_filter = the_loci,
                font.size = chr_font_size, 
                colnames_offset_y = -40, 
                colAnnotMat = colAnnotMat, 
                colAnnotWidthMult = max(2, round(.1*ncol(mat))))
  p <- p + scale_fill_manual(breaks=breaks, 
                             values=cn_colours, 
                             labels = labels, 
                             guide = guide_legend(nrow = 1, 
                                                  byrow = TRUE, 
                                                  direction='horizontal', 
                                                  label.position = 'bottom', 
                                                  label.theme = element_text(angle = 90)))
  
  # Add column-annotation annotations
  pdat <- p$layers
  #p + geom_tile()
  
  
  p <- p + theme(legend.position="top")
  # Add the cluster annotations
  p <- annoate_tree_aug(p = p, g = g, candiate_edges = NULL, aug_cut = aug_cut, the_height = the_height, show_guide = F, show_label = F)
  #p2 <-  annoate_tree_aug(p = NULL, g = g, candiate_edges = NULL, aug_cut = aug_cut, the_height = the_height, show_guide = F, show_label = F)
  p    
}





### Matrix functions
################################################################################################################################################
sort_mat_by_bins <- function(the_mat) {
  # prevent scientific notation
  options(scipen=999)
  options(stringsAsFactors=FALSE) 
  cnv_txt <- parse_bin_names(rownames(the_mat), as_factors = F)
  the_mat <- cbind(cnv_txt, the_mat)
  
  # Sort the matrix by their chromosome (just inside the chromosome)
  the_mat$chr[the_mat$chr == 'X'] <- '40' 
  the_mat$chr[the_mat$chr == 'Y'] <- '55' 
  the_mat$chr <- as.numeric(the_mat$chr)
  the_mat <- the_mat[order(the_mat$chr, the_mat$start), ]
  the_mat$chr[the_mat$chr == '40'] <- 'X' 
  the_mat$chr[the_mat$chr == '55'] <- 'Y' 
  
  # Remove chr, start, end
  the_mat$chr <- NULL
  the_mat$start <- NULL
  the_mat$end <- NULL
  
  the_mat
}

remove_pad_bins <- function(the_mat) {
  pad_index <- get_pad_bin_index_for_mat(the_mat = the_mat)
  if (length(pad_index) > 0)
    the_mat <- the_mat[-pad_index, ]
  
  the_mat
}

get_pad_bin_index_for_mat <- function(the_mat) {
  dd <- dim(the_mat)
  if (max(dd)/min(dd) > 1000) {
    bin_names <- unique(the_mat$loci)
  } else {
    bin_names <- rownames(the_mat)
  }
  
  get_pad_bin_index(bin_names)
}

get_pad_bin_index <- function(bin_names) {
  bin_dat <- parse_bin_names(bin_names)
  which(abs(bin_dat$start - bin_dat$end) == 1)
}

### Tree
################################################################################################################################################
ltm_clust_to_cuttree <- function(ltm_clust, clone_ids) {
  df <- ltm_clust_to_df(ltm_clust, clone_ids)
  res <- df$genotype
  names(res) <- df$single_cell_id
  res
}

ltm_clust_to_df <- function(all_clones, clone_ids) {
  df <- NULL
  for (i in seq(length(all_clones))) {
    clone_name <- names(all_clones)[[i]]
    temp <- data.frame(single_cell_id=all_clones[[clone_name]], genotype=clone_ids[[i]])  
    if (is.null(df))
      df <- temp
    else
      df <- rbind(df, temp)
  }
  df
}

data_frame_to_list <- function(clustering) {
  clusts <- unique(clustering$clone_id)
  res <- list()
  for (cl in clusts) {
    res[[cl]] <- clustering$cell_id[clustering$clone_id == cl]
  }
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
  tr$edge.length <- rep(4, Nedge)  # replace 1 by 4
  
  list(tree=tr, node_names=d)
}

find_tail_height <- function(g = NULL) {
  h <- get_height_dat(g)
  res <- convert_edge_list_to_ape(edge_list_mat = igraph::as_edgelist(g))
  
  tree <- res$tree
  node_names <- res$node_names
  p <- ggtree(tree) + theme_tree2()
  qq <- p$data
  qq <- qq[qq$isTip == TRUE, ]
  
  # Ignore SA928 gm cells
  qq <- qq[!grepl('SA928', qq$label), ]
  uml <- qq$label[which.max(qq$y)]
  lml <- qq$label[which.min(qq$y)]
  mrca <- find_MRCA(the_g = g, the_nodes = c(uml, lml), height_dat = h)
  mrca_node <- node_names$node_number[node_names$node_name == mrca]
  tail_length <- h$dist[h$id == mrca]
  
  result <- list(h_cut = NULL, tail_length = tail_length, nodes = c(uml, lml, mrca))
  
  if (tail_length > 15)
    result$h_cut <- tail_length - 5
  
  result  
}

graph_to_tree_dic <- function(the_graph) {
  edge_list <- as.data.frame(igraph::as_edgelist(the_graph), stringsAsFactors = F)
  colnames(edge_list) <- c('source', 'target')
  res <- convert_edge_list_to_ape(edge_list_mat = edge_list)
  res
}

find_MRCA <- function(the_g, the_nodes, height_dat) {
  if (the_nodes[1] == the_nodes[2]) return(the_nodes[1])
  if (height_dat$id[1] == 'root' | (height_dat$dist[1] < height_dat$dist[nrow(height_dat)])) {
    print("Warning! Sorting height dat in decreasing order...")
    height_dat <- height_dat[order(height_dat$dist, decreasing = T), ]
  }
  parents <- list()
  for (node in the_nodes) {
    print(node)
    i = which(height_dat$id == node)
    # Append the node itself, so if it's a parent of the others, it's taken care of
    parents[[node]] <- c(find_all_parents(the_g, height_dat$id[i], height_dat$dist[i]), height_dat$id[i])
  }
  # TODO: support multiple subgroups
  cp <- intersect(parents[[1]], parents[[2]])
  
  height_dat$id[height_dat$id %in% cp][1]
}

trim_tree_before_height <- function(tree, node_names, g, the_height, after_height=NULL) {
  h <- get_height_dat(g)
  max_height <- max(h$dist)
  h <- h[h$dist < the_height, ]
  
  # Find edge row_id in ape format
  qq = tree$edge
  pp = node_names[order(node_names$node_number), ]
  q1 = pp$node_name[qq[, 1]]
  q2 = pp$node_name[qq[, 2]]
  named_edges = matrix(c(q1, q2), ncol = 2, byrow = FALSE)
  named_edges = as.data.frame(named_edges)
  colnames(named_edges) = c('source', 'target')
  tree_edge_ids = which((named_edges$source %in% h$id) & (named_edges$source %in% h$id))
  
  # tree$edge.length[tree_edge_ids] = .01
  tree$edge.length[tree_edge_ids] = .09
  if (!is.null(after_height)) {
    tree$edge.length[setdiff(seq(tree$edge.length), tree_edge_ids)] = after_height
  }
  
  list(tree=tree, height=max_height, node_names=node_names)
}

### Heatmap
################################################################################################################################################
heatmap_mat_from_bin_cellID <- function(res_pick, update_colnames=TRUE) {
  the_mat = t(res_pick)
  the_mat <- prepare_mat_by_chr(the_mat, update_colnames)
  if (length(which(unname(is.na(the_mat[2, ])))) > 0)
    the_mat = the_mat[, -c(which(unname(is.na(the_mat[2, ]))))]
  the_mat
}

prepare_mat_by_chr <- function(the_mat, update_col_names=TRUE) {
  the_array = colnames(the_mat)
  chunks <- unlist(strsplit(the_array, '_'))
  chrs <- matrix(chunks, nrow=length(the_array), ncol=3, byrow = T)[, 1]
  unique(chrs)
  the_mat <- the_mat[, order_by_chr_name(chrs)]
  if (update_col_names) {
    colnames(the_mat) <- strip_chr_names(the_array = colnames(the_mat))
  } 
  
  the_mat
}

legacy_colours_for_CN_heatmap <- function() {
  # From CN = 0 to CN = 11
  c('#4880B8', '#A7C9DF','#CCCCCC','#F5CE93','#ED9364','#D2553E','#A42116','#8B1A43','#CB3576','#D06CAD','#C196C4','#D0BAD8')
}

sos_heat <- function (p, data, offset = 0, width = 1, low = "green", high = "red", 
                      color = "white", colnames = TRUE, colnames_position = "bottom", 
                      colnames_angle = 0, colnames_level = NULL, colnames_offset_x = 0, 
                      colnames_offset_y = 0, font.size = 4, hjust = 0.5, colnames_filter = NULL, loci_filter = NULL,
                      colAnnotMat = NULL, colAnnotWidthMult = 1) 
{
  colnames_position <- colnames_position %>% match.arg(c("bottom", "top"))
  variable <- value <- lab <- y <- NULL
  width <- width * (p$data$x %>% range(na.rm = TRUE) %>% diff)/ncol(data)
  isTip <- x <- y <- variable <- value <- from <- to <- NULL
  df <- p$data
  df <- df[df$isTip, ] # for the tree
  start <- max(df$x, na.rm = TRUE) + offset
  dd <- as.data.frame(data) # cnv matrix
  i <- order(df$y)
  i <- i[!is.na(df$y[i])]
  lab <- df$label[i]
  dd <- dd[lab, , drop = FALSE]
  dd$y <- sort(df$y)
  dd$lab <- lab
  dd <- gather(dd, variable, value, -c(lab, y))
  i <- which(dd$value == "")
  if (length(i) > 0) {
    dd$value[i] <- NA
  }
  if (is.null(colnames_level)) {
    dd$variable <- factor(dd$variable, levels = colnames(data))
  } else {
    dd$variable <- factor(dd$variable, levels = colnames_level)
  }
  
  # Manually adding x coordinates for the heatmap
  V2 <- start + as.numeric(dd$variable) * width
  mapping <- data.frame(from = dd$variable, to = V2)
  mapping <- unique(mapping)
  dd$x <- V2
  dd$width <- width
  
  
  function() {
    # dummy colAnnotMat
    colAnnotMat <- matrix(NA, nrow = lenu(dd$y), ncol = 2)
    rownames(colAnnotMat) <- unique(dd$lab)
    colAnnotMat[, 1] <- '0'
    colAnnotMat[, 2] <- '1'
  }
  
  # Push all x by the width of the annotations
  nColumnAnnotations <- 0
  if (!is.null(colAnnotMat)) {
    nColumnAnnotations <- ncol(colAnnotMat)
  }
  dd$x <- dd$x + (2*nColumnAnnotations+1)*colAnnotWidthMult*width
  # Update column names (chr)
  mapping$to <- mapping$to + (2*nColumnAnnotations+1)*colAnnotWidthMult*width
  
  function() {
    # They have to line up with the loci, i.e., the columns of the heatmap
    # dummy rowAnnotMat
    rowAnnotMat <- matrix(NA, nrow = 2, ncol = lenu(dd$x))
    colnames(rowAnnotMat) <- unique(dd$lab)
    rowAnnotMat[, 1] <- '0'
    rowAnnotMat[, 2] <- '1'
  }
  
  
  if (is.null(color)) {
    p2 <- p + geom_tile(data = dd, aes(x, y, fill = value), width = width, inherit.aes = FALSE)
  } else {
    p2 <- p + geom_tile(data = dd, aes(x, y, fill = value), width = width, color = color, inherit.aes = FALSE)
  }
  if (is(dd$value, "numeric")) {
    p2 <- p2 + scale_fill_gradient(low = low, high = high, na.value = NA)
  } else {
    p2 <- p2 + scale_fill_discrete(na.value = NA)
  }
  
  if (colnames_offset_x == 0.0) {
    colnames_offset_x <- -width*.5
  }
  
  if (colnames) {
    if (colnames_position == "bottom") {
      y <- 0
    } else {
      y <- max(p$data$y) + 1
    }
    # save a temp for the filters
    mapping_l <- mapping
    mapping$y <- y
    if (!is.null(colnames_filter)) {
      mapping <- mapping[mapping$from %in% colnames_filter$from, ]
      mapping <- dplyr::right_join(mapping, colnames_filter)
      mapping$from <- mapping$new_label
    }
    
    p2 <- p2 + geom_text(data = mapping, aes(x = to, y = y, label = from), size = font.size, inherit.aes = FALSE, 
                         angle = colnames_angle, nudge_x = colnames_offset_x, 
                         nudge_y = colnames_offset_y, hjust = hjust)
    
    # @Sohrab: add the loci name here
    if (!is.null(loci_filter)) {
      mapping <- mapping_l
      mapping$y <- max(p$data$y) + 1 - y
      
      mapping <- mapping[mapping$from %in% loci_filter$from, ]
      mapping <- dplyr::right_join(mapping, loci_filter)
      mapping$from <- get_pretty_names_for_loci(mapping$new_label)
      # p2 <- p2 + geom_text(data = mapping, aes(x = to, y = y, label = from), size = font.size, inherit.aes = FALSE, 
      #                      angle = 45, nudge_x = colnames_offset_x, 
      #                      nudge_y = -colnames_offset_y, hjust = hjust)
      p2 <- p2 + geom_text_repel(data = mapping, 
                                 mapping = aes(x = to, y = y, label = from), 
                                 colour = 'darkblue', 
                                 size = font.size, 
                                 nudge_y = -colnames_offset_y,
                                 vjust = 1,
                                 nudge_x = 0, angle = 90, 
                                 fontface = 'bold',
                                 segment.color = 'black',
                                 segment.size = .2, 
                                 inherit.aes = FALSE)
      p2 <- p2 + scale_y_continuous(expand = c(.2,0))
    }
  }
  
  # q2 = p2
  p2 <- p2 + theme(legend.position = "right", legend.title = element_blank())
  attr(p2, "mapping") <- mapping
  
  
  # Add column annotation (for libs and samples, etc)
  p2 <- 
    add_column_annot(p = p2, dd = dd, colAnnotMat = colAnnotMat, colAnnotWidthMult = colAnnotWidthMult)
  
  p3 <- p2
  #add_row_annot(p = p2, dd = dd, rowAnnotMat = rowAnnotMat, rowAnnotHeightMult = rowAnnotHeightMult)
  
  return(p3)
}

# colAnnotMat: matrix of nColumnAnnotations columns & rownames set to cells: values and 
add_column_annot <- function(p, dd, colAnnotMat, colAnnotWidthMult) {
  nColumnAnnotations <- ncol(colAnnotMat)
  offset <- (2*nColumnAnnotations + 1)
  width <- max(dd$width)
  for (cl in 1:nColumnAnnotations) {
    colIndex <- cl*2
    ds <- dd
    ds$value <- NULL
    ds <- dplyr::left_join(ds, data.frame(lab = rownames(colAnnotMat), value = colAnnotMat[, cl]))
    ds$x <- ds$x - (offset - colIndex + 1)*colAnnotWidthMult*width
    # Keep one variable for each label
    ds <- ds[ds$variable == levels(ds$variable)[1], ] 
    
    p <- p + geom_tile(data = ds, aes(x, y, fill = value), width = width*colAnnotWidthMult, inherit.aes = FALSE)
  }
  
  p
}


add_row_annot <- function(p, dd, rowAnnotMat, rowAnnotHeightMult) {
  nRowAnnotations <- ncol(rowAnnotMat)
  offset <- (2*nRowAnnotations + 1)
  height <- max(dd$height)
  for (cl in 1:nRowAnnotations) {
    rowIndex <- cl*2
    ds <- dd
    ds$value <- NULL
    ds <- dplyr::left_join(ds, data.frame(lab = rownames(rowAnnotMat), value = rowAnnotMat[, cl]))
    ds$x <- ds$x - (offset - colIndex + 1)*rowAnnotHeightMult*height
    # Keep one variable for each label
    ds <- ds[ds$variable == levels(ds$variable)[1], ] 
    
    p <- p + geom_tile(data = ds, aes(x, y, fill = value), height = height*rowAnnotHeightMult, inherit.aes = FALSE)
  }
  
  p
}


collapse_tree_by_edge <- function(the_g, nodes, do_plot=FALSE, annotate_nodes=NULL) {
  h <- get_height_dat(the_g)
  h <- h[h$id %in% nodes, ]
  h <- h[order(h$dist, decreasing = T), ]
  
  # The neighbour with the smallest height is the parent
  edge_list <- matrix(NA, nrow=nrow(h), ncol=2, dimnames = list(c(), c('source', 'target')))
  for (i in seq(nrow(h))) {
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
  if (do_plot) {
    out_path = file.path('~/projects/fitclone/presentations/August28/edge_cut/', paste0(datatag, '_clone_tree.png'))
    plot_graph(ft = as.matrix(edge_list), out_path = out_path, annotate_nodes = get_pretty_names_for_loci(annotate_nodes))
    print(sprintf('Results in %s', out_path))
  }
  
  gg <- igraph::graph_from_edgelist(el = orig_edge_list)
  V(gg)$id <- seq(vcount(gg))
  list(edge_list=orig_edge_list, graph=gg)
}


# TODO: rename these
fitclone_get_theme_no_grid <- function() {
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  # axis.line = element_line(colour = "black"), 
  # panel.background = element_rect(fill = NA, color = "black")) 
}

fitclone_get_theme_no_legend <- function() {
  theme(legend.position = 'none')
}

fitclone_get_theme_no_axis <- function() {
  theme(axis.line = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.title = element_blank(), 
        axis.ticks = element_blank())
}


legend_plot <- function(clone_dic, vertical = FALSE, font.size = 7, point_size = 5, main = NULL) {
  if (vertical) {
    ll.dat <- data.frame(x = c(rep(1, length(clone_dic$letters))), y = rev(seq_along(clone_dic$letters)), label = clone_dic$letters, stringsAsFactors = F)
  } else {
    ll.dat <- data.frame(x = seq_along(clone_dic$letters), y = c(rep(1, length(clone_dic$letters))), label = clone_dic$letters, stringsAsFactors = F)
  }
  
  # Better results with spaces
  ll.dat$label <- paste0('   ', ll.dat$label)
  clone_dic$letters <- paste0('   ', clone_dic$letters)
  
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
    geom_text(aes(label = label), colour = 'black', nudge_x = nudge_x, nudge_y = nudge_y, fontface = "bold", size = font.size, hjust = 0) + 
    fitclone_get_theme_no_grid() + 
    fitclone_get_theme_no_axis() + 
    fitclone_get_theme_no_legend() + 
    scale_color_manual(values = myColors) + 
    theme(panel.background = element_blank())
  if (point_size > 0) {
    ll.dat.point <- ll.dat[-c(1), ]
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
  }
  
  pg
}




fast_tree_aug <- function(g, aug_cut, show_guide = TRUE) {
  label_offset = .8; angle=0; fontsize=8
  candiate_edges <- get_dummy_candid_edge(g, aug_cut)
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
  p
}


annoate_tree_aug <- function(p=NULL, g=NULL, candiate_edges, aug_cut, clone_dic=NULL, the_height, label_offset = .8, angle=0, fontsize=8, draw_bar=FALSE, show_guide=FALSE, show_label=TRUE, branch_width = 1) {
  if (is.null(candiate_edges)) {
    candiate_edges <- get_dummy_candid_edge(g, aug_cut)
    print(candiate_edges)
  }
  
  # A hack to make colouring work
  orig_p <- p
  
  # TODO: Sohrab, sort this out... - it does not convert to newick with loci name intact
  res <- graph_to_tree_dic(g)
  
  tree <- res$tree
  node_names <- res$node_names
  
  
  # Cut the tail
  if (!is.null(the_height)) {
    tree_res <- trim_tree_before_height(tree, node_names, g, the_height)
    tree <- tree_res$tree
    node_names <- tree_res$node_names
  }
  
  myColors <- get_cluster_colours(length(aug_cut))
  names(myColors) <- names(aug_cut)
  
  # TODO: Since 0 is shared with the heatmap CNV state, use something else for the Un-assigned ones
  myColors <- c(myColors, '0'='black')
  
  cdat <- data.frame(clone_id=names(myColors), colours=myColors, stringsAsFactors = F)
  
  # Group 0, the default of the tree should be black
  if (is.null(p)) {
    xtree <- ggtree::groupOTU(tree, aug_cut)
    p <- ggtree(xtree, aes(color=group), size = branch_width)
  }
  
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
    if (is.null(orig_p)) {
      clone_name_label <- paste0(clone_name_label, ' (', format(candiate_edges$frac[[i]], digits = 2), ')')
    }
    
    if (!draw_bar) {
      if (show_label) {
        p <- p + geom_strip(max_node, min_node, barsize=0, 
                            label = candiate_edges$letters[[i]], 
                            color = candiate_edges$colours[[i]], 
                            offset = (label_offset + the_x_offset), fontsize = fontsize, angle = angle, hjust = 0, offset.text = -.5)
        
        p <- p + geom_strip(max_node, min_node, barsize=0, 
                            label = sprintf('(%.2f)', candiate_edges$frac[[i]]),
                            color = candiate_edges$colours[[i]], 
                            offset = (label_offset + the_x_offset + 5), fontsize = fontsize-3, angle = angle, hjust = 0, offset.text = -.5)
      }
    } else {
      p <- p + geom_strip(max_node, min_node, barsize=2, 
                          label = clone_name_label, 
                          color=candiate_edges$colours[[i]], 
                          offset=(label_offset + the_x_offset), fontsize = fontsize, angle = angle, hjust = 0, offset.text = 1)
    }
    
  }
  
  
  if (!show_label) {
    # Add clone counts
    ntotal <- g.count.cells(g)
    counts <- unlist(lapply(names(aug_cut), function(x) length(aug_cut[[x]])))
    clone_counts <- data.frame(clone_id = names(aug_cut), counts = counts, frac = counts/ntotal, stringsAsFactors = F)
    
    # Add a legend plot
    mode_clone_dic <- dplyr::left_join(clone_counts, cdat)
    mode_clone_dic$color <- mode_clone_dic$colours
    
    mode_clone_dic$letters <- sprintf('%s (n = %d, %.1f%%)', mode_clone_dic$clone_id, mode_clone_dic$counts, mode_clone_dic$frac*100)
    main <- sprintf('Total cells = %d', sum(mode_clone_dic$counts))
    
    # Add as another row
    mode_clone_dic <- dplyr::bind_rows(data.frame(letters = main, color = 'black'), mode_clone_dic)
    
    pg <- legend_plot(clone_dic = mode_clone_dic, vertical = TRUE, font.size = 3, point_size = 2)
    pg <- pg + theme(plot.subtitle = element_text(size = 10), 
                     plot.margin = margin(0,0,0,0), 
                     #rect = element_rect(fill = "transparent", colour = "transparent"), 
                     rect = element_blank(),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.grid = element_blank(),
                     axis.line = element_blank(),
                     panel.border = element_blank())
    
    p <- 
      ggdraw() +
      draw_plot(p) +
      draw_plot(pg, x = -1.2, y = .3, width = 2.4, scale = .45)
  }
}



### Util_utils
################################################################################################################################################


get_ploidy_for_mat <- function(the_mat) {
  dmat <- wide_to_tidy_for_mat(the_mat)
  pa <- dmat %>% dplyr::count(cell_names, copy_number) %>% dplyr::group_by(cell_names) %>% dplyr::filter(n == max(n)) %>% dplyr::rename(ploidy = copy_number) %>% dplyr::select(cell_names, ploidy) %>% as.data.frame()
  nrow(pa)
  dim(the_mat)
  pa$cell_names <- as.character(pa$cell_names)
  pa
}

wide_to_tidy_for_mat <- function(mat, chr_filter = NULL) {
  cnv_txt <- parse_bin_names(rownames(mat))
  mat <- cbind(cnv_txt, mat)
  if (!is.null(chr_filter))
    mat <- mat[mat$chr %in% chr_filter, ]
  reshape2::melt(mat, id.vars = c('chr', 'start', 'end'), value.name='copy_number', variable.name = 'cell_names')
}


get_candi_edge <- function(the_cut, aug_cut, ts, tcc, new_clades) {
  if (is.null(aug_cut) & !is.null(the_cut)) {
    candi_edge <- ts[ts$clone_id %in% unlist(the_cut), ]
  } else {
    candi_edge <- ts[ts$clone_id %in% names(aug_cut), ]
    candi_edge <- candi_edge[!(candi_edge$clone_id %in% new_clades), ]
    for (nc in new_clades) {
      frac <- length(aug_cut[[nc]]) / sum(tcc$Freq)
      temp <- data.frame(clone_id = nc, frac = frac)
      candi_edge <- bind_rows(candi_edge, temp)
    }
  }
  candi_edge
}


get_clone_dic_for_seq <- function(a_seq) {
  data.frame(old_K = a_seq, 
             is_ref = FALSE, letters = LETTERS[seq(length(a_seq))], color = get_cluster_colours(length(a_seq))[1:length(a_seq)],
             pretty_names = get_pretty_names_for_loci(a_seq), stringsAsFactors = F)
}

median_genotype_from_list <- function(clone_list, mat) {
  clones <- names(clone_list)
  nClones <- length(clones)
  res <- matrix(NA, nrow = nClones, ncol = nrow(mat))
  rownames(res) <- clones
  
  for (clone_i in seq_along(clones)) {
    clone <- clones[[clone_i]]
    sub_mat <- mat[, colnames(mat) %in% clone_list[[clone]]]
    res[clone_i, ] <- apply(as.matrix(sub_mat), 1, median)
  }
  colnames(res) <- rownames(mat)
  res
}


get_leaves_names_from_graph <- function(the_g, only_loci=FALSE) {
  edge_list <- as.data.frame(igraph::as_edgelist(the_g), stringsAsFactors = F)
  colnames(edge_list) <- c('source', 'target')
  res <- get_leaves_names(edge_list)
  if (only_loci) {
    res <- grep('locus_', res, value = T)
  }
  res
}

sugar_dummy_clone_dic <- function(ts, aug_cut, new_clades, tcc) {
  candi_edge <- get_candi_edge(the_cut = NULL, aug_cut = aug_cut, ts = ts, tcc = tcc, new_clades = new_clades)
  get_dummy_clone_dic(candi_edge)
}

get_dummy_clone_dic <- function(candi_edge) {
  data.frame(old_K = candi_edge$clone_id, 
             is_ref = FALSE, letters = LETTERS[seq(nrow(candi_edge))],
             pretty_names = get_pretty_names_for_loci(candi_edge$clone_id), stringsAsFactors = F)
}





get_dummy_candid_edge <- function(g, aug_cut) {
  ntotal <- g.count.cells(g)
  frac <- unlist(lapply(names(aug_cut), function(x) length(aug_cut[[x]])/ntotal))
  data.frame(clone_id = names(aug_cut), frac = frac, stringsAsFactors = F)
}

get_exp_path_for_datatag <- function(dt) {
  sa609_config <- load_main_configs(datatag = dt)
  dt_batch_path <- file.path('~/Desktop/SC-1311', dt, 'batch_runs', sa609_config$all_cells$batch_name)
  file.path(dt_batch_path, 'outputs', sa609_config$all_cells$exp_dir)
}


### Parsing inputs
################################################################################################################################################
read_ltm_tree <- function(edge_list) {
  # Find the root
  g <- igraph::graph_from_edgelist(as.matrix(edge_list))
  V(g)$id <- seq(vcount(g))
  return(g)
}


get_libid_from_cell_names <- function(cell_names) {
  # TODO: add the condition for SA906 and SA666
  gsub('SA([0-9]|[A-Z]|[a-z])+-(A([0-9]|[A-Z])+)-.*', '\\2', cell_names)
}



collapse_bin_names <- function(bin_dat) {
  paste0(bin_dat$chr, '_', bin_dat$start, '_', bin_dat$end)
}


# For Tyler version
# parse_bin_names <- function(bin_names, as_factors = FALSE) {
#   # Remove corrupt_tree locus tag if it's there
#   bin_names <- gsub('locus_', '', bin_names)
#   # ex: bin_names = "7_-1_0"  
#   # ex: bin_names = 7_150500001_151000000"
#   chr <- gsub('([0-9]+|X|Y)_(-1|[0-9]+)_[0-9]+', '\\1', bin_names)
#   start <- as.numeric(gsub('([0-9]+|X|Y)_(-1|[0-9]+)_[0-9]+', '\\2', bin_names))
#   end <- as.numeric(gsub('([0-9]+|X|Y)_(-1|[0-9]+)_([0-9]+)', '\\3', bin_names))
#   # print(paste0(chr,'_',start,'_',end))
#   data.frame(chr = chr, start = start, end = end, stringsAsFactors = as_factors)
# }

parse_bin_names <- function(bin_names, as_factors = FALSE) {
  # Remove corrupt_tree locus tag if it's there
  bin_names <- gsub('locus_', '', bin_names)
  chr <- gsub('([0-9]+|X|Y)_[0-9]+_[0-9]+', '\\1', bin_names)
  start <- as.numeric(gsub('([0-9]+|X|Y)_([0-9]+)_[0-9]+', '\\2', bin_names))
  end <- as.numeric(gsub('([0-9]+|X|Y)_([0-9]+)_([0-9]+)', '\\3', bin_names))
  data.frame(chr = chr, start = start, end = end, stringsAsFactors = as_factors)
}

### Graph general
################################################################################################################################################
get_height_dat <- function(g) {
  the_root <- find_root(g)
  height_search <- bfs(g, root=c(the_root), order=TRUE, neimode = 'out', unreachable=FALSE, dist = TRUE, rank=TRUE, succ=TRUE, father=TRUE, pred=TRUE)
  res <- data.frame(id=names(height_search$dist), dist=as.numeric(height_search$dist), stringsAsFactors = F)
  res
}

find_root <- function(g) {
  the_root <- NULL
  n_nodes <- vcount(g)
  for (i in seq(n_nodes)) {
    s1 = bfs(g, root=i, order=TRUE, neimode = 'out', unreachable=FALSE)
    if (any(is.na(s1$order)) == FALSE)  {
      the_root <- i
    }
  }
  if (is.null(the_root)) {
    print('Warning! No root was found...')
  }
  the_root
}

get_decendents <- function(clone_roots, the_graph, min_cell_per_clone=1) {
  decendents <- list()
  for (clone_root in clone_roots) {
    decend_node_ids <- decendents_ids(clone_root, the_graph)
    cells_in_clone <- length(decend_node_ids)
    if (cells_in_clone > min_cell_per_clone) {
      decendents[[clone_root]] <- decend_node_ids
    } 
  }
  decendents
}

decendents_ids <- function(root_id, the_graph) {
  decend_node_ids = bfs(the_graph, root=c(root_id), order = TRUE, neimode = 'out', unreachable=FALSE)
  decend_node_ids <- as.numeric(decend_node_ids$order)
  decend_node_ids[!is.na(decend_node_ids)]
}

find_parent <- function(the_g, node) {
  igraph::neighborhood(the_g, nodes=node, mode='in')[[1]]$name[[2]]
}

find_all_parents <- function(the_g, node, max_height) {
  igraph::neighborhood(the_g, nodes=node, mode='in', order = max_height)[[1]]$name[-c(1)]
}

find_all_descendents <- function(the_g, node, max_height) {
  igraph::neighborhood(the_g, nodes=node, mode='out', order = max_height)[[1]]$name[-c(1)]
}

find_immediate_children <- function(the_g, node) {
  children <- igraph::neighborhood(the_g, nodes=node, mode='out', order = 1)[[1]]$name
  children[!grepl(node, children)]
}

get_desc_names <- function(datatag, the_graph, the_node) {
  desc <- get_decendents(the_node, the_graph = the_graph)
  remove_na_from_list(node_number_to_label(the_graph, desc))
}

get_non_overlapping_desc <- function(the_graph, parent, node, cell_in_cut) {
  children <- find_immediate_children(the_g = the_graph, node = parent)
  children <- children[children != node & !is.na(children)]
  new_desc <- c()
  for (cc in children) {
    hh <- get_desc_names(datatag = datatag, the_graph = the_graph, the_node = cc)
    if (length(hh) >0) {
      hh <- hh[[1]]
    }
    hh <- hh[!is.na(hh)]
    
    if (length(hh) == 0) {
      new_desc <- c(new_desc, cc)
    } else if (!any(cell_in_cut %in% hh)) {
      new_desc <- c(new_desc, cc, hh)
    }
  }
  new_desc
}

get_children_from_parents <- function(g, parents) {
  h <- data.frame(id=V(g)$name, dist=unlist(V(g)$height), stringsAsFactors = F)
  h <- h[h$id %in% names(parents), ]
  h <- h[order(h$dist, decreasing = F), ]
  children <- c()
  for (ii in seq(nrow(h)-1)) {
    node <- h$id[[ii]]
    higer_nodes <- h$id[(ii+1):nrow(h)]
    if (any(parents[[node]] %in% higer_nodes)) {
      children <- c(children, node) 
    }
  }
  children
}

pick_elder_parent <- function(g, nodes) {
  h <- data.frame(id=V(g)$name, dist=unlist(V(g)$height), stringsAsFactors = F)
  h <- h[h$id %in% nodes, ]
  h <- h[order(h$dist, decreasing = F), ]
  h$id[1]
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

setup_graph <- function(g=NULL) {
  V(g)$id <- seq(vcount(g))
  h <- get_height_dat(g)
  stopifnot(all(h$id == V(g)$name))
  V(g)$height <- as.integer(h$dist)
  g
}

### Time points
################################################################################################################################################
node_number_to_label <- function(the_graph, groups) {
  all_dat <- data.frame(id=V(the_graph)$id, label=V(the_graph)$name, stringsAsFactors = F)
  res <- list()
  for (clade_root in names(groups)) {
    cells_in_clone <- groups[[clade_root]]
    cells <- all_dat$label[match(c(-1, cells_in_clone), all_dat$id)]
    res[[clade_root]] <- cells
  }
  res
}

### Misc
################################################################################################################################################
remove_na_from_list <- function(some_list) {
  # Why doesn't it keep names?
  tmp <- lapply(names(some_list), function(x) {some_list[[x]][!is.na(some_list[[x]])] })
  names(tmp) <- names(some_list)
  tmp
}

drop_non_cells <- function(aug_cut) {
  res <- lapply(names(aug_cut), function(x) aug_cut[[x]][!is.na(aug_cut[[x]]) & !grepl('(SA928|locus_)', aug_cut[[x]])]  )
  names(res) <- names(aug_cut)
  res
}

convert_to_freq <- function(time_series_dat, total_cell_counts_breakdown) {
  if (is.null(total_cell_counts_breakdown)) {
    print('Normalising to themselves')
    for (i in 2:ncol(time_series_dat)) {
      the_sum <- sum(time_series_dat[, i])
      if (the_sum > 0)
        time_series_dat[, i] <- time_series_dat[, i]/the_sum
    }
  } else {
    print('Normalising to Total count')
    for (i in 2:ncol(time_series_dat)) { 
      time_series_dat[, i] <- time_series_dat[, i]/total_cell_counts_breakdown$Freq[i-1]
    }
  }
  
  time_series_dat
}

order_cols <- function(ts) {
  x_cols = grep('X[0-9][0-9]*', colnames(ts))
  temp = ts[, x_cols, drop =F]
  temp = temp[, order(as.numeric(gsub('X', '', colnames(temp)))), drop=F]
  temp$clone_id <- ts$clone_id
  temp = temp[, c( ncol(temp), seq(ncol(temp)-1)), drop=F]
  temp
}


order_by_chr_name <- function(the_array) {
  the_array[the_array == 'X'] <- 44
  the_array[the_array == 'Y'] <- 100
  
  order(as.numeric(the_array))
}

strip_chr_names <- function(the_array) {
  print('In strip_chr_names')
  print(the_array)
  chunks <- unlist(strsplit(the_array, '_'))
  chrs <- matrix(chunks, nrow=length(the_array), ncol=3, byrow = T)[, 1]
  uchr <- unique(chrs)
  
  for (cc in uchr) {
    cc_idx <- which(chrs == cc)
    the_array[cc_idx] <- ""
    the_array[cc_idx[1]] <- cc
  }
  the_array
}

g.count.cells <- function(g=NULL) {
  length(g.get.cells(g))
}

g.get.cells <- function(g=NULL) {
  get_cells(V(g)$name)
}


g.get.loci <- function(g=NULL, edge_list_path=NULL) {
  if (is.null(g)) g <- read_ltm_tree(edge_list_path = edge_list_path)
  get_loci(V(g)$name)
}

get_loci <- function(str_arry) grep('locus_', str_arry, value = T)

get_cells <- function(str_arry) str_arry[!grepl('locus_', str_arry)]

get_cluster_colours <- function(nClusters) {
  if (nClusters > 8) {
    clust_colours <- colorRampPalette(brewer.pal(8, "Set2"))(nClusters)
  } else {
    clust_colours <- brewer.pal(nClusters, "Set2")
  }
  clust_colours
}


find_closest_index <- function(an_array, a_value) {
  tt <- an_array-a_value
  tt[tt <  0] <- -tt[tt <  0]
  which.min(tt)
}


get_colour_for_base_colour <- function(values, base_colours) {
  # Values are in (0, 11) inclusive
  # base_colours = legacy_colours_for_CN_heatmap()
  # values = unique(dat$copy_number)
  
  # Set values greater than 11 to 12 and set the colour to black
  base_colours <- c(base_colours, 'black')
  
  # Generate enhanced colours 
  shade_n <- 10
  c_mat <- matrix(NA, nrow = length(base_colours), ncol = shade_n)
  for (i in seq(length(base_colours))) {
    if (i < length(base_colours)) {
      c_mat[i, ] <- colorRampPalette(c(base_colours[i], base_colours[i+1]))(shade_n)
    } else {
      c_mat[i, ] <- colorRampPalette(c(base_colours[i], 'black'))(shade_n)
    }
  }
  
  res <- list()
  for (i in seq(length(values))) {
    #print(i)
    xx <- find_closest_index( ((1:10)-1)/shade_n, values[i] - round(values[[i]]))
    #print(sprintf('x %d', xx))
    res[[i]] <- c_mat[round(values[[i]]) + 1, xx ]
  }
  
  stopifnot(length(unlist(res)) == length(values))
  
  unlist(res)
}


get_corrupt_posterior_values <- function(offset = 50) {
  offset + (seq(1,11)-1)/10
}


library(tools)


sample_run_fast <- function(ctdir = NULL, just_normal = FALSE, use_greedy_tree = FALSE, ...) {
  if (!is.null(ctdir)) {
    newick <- file.path(ctdir, 'tree.newick')
    corrupt_input <- file.path(ctdir, 'filtered.csv')  # corrupt_filter func, input: straighten_output.csv, output: filtered.csv
    posterior_mat <- file.path(ctdir, 'average.csv')
    if (use_greedy_tree) {
      newick <- file.path(ctdir, 'consensus.newick')
    }
  }
  print('Generating the ordinary one')
  if (just_normal) {
    sample_run(newick = newick, corrupt_input = NULL, posterior_mat = NULL, ...)
  } else {
    sample_run(newick = newick, corrupt_input = corrupt_input, posterior_mat = posterior_mat, ...)
  }
}


sample_run <- function(newick, copy_number, cell_clones = NULL, output = NULL, use_all_loci = TRUE, filter_cells = NULL, 
                       corrupt_input = NULL, posterior_mat = NULL, chr_filter = NULL, dev = 'png', g = NULL, drop_loci_names = FALSE) {
  if (is.null(output)) {
    tag <- ifelse(use_all_loci, 'all_loci', 'just_tree_loci')
    #output = sprintf('%s/%s_heatmap_original_%s_%s_.png', dirname(newick), file_path_sans_ext(basename(newick)), generate_random_str(), tag)
    if (dev == 'png' & !is.null(corrupt_input) & !is.null(posterior_mat)) {
      # split by chromosome
      # output <- sprintf('%s/%s/%s_heatmap_original_%s_%s.%s', dirname(newick), generate_random_str(), file_path_sans_ext(basename(newick)), generate_random_str(), tag, dev)
      output <- sprintf('%s/%s/%s_heatmap_original_%s_%s.%s', dirname(newick), 'debug_tree', file_path_sans_ext(basename(newick)), generate_random_str(), tag, dev)
      dir.create(dirname(output), recursive = T)
    } else {
      output <- sprintf('%s/%s_heatmap_original_%s_%s.%s', dirname(newick), file_path_sans_ext(basename(newick)), generate_random_str(), tag, dev)
    }
    
  }
  
  if (is.null(g)) {
    tree <- read.tree(newick)
    g <- read_ltm_tree(tree_2_edge_list(tree))
  }
  
  copynumber <- read.csv(copy_number, header=T, row.names=1, check.names = F,stringsAsFactors = FALSE)
  # copynumber <- read.delim(copy_number, check.names=FALSE, stringsAsFactors=FALSE, sep=",", row.names = 1)
  copynumber <- format_copynumber_matrix(copynumber)
  
  # JUST FOR 80 percent
  if (!is.null(filter_cells)) {
    copynumber <- copynumber[, (colnames(copynumber) %in% filter_cells)]
    bad_cells <- g.get.cells(g)[!(g.get.cells(g) %in% filter_cells)]
    bad_cells <- bad_cells[bad_cells != 'root']
    g <- igraph::delete.vertices(g, bad_cells)
  }
  
  copynumber <- copynumber[, colnames(copynumber) %in% g.get.cells(g)]
  
  stopifnot(all(colnames(copynumber) %in% g.get.cells(g)))
  
  # Interleave the binary space and posterior_genotype into the output
  if (!is.null(posterior_mat)) {
    print('Interleaving posterior mat (average.csv)...')
    pmat <- fasttable(posterior_mat)
    #fmat <- fasttable(corrupt_input)
    #all(pmat$loci %in% fmat$loci)
    head(pmat)
    lenu(pmat$tipInclusionProbabilities)
    pos_offset <- 50
    pmat$tipInclusionProbabilities <- pmat$tipInclusionProbabilities + pos_offset
    # Bin the posterior values into 10 shades to make manual colour assignment easier
    values <- get_corrupt_posterior_values(pos_offset)
    binned <- values[unlist(lapply(pmat$tipInclusionProbabilities, function(x) find_closest_index(values, x)))]
    pmat$tipInclusionProbabilities <- binned
    #tt = get_colour_for_base_colour(values = unique(binned), base_colours = brewer.pal(n = 9, name = "Greys"))
    #colorRampPalette(brewer.pal(8, "Greys"))(lenu(get_corrupt_posterior_values()))
    pmat <- pmat %>% tidyr::spread(key = cells, value = tipInclusionProbabilities)
    pmat$loci <- gsub('locus_', '', pmat$loci)
    rownames(pmat) <- pmat$loci
    pmat$loci <- NULL
    headm(pmat)
    colnames(pmat) <- gsub('cell_', '', colnames(pmat))
    stopifnot(all(colnames(pmat) %in% colnames(copynumber)))
    
    # Add +4 to the start value of the loci (just used for sorting)
    options(scipen=999) # prevent scientific notation
    bins <- parse_bin_names(rownames(pmat), as_factors = F) 
    bins$start <- bins$start + 3
    rownames(pmat) <- collapse_bin_names(bins)
    
    # Add to the copynumber and resort
    s1 <- as.data.frame(copynumber)
    s1$bin_name <- rownames(s1)
    pmat$bin_name <- rownames(pmat)
    copynumber1 <- dplyr::bind_rows(s1, pmat)
    rownames(copynumber1) <- copynumber1$bin_name
    copynumber1$bin_name <- NULL
    copynumber1 <- as.matrix(copynumber1)
    str(copynumber1)
    
    copynumber <- sort_mat_by_bins(the_mat = copynumber1)
  }
  
  if (!is.null(corrupt_input)) {
    print('Interleaving corrupt_input mat (filtered.csv)...')
    # Change values in fdat, start from 40...
    fmat <- fasttable(corrupt_input)
    fmat$tipInclusionProbabilities <- fmat$tipInclusionProbabilities + 40
    # update cell names
    fmat$cells <- gsub('cell_', '', fmat$cells)
    # sanity check
    stopifnot(all(unique(fmat$cells) %in% colnames(copynumber)))
    which(unique(fmat$cells) %ni% colnames(copynumber))
    head(fmat)
    
    # Add binary loci after the two loci
    # Padding?
    fdat <- fmat %>% tidyr::spread(key = cells, value = tipInclusionProbabilities)
    fdat$loci <- gsub('locus_', '', fdat$loci)
    rownames(fdat) <- fdat$loci
    fdat$loci <- NULL
    headm(fdat)
    
    # Add +2 to the start value of the loci
    options(scipen=999) # prevent scientific notation
    bins <- parse_bin_names(rownames(fdat), as_factors = F) 
    bins$start <- bins$start + 2
    rownames(fdat) <- collapse_bin_names(bins)
    
    # Add to the copynumber and resort
    s1 <- as.data.frame(copynumber)
    s1$bin_name <- rownames(s1)
    fdat$bin_name <- rownames(fdat)
    copynumber1 <- dplyr::bind_rows(s1, fdat)
    rownames(copynumber1) <- copynumber1$bin_name
    copynumber1$bin_name <- NULL
    copynumber1 <- as.matrix(copynumber1)
    str(copynumber1)
    
    copynumber <- sort_mat_by_bins(the_mat = copynumber1)
  }
  
  
  function() {
    qq = rownames(copynumber)
    ff = parse_bin_names(qq)
    ff
  }
  
  if (!is.null(cell_clones)) {
    cellclones <- read.delim(cell_clones, stringsAsFactors = FALSE, sep=",")
  } else {
    # Assign all cells to one clone
    cellclones <- data.frame(cell_id = colnames(copynumber), clone_id = 'A', stringsAsFactors = F)
  }
  
  # TODO: check inputs
  # 1. Are all the cells in the matrix?
  # 2. Are there double assignments?
  
  # TODO: input clone colours
  #fitclone_plot_tree_heatmap(g = g, cell_clones = cellclones, output =  output,  mat = copynumber, just_tree_loci = TRUE, chr_filter = c(1,2,3,4))
  
  #fitclone_plot_tree_heatmap(g = g, cell_clones = cellclones, output =  '/Users/sohrabsalehi/projects/tree_viz/ss_viz/outputs/test.png',  mat = copynumber, just_tree_loci = TRUE, chr_filter = c(1,2,3,4))
  #fitclone_plot_tree_heatmap(g = g, cell_clones = cellclones, output = argv$output,  mat = copynumber, just_tree_loci = FALSE, chr_filter = NULL)
  # aug_cut <- data_frame_to_list(cell_clones)  
  # fast_tree_aug(g = g, aug_cut = aug_cut)
  # 
  #fitclone_plot_tree_heatmap(g = g, cell_clones = cellclones, output = output,  mat = copynumber, just_tree_loci = FALSE, chr_filter = NULL)
  if (dev == 'png' & !is.null(corrupt_input) & !is.null(posterior_mat)) {
    if (is.null(chr_filter)) {
      # break it down to multiple ones 
      chr_filter <- list(c(1,2), c(3,4), c(5,6,7), c(8,9,10,11), c(12, 13, 14, 15), c(16, 17, 18, 19, 20), c(21, 22, 'X', 'Y')  )
    }
    for (cf in chr_filter) {
      print(sprintf('chr_filter is %s', cf))
      output <- gsub('.png', paste0('_', paste0(cf, collapse = ''), '.png'), output)
      fitclone_plot_tree_heatmap(g = g, cell_clones = cellclones, output = output,  mat = copynumber, just_tree_loci = !use_all_loci, chr_filter = cf)
    }
    
    system(sprintf('open %s', dirname(output)))
    
    # do the conversion
    # cmd_str <- sprintf('convert $(ls -v1 *.png) collection_%s.pdf', tolower(format(Sys.time(), "%b_%d")))
    # ss <- getwd()
    # setwd(dirname(output))
    # #system2(command = '/Users/sohrabsalehi/Desktop/SC-1311/plots/may_15/timescape/run.sh', args = dt)
    # system(cmd_str)
    # setwd(ss)
    
    
  } else {
    # cell_clones = cellclones; mat = copynumber; just_tree_loci = !use_all_loci
    fitclone_plot_tree_heatmap(g = g, cell_clones = cellclones, output = output,  mat = copynumber, just_tree_loci = !use_all_loci, chr_filter = chr_filter, drop_loci_names = drop_loci_names)
    system(sprintf('open %s', output))
  }
  #fitclone_plot_tree_heatmap(g = g, cell_clones = cellclones, output = output,  mat = copynumber, just_tree_loci = !use_all_loci, chr_filter = c(1,2))
  
  print(output)
}


annotate_nodes_on_tree <- function(p = NULL, tree_node_dic = NULL, edge_list_path, node_names, show_label=TRUE, leaf_shape=21, leaf_colour='steelblue', leaf_size=3, equal_branch=FALSE, use_pretty_names=T) {
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
  
  print(internal_numbers)
  
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
  p <- p + geom_point2(aes(subset=(node %in% internal_numbers)), size=1, shape=23, fill="red") 
  p$data <- dplyr::left_join(p$data, vertex_names, by=c('node'='node_number'))
  #p <- p + geom_text2(aes(subset=(node %in% internal_numbers), label=node_name), hjust=-.3, color = 'red', position = position_jitter(height=150)) 
  if (use_pretty_names) {
    p$data$node_name <- get_pretty_names_for_loci(p$data$node_name)
  }
  
  p <- 
    p + geom_text_repel(data = p$data[p$data$node %in% internal_numbers, ], 
                        mapping = aes(x = x, y = y, label = node_name), colour = '#ED7953FF', size = 3, nudge_y = 0, nudge_x = 0, angle = 45, 
                        #arrow = arrow(length = unit(0.03, "npc"), type = "open", ends = "last"),
                        fontface = 'bold',
                        segment.color = 'black',
                        segment.size = .2)
  
  #p <- p + geom_text2(aes(subset=(node %in% internal_numbers), label=node_name), hjust=-.3, color = '#ED7953FF', angle = 45, face = "bold") 
  p  
}




complete_run <- function(newick = NULL, copy_number = NULL, cell_clones = NULL, use_all_loci = TRUE, dev = 'png', chr_filter = NULL, corrupt_input = NULL, posterior_mat = NULL) {
  # No augmentation 
  sample_run(newick = newick, 
             copy_number = copy_number,
             cell_clones = cell_clones,
             use_all_loci = use_all_loci,
             corrupt_input = NULL,
             posterior_mat = NULL,
             dev = dev,
             chr_filter = chr_filter)
  
  # With augmentation
  if (!is.null(corrupt_input) | !is.null(posterior_mat)) {
    sample_run(newick = newick, 
               copy_number = copy_number,
               cell_clones = cell_clones,
               use_all_loci = use_all_loci,
               corrupt_input = corrupt_input,
               posterior_mat = posterior_mat,
               dev = dev,
               chr_filter = chr_filter)
  }
  
  # TODO: add teh blue one, with just corrupt tree input here
  
  # Add libraries
}




'%ni%' <- function(x,y)!('%in%'(x,y))

fasttable <- function(...) as.data.frame(fread(..., stringsAsFactors = F))

len <- function(...) length(...)

lenu <- function(...) length(unique(...))

headm <- function(amat, n = 10, ncol= 10) head(amat[, 1:ncol], n = n)

generate_random_str <- function(n = 1) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}





# option_list <- list(make_option(c("-t", "--newick"), type="character", default=NULL, help="output from corrupt_tree", metavar="character"),
#                     make_option(c("-n", "--copynumber"), type="character", default=NULL, help="filtered_cell_cn", metavar="character"),
#                     make_option(c("-c", "--cellclones"), type="character", default=NULL, help="output from tree_cutting", metavar="character"),
#                     make_option(c("-o", "--output"), type="character", default=NULL, help="tree_heatmap.png", metavar="character")
# )

#opt_parser <- OptionParser(option_list=option_list)
#opt <- parse_args(opt_parser)
#sample_run(opt$newick, opt$copynumber, opt$cellclones, opt$output)


# function() {
#   sample_run(newick = '~/Desktop/SC-1311/SA922/processed_data/bincount/10k_10chains_filtered_fewer_loci/np_pt_work_filtered_cycle_0.85_SA922_fewer_loci_05/no_pt_.05/tree.newick', 
#              copy_number = '~/Desktop/SC-1311/SA922/processed_data/SA922_cnv_data.csv',
#              cell_clones = NULL,
#              use_all_loci = FALSE,
#              corrupt_input = NULL,
#              posterior_mat = NULL,
#              dev = 'png',
#              output = NULL, 
#              #chr_filter = c(1,2))
#   chr_filter = NULL)
#   
# }





# Account for synchronized cnvs
pad_mod_cnv_mat <- function(Y, output_file, pad_left = TRUE) {
  # TODO: We're missing synch events across the genome, i.e., those that happen in multiple cells, in more than one bin. They stand out in say SA919_Hoa
  # Binary matrix B, CNV matrix Y, cell c, locus l, set of cells C
  # 1. For each locus l, after removal of multi-chr-spanning loci, 
  # 1.a. find set of cells C for which B_l,[C] == 1 & |unique(Y_l,[C])| > 1
  # Note: this will fail at chr ends, handle via chr-mod-padding
  # 2. Add a new bin such that
  # B_l',i =  1, if Y_l,i == j, where j \prop Cat(unique(Y_l,[C]); proportions of diff. CNVs in [C])
  #           0, otherwise
  # Note: consider using B after jitterfix
  # Note: this will miss whole chr events -> Add a locus to B at either end of each chr where it's 1 for the mod and 0 o.w. (consider setting the second most frequent CNV bins to 1)
  # 3. Only use bins with |max(j)| < p%
  # Only consider loci that have |B_l,. == 1| > 90%
  # 4. Set l' (the name of the new locus in B) to chr_{(l+1)_(start) - l_start}_{0}
  
  # Inputs the_mat: Y; 
  # out_dir <- file.path(dirname(cnv_path), 'results'); dir.create(out_dir, showWarnings = F, recursive = T)
  # Y <- fasttable(cnv_path)
  # rownames(Y) <- Y$V1; Y$V1 <- NULL
  
  
  Y <- sort_mat_by_bins(Y)
  bin_dat <- parse_bin_names(rownames(Y))
  
  ########################################
  B <- NULL
  uchrs <- unique(bin_dat$chr)
  for (chr in uchrs) {
    print(sprintf('chr == %s', chr))
    sub_Y <- Y[bin_dat$chr == chr, ]
    tmp <- abs(sub_Y[-c(nrow(sub_Y)), ] - sub_Y[-c(1), ])
    B <- addrow(B, tmp)
  }
  
  B[B > 1] <- 1
  
  ## Now has the right inputs, Y and B, both sorted
  # Note: For B[l, ] == 0, chromosome padding will catch them, otherwise they will be caught on a different bin
  # TODO: jitter fix B
  # TODO: DO THIS AFTER JITTERFIX in nextflow pipeline so as to not allow jitterfix to kill the signal
  # 1. Pick loci s.t. |B_l,. == 1| > 85%
  freq <- rowMeans(B)
  #L  <- which(unname(freq > .85))
  L  <- which(unname(freq > .63))
  freq[L]
  new_B <- B
  new_bins <- c()
  for (l in L) {
    # l = L[1]
    indx <- which(B[l, ] == 1)
    # 2. Add a new bin such that everything is great
    l_y <- which(rownames(Y) == rownames(B)[l])
    counts <- as.data.frame(table(as.vector(as.matrix(Y[l_y, indx]))), stringsAsFactors = F)
    counts$Freq <- counts$Freq/sum(counts$Freq)
    counts$Var1 <- as.numeric(counts$Var1)
    counts <- counts[order(counts$Freq, decreasing = T), ]
    if (nrow(counts) > 1) {
      if (counts$Freq[2] > .05) {
        print(sprintf('%s: L_B = %d, L_Y = %d',rownames(B)[l], l, l_y))
        tmp <- B[l, ]
        tmp[] <- 0
        tmp[which(Y[l_y, ] == counts$Var1[2])] <- 1
        fake_start <- floor((bin_dat$start[l_y] + bin_dat$end[l_y])/3)
        rownames(tmp) <- sprintf('%s_%d_%d', bin_dat$chr[l_y], fake_start, fake_start + 1)
        print(rownames(tmp))
        new_B <- rbind(new_B, tmp)
        new_bins <- c(new_bins, rownames(tmp))
      }
    }
  }
  
  
  # Now add the mode chromosome padding: use the categorical here
  new_bins <- c()
  for (chr in uchrs) {
    print(sprintf('Adding padding for chr == %s', chr))
    # chr = "10"
    last_Y <- Y[bin_dat$chr == chr, ]
    last_Y <- last_Y[nrow(last_Y), ]
    last_bin_name <- rownames(last_Y)
    l <- which(rownames(Y) == last_bin_name) # index in bin_names
    counts <- as.data.frame(table(as.vector(as.matrix(last_Y))), stringsAsFactors = F) %>% 
      dplyr::mutate(freq = Freq/sum(Freq), j = as.numeric(Var1)) %>% 
      dplyr::arrange(desc(freq)) %>% 
      dplyr::select(j, freq) 
    
    if (nrow(counts) > 1) {
      print(sprintf('%s: L_Y = %d', last_bin_name, l))
      tmp <- B[1, ]
      tmp[] <- 0
      tmp[which(sub_Y[l_y, ] == counts$j[2])] <- 1
      rownames(tmp) <- sprintf('%s_%d_%d', bin_dat$chr[l], bin_dat$end[l] + 1, bin_dat$end[l] + 2) # It is always larger
      print(rownames(tmp))
      new_B <- rbind(new_B, tmp)
      new_bins <- c(new_bins, rownames(tmp))
    }
    
    
  }
  
  # Do we want to add padding on the left too? Technically they will reinforce the whole chromosome and initila chromosome events...
  # TODO: where are the bins that we're removing??
  if (pad_left) {
    new_bins <- c()
    for (chr in uchrs) {
      print(sprintf('Adding left-padding for chr  %s', chr))
      # chr = "10"
      first_Y <- Y[bin_dat$chr == chr, ][1, ]
      first_bin_name <- rownames(first_Y)
      l <- which(rownames(Y) == first_bin_name) # index in bin_names
      counts <- as.data.frame(table(as.vector(as.matrix(first_Y))), stringsAsFactors = F) %>% 
        dplyr::mutate(freq = Freq/sum(Freq), j = as.numeric(Var1)) %>% 
        dplyr::arrange(desc(freq)) %>% 
        dplyr::select(j, freq) 
      
      # Found multi cnvs?
      if (nrow(counts) > 1) {
        print(sprintf('%s: L_Y = %d', first_bin_name, l))
        tmp <- B[1, ]
        tmp[] <- 0
        tmp[which(sub_Y[l_y, ] == counts$j[2])] <- 1
        rownames(tmp) <- sprintf('%s_%d_%d', chr, 0, 1) # It is always smaller
        print(rownames(tmp))
        new_B <- rbind(new_B, tmp)
        new_bins <- c(new_bins, rownames(tmp))
      }
    }
  }
  
  
  # sort new_B by rows again to put the new loci in their place...
  new_B <- sort_mat_by_bins(the_mat = new_B)
  nrow(new_B) - nrow(B) 
  
  
  new_B$loci <- rownames(new_B)
  rownames(new_B) <- NULL
  
  # cells,loci,tipInclusionProbabilities
  new_B <- new_B %>% tidyr::gather(key = 'cells', value = 'tipInclusionProbabilities', -loci)
  
  new_B <- new_B[, c(2,1,3)]
  head(new_B); summary(new_B$tipInclusionProbabilities); boxplot(new_B$tipInclusionProbabilities); table(new_B$tipInclusionProbabilities)
  
  # output_file <- file.path(out_dir, 'bin_cnvs_corrupt_double_padding.csv')
  data.table::fwrite(new_B, output_file, row.names = F, quote = F)  
  
  # ./nextflow run binary-infer-pipeline.nf -resume --tipInclusionProbabilities ~/Desktop/SC-1311/SA609/processed_data/bincount/sa609_cnvs_corrupt.csv
  # pipeline <- 'no_pt_fast'
  # print(sprintf("./nextflow run %s.nf -resume --tipInclusionProbabilities %s", pipeline, file_path))
}


# remove odd cells
# Does not count whole-chromosome events
compute_jump_cells <- function(mat) {
  # Cells that have too many CN up and downs may represent douplet? or s-phases
  # Exclude them after a threshold
  # dat <- load_new_cn_data(datatag)
  dat <- sort_mat_by_bins(mat)
  
  # Bins X cells
  mat_delta <- abs(dat[1:(nrow(dat)-1), ] - dat[2:nrow(dat), ])
  mat_delta[mat_delta > 1] <- 1
  njumps <- colSums(mat_delta)
  avgCNA <- colMeans(dat)
  # res <- sort(njumps, decreasing = T)
  stopifnot(colnames(mat_delta) == colnames(dat))
  data.frame(cell_id = colnames(dat), njumps = njumps, avgCNA = avgCNA, stringsAsFactors = F) %>% dplyr::arrange(desc(njumps))
}

cn2binary_newencoding <- function(cnv_path, output_file, prob = .90) {
  # TODO: compute average ploidy too for each cell
  
  # cnv_path <- paste0(data_dir,'total_merged_filtered_states_original.csv')
  mat <- read.delim(cnv_path, check.names = FALSE, stringsAsFactors = FALSE, sep = ",")
  if(prob > 0){
    jump_rank <- compute_jump_cells(mat)
    
    
    # # Filter cells by jump/stuff
    bad_cells <- jump_rank %>% dplyr::filter( njumps > quantile(x=njumps, probs = prob) | avgCNA > quantile(x=avgCNA, probs = prob)) %>% dplyr::select(cell_id)
    
    # bad_cells <- jump_rank %>% dplyr::filter( njumps > quantile(x=njumps, probs = prob)) %>% dplyr::select(cell_id)
    
    bad_cells <- bad_cells$cell_id
    print("Number of bad cells: ")
    print(length(bad_cells))
    cc <- mat[, colnames(mat) %ni% bad_cells]
    print(paste0("Filtered output: ", dim(cc)[1]," ",dim(cc)[2]))
    write.csv(cc, opt$input,  row.names = T, quote=F)
    pad_mod_cnv_mat(cc, output_file, pad_left = TRUE)
    
  } else{
    print("Convert data to binary without filtering data: ")
    pad_mod_cnv_mat(mat, output_file, pad_left = TRUE)
  }
  
  # The results would be in:
  # bin_cnvs_corrupt_double_padding.csv
  
}

# Convert copy number to binary using new encoding scheme
# print(paste0("Filtered outlier threshold is: ", opt$filterthrs))
# cn2binary_newencoding(opt$input, opt$outfile, as.numeric(opt$filterthrs))

# test_run <- function() {
#   # TODO: compute average ploidy too for each cell
#   cnv_path <- 'data/total_merged_filtered_states.csv'
#   mat <- fasttable(cnv_path)
#   rownames(mat) <- mat$V1; mat$V1 <- NULL
#   jump_rank <- compute_jump_cells(mat = mat)
#   
#   # Filter cells by jump/stuff
#   prob <- .80
#   #bad_cells <- jump_rank %>% dplyr::filter( njumps > quantile(x=njumps, probs = prob) | avgCNA > quantile(x=avgCNA, probs = prob)) %>% dplyr::select(cell_id)
#   bad_cells <- jump_rank %>% dplyr::filter( njumps > quantile(x=njumps, probs = prob)) %>% dplyr::select(cell_id)
#   bad_cells <- bad_cells$cell_id
#   
#   # uncomment for enhanced viz
#   # sample_run_fast(ctdir = '/Users/sohrabsalehi/projects/corrupt-nextflow/deliverables/no_pt_fast_default_error_global_cell_filter_full_padding_tyler_more_loci/',
#   #                 copy_number = cnv_path,
#   #                 cell_clones = NULL,
#   #                 use_all_loci = TRUE,
#   #                 #filter_cells = bad_cells,
#   #                 filter_cells = setdiff(jump_rank$cell_id, bad_cells),
#   #                 dev = 'png',
#   #                 #chr_filter = c(1),
#   #                 chr_filter = NULL,
#   #                 just_normal = TRUE, 
#   #                 drop_loci_names = TRUE)
#   
#   
#   # Filter cells
#   cc <- mat[, colnames(mat) %ni% bad_cells]; dim(cc)
#   write.table(cc, 'data/total_merged_filtered_states_filtered_80.csv', sep = ',', row.names = T)
#   pad_mod_cnv_mat('data/total_merged_filtered_states_filtered_80.csv')
#   
#   # The results would be in:
#   # ./data/results/bin_cnvs_corrupt_double_padding.csv
#   
#   # Viz
#   # sample_run_fast(ctdir = '??',
#   #                 copy_number = cnv_path,
#   #                 cell_clones = NULL,
#   #                 use_all_loci = TRUE,
#   #                 dev = 'png',
#   #                 #chr_filter = c(1),
#   #                 chr_filter = NULL,
#   #                 just_normal = TRUE, 
#   #                 drop_loci_names = TRUE)
#   
# }


