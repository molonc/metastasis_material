# print(.libPaths())


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

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

source(paste0(script.basename, "/utils.R"))

fitclone_plot_tree_heatmap <- function(the_graph,
 cell_clones = NULL,
 copy_number_matrix = copy_number_matrix,
 output = NULL,
 just_tree_loci = TRUE,
 chr_filter = c(1,2)
 ) {
  fitclone_plot_tree_heatmap_more(the_graph = the_graph,
    cell_clones = cell_clones,
    copy_number_matrix = copy_number_matrix,
    output = output,
    just_tree_loci = just_tree_loci,
    chr_filter = chr_filter
    )
}

fitclone_plot_tree_heatmap_more <- function(the_graph,
  cell_clones,
  copy_number_matrix = copy_number_matrix,
  output = NULL,
  aug_cut = NULL,
  chr_filter = NULL,
  just_tree_loci = NULL,
  chr_font_size = 3.8
  ) {
  all_cells <- g.count.cells(the_graph)
  covered_cells <- nrow(cell_clones)
  title_str <- sprintf("%d/%d", covered_cells, all_cells)

  print("Generating heatmap object...")
  pp <- tree_and_heatmap(the_graph = the_graph,
   cell_clones = cell_clones,
   copy_number_matrix = copy_number_matrix,
   chr_filter = chr_filter,
   just_tree_loci = just_tree_loci,
   chr_font_size = chr_font_size
   )
  pp <- pp +
  ggtitle(title_str) +
  theme(legend.text = element_text(size = 15),
    plot.title = element_text(size = 21),
    plot.subtitle = element_text(size = 20, face = "italic", color = "black")
    )

  print("Saving the heatmap object...")
  coef <- 1
  height <- 8 * coef
  width <- 14 * coef
  ggsave(plot = pp, filename = output,  width = width, height = height, units = "in", limitsize = FALSE)
}

tree_and_heatmap <- function(the_graph, cell_clones, copy_number_matrix, just_tree_loci = FALSE, chr_filter = NULL, chr_font_size = 3.8) {
  if (!is.null(chr_filter)) {
    chr <- gsub("([0-9]+|X|Y)_[0-9]+_[0-9]+", "\\1", rownames(copy_number_matrix))
    copy_number_matrix <- copy_number_matrix[chr %in% chr_filter, ]
  }

  copy_number_matrix <- sort_mat_by_bins(copy_number_matrix = copy_number_matrix)
  copy_number_matrix <- remove_pad_bins(copy_number_matrix)

  copy_number_matrix <- copy_number_matrix[, colnames(copy_number_matrix) %in% cell_clones$cell_id]
  copy_number_matrix <- copy_number_matrix[, match(cell_clones$cell_id, colnames(copy_number_matrix))]


  # The tree
  res <- convert_edge_list_to_ape(edge_list_mat = igraph::as_edgelist(the_graph))

  tree <- res$tree
  node_names <- res$node_names
  the_height <- find_tail_height(the_graph = the_graph)$h_cut

  if (!is.null(the_height)) {
    tree_res <- trim_tree_before_height(tree, node_names, the_graph, the_height, 20)
    tree <- tree_res$tree
  }

  # Add the heatmap
  sub_mat <- heatmap_mat_from_bin_cellID(copy_number_matrix, update_colnames = FALSE)

  if (just_tree_loci) {
    the_loci <- grep("locus", node_names$node_name, value = TRUE)
    the_loci <- gsub("locus_", "", the_loci)
    sub_mat <- sub_mat[, colnames(sub_mat) %in% the_loci]
  }

  colnames_filter <- strip_chr_names(colnames(sub_mat))
  ss <- which(colnames_filter != "")
  colnames_filter <- data.frame(from = colnames(sub_mat)[ss], new_label = colnames_filter[ss])

  # Add white borders between chromosomes
  genotype <- as_tibble(sub_mat)
  white_border_size <- 0
  if (!is.null(just_tree_loci)) {
    white_border_size <- 0
  }

  if (white_border_size > 0) {
    for (cname in colnames_filter$from[-c(1)]) {
      for (i in 1:white_border_size) {
        new_name <- gsub("([0-9]+|X|Y)_(.*)", paste0("\\1_", i - 1, "_\\2"), cname)
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
  genotype <- matrix(data = genotype, nrow = d[1], ncol = d[2])
  colnames(genotype) <- cn
  rownames(genotype) <- rn

  colnames_level <- unique(colnames(genotype))

  # Add colours for CNV value
  cn_colours <- legacy_colours_for_CN_heatmap()
  names(cn_colours) <- paste0(seq_along(cn_colours) - 1)

  tip.loci <- grep("locus", tree$tip.label, value = TRUE)
  tree <- drop.tip(tree, tip.loci)

  # Plot
  aug_cut <- data_frame_to_list(cell_clones)

  #browser()

  xtree <- ggtree::groupOTU(tree, aug_cut)
  p <- ggtree(xtree, aes(color = group))
  offset <- 10
  breaks <- factor(x = names(cn_colours), levels = names(cn_colours), ordered = TRUE)

  p <- sos_heat(p = p,
    data = genotype,
    offset = offset,
    width = 8,
    colnames = TRUE,
    color = NULL,
    colnames_level = colnames_level,
    colnames_filter = colnames_filter,
    font.size = chr_font_size,
    colnames_offset_y = -40
    ) +
  scale_fill_manual(breaks = breaks,
   values = cn_colours,
   guide = guide_legend(nrow = 1, direction = "horizontal", label.position = "bottom")
   )

  p <- p + theme(legend.position = "top")
  # Add the cluster annotations
  p <- annoate_tree_aug(p = p,
    the_graph = the_graph,
    candiate_edges = NULL,
    aug_cut = aug_cut,
    the_height = the_height,
    show_guide = FALSE,
    show_label = FALSE
    )
  return(p)
}

sos_heat <- function (p,
  data,
  offset = 0,
  width_mut = 1,
  low = "green",
  high = "red",
  color = "white",
  colnames = TRUE,
  colnames_position = "bottom",
  colnames_angle = 0,
  colnames_level = NULL,
  colnames_offset_x = 0,
  colnames_offset_y = 0,
  font.size = 4,
  hjust = 0.5,
  colnames_filter = NULL
  ) {
  colnames_position <- colnames_position %>% match.arg(c("bottom", "top"))
  variable <- value <- lab <- y <- NULL
  width <- width_mut * (p$data$x %>% range(na.rm = TRUE) %>% diff) / ncol(data)
  isTip <- x <- y <- variable <- value <- from <- to <- NULL
  df <- p$data
  df <- df[df$isTip, ]
  start <- max(df$x, na.rm = TRUE) + offset
  i <- order(df$y)
  i <- i[!is.na(df$y[i])]
  lab <- df$label[i]
  dd <- as.data.frame(data)
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

  V2 <- start + as.numeric(dd$variable) * width
  mapping <- data.frame(from = dd$variable, to = V2)
  mapping <- unique(mapping)
  dd$x <- V2
  dd$width <- width
  dd$x_lib <- min(dd$x) - width_mut
  dd$value <- sapply(dd$value, trimws)

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

  if (colnames) {
    if (colnames_position == "bottom") {
      y <- 0
    } else {
      y <- max(p$data$y) + 1
    }
    mapping$y <- y
    if (!is.null(colnames_filter)) {
      mapping <- mapping[mapping$from %in% colnames_filter$from, ]
      mapping <- dplyr::right_join(mapping, colnames_filter)
      mapping$from <- mapping$new_label
    }

    p2 <- p2 + geom_text(data = mapping,
      aes(x = to, y = y, label = from),
      size = font.size,
      inherit.aes = FALSE,
      angle = colnames_angle,
      nudge_x = colnames_offset_x,
      nudge_y = colnames_offset_y,
      hjust = hjust
    )
  }

  p2 <- p2 + theme(legend.position = "right", legend.title = element_blank())
  attr(p2, "mapping") <- mapping
  return(p2)
}

fast_tree_aug <- function(the_graph, aug_cut, show_guide = TRUE) {
  label_offset <- .8
  angle <- 0
  fontsize <- 8
  candiate_edges <- get_dummy_candid_edge(the_graph, aug_cut)
  res <- graph_to_tree_dic(the_graph)

  tree <- res$tree
  node_names <- res$node_names

  myColors <- get_cluster_colours(length(aug_cut))
  candiate_edges$colours <- myColors
  candiate_edges$letters <- candiate_edges$clone_id
  names(myColors) <- names(aug_cut)

  # TODO: Since 0 is shared with the heatmap CNV state, use something else for the Un-assigned ones
  myColors <- c(myColors, "0" = "black")

  cdat <- data.frame(clone_id = names(myColors), colours = myColors, stringsAsFactors = FALSE)

  # Group 0, the default of the tree should be black
  xtree <- groupOTU(tree, .node = aug_cut)
  p <- ggtree(xtree, aes(color = group), size = 1)

  # TODO: where are the f..ing clone names? Why aren't they shown
  p <- p + scale_color_manual(values = myColors, guide = show_guide)

  qq <- p$data
  qq <- qq[qq$isTip == TRUE, ]

  compute_mean_clade_size <- function(candiate_edges, qq) {
    clade_size <- c()
    for (i in seq_along(candiate_edges$clone_id)) {
      the_key <- candiate_edges$clone_id[[i]]
      sq <- qq[qq$label %in% aug_cut[[the_key]], ]
      clade_size <- c(clade_size, max(sq$y) - min(sq$y))
    }
    min(clade_size)
  }

  ave_clade_size <- compute_mean_clade_size(candiate_edges, qq)

  for (i in seq_along(candiate_edges$clone_id)) {
    # Find the two ends of this pseudo-clade
    # TODO: handle the case where the two ends of the pseudo-clade run over another clade
    the_key <- candiate_edges$clone_id[[i]]
    sq <- qq[qq$label %in% aug_cut[[the_key]], ]
    max_node <- sq$node[which.max(sq$y)]
    min_node <- sq$node[which.min(sq$y)]

    # Find the edge of the left-most side of the tree-box
    y_cent <- min(sq$y) + ((max(sq$y) - min(sq$y)) / 2)
    if (ave_clade_size > .05 * max(qq$y)) {
      ave_clade_size <- .05 * max(qq$y)
    }

    # where we want the label to be, i.e., close to its clade, but not overlapping the tree
    max_local_x <- max(qq$x[qq$y < (y_cent + (ave_clade_size / 2)) & qq$y > (y_cent - (ave_clade_size / 2))])
    max_x <- max(qq$x) # Where the label will be shown by default

    #the_x_offset <- max_local_x - .9 * max_x
    the_x_offset <- max_local_x - max_x

    # Don't show the label with the heatmap
    clone_name_label <- paste0(candiate_edges$letters[[i]])
    clone_name_label <- paste0(clone_name_label, " (", format(candiate_edges$frac[[i]], digits = 2), ")")


    p <- p + geom_strip(max_node,
      min_node,
      barsize = 0,
      label = clone_name_label,
      color = candiate_edges$colours[[i]],
      offset = (label_offset + the_x_offset),
      fontsize = fontsize,
      angle = angle,
      hjust = 0,
      offset.text = 1
    )
  }
  return(p)
}

annoate_tree_aug <- function(p = NULL,
  the_graph = NULL,
  candiate_edges,
  aug_cut,
  clone_dic = NULL,
  the_height,
  label_offset = .8,
  angle = 0,
  fontsize = 8,
  draw_bar = FALSE,
  show_guide = FALSE,
  show_label = TRUE,
  branch_width = 1
  ) {

  if (is.null(candiate_edges)) {
    candiate_edges <- get_dummy_candid_edge(the_graph, aug_cut)
  }

  # A hack to make colouring work
  orig_p <- p

  # TODO: Sohrab, sort this out...
  res <- graph_to_tree_dic(the_graph)

  tree <- res$tree
  node_names <- res$node_names


  # Cut the tail
  if (!is.null(the_height)) {
    tree_res <- trim_tree_before_height(tree, node_names, the_graph, the_height)
    tree <- tree_res$tree
  }

  myColors <- get_cluster_colours(length(aug_cut))
  if(length(myColors)>length(aug_cut)){  # in case there are only 2 levels but brewer.pal(nClusters, "Set2") generate at least 3 levels
    myColors <- myColors[1:length(aug_cut)]
  }
  names(myColors) <- names(aug_cut)

  # TODO: Since 0 is shared with the heatmap CNV state, use something else for the Un-assigned ones
  myColors <- c(myColors, "0" = "black")

  cdat <- data.frame(clone_id = names(myColors), colours = myColors, stringsAsFactors = FALSE)

  # Group 0, the default of the tree should be black
  if (is.null(p)) {
    xtree <- ggtree::groupOTU(tree, aug_cut)
    p <- ggtree(xtree, aes(color = group), size = branch_width)
  }


  # TODO: where are the f..ing clone names? Why aren't they shown
  p <- p + scale_color_manual(values = myColors, guide = show_guide)

  qq <- p$data
  qq <- qq[qq$isTip == TRUE, ]

  compute_mean_clade_size <- function(candiate_edges, qq) {
    clade_size <- c()
    for (i in seq_along(candiate_edges$clone_id)) {
      the_key <- candiate_edges$clone_id[[i]]
      sq <- qq[qq$label %in% aug_cut[[the_key]], ]
      clade_size <- c(clade_size, max(sq$y) - min(sq$y))
    }
    min(clade_size)
  }

  ave_clade_size <- compute_mean_clade_size(candiate_edges, qq)

  for (i in seq_along(candiate_edges$clone_id)) {
    # Find the two ends of this pseudo-clade
    # TODO: handle the case where the two ends of the pseudo-clade run over another clade
    the_key <- candiate_edges$clone_id[[i]]
    sq <- qq[qq$label %in% aug_cut[[the_key]], ]
    max_node <- sq$node[which.max(sq$y)]
    min_node <- sq$node[which.min(sq$y)]

    # Find the edge of the left-most side of the tree-box
    y_cent <- min(sq$y) + ((max(sq$y) - min(sq$y)) / 2)
    if (ave_clade_size > .05 * max(qq$y)) {
      ave_clade_size <- .05 * max(qq$y)
    }

    # where we want the label to be, i.e., close to its clade, but not overlapping the tree
    max_local_x <- max(qq$x[qq$y < (y_cent + (ave_clade_size / 2)) & qq$y > (y_cent - (ave_clade_size / 2))])
    max_x <- max(qq$x) # Where the label will be shown by default

    #the_x_offset <- max_local_x - .9 * max_x
    the_x_offset <- max_local_x - max_x

    # Don't show the label with the heatmap
    clone_name_label <- paste0(candiate_edges$letters[[i]])
    if (is.null(orig_p)) {
      clone_name_label <- paste0(clone_name_label, " (", format(candiate_edges$frac[[i]], digits = 2), ")")
    }

    if (!draw_bar) {
      if (show_label) {
        p <- p + geom_strip(max_node,
          min_node,
          barsize = 0,
          label = candiate_edges$letters[[i]],
          color = candiate_edges$colours[[i]],
          offset = (label_offset + the_x_offset),
          fontsize = fontsize,
          angle = angle,
          hjust = 0,
          offset.text = -.5
          )

        p <- p + geom_strip(max_node,
          min_node,
          barsize = 0,
          label = sprintf("(%.2f)", candiate_edges$frac[[i]]),
          color = candiate_edges$colours[[i]],
          offset = (label_offset + the_x_offset + 5),
          fontsize = (fontsize - 3),
          angle = angle,
          hjust = 0,
          offset.text = -.5
          )
      }
    } else {
      p <- p + geom_strip(max_node,
        min_node,
        barsize = 2,
        label = clone_name_label,
        color = candiate_edges$colours[[i]],
        offset = (label_offset + the_x_offset),
        fontsize = fontsize,
        angle = angle,
        hjust = 0,
        offset.text = 1
      )
    }
  }


  if (!show_label) {
    # Add clone counts
    ntotal <- g.count.cells(the_graph)
    counts <- unlist(lapply(names(aug_cut), function(x) length(aug_cut[[x]])))
    clone_counts <- data.frame(clone_id = names(aug_cut),
      counts = counts,
      frac = (counts / ntotal),
      stringsAsFactors = FALSE
    )

    # Add a legend plot
    mode_clone_dic <- dplyr::left_join(clone_counts, cdat)
    mode_clone_dic$color <- mode_clone_dic$colours

    mode_clone_dic$letters <- sprintf("%s (n = %d, %.1f%%)",
      mode_clone_dic$clone_id,
      mode_clone_dic$counts,
      (mode_clone_dic$frac * 100)
    )

    main <- sprintf("Total cells = %d", sum(mode_clone_dic$counts))

    # Add as another row
    mode_clone_dic <- dplyr::bind_rows(data.frame(letters = main, color = "black"), mode_clone_dic)

    pg <- legend_plot(clone_dic = mode_clone_dic, vertical = TRUE, font.size = 3, point_size = 2)
    pg <- pg + theme(plot.subtitle = element_text(size = 10), plot.margin = margin(0,0,0,0))

    p <- ggdraw() +
    draw_plot(p) +
    draw_plot(pg, x = -1.2, y = .3, width = 2.4, scale = .45)
  }
}

launch_run <- function(newick, copy_number, cell_clones, output) {

  tree <- read.tree(newick)
  the_graph <- read_ltm_tree(tree_2_edge_list(tree))
  
  # Hoa, modification to new encoding version
  copynumber <- read.csv(copy_number, header=T, row.names=1, check.names = F,stringsAsFactors = FALSE)
  # copynumber <- read.delim(copy_number, check.names = FALSE, stringsAsFactors = FALSE, sep = ",")
  
  copynumber <- format_copynumber_matrix(copynumber)
  cellclones <- read.delim(cell_clones, stringsAsFactors = FALSE, sep = ",")

  # TODO: check inputs
  # 1. Are all the cells in the matrix?
  # 2. Are there double assignments?

  # TODO: input clone colours

  # What does this section do?
  # aug_cut <- data_frame_to_list(cell_clones)
  # fast_tree_aug(g = g, aug_cut = aug_cut)

  fitclone_plot_tree_heatmap(the_graph = the_graph,
    cell_clones = cellclones,
    output = output,
    copy_number_matrix = copynumber,
    just_tree_loci = FALSE,
    chr_filter = NULL
  )
  print("DONE!!")
}

option_list <- list(make_option(c("-t", "--newick"),
  type = "character",
  default = NULL,
  help = "output from corrupt_tree",
  metavar = "character"
  ),
  make_option(c("-n", "--copynumber"),
    type = "character",
    default = NULL,
    help = "filtered_cell_cn",
    metavar = "character"
  ),
  make_option(c("-c", "--cellclones"),
    type = "character",
    default = NULL,
    help = "output from tree_cutting",
    metavar = "character"
  ),
  make_option(c("-o", "--output"),
    type = "character",
    default = NULL,
    help = "tree_heatmap.png",
    metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
launch_run(opt$newick, opt$copynumber, opt$cellclones, opt$output)
