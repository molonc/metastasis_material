# Parsing inputs ------------------------------------------------------------------------------------------------------
# read_ltm_tree <- function(edge_list) {
#     # Find the root
#     g <- igraph::graph_from_edgelist(as.matrix(edge_list))
#     V(g)$id <- seq(vcount(g))
#     return(g)
# }

tree_2_edge_list <- function(tree) {
    tree$node.label[1] <- "root"
    node_names <- c(tree$tip.label, tree$node.label)

    edges <- tree$edge
    edges <- data.frame(source = node_names[edges[, 1]],
                        target = node_names[edges[, 2]],
                        stringsAsFactors = FALSE
    )
    edges$source <- gsub("cell_", "", edges$source)
    edges$target <- gsub("cell_", "", edges$target)

    return(edges)
}

format_copynumber_matrix <- function(copynumber) {
    # Uncomment if using Tyler's format, Comment if using SS's format
    # bin_ids <- paste(copynumber$chr, copynumber$start, copynumber$end, sep="_")
    # rownames(copynumber) <- bin_ids
    # copynumber <- subset(copynumber, select=-c(chr, start, end, width))
    copynumber <- as.matrix(copynumber)
    return(copynumber)
}

parse_bin_names <- function(bin_names) {
    # Remove corrupt_tree locus tag if it's there
    bin_names <- gsub("locus_", "", bin_names)
    chr <- gsub("([0-9]+|X|Y)_[0-9]+_[0-9]+", "\\1", bin_names)
    start <- as.numeric(gsub("([0-9]+|X|Y)_([0-9]+)_[0-9]+", "\\2", bin_names))
    end <- as.numeric(gsub("([0-9]+|X|Y)_([0-9]+)_([0-9]+)", "\\3", bin_names))
    data.frame(chr = chr, start = start, end = end)
}



get_pad_bin_index <- function(bin_names) {
    bin_dat <- parse_bin_names(bin_names)
    which(abs(bin_dat$start - bin_dat$end) == 1)
}

data_frame_to_list <- function(clustering) {
    clusts <- unique(clustering$clone_id)
    res <- list()
    for (cl in clusts) {
        res[[cl]] <- clustering$cell_id[clustering$clone_id == cl]
    }
    return(res)
}

order_by_chr_name <- function(the_array) {
    the_array[the_array == "X"] <- 44
    the_array[the_array == "Y"] <- 100

    order(as.numeric(the_array))
}

strip_chr_names <- function(the_array) {
    chunks <- unlist(strsplit(the_array, "_"))
    chrs <- matrix(chunks, nrow = length(the_array), ncol = 3, byrow = TRUE)[, 1]
    uchr <- unique(chrs)

    for (cc in uchr) {
        cc_idx <- which(chrs == cc)
        the_array[cc_idx] <- ""
        the_array[cc_idx[1]] <- cc
    }
    return(the_array)
}

# Matrix functions ------------------------------------------------------------------------------------------------------
sort_mat_by_bins <- function(copy_number_matrix) {
    # prevent scientific notation
    options(scipen = 999)
    options(stringsAsFactors = FALSE)
    cnv_txt <- parse_bin_names(rownames(copy_number_matrix))
    copy_number_matrix <- cbind(cnv_txt, copy_number_matrix)

    # Sort the matrix by their chromosome (just inside the chromosome)
    copy_number_matrix$chr[copy_number_matrix$chr == "X"] <- "40"
    copy_number_matrix$chr[copy_number_matrix$chr == "Y"] <- "55"
    copy_number_matrix$chr <- as.numeric(copy_number_matrix$chr)
    copy_number_matrix <- copy_number_matrix[order(copy_number_matrix$chr, copy_number_matrix$start), ]
    copy_number_matrix$chr[copy_number_matrix$chr == "40"] <- "X"
    copy_number_matrix$chr[copy_number_matrix$chr == "55"] <- "Y"

    # Remove chr, start, end
    copy_number_matrix$chr <- NULL
    copy_number_matrix$start <- NULL
    copy_number_matrix$end <- NULL

    return(copy_number_matrix)
}

remove_pad_bins <- function(copy_number_matrix) {
    pad_index <- get_pad_bin_index_for_mat(copy_number_matrix = copy_number_matrix)
    if (length(pad_index) > 0) {
        copy_number_matrix <- copy_number_matrix[-pad_index, ]
    }

    return(copy_number_matrix)
}

get_pad_bin_index_for_mat <- function(copy_number_matrix) {
    dd <- dim(copy_number_matrix)
    if (max(dd) / min(dd) > 1000) {
        bin_names <- unique(copy_number_matrix$loci)
    } else {
        bin_names <- rownames(copy_number_matrix)
    }

    get_pad_bin_index(bin_names)
}

prepare_mat_by_chr <- function(copy_number_matrix, update_col_names = TRUE) {
    the_array <- colnames(copy_number_matrix)
    chunks <- unlist(strsplit(the_array, "_"))
    chrs <- matrix(chunks, nrow = length(the_array), ncol = 3, byrow = TRUE)[, 1]
    unique(chrs)
    copy_number_matrix <- copy_number_matrix[, order_by_chr_name(chrs)]
    if (update_col_names) {
        colnames(copy_number_matrix) <- strip_chr_names(the_array = colnames(copy_number_matrix))
    }

    return(copy_number_matrix)
}

# Tree functions ------------------------------------------------------------------------------------------------------
convert_edge_list_to_ape <- function(edge_list_mat = NULL) {
    options(stringsAsFactors = FALSE)
    edge_list <- as.data.frame(edge_list_mat)
    colnames(edge_list) <- c("source", "target")

    leaf_names <- get_leaves_names(edge_list)
    root_name <- get_root_name(edge_list)
    all_nodes <- unique(c(edge_list$source, edge_list$target))

    internal_nodes <- setdiff(all_nodes, leaf_names)
    internal_nodes <- setdiff(internal_nodes, root_name)

    n_leaves <- length(leaf_names)
    n_internal <- length(internal_nodes) + 1

    d1 <- data.frame(node_name = leaf_names, node_number = seq(n_leaves))
    d2 <- data.frame(node_name = root_name, node_number = (n_leaves + 1))

    # Account for a start where the only loci is the root
    if (length(internal_nodes) == 0) {
        d <- rbind(d1, d2)
    } else {
        d3 <- data.frame(node_name = internal_nodes, node_number = seq(n_internal - 1))
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

    edges <- as.matrix(edges)
    colnames(edges) <- NULL

    d <- d[order(d$node_number), ]
    tr <- list(edge = edges, tip.label = d$node_name[1:n_leaves], Nnode = n_internal)
    class(tr) <- "phylo"

    Nedge <- nrow(tr$edge)
    tr$edge.length <- rep(1, Nedge)
    list(tree = tr, node_names = d)
}

trim_tree_before_height <- function(tree, node_names, g, the_height, after_height = NULL) {
    h <- get_height_dat(g)
    max_height <- max(h$dist)
    h <- h[h$dist < the_height, ]

    # Find edge row_id in ape format
    qq <- tree$edge
    pp <- node_names[order(node_names$node_number), ]
    q1 <- pp$node_name[qq[, 1]]
    q2 <- pp$node_name[qq[, 2]]
    named_edges <- matrix(c(q1, q2), ncol = 2, byrow = FALSE)
    named_edges <- as.data.frame(named_edges)
    colnames(named_edges) <- c("source", "target")
    tree_edge_ids <- which((named_edges$source %in% h$id) & (named_edges$source %in% h$id))

    tree$edge.length[tree_edge_ids] <- .01
    if (!is.null(after_height)) {
        tree$edge.length[setdiff(seq(tree$edge.length), tree_edge_ids)] <- after_height
    }

    list(tree = tree, height = max_height, node_names = node_names)
}

get_leaves_names <- function(edge_list) {
    leaves <- unique(edge_list$target)[!(unique(edge_list$target) %in%  unique(edge_list$source))]
    target_freq <- as.data.frame(table(edge_list$target))
    stopifnot(all(target_freq$Freq[target_freq$Var1 %in% leaves] == 1))
    return(leaves)
}

get_root_name <- function(edge_list) {
    unique(edge_list$source)[!(unique(edge_list$source) %in% unique(edge_list$target))]
}

# Graph functions ------------------------------------------------------------------------------------------------------
get_height_dat <- function(the_graph) {
    the_root <- find_root(the_graph)
    height_search <- bfs(the_graph,
                         root = c(the_root),
                         order = TRUE,
                         neimode = "out",
                         unreachable = FALSE,
                         dist = TRUE,
                         rank = TRUE,
                         succ = TRUE,
                         father = TRUE,
                         pred = TRUE
                         )
    res <- data.frame(id = names(height_search$dist), dist = as.numeric(height_search$dist), stringsAsFactors = FALSE)
    return(res)
}

find_root <- function(the_graph) {
    the_root <- NULL
    n_nodes <- vcount(the_graph)
    for (i in seq(n_nodes)) {
        s1 <- bfs(the_graph, root = i, order = TRUE, neimode = "out", unreachable = FALSE)
        if (any(is.na(s1$order)) == FALSE)  {
            the_root <- i
        }
    }
    if (is.null(the_root)) {
        print("Warning! No root was found...")
    }
    return(the_root)
}

fast_read_ltm_tree <- function(edge_list) {
  # Find the root
  g <- igraph::graph_from_edgelist(as.matrix(edge_list))
  V(g)$id <- seq(vcount(g))
  return(g)
}
find_all_parents <- function(the_graph, node, max_height) {
    igraph::neighborhood(the_graph, nodes = node, mode = "in", order = max_height)[[1]]$name[-c(1)]
}

find_tail_height <- function(the_graph = NULL) {
    h <- get_height_dat(the_graph)
    res <- convert_edge_list_to_ape(edge_list_mat = igraph::as_edgelist(the_graph))

    tree <- res$tree
    node_names <- res$node_names
    p <- ggtree(tree) + theme_tree2()
    qq <- p$data
    qq <- qq[qq$isTip == TRUE, ]

    # Ignore SA928 gm cells
    qq <- qq[!grepl("SA928", qq$label), ]
    uml <- qq$label[which.max(qq$y)]
    lml <- qq$label[which.min(qq$y)]
    mrca <- find_MRCA(the_graph = the_graph, the_nodes = c(uml, lml), height_dat = h)
    mrca_node <- node_names$node_number[node_names$node_name == mrca]
    tail_length <- h$dist[h$id == mrca]

    result <- list(h_cut = NULL, tail_length = tail_length, nodes = c(uml, lml, mrca))

    if (tail_length > 15) {
        result$h_cut <- tail_length - 5
    }

    return(result)
}

graph_to_tree_dic <- function(the_graph) {
    edge_list <- as.data.frame(igraph::as_edgelist(the_graph), stringsAsFactors = FALSE)
    colnames(edge_list) <- c("source", "target")
    res <- convert_edge_list_to_ape(edge_list_mat = edge_list)
    return(res)
}

find_MRCA <- function(the_graph, the_nodes, height_dat) {
    if (the_nodes[1] == the_nodes[2]) {
        return(the_nodes[1])
    }
    if (height_dat$id[1] == "root" | (height_dat$dist[1] < height_dat$dist[nrow(height_dat)])) {
        print("Warning! Sorting height dat in decreasing order...")
        height_dat <- height_dat[order(height_dat$dist, decreasing = TRUE), ]
    }
    parents <- list()
    for (node in the_nodes) {
        print(node)
        i <- which(height_dat$id == node)
        # Append the node itself, so if it's a parent of the others, it's taken care of
        parents[[node]] <- c(find_all_parents(the_graph, height_dat$id[i], height_dat$dist[i]), height_dat$id[i])
    }
    # TODO: support multiple subgroups
    cp <- intersect(parents[[1]], parents[[2]])

    height_dat$id[height_dat$id %in% cp][1]
}

get_dummy_candid_edge <- function(the_graph, aug_cut) {
    ntotal <- g.count.cells(the_graph)
    frac <- unlist(lapply(names(aug_cut), function(x) (length(aug_cut[[x]]) / ntotal)))
    data.frame(clone_id = names(aug_cut), frac = frac, stringsAsFactors = FALSE)
}

g.count.cells <- function(the_graph = NULL) {
    length(g.get.cells(the_graph))
}

g.get.cells <- function(the_graph = NULL) {
    get_cells(V(the_graph)$name)
}

get_cells <- function(str_arry) {
    str_arry[!grepl("locus_", str_arry) & !grepl("root", str_arry)]
}

# Plotting ------------------------------------------------------------------------------------------------------
heatmap_mat_from_bin_cellID <- function(res_pick, update_colnames = TRUE) {
    the_mat <- t(res_pick)
    the_mat <- prepare_mat_by_chr(the_mat, update_colnames)
    if (length(which(unname(is.na(the_mat[2, ])))) > 0) {
        the_mat <- the_mat[, -c(which(unname(is.na(the_mat[2, ]))))]
    }
    return(the_mat)
}

legacy_colours_for_CN_heatmap <- function() {
    # From CN = 0 to CN = 11
    c("#4880B8",
      "#A7C9DF",
      "#CCCCCC",
      "#F5CE93",
      "#ED9364",
      "#D2553E",
      "#A42116",
      "#8B1A43",
      "#CB3576",
      "#D06CAD",
      "#C196C4",
      "#D0BAD8"
      )
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
    theme(legend.position = "none")
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
        ll.dat <- data.frame(x = c(rep(1, length(clone_dic$letters))),
                             y = rev(seq_along(clone_dic$letters)),
                             label = clone_dic$letters,
                             stringsAsFactors = FALSE
                             )
    } else {
        ll.dat <- data.frame(x = seq_along(clone_dic$letters),
                             y = c(rep(1, length(clone_dic$letters))),
                             label = clone_dic$letters,
                             stringsAsFactors = FALSE
                             )
    }

    # Better results with spaces
    ll.dat$label <- paste0("   ", ll.dat$label)
    clone_dic$letters <- paste0("   ", clone_dic$letters)

    nudge_x <- 0
    nudge_y <- 0
    if (vertical) {
        #nudge_y = .5
        nudge_x <- 0
    } else {
        nudge_x <- .5
    }

    myColors <- clone_dic$color
    names(myColors) <- clone_dic$letters
    pg <- ggplot(data = ll.dat, mapping = aes(x = x, y = y, color = label)) +
          geom_text(aes(label = label),
                    colour = "black",
                    nudge_x = nudge_x,
                    nudge_y = nudge_y,
                    fontface = "bold",
                    size = font.size,
                    hjust = 0
                    ) +
          fitclone_get_theme_no_grid() +
          fitclone_get_theme_no_axis() +
          fitclone_get_theme_no_legend() +
          scale_color_manual(values = myColors) +
          theme(panel.background = element_blank(), rect = element_rect(fill = "transparent"))

    if (point_size > 0) {
        ll.dat.point <- ll.dat[-c(1), ]
        pg <- pg + geom_point(data = ll.dat.point, size = point_size, shape = 15)
    }

    if (!is.null(main)) {
        pg <- pg + ggtitle(label = "", subtitle = main)
    }

    if (vertical) {
        #pg <- pg + xlim(c(1-.1, 1+.1))
        pg <- pg + xlim(c(1 - nudge_x, 1 + nudge_x * 2))
    } else {
        pg <- pg + ylim(c(1 - .1, 1 + .1))
    }

    return(pg)
}

get_cluster_colours <- function(nClusters) {
    if (nClusters > 8) {
        clust_colours <- colorRampPalette(brewer.pal(8, "Set2"))(nClusters)
    } else {
        clust_colours <- brewer.pal(nClusters, "Set2")
    }
    return(clust_colours)
}


