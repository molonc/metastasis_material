# source("/datadrive/tree_viz/ss_viz/general_utils.R")
# source("/datadrive/tree_viz/ss_viz/tree_and_heatmap_utils_utils.R")

### Matrix functions
################################################################################################################################################
sort_mat_by_bins <- function(the_mat) {
    # prevent scientific notation
    options(scipen=999) 
    chr = gsub('([0-9]+|X|Y)_[0-9]+_[0-9]+', '\\1', rownames(the_mat))
    start = as.numeric(gsub('([0-9]+|X|Y)_([0-9]+)_[0-9]+', '\\2', rownames(the_mat)))
    cnv_txt = data.frame(chr=chr, start=start, stringsAsFactors = F)
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

convert_edge_list_to_ape <- function(edge_list_path, edge_list_mat=NULL) {
    options(stringsAsFactors = F)
    if (!is.null(edge_list_mat)) {
        edge_list = as.data.frame(edge_list_mat)
    } else {
        edge_list <- read.csv(edge_list_path, header = T)
    }
    
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
    
    edges = edge_list
    
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

find_tail_height <- function(g = NULL, edge_list_path) {
    if (!is.null(edge_list_path)) {
        outpath <- sprintf('%s_tail_height.rds', edge_list_path) 
        if (file.exists(outpath))
            return(readRDS(outpath))
    }
    
    if (is.null(g)) {
        g <- read_ltm_tree(edge_list_path)
    }

    h <- get_height_dat(g, edge_list_path)
    if (!is.null(edge_list_path)) {
        res <- convert_edge_list_to_ape(edge_list_path)
    } else {
        res <- graph_to_tree_dic(g)
    }
    
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
    if (!is.null(edge_list_path)) {
        saveRDS(result, outpath)
    }
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
    
    tree$edge.length[tree_edge_ids] = .01
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

sos_heat <- function (p, data, offset = 0, width = 1, low = "green", high = "red", 
                                            color = "white", colnames = TRUE, colnames_position = "bottom", 
                                            colnames_angle = 0, colnames_level = NULL, colnames_offset_x = 0, 
                                            colnames_offset_y = 0, font.size = 4, hjust = 0.5, colnames_filter = NULL) 
{
    colnames_position <- colnames_position %>% match.arg(c("bottom", "top"))
    variable <- value <- lab <- y <- NULL
    width <- width * (p$data$x %>% range(na.rm = TRUE) %>% diff)/ncol(data)
    isTip <- x <- y <- variable <- value <- from <- to <- NULL
    df <- p$data
    df <- df[df$isTip, ]
    start <- max(df$x, na.rm = TRUE) + offset
    dd <- as.data.frame(data)
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
    }
    else {
        dd$variable <- factor(dd$variable, levels = colnames_level)
    }

    V2 <- start + as.numeric(dd$variable) * width
    mapping <- data.frame(from = dd$variable, to = V2)
    mapping <- unique(mapping)
    dd$x <- V2
    dd$width <- width

    if (is.null(color)) {
        p2 <- p + geom_tile(data = dd, aes(x, y, fill = value), width = width, inherit.aes = FALSE)
    }
    else {
        p2 <- p + geom_tile(data = dd, aes(x, y, fill = value), width = width, color = color, inherit.aes = FALSE)
    }
    if (is(dd$value, "numeric")) {
        p2 <- p2 + scale_fill_gradient(low = low, high = high, na.value = NA)
    }
    else {
        p2 <- p2 + scale_fill_discrete(na.value = NA)
    }

    if (colnames) {
        if (colnames_position == "bottom") {
            y <- 0
        }
        else {
            y <- max(p$data$y) + 1
        }
        mapping$y <- y
        if (!is.null(colnames_filter)) {
            mapping = mapping[mapping$from %in% colnames_filter$from, ]
            mapping <- dplyr::right_join(mapping, colnames_filter)
            mapping$from <- mapping$new_label
        }
        
        p2 <- p2 + geom_text(data = mapping, aes(x = to, y = y, label = from), size = font.size, inherit.aes = FALSE, 
                                                 angle = colnames_angle, nudge_x = colnames_offset_x, 
                                                 nudge_y = colnames_offset_y, hjust = hjust)
    }

    p2 <- p2 + theme(legend.position = "right", legend.title = element_blank())
    attr(p2, "mapping") <- mapping
    return(p2)
}

annoate_tree_aug <- function(p=NULL, g=NULL, edge_list_path=NULL, candiate_edges, aug_cut, clone_dic=NULL, the_height, label_offset = .8, angle=0, fontsize=8, draw_bar=FALSE, show_guide=FALSE, show_label=TRUE, branch_width = 1) {
    if (is.null(clone_dic)) {
        clone_dic <- get_clone_dic_for_seq(names(aug_cut))
    }
    
    if (is.null(candiate_edges)) {
        candiate_edges <- get_dummy_candid_edge(g, aug_cut)
    }
    
    # A hack to make colouring work
    orig_p <- p
    
    if (is.null(g)) {
        res <- convert_edge_list_to_ape(edge_list_path)
    } else  {
        res <- graph_to_tree_dic(g)
    }
    tree <- res$tree
    node_names <- res$node_names
    
    
    # Cut the tail
    if (!is.null(the_height)) {
        tree_res <- trim_tree_before_height(tree, node_names, g, the_height)
        tree <- tree_res$tree
    }
    
    if (is.null(clone_dic$color)) {
        myColors <- get_cluster_colours(length(clone_dic$letters))
        names(myColors) <- clone_dic$old_K

        # Drop the one not in candi_edges
        myColors <- myColors[names(myColors) %in% unique(candiate_edges$clone_id)]
    } else {
        myColors <- clone_dic$color
        names(myColors) <- clone_dic$old_K
        clone_dic$color <- NULL
        clone_dic$colours <- NULL
        clone_dic$colour <- NULL
    }
        
    # TODO: Since 0 is shared with the heatmap CNV state, use something else for the Un-assigned ones
    myColors <- c(myColors, '0'='black')
    
    cdat <- data.frame(clone_id=names(myColors), colours=myColors, stringsAsFactors = F)
    cdat <- dplyr::right_join(cdat, clone_dic, by=c('clone_id'='old_K'))
    candiate_edges <- dplyr::right_join(candiate_edges, cdat)

    # Group 0, the default of the tree should be black
    if (is.null(p)) {
        xtree <- groupOTU(tree, .node=aug_cut)
        p <- ggtree(xtree, aes(color=group), size = branch_width)
    }

    
    # TODO: where are the f..ing clone names? Why aren't they shown 
    # TODO: this maybe inconsistent with the clone_dic... ; removing the confusing guide
    if (show_guide) {
        p <- p + scale_color_manual(values = myColors)
    } else {
        p <- p + scale_color_manual(values = myColors, guide=FALSE)
    }
    
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
        ntotal <- sum(get_tcc(g)$Freq)
        counts <- unlist(lapply(names(aug_cut), function(x) length(aug_cut[[x]])))
        clone_counts <- data.frame(clone_id = names(aug_cut), counts = counts, stringsAsFactors = F)
        
        # Add a legend plot
        mode_clone_dic <- dplyr::left_join(candiate_edges[, c('frac', 'letters', 'colours')], clone_dic)
        mode_clone_dic$color <- mode_clone_dic$colours
        mode_clone_dic <- dplyr::left_join(mode_clone_dic, clone_counts, by=c('old_K'='clone_id'))
        
        mode_clone_dic$letters <- sprintf('%s (n = %d, %.1f%%)', mode_clone_dic$letters, mode_clone_dic$counts, mode_clone_dic$frac*100)
        main <- sprintf('Total cells = %d', sum(mode_clone_dic$counts))
        
        # Add as another row
        mode_clone_dic <- dplyr::bind_rows(data.frame(letters = main, color = 'black'), mode_clone_dic)
        
        pg <- legend_plot(clone_dic = mode_clone_dic, vertical = TRUE, font.size = 3, point_size = 2)
        pg <- pg + theme(plot.subtitle = element_text(size = 10), plot.margin = margin(0,0,0,0))
        
        if ( !(datatag %in% c('SA609', 'SA532', 'SA609X3X8a'))) {
            if (datatag == 'SA906b') {
                p <- ggdraw() +
                    draw_plot(p) +
                    draw_plot(pg, x = -.7, y = .3, width = 2.4, scale = .45)
            } else {
                p <- ggdraw() +
                    draw_plot(p) +
                    draw_plot(pg, x = -.7, y = -.12, width = 2.4, scale = .5)
            }
        } else {
            if (datatag == 'SA609X3X8a') {
                p <- ggdraw() +
                    draw_plot(p) +
                    draw_plot(pg, x = -1.1, y = .25, width = 2.4, scale = .5)
            } else {
                p <- ggdraw() +
                    draw_plot(p) +
                    draw_plot(pg, x = -1.1, y = .3, width = 2.4, scale = .45)
            }
        }
    }
    
    # TODO: use height of the tree and the_height argument to fix this 
    # if (is.null(orig_p))  {
    #     if (!is.null(the_height)) {
    #         # new_x_lim <- (the_height * .1) + (tree_res$height - the_height) + 5
    #         # p + theme(plot.margin=unit(c(1,.5,1,1),"cm")) + ggplot2::xlim(0, new_x_lim)
    #         p
    #     }
    #     else  {
    #     p
    #     }
    # } else {
    #     # Add some room for the labels
    #     # coef <- .2
    #     # plot_width <- max(p$dat$x) - min(p$dat$x) 
    #     # lim_bottom_x <- max(p$dat$x) + coef*plot_width
    #     # 
    #     # p + xlim(c(0, lim_bottom_x))
    #     p
    # }
}


### Fitclone_plot_tree_heatmap_more dependancies
################################################################################################################################################
get_title_str <- function(datatag, n_cells, clone_name = NULL) {
    if (is.null(clone_name)) 
        sprintf('%s (n = %d)', datatag, n_cells)
    else 
        sprintf('%s - Clone %s (n = %d)', datatag, clone_name, n_cells)
}

get_ID_from_edge_list_path <- function(edge_list_path, keep_timestamp = TRUE) {
    if (!is.null(edge_list_path)) {
        if (edge_list_path == '') return('')
        tokens <- strsplit(file_path_sans_ext(basename(edge_list_path)), '__')[[1]]
        temp_str <- tokens[[2]]
        timestamp <- tokens[[3]]
        get_pretty_str_from_funny_str(input_str = paste0('__', temp_str, '__', timestamp), keep_timestamp)
    } else {
        get_pretty_str_from_funny_str(generate_random_funny_string(), keep_timestamp)
    }
}

get_file_name_postfix <- function(datatag, edge_list_path, tag = '', file_extension = '', add_rnd_str = TRUE) {
    if (edge_list_path == '') {
        id = ''
    } else {
        id <- strsplit(file_path_sans_ext(basename(edge_list_path)), '__')[[1]]
        id <- id[[2]]
    }
    rnd_str <-ifelse(add_rnd_str,  paste0('_', generate_random_str()), '')
    sprintf('%s_%s%s__%s%s', tag, datatag, rnd_str, id, file_extension)
}
