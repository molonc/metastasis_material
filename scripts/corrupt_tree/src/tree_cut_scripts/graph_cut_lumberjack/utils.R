## Installation
## See notes at run_treecut_parameters.R 
## How to run program, see main() function at run_treecut_parameters.R
suppressPackageStartupMessages({
    require("dplyr")
    require("ggtree")
    require("ape")
    require("gtools")
    require("igraph")
    require("data.table")
})


# Parsing inputs ------------------------------------------------------------------------------------------------------
tree_2_edge_list <- function(tree) {
    tree$node.label[1] <- "root"
    node_names <- c(tree$tip.label, tree$node.label)

    edges <- tree$edge
    edges <- data.frame(source = node_names[edges[, 1]], target = node_names[edges[, 2]], stringsAsFactors = FALSE)
    edges$source <- gsub("cell_", "", edges$source)
    edges$target <- gsub("cell_", "", edges$target)

    return(edges)
}

format_copynumber_matrix <- function(copynumber) {
	# Leave commented out if using SS's format, uncomment if using Tyler's format 
    # bin_ids <- paste(copynumber$chr, copynumber$start, copynumber$end, sep="_")
    # rownames(copynumber) <- bin_ids
    # copynumber <- subset(copynumber, select=-c(chr, start, end, width))
    copynumber <- as.matrix(copynumber)
    return(copynumber)
}

list_to_data_frame <- function(clust_list, clone_dic = NULL) {
    # TODO: incorporate the clone_dic
    dat <- data.frame(cell_id = unlist(clust_list, use.names = FALSE), clone_id = NA, stringsAsFactors = FALSE)
    for (nl in names(clust_list)) {
        dat$clone_id[dat$cell_id %in% clust_list[[nl]]] <- nl
    }

    if (!is.null(clone_dic)) {
        dat <- dplyr::left_join(dat, clone_dic, by = c("clone_id"="old_K"))
        dat$clone_id <- dat$letters
        dat <- dat[, c("cell_id", "clone_id")]
    }
    return(dat)
}

remove_na_from_list <- function(some_list) {
    # Why doesn't it keep names?
    tmp <- lapply(names(some_list), function(x) {some_list[[x]][!is.na(some_list[[x]])]})
    names(tmp) <- names(some_list)
    return(tmp)
}

# Tree functions ------------------------------------------------------------------------------------------------------
read_ltm_tree <- function(edge_list) {
    # Find the root
    g <- igraph::graph_from_edgelist(as.matrix(edge_list))
    V(g)$id <- seq(vcount(g))
    return(g)
}

get_leaves_names <- function(edge_list) {
    leaves <- unique(edge_list$target)[!(unique(edge_list$target) %in% unique(edge_list$source))]
    target_freq <- as.data.frame(table(edge_list$target))
    stopifnot(all(target_freq$Freq[target_freq$Var1 %in% leaves] == 1))
    return(leaves)
}

get_decendents <- function(clone_roots, the_graph, min_cell_per_clone = 1) {
    decendents <- list()
    for (clone_root in clone_roots) {
        decend_node_ids <- decendents_ids(clone_root, the_graph)
        cells_in_clone <- length(decend_node_ids)
        if (cells_in_clone > min_cell_per_clone) {
            decendents[[clone_root]] <- decend_node_ids
        }
    }

    return(decendents)
}

decendents_ids <- function(root_id, the_graph) {
    decend_node_ids <- bfs(the_graph,
                           root = c(root_id), 
                           order = TRUE, 
                           neimode = "out",
                           unreachable = FALSE
                           )

    decend_node_ids <- as.numeric(decend_node_ids$order)
    decend_node_ids <- decend_node_ids[!is.na(decend_node_ids)]
    return(decend_node_ids)
}

# Graph functions ------------------------------------------------------------------------------------------------------
g.count.cells <- function(the_graph = NULL, edge_list_path = NULL) {
    length(g.get.cells(the_graph, edge_list_path))
}

g.get.cells <- function(the_graph = NULL, edge_list_path = NULL) {
    if (is.null(the_graph)) {
        the_graph <- read_ltm_tree(edge_list_path = edge_list_path)
    }
    get_cells(V(the_graph)$name)
}

get_cells <- function(str_arry) {
    str_arry[!grepl("locus_", str_arry)]
}

node_number_to_label <- function(the_graph, groups) {
    all_dat <- data.frame(id = V(the_graph)$id, label = V(the_graph)$name, stringsAsFactors = FALSE)

    res <- list()
    for (clade_root in names(groups)) {
        cells_in_clone <- groups[[clade_root]]
        # find the labels
        cells <- all_dat$label[match(c(-1, cells_in_clone), all_dat$id)]
        res[[clade_root]] <- cells
    }

    return(res)
}

retrieve_time_point_per_clade <- function(the_graph, decendents, trim_cell_name = TRUE) {
    if (is.numeric(unlist(decendents))) {
        decendents <- node_number_to_label(the_graph, decendents)
    }

    res <- list()
    for (clade_root in names(decendents)) {
        cells <- decendents[[clade_root]]
        if (trim_cell_name) {
            cells <- cells[!grepl("locus_", cells)]
            cells <- cells[cells != clade_root]
        }
        res[[clade_root]] <- cells
    }

    return(res)
}

find_all_parents <- function(the_graph, node, max_height) {
    igraph::neighborhood(the_graph, nodes = node, mode = "in", order = max_height)[[1]]$name[-c(1)]
}

find_root <- function(the_graph) {
    the_root <- NULL
    n_nodes <- vcount(the_graph)
    for (i in seq(n_nodes)) {
        s1 <- bfs(the_graph, root = i, order = TRUE, neimode = "out", unreachable = FALSE)
        if (any(is.na(s1$order)) == FALSE) {
            the_root <- i
        }
    }
    if (is.null(the_root)) {
        print("Warning! No root was found...")
    }
    return(the_root)
}

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

find_parent <- function(the_graph, node) {
    igraph::neighborhood(the_graph, nodes = node, mode = "in")[[1]]$name[[2]]
}

find_immediate_children <- function(the_graph, node) {
    children <- igraph::neighborhood(the_graph, nodes = node, mode = "out", order = 1)[[1]]$name
    children[!grepl(node, children)]
}

get_desc_names <- function(the_graph, the_node) {
    desc <- get_decendents(the_node, the_graph = the_graph)
    remove_na_from_list(node_number_to_label(the_graph, desc))
}

get_non_overlapping_desc <- function(the_graph, parent, node, cell_in_cut) {
    children <- find_immediate_children(the_g = the_graph, node = parent)
    children <- children[children != node & !is.na(children)]
    new_desc <- c()
    for (cc in children) {
        hh <- get_desc_names(the_graph = the_graph, the_node = cc)
        if (length(hh) > 0) {
            hh <- hh[[1]]
        }
        hh <- hh[!is.na(hh)]

        if (length(hh) == 0) {
            new_desc <- c(new_desc, cc)
        } else if (!any(cell_in_cut %in% hh)) {
            new_desc <- c(new_desc, cc, hh)
        }
    }
    return(new_desc)
}

get_dfs_queue <- function(the_graph, edge) {
    s1 <- dfs(the_graph, root = edge, order = TRUE, neimode = "out", unreachable = FALSE)
    qq <- na.omit(s1$order)
    unlist(lapply(seq(length(qq)), function(x) qq[x]$name))
}

setup_graph <- function(the_graph = NULL) {
    V(the_graph)$id <- seq(vcount(the_graph))
    h <- get_height_dat(the_graph)
    stopifnot(all(h$id == V(the_graph)$name))
    V(the_graph)$height <- as.integer(h$dist)
    return(the_graph)
}

get_leaves_names_from_graph <- function(the_graph, only_loci = FALSE) {
    edge_list <- as.data.frame(igraph::as_edgelist(the_graph), stringsAsFactors = FALSE)
    colnames(edge_list) <- c("source", "target")
    res <- get_leaves_names(edge_list)
    if (only_loci) {
        res <- grep("locus_", res, value = TRUE)
    }
    return(res)
}

find_all_descendents <- function(the_graph, node, max_height) {
    igraph::neighborhood(the_graph, nodes = node, mode = "out", order = max_height)[[1]]$name[-c(1)]
}

# Naming functions ------------------------------------------------------------------------------------------------------
get_clone_dic_for_seq <- function(a_seq) {
    data.frame(old_K = a_seq,
               is_ref = FALSE, 
               letters = LETTERS[seq(length(a_seq))],
               pretty_names = get_pretty_names_for_loci(a_seq), 
               stringsAsFactors = FALSE
               )
}

get_pretty_names_for_loci <- function(str_array) {
    str_array <- gsub("locus_", "", str_array)
    chr <- gsub("([0-9]+|X|Y)_.*", "\\1", str_array)
    str_array <- gsub(".*_([0-9]+)_.*", "\\1", str_array)
    str_array <- substr(str_array, 1, 4)
    str_array <- paste0("chr", chr, "_", str_array, "")
    gsub("chrroot_root", "root", str_array)
}

rename_clones <- function(cell_clones) {
    cell_clones$clone_id <- as.integer(factor(cell_clones$clone_id))

    if (length(unique(cell_clones$clone_id)) <= 26) {
        cell_clones$clone_id <- LETTERS[cell_clones$clone_id]
    }
    return(cell_clones)
}

# Misc functions ------------------------------------------------------------------------------------------------------
median_genotype_from_list <- function(clone_list, copy_number_matrix) {
    clones <- names(clone_list)
    nClones <- length(clones)
    res <- matrix(NA, nrow = nClones, ncol = nrow(copy_number_matrix))
    rownames(res) <- clones

    for (clone_i in seq_along(clones)) {
        clone <- clones[[clone_i]]
        sub_mat <- copy_number_matrix[, colnames(copy_number_matrix) %in% clone_list[[clone]]]
        res[clone_i, ] <- apply(as.matrix(sub_mat), 1, median)
    }
    colnames(res) <- rownames(copy_number_matrix)
    return(res)
}

compute_dist_mat <- function(mg_mat, use_hamming = FALSE) {
    # print("Testing")
    if(is.null(mg_mat)){
        # print("Testing 1")
        return(NULL)
    } else if(nrow(mg_mat)==1){
        # print("Testing 2")
        return(matrix(0, nrow = 1, ncol = 1))
    } else{
        # print("Testing 3")
        n_clones <- nrow(mg_mat)
        dist_to_median <- matrix(0, nrow = n_clones, ncol = n_clones)
        
        for (i in seq(n_clones)) {
            for (j in seq(n_clones)) {
                if (i < j) {
                    if (use_hamming) {
                        dist_to_median[i, j] <- mean(mg_mat[i, ] != mg_mat[j, ])
                    } else {
                        dist_to_median[i, j] <- mean(abs(mg_mat[i, ] - mg_mat[j, ]))
                    }
                }
            }
        }
        return(dist_to_median)
    }
    
}

is_outlier <- function(x, prob, type = "t") {
    tt <- outliers::scores(x = x, type = type, prob = prob)
    ff <- outliers::scores(x = x, type = type) < 0
    tt & ff
}

merge_clones <- function(cdat) {
    list_of_clones <- list()
    for (i_cn in seq(nrow(cdat))) {
        xx <- unname(unlist(cdat[i_cn, , drop = TRUE]))
        p1 <- who_has_it(xx[1], list_of_clones)
        p2 <- who_has_it(xx[2], list_of_clones)
        if (is.na(p1) & is.na(p2)) {
            list_of_clones[[length(list_of_clones) + 1]] <- xx
        } else if (is.na(p1)) {
            list_of_clones[[p2]] <- unique(c(list_of_clones[[p2]], xx))
        } else if (is.na(p2)) {
            list_of_clones[[p1]] <- unique(c(list_of_clones[[p1]], xx))
        } else {
            if (p1 != p2) {
                list_of_clones[[p1]] <- unique(c(list_of_clones[[p1]], list_of_clones[[p2]]))
                list_of_clones[[p2]] <- NULL
            }
            list_of_clones[[p1]] <- unique(c(list_of_clones[[p1]], xx))
        }
    }
    return(list_of_clones)
}

who_has_it <- function(it, list_of_clones) {
    for (i in seq_along(list_of_clones)) {
        if (it %in% list_of_clones[[i]]) {
            return(i)
        }
    }
    return(NA)
}

is_consecuitive <- function(a_list) {
    # TODO: check for B_p...

    if (length(a_list) < 2) {
        return(TRUE)
    }

    a_list <- sort(a_list)
    i1 <- a_list[[1]]
    i2 <- a_list[[length(a_list)]]

    ideal_string <- gsub(sprintf(".*(%s.*%s).*", i1, i2), "\\1", paste0(LETTERS, collapse = ""))
    if (nchar(ideal_string) != length(a_list)) {
        return(FALSE)
    }

    return(TRUE)
}

find_non_consecutive_ones <- function(a_list, orig_name) {
    if (length(a_list) < 2) {
        return(NULL)
    }
    a_list <- sort(a_list)
    i1 <- a_list[[1]] # start
    i2 <- a_list[[length(a_list)]] # end

    ideal_string <- gsub(sprintf(".*(%s.*%s).*", i1, i2), "\\1", paste0(LETTERS, collapse = ""))
    ideal_list <- strsplit(ideal_string, "")[[1]]

    # Use the non-existing ones as delimiter to divide
    # -- all after should not appear
    del <- ideal_list[which(!(ideal_list %in% a_list))]
    if (length(del) > 1) {
        print("WARNING! NOT IMPLEMENTED!!!")
    }

    del <- sort(del)
    new_list <- list()
    remaining_list <- ideal_list
    for (ii in del) {
        ii_i <- which(remaining_list == ii)
        if (length(ii_i) > 0) {
            if (ii_i == 1) {
                # don't add this one
            } else {
                new_name <- sprintf("%s_(%s)", orig_name, paste0(remaining_list[1:(ii_i - 1)], collapse = "+"))
                new_list[[new_name]] <- remaining_list[1:(ii_i - 1)]
            }
            remaining_list <- remaining_list[-c(1:ii_i)]
        } else {
            print("ERROR!")
        }
    }
    if (length(remaining_list) > 0) {
        new_name <- sprintf("%s_(%s)", orig_name, paste0(remaining_list, collapse = "+"))
        new_list[[new_name]] <- remaining_list
    }

    return(new_list)
}


