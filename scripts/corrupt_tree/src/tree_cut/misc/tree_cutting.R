#!/usr/bin/env Rscript

library(ape)
library(argparse)
library(gtools)
library(igraph)

get_args <- function() {
    descr <- "Cut a corrupt tree to produce cell-clone assignments tsv file."
    p <- ArgumentParser(description=descr)

    p$add_argument("newick", help="Corrupt tree input newick file")
    p$add_argument(
        "copynumber",
        help=paste(
            "Copynumber input tsv file with columns for bins",
            "[chrom, start, end, width] and the rest should be copynumber",
            "for each cell."
        )
    )
    p$add_argument("output", help="Cell clone output tsv file")

    p$add_argument(
        "--minimum-fraction", type="double", default=0.02,
        help="Minimum fraction of cells in initial clone selection"
    )
    p$add_argument(
        "--maximum-fraction", type="double", default=0.38,
        help="Maximum fraction of cells in initial clone selection"
    )

    return(p$parse_args())
}

format_copynumber_matrix <- function(copynumber) {
    bin_ids <- paste(
        copynumber$chr, copynumber$start, copynumber$end, sep="_")

    rownames(copynumber) <- bin_ids
    copynumber <- subset(copynumber, select=-c(chr, start, end, width))
    copynumber <- as.matrix(copynumber)

    return(copynumber)
}

split_merge_drop_clones <- function(internal_dist, sub_subclone_list,
                                    prob=0.75, min_clone_fraction=0.01,
                                    min_merge_back_fraction=0.05) {
    # Check for duplicates
    stopifnot(!any(duplicated(unlist(unname(sub_subclone_list)))))

    # 1.a. Use t-test scores to merge subclades within each clade
    new_clades <- list()
    for (cn in unique(internal_dist$original_clone)) {
        sub_internal_dist <- internal_dist[internal_dist$original_clone == cn, ]
        if (nrow(sub_internal_dist) > 1) {
            tt <- is_outlier(x=sub_internal_dist$sub_dist, prob=prob)
        } else {
            tt <- is_outlier(x=internal_dist$sub_dist, prob=prob)
            if (cn %in% internal_dist$original_clone[tt]) {
                tt <- c(TRUE)
            } else {
                tt <- c(FALSE)
            }
        }

        sub_rows <- sub_internal_dist[tt, c(1, 2)]
        all_subclones <- unique(c(
            sub_internal_dist$src_clone, sub_internal_dist$trg_clone))

        if (nrow(sub_rows) > 0) {
            tmp_clones <- unique(c(sub_rows$src_clone, sub_rows$trg_clone))
            non_merge_subclones <- setdiff(all_subclones, tmp_clones)
            new_clades[[cn]] <- list(
                merged_clones=list(merge_clones(sub_rows)),
                indiv_clones=non_merge_subclones
            )
        } else {
            new_clades[[cn]] <- list(indiv_clones=all_subclones)
        }
    }
    # str(new_clades)

    new_new_clades <- new_clades
    # 1.a.0 Do not merge subclones that envelope a third one
    for (cn in names(new_clades)) {
        for (ii in new_clades[[cn]]$merged_clones) {
            i_index <- 0
            for (jj in ii) {
                i_index <- i_index + 1
                if (!is_consecuitive(a_list=jj)) {
                    # Remove the offending ones
                    the_new_split <- find_non_consecutive_ones(jj, cn)
                    # print('the_new_split')
                    # print(the_new_split)
                    # remove this group from merged groups
                    new_new_clades[[cn]]$merged_clones[[1]][[i_index]] <- NA
                    # for each one in the new group, either add them as a new
                    # merged or a new nonmerged
                    for (qq in names(the_new_split)) {
                        if (nchar(qq) > 5) { # e.g., B_(E)
                            ff <- unname(the_new_split[qq])
                            new_new_clades[[cn]]$merged_clones[[1]] <- append(
                                new_new_clades[[cn]]$merged_clones[[1]], ff)
                        } else {
                            # print(sprintf('Adding to indiv %s', qq))
                            # print(the_new_split[[qq]])
                            new_new_clades[[cn]]$indiv_clones <- append(
                                new_new_clades[[cn]]$indiv_clones,
                                the_new_split[[qq]]
                            )
                        }
                    }
                    # print(sprintf('----------------------- begin{State at the end of %s} --------------------', cn))
                    # print(new_new_clades[[cn]])
                    # print(sprintf('----------------------- end{State at the end of %s} --------------------', cn))
                }
            }
        }

        # print(str(new_new_clades[[cn]]))
    }

    new_clades <- new_new_clades


    # 1.b. Sow the merged and unmerged clones together (since the output of 1.a. is involved)
    new_aug <- list()
    for (cn in names(new_clades)) {
        # Add the non_merged clones
        for (ii in new_clades[[cn]]$indiv_clones) {
            new_name <- paste0(cn, '_', ii)
            new_aug[[new_name]] <- sub_subclone_list[[cn]][[ii]]
        }

        # Now add the merged ones
        for (ii in new_clades[[cn]]$merged_clones) {
            # There may be multiple merged subgroups
            for (jj in ii) {
                new_name <- paste0(cn, '_(', paste0(jj, collapse='+'), ')')
                for (kk in jj) {
                    new_aug[[new_name]] <- c(new_aug[[new_name]], sub_subclone_list[[cn]][[kk]])
                }
            }
        }
    }
    # str(new_aug)

    # Check for duplicates
    stopifnot(!any(duplicated(unlist(unname(new_aug)))))

    # 2. Drop stand-alone clades that are less than min_clone_fraction
    n_cells <- length(unname(unlist(new_aug)))
    lens <- unlist(lapply(new_aug, function(x) length(x))) / n_cells
    drop_clones <- names(lens)[lens < min_clone_fraction]
    if (length(drop_clones) > 0) {
        print('Dropping '); print(drop_clones)
        new_aug[drop_clones] <- NULL
    }

    # Check for duplicates
    stopifnot(!any(duplicated(unlist(unname(new_aug)))))

    new_new_aug <- new_aug
    # 3. Merge the too small ones back to the main clone
    print('ERROR! NOT MERGING TINY CLONES...')
    # tiny_clones <- names(lens)[lens < min_merge_back_fraction]
    # for (tc in tiny_clones) {
    #     print(sprintf('Merging tiny %s', tc))
    #     tokens <- strsplit(tc, '_')[[1]]
    #     orig_clone_name <- lens[grep(paste0(tokens[[1]], '_'), names(lens))]
    #     orig_clone_name <- names(orig_clone_name)[which.max(orig_clone_name)]
    #     if (orig_clone_name != tc) {
    #         new_new_aug[[orig_clone_name]] <- c(
    #             new_new_aug[[orig_clone_name]], new_aug[[tc]])
    #         new_new_aug[[tc]] <- NULL
    #     }
    # }

    # str(new_new_aug)

    # Check for double assignment
    stopifnot(!any(duplicated(unlist(unname(new_new_aug)))))

    new_new_aug
}

read_ltm_tree <- function(edge_list) {
    # Find the root
    g <- igraph::graph_from_edgelist(as.matrix(edge_list))
    V(g)$id <- seq(vcount(g))
    return(g)
}

g.count.cells <- function(g=NULL, edge_list_path=NULL) {
    length(g.get.cells(g, edge_list_path))
}

g.get.cells <- function(g=NULL, edge_list_path=NULL) {
    if (is.null(g)) g <- read_ltm_tree(edge_list_path=edge_list_path)
    get_cells(V(g)$name)
}

get_cells <- function(str_arry) str_arry[!grepl('locus_', str_arry)]

parse_cell_names <- function(cell_names) {
    sample_library_ids <- sapply(strsplit(cell_names, "-"), function(x) {
        if (length(x) >= 2) return(paste0(x[[1]], "-", x[[2]]))
        return(x)
    })
    return(sample_library_ids)
}

decendents_ids <- function(root_id, the_graph) {
    decend_node_ids <- bfs(
        the_graph,
        root=c(root_id), order=TRUE, neimode='out',
        unreachable=FALSE
    )
    decend_node_ids <- as.numeric(decend_node_ids$order)
    decend_node_ids <- decend_node_ids[!is.na(decend_node_ids)]
    return(decend_node_ids)
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

    return(decendents)
}

node_number_to_label <- function(the_graph, groups) {
    all_dat <- data.frame(
        id=V(the_graph)$id, label=V(the_graph)$name, stringsAsFactors=FALSE
    )

    res <- list()
    for (clade_root in names(groups)) {
        cells_in_clone <- groups[[clade_root]]
        # find the labels
        cells <- all_dat$label[match(c(-1, cells_in_clone), all_dat$id)]
        res[[clade_root]] <- cells
    }

    return(res)
}

retrieve_time_point_per_clade <- function(the_graph, decendents,
                                          trim_cell_name=TRUE) {
    if (is.numeric(unlist(decendents))) {
        decendents <- node_number_to_label(the_graph, decendents)
    }

    res <- list()
    for (clade_root in names(decendents)) {
        cells <- decendents[[clade_root]]
        if (trim_cell_name) {
            cells <- cells[!grepl('locus_', cells)]
            cells <- cells[cells != clade_root]
        }
        res[[clade_root]] <- cells
    }

    return(res)
}

order_cols <- function(ts) {
    # TODO fix this
    x_cols <- grep('X[0-9][0-9]*', colnames(ts))
    temp <- ts[, x_cols, drop=FALSE]
    temp <- temp[, mixedorder(colnames(temp)), drop=FALSE]
    temp$clone_id <- ts$clone_id
    temp <- temp[, c(ncol(temp), seq(ncol(temp) - 1)), drop=FALSE]
    return(temp)
}

get_time_series <- function(clade_desc_list, timepoints=NULL) {
    # Also take the clade root into account
    time_series_dat <- NULL
    for (clade_root in names(clade_desc_list)) {
        break_down <- table(clade_desc_list[[clade_root]])
        temp <- data.frame(
            clone_id=clade_root, as.list(break_down), stringsAsFactors=FALSE)

        if (is.null(time_series_dat)) {
            time_series_dat <- temp
        } else {
            time_series_dat <- plyr::rbind.fill(time_series_dat, temp)
        }
    }

    for (i in seq(nrow(time_series_dat))) {
        time_series_dat[i, is.na(time_series_dat[i, ])] <- 0
    }

    # Ensure all time points exits
    if (!is.null(timepoints)) {
        for (xx in timepoints) {
            if (!(xx %in% colnames(time_series_dat))) {
                time_series_dat[[xx]] <- 0
            }
        }
        time_series_dat_cols <- c(
            1, 1 + order(colnames(time_series_dat)[2:ncol(time_series_dat)]))
        time_series_dat <- time_series_dat[, time_series_dat_cols]
    }

    for (clade_root in names(clade_desc_list)) {
        root_label <- parse_cell_names(clade_root)
        # TODO fix this
        if (!grepl('SA', root_label)) {
            if (root_label %in% colnames(time_series_dat)) {
                clade_root_mask <- time_series_dat$clone_id == clade_root
                time_series_dat[
                    clade_root_mask, root_label
                ] <- time_series_dat[clade_root_mask, root_label] + 1
            }
        }
    }

    order_cols(ts=time_series_dat)
}

convert_to_freq <- function(time_series_dat, total_cell_counts_breakdown) {
    if (is.null(total_cell_counts_breakdown)) {
        # print('Normalising to themselves')
        for (i in 2:ncol(time_series_dat)) {
            the_sum <- sum(time_series_dat[, i])
            if (the_sum > 0) {
                time_series_dat[, i] <- time_series_dat[, i] / the_sum
            }
        }
    } else {
        # print('Normalising to Total count')
        # the first column is clone_id
        for (i in 2:ncol(time_series_dat)) {
            time_series_dat[, i] <- time_series_dat[, i] /
                total_cell_counts_breakdown$Freq[i - 1]
        }
    }

    time_series_dat
}

get_time_series_for_cluster_cut <- function(res, total_cell_counts_break_down,
                                            convert_to_freq=TRUE) {
    time_series_dat <- get_time_series(
        clade_desc_list=res,
        timepoints=as.character(total_cell_counts_break_down$time)
    )

    if (convert_to_freq) {
        time_series_dat <- convert_to_freq(
            time_series_dat=time_series_dat,
            total_cell_counts_breakdown=total_cell_counts_break_down
        )
    }

    return(order_cols(time_series_dat))
}


find_all_parents <- function(the_g, node, max_height) {
    igraph::neighborhood(
        the_g, nodes=node, mode='in', order=max_height
    )[[1]]$name[-c(1)]
}


find_root <- function(g) {
    the_root <- NULL
    n_nodes <- vcount(g)
    for (i in seq(n_nodes)) {
        s1 <- bfs(g, root=i, order=TRUE, neimode='out', unreachable=FALSE)
        if (any(is.na(s1$order)) == FALSE) {
            the_root <- i
        }
    }
    if (is.null(the_root)) {
        print('Warning! No root was found...')
    }
    the_root
}

get_height_dat <- function(g) {
    the_root <- find_root(g)
    height_search <- bfs(
        g,
        root=c(the_root), order=TRUE, neimode='out', unreachable=FALSE,
        dist=TRUE, rank=TRUE, succ=TRUE, father=TRUE, pred=TRUE
    )
    res <- data.frame(
        id=names(height_search$dist),
        dist=as.numeric(height_search$dist),
        stringsAsFactors=FALSE
    )
    return(res)
}

get_edge_data <- function(g) {
    # Number of edges=Number of nodes - number of leaves
    # TODO: generalise this...

    # TODO: @Sohrab: remove this
    # out_ts <- sprintf('%s_%s_ts_dat.rds', './outputs/', edge_density(g))
    # if (file.exists(out_ts)) return(readRDS(out_ts))

    n_nodes <- length(V(g)$name)
    n_leaves <- length(V(g)$name[
        !grepl('root', V(g)$name) & !grepl('locus', V(g)$name)
    ])
    n_edges <- n_nodes - n_leaves

    # internal_nodes
    i_nodes <- V(g)$name[grepl('root', V(g)$name) | grepl('locus', V(g)$name)]

    # Account for GM cells
    # TODO REMOVE
    GM_nodes <- grep('SA928', V(g)$name, value=TRUE)

    # Sanity check
    # stopifnot(sum(tcc$Freq) == (n_leaves - length(GM_nodes)))

    # Find all descendents
    desc <- get_decendents(
        clone_roots=i_nodes, the_graph=g, min_cell_per_clone=1)
    # Find out which one is a cells and which is a locus (will be NA)
    res <- retrieve_time_point_per_clade(the_graph=g, decendents=desc)

    # Remove all NAs entries in res
    for (cn in names(res)) {
        if (all(is.na(res[[cn]]))) {
            res[[cn]] <- NULL
        } else {
            res[[cn]] <- res[[cn]][!is.na(res[[cn]])]
        }
    }

    # Count the number of cells in each subgroup
    ts <- lapply(names(res), function(x) {
        length(res[[x]])
    })
    ts <- data.frame(clone_id=names(res), N=unlist(ts), stringsAsFactors=FALSE)
    stopifnot(ts$N[ts$clone_id == 'root'] == n_leaves)

    ts$frac <- ts$N / n_leaves
    stopifnot(max(ts$frac) <= 1 & min(ts$frac) >= 0)

    # Add height that
    height_dat <- get_height_dat(g)
    height_dat <- height_dat[height_dat$id %in% ts$clone_id, ]
    colnames(height_dat) <- c('clone_id', 'height')

    ts <- dplyr::right_join(ts, height_dat, by=c('clone_id'))

    ts <- ts[order(ts$height), ]
    ts$id <- paste0(
        '(height=', ts$height, ', frac=', format(ts$frac, digits=2), ')'
    )

    # TODO:@Sohrab remove this
    # saveRDS(ts, out_ts)

    return(ts)
}

collapse_tree_by_edge <- function(the_g, nodes, do_plot=FALSE,
                                  annotate_nodes=NULL) {
    h <- get_height_dat(the_g)
    h <- h[h$id %in% nodes, ]
    h <- h[order(h$dist, decreasing=TRUE), ]

    # The neighbour with the smallest height is the parent
    edge_list <- matrix(
        NA, nrow=nrow(h), ncol=2, dimnames=list(c(), c('source', 'target')))
    for (i in seq(nrow(h))) {
        p <- find_all_parents(the_g, h$id[i], h$dist[i])
        parent <- 'root'
        if (nrow(h[h$id %in% p, ]) > 0)
            parent <- h$id[h$id %in% p][1]

        edge_list[i, ] <- c(parent, h$id[i])
    }

    orig_edge_list <- edge_list

    # Convert cell name to clone number
    clone_g <- igraph::graph_from_edgelist(el=edge_list)
    gg <- igraph::graph_from_edgelist(el=orig_edge_list)
    V(gg)$id <- seq(vcount(gg))
    list(edge_list=orig_edge_list, graph=gg)
}

h1_auto_cut <- function(the_g, candiate_edges, minimum_fraction,
                        maximum_fraction=0.5) {

    # Sort edges by height
    candiate_edges <- candiate_edges[
        order(candiate_edges$height, decreasing=FALSE), ]
    tau <- collapse_tree_by_edge(the_g=the_g, nodes=candiate_edges$clone_id)

    qeueu <- get_dfs_queue(tau$graph, 'root')
    # Remove the root
    qeueu <- qeueu[-c(1)]
    covered_nodes <- list()
    the_cut <- list()
    for (edge in qeueu) {
        if (edge %in% covered_nodes) next
        tmp <- candiate_edges[candiate_edges$clone_id == edge, ]

        # If passed
        filter <- (tmp$frac < maximum_fraction)

        if (filter) {
            the_cut <- append(the_cut, edge)
            # Remove all children from the queue
            covered_nodes <- append(
                covered_nodes,
                find_all_descendents(
                    tau$graph, edge, max(candiate_edges$height)
                )
            )
        } else {
            # If failed at a leaf, pick the first node that is a child of a
            # parent with 2 or more children
            if (length(find_all_descendents(tau$graph, edge, 1)) == 0) {
                t1 <- find_all_parents(tau$graph, edge, tmp$height)
                t2 <- candiate_edges[candiate_edges$clone_id %in% t1, ]
                degrees <- igraph::degree(graph=tau$graph, v=t1, mode='out')
                degrees <- data.frame(
                    clone_id=names(degrees), degree=as.vector(degrees))
                t2 <- dplyr::right_join(t2, degrees)
                t2 <- t2[order(t2$height, decreasing=TRUE), ]
                branch_parent_index <- which(t2$degree > 1)[1]
                if (branch_parent_index > 1) {
                    edge <- t2[branch_parent_index - 1, ]$clone_id
                    the_cut <- append(the_cut, edge)
                    # Remove all children from the queue
                    covered_nodes <- append(
                        covered_nodes,
                        find_all_descendents(
                            tau$graph, edge, max(candiate_edges$height)
                        )
                    )
                }
                # Otherwise ignore this part of the tree
            }
        }
    }
    return(the_cut)
}

get_the_cut <- function(g, minimum_fraction, maximum_fraction) {
    ts <- get_edge_data(g)

    candiate_edges <- ts[
        ts$frac > minimum_fraction,
    ]

    the_cut <- h1_auto_cut(
        the_g=g, candiate_edges=candiate_edges,
        minimum_fraction=minimum_fraction, maximum_fraction=maximum_fraction
    )

    return(the_cut)
}

find_parent <- function(the_g, node) {
    igraph::neighborhood(the_g, nodes=node, mode='in')[[1]]$name[[2]]
}

find_immediate_children <- function(the_g, node) {
    children <- igraph::neighborhood(
        the_g, nodes=node, mode='out', order=1)[[1]]$name
    children[!grepl(node, children)]
}

remove_na_from_list <- function(some_list) {
    # Why doesn't it keep names?
    tmp <- lapply(names(some_list), function(x) {
        some_list[[x]][!is.na(some_list[[x]])]
    })
    names(tmp) <- names(some_list)
    return(tmp)
}

get_desc_names <- function(the_graph, the_node) {
    desc <- get_decendents(the_node, the_graph=the_graph)
    remove_na_from_list(node_number_to_label(the_graph, desc))
}

get_non_overlapping_desc <- function(the_graph, parent, node, cell_in_cut) {
    children <- find_immediate_children(the_g=the_graph, node=parent)
    children <- children[children != node & !is.na(children)]
    new_desc <- c()
    for (cc in children) {
        hh <- get_desc_names(the_graph=the_graph, the_node=cc)
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

add_siblings_cut <- function(the_cut, the_g, minimum_fraction,
                             maximum_fraction, oscilation) {
    # Now use the manually curated clones. Drop the subclades that threaten
    # the purity of the clade
    ts <- get_edge_data(the_g)
    # heights
    q <- ts %>%
        dplyr::filter(clone_id %in% unlist(the_cut)) %>%
        dplyr::select(clone_id, height) %>%
        dplyr::arrange(desc(height)) %>%
        as.data.frame()

    # Add the descendents
    new_cut <- get_desc_names(the_g, the_cut)

    # The clones should not include anyone else or any children that are a
    # parent to another existing clone
    cell_in_cut <- unlist(unname(get_desc_names(
        the_graph=the_g, the_node=the_cut)))
    cell_in_cut <- cell_in_cut[!is.na(cell_in_cut)]

    # TODO: Check and remove all immediate children that are a
    # parent of someone else
    for (i in seq(nrow(q))) {
        node <- q$clone_id[i]; node_height <- q$height[i]
        parent <- find_parent(the_g, node=node)

        # Don't overright a parent that has already been processed
        # TODO: should we generlaly break here?
        if (parent %in% names(new_cut)) next

        p_desc <- get_non_overlapping_desc(
            the_graph=the_g, parent=parent, node=node, cell_in_cut=cell_in_cut
        )
        if (length(p_desc) != 0) {
            new_cut[[parent]] <- p_desc
            cell_in_cut <- c(cell_in_cut, p_desc)
        }
    }

    # TODO: replace with a black list
    # Remove GM cells
    for (nc in names(new_cut)) {
        new_cut[[nc]] <- new_cut[[nc]][!grepl('SA928', new_cut[[nc]])]
    }

    # Remove loci
    for (nc in names(new_cut)) {
        new_cut[[nc]] <- new_cut[[nc]][!grepl('locus_', new_cut[[nc]])]
    }

    # Remove NAs
    for (nc in names(new_cut)) {
        new_cut[[nc]] <- new_cut[[nc]][!is.na(new_cut[[nc]])]
    }

    # Impose min and max fraction
    new_clades <- setdiff(names(new_cut), unlist(the_cut))
    bad_clades <- c()
    n_cells <- ts$N[ts$clone_id == 'root']
    for (nc in new_clades) {
        frac <- length(new_cut[[nc]]) / n_cells
        if (frac < minimum_fraction) {
            bad_clades <- c(bad_clades, nc)
        }
    }

    if (length(bad_clades) > 0) {
        for (nc in bad_clades) {
            new_cut[nc] <- NULL
        }
    }


    return(new_cut)
}

get_aug_cut <- function(g, minimum_fraction, maximum_fraction) {
    the_cut <- get_the_cut(
        g=g, minimum_fraction=minimum_fraction,
        maximum_fraction=maximum_fraction
    )
    add_siblings_cut(
        the_cut=the_cut, the_g=g, minimum_fraction=minimum_fraction,
        maximum_fraction=maximum_fraction
    )
}

get_full_subgraphs_for_aug_cut <- function(g, aug_cut) {
    # For each clone, find the edge'
    # Gather all loci below immediate children
    children <- list()
    for (ii in names(aug_cut)) {
        children[[ii]] <- find_immediate_children(the_g=g, node=ii)
        children[[ii]] <- children[[ii]][!is.na(children[[ii]])]
        children[[ii]] <- children[[ii]][!(children[[ii]] %in% names(aug_cut))]
    }

    # Purge level 2
    bad_children <- list()
    for (ii in names(children)) {
        just_loci <- children[[ii]][grepl('locus_', children[[ii]])]
        temp <- get_desc_names(the_graph=g, the_node=just_loci)
        bad_children[[ii]] <- c()
        for (tmp in names(temp)) {
            if (any(names(aug_cut) %in% temp[[tmp]])) {
                bad_children[[ii]] <- c(bad_children[[ii]], tmp)
            }
        }
    }
    for (ii in names(bad_children)) {
        children[[ii]] <- children[[ii]][
            !(children[[ii]] %in% bad_children[[ii]])]
    }

    # All loci: children of immediate loci + immediate children that are cells
    desc_loci <- list()
    for (ii in names(children)) {
        loci_children <- grep('locus_', children[[ii]], value=TRUE)
        desc <- unname(unlist(get_desc_names(
            the_graph=g, the_node=loci_children)))
        desc <- desc[!(desc %in% loci_children)]

        # Add all immediate children
        desc <- c(desc, children[[ii]])
        desc_loci[[ii]] <- desc
    }

    # Remove tip loci
    tip_loci <- get_leaves_names_from_graph(the_g=g, only_loci=TRUE)

    if (length(tip_loci) > 0) {
        for (ii in names(desc_loci)) {
            desc_loci[[ii]] <- desc_loci[[ii]][
                !(desc_loci[[ii]] %in% tip_loci)]
        }
    }

    # TODO fix this
    # Remove GM cells
    if (length(tip_loci) > 0) {
        for (ii in names(desc_loci)) {
            desc_loci[[ii]] <- desc_loci[[ii]][
                !(grepl('SA928', desc_loci[[ii]]))]
        }
    }

    # Create subgraphs
    subgraphs <- list()
    for (ii in names(desc_loci)) {
        vids <- c(desc_loci[[ii]], ii)
        subgraphs[[ii]] <- igraph::induced.subgraph(graph=g, vids=vids)
        subgraphs[[ii]] <- setup_graph(g=subgraphs[[ii]])
    }

    subgraphs
}

split_by_genotype <- function(g, minimum_fraction, maximum_fraction,
                              clone_size_threshold, edge_diff_threshold,
                              use_hamming=FALSE) {
    # get the subgraphs for each clade
    aug_cut <- get_aug_cut(g=g, minimum_fraction, maximum_fraction)
    subgraphs <- get_full_subgraphs_for_aug_cut(g, aug_cut)

    cn <- as.list(names(subgraphs))
    names(cn) <- LETTERS[1:length(cn)]

    # check each clone for split
    internal_dist <- NULL
    sub_subclone_list <- NULL
    clone_dic_list <- NULL
    sub_median_mat_list <- NULL

    for (ii in names(subgraphs)) {
        # ii <- cn$D
        gg <- subgraphs[[ii]]

        # Edges that are at least 20% of the clade
        ts <- get_edge_data(g=gg)
        candid_edges <- ts[ts$frac > clone_size_threshold, ]

        # Further filter edges to ignore too close parental edges
        # ss is ordered so that closests to the root is first
        ss <- get_desc_names(gg, as.list(candid_edges$clone_id))

        ss_temp <- ss # Keep the updated fractions here
        total_cells <- g.count.cells(gg)
        n_test_edges <- nrow(candid_edges)
        delta_frac <- matrix(0, nrow=n_test_edges, ncol=n_test_edges)
        for (i in seq(n_test_edges)) {
            for (j in seq(n_test_edges)) {
                if (i < j) {
                    delta_frac[i, j] <- length(setdiff(ss[[i]], ss[[j]])) /
                        total_cells
                }
            }
        }

        # Start from all the hanging-leaf-like loci
        test_candids <- candid_edges$clone_id
        
        # @Debug
        #annotate_nodes_on_tree(p = NULL, edge_list_path = NULL, node_names = test_candids, tree_node_dic = graph_to_tree_dic(gg))

        cc <- candid_edges[, c('clone_id', 'frac', 'height')]

        # While moving up the tree, remove the cells from the frac of a
        # parent that belong to the already seen ones
        bad_cadidates <- c()
        for (ii in rev(test_candids)) {
            p1 <- find_all_parents(gg, ii, max_height=max(candid_edges$height))
            if (length(p1) > 0) {
                p2 <- test_candids[test_candids %in% p1]
                # Find immediate parent
                immediate_parent <- rev(
                    candid_edges$clone_id[candid_edges$clone_id %in% p2])[1]
                i_index <- which(candid_edges$clone_id == ii)
                j_index <- which(candid_edges$clone_id == immediate_parent)
                if (delta_frac[j_index, i_index] < edge_diff_threshold) {
                    bad_cadidates <- c(bad_cadidates, ii)
                } else {
                    # For all of j's children except i, update the difference
                    # by the negative of size of i
                    ch_in_list <- ss[[j_index]][
                        ss[[j_index]] %in% test_candids[-c(i_index)]]
                    if (length(ch_in_list)) {
                        for (x in ch_in_list) {
                            x_index <- which(candid_edges$clone_id == x)
                            delta_frac[j_index, x_index] <-
                                delta_frac[j_index, x_index] -
                                length(ss[[i_index]])
                        }
                    }
                }
            }
        }

        passed_candids <- setdiff(test_candids, bad_cadidates)
        # Sort by height
        passed_candids <- rev(
            candid_edges$clone_id[candid_edges$clone_id %in% passed_candids])
        private_desc <- list()
        i_index <- 0
        for (ii in passed_candids) {
            i_index <- i_index + 1
            if (i_index > 1) {
                private_desc[[ii]] <- setdiff(
                    ss[[ii]], unlist(ss[passed_candids[seq(i_index - 1)]]))
            } else {
                private_desc[[ii]] <- ss[[ii]]
            }
        }

        # Drop all loci
        for (ii in names(private_desc)) {
            private_desc[[ii]] <- private_desc[[ii]][
                !grepl('locus', private_desc[[ii]])]
            if (length(private_desc[[ii]]) == 0) private_desc[[ii]] <- NULL
        }

        # Add subclone letters ...
        clone_dic <- get_clone_dic_for_seq(a_seq=names(private_desc))
        clone_dic_list[[names(cn)[cn == ii]]] <- clone_dic

        # @Debug
        #p <- fast_tree_aug(g = gg, aug_cut = private_desc, outdir = '~/projects/tree_cutting/outputs/debug/July11/', title = names(cn)[cn == ii])
        
        # Keep track of sub_subclones for merging, dividing later on
        tmp_private_desc <- private_desc
        names(tmp_private_desc) <- clone_dic$letters
        sub_subclone_list[[names(cn)[cn == ii]]] <- tmp_private_desc
    }

    list(sub_subclone_list=sub_subclone_list, clone_dic_list=clone_dic_list)
}

split_by_ploidy <- function(sub_subclone_list, clone_dic_list) {
    mat <- load_new_cn_data(datatag)

    new_sub_subclone_list <- list()
    new_clone_dic_list <- list()
    for (cn in names(sub_subclone_list)) {
        private_desc <- sub_subclone_list[[cn]]
        clone_dic <- clone_dic_list[[cn]]
        passed_mg <- median_genotype_from_list(private_desc, mat)
        pa <- get_ploidy_for_mat(t(passed_mg))
        ploidies <- unique(pa$ploidy)

        if (length(ploidies) > 1) {
            for (pl in ploidies) {
                new_sub_subclone_list[[paste0(cn, '-', 'p', pl)]] <-
                    private_desc[which(pa$ploidy == pl)]
                new_clone_dic_list[[paste0(cn, '-', 'p', pl)]] <-
                    clone_dic[which(pa$ploidy == pl), ]
            }
        } else {
            new_sub_subclone_list[[cn]] <- private_desc
            new_clone_dic_list[[cn]] <- clone_dic
        }
    }

    list(
        sub_subclone_list=new_sub_subclone_list,
        clone_dic_list=new_clone_dic_list
    )
}


split_by_genotype_driver <- function(
         g,
         minimum_fraction,
         maximum_fraction,
         mat,
         clone_size_threshold=0.13,
         edge_diff_threshold=0.04,
         adjust_ploidy=FALSE,
         prob=0.75,
         min_merge_back_fraction=0.03) {

    # TODO remove... but removing changes results
    options(stringsAsFactors=FALSE)

    # 1. Split clones by their fractions
    res <- split_by_genotype(
        g,
        minimum_fraction,
        maximum_fraction,
        clone_size_threshold=clone_size_threshold,
        edge_diff_threshold=edge_diff_threshold
    )
    sub_subclone_list <- res$sub_subclone_list
    clone_dic_list <- res$clone_dic_list

    # 2. Split subclone by their ploidy
    if (adjust_ploidy) {
        res <- split_by_ploidy(sub_subclone_list, clone_dic_list)
        sub_subclone_list <- res$sub_subclone_list
        clone_dic_list <- res$clone_dic_list
    }

    # 3. Compute distance between subclones in each clone
    res <- split_compute_intra_dist(
        sub_subclone_list=sub_subclone_list,
        clone_dic_list=clone_dic_list,
        mat=mat
    )
    internal_dist <- res$internal_dist
    sub_median_mat_list <- res$sub_median_mat_list

    # 4. Merge similar clones and drop too small ones
    new_aug <- split_merge_drop_clones(
        internal_dist=internal_dist,
        sub_subclone_list=sub_subclone_list,
        prob=prob,
        min_merge_back_fraction=min_merge_back_fraction
    )

    res <- list(
        new_aug=new_aug,
        internal_dist=internal_dist,
        clone_dic_list=clone_dic_list,
        sub_subclone_list=sub_subclone_list
    )
    return(res)
}

get_dfs_queue <- function(the_g, edge) {
    s1 <- dfs(the_g, root=edge, order=TRUE, neimode='out', unreachable=FALSE)
    qq <- na.omit(s1$order)
    unlist(lapply(seq(length(qq)), function(x) qq[x]$name))
}

get_leaves_names <- function(edge_list) {
    leaves <- unique(edge_list$target)[
        !(unique(edge_list$target) %in% unique(edge_list$source))]
    target_freq <- as.data.frame(table(edge_list$target))
    stopifnot(all(target_freq$Freq[target_freq$Var1 %in% leaves] == 1))
    leaves
}

setup_graph <- function(g=NULL) {
    V(g)$id <- seq(vcount(g))
    h <- get_height_dat(g)
    stopifnot(all(h$id == V(g)$name))
    V(g)$height <- as.integer(h$dist)
    g
}

get_leaves_names_from_graph <- function(the_g, only_loci=FALSE) {
    edge_list <- as.data.frame(
        igraph::as_edgelist(the_g), stringsAsFactors=FALSE)
    colnames(edge_list) <- c('source', 'target')
    res <- get_leaves_names(edge_list)
    if (only_loci) {
        res <- grep('locus_', res, value=TRUE)
    }
    res
}

find_all_descendents <- function(the_g, node, max_height) {
    igraph::neighborhood(
        the_g, nodes=node, mode='out', order=max_height
    )[[1]]$name[-c(1)]
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

split_compute_intra_dist <- function(sub_subclone_list, clone_dic_list, mat,
                                     use_hamming=FALSE) {
    sub_median_mat_list <- list()
    internal_dist <- list()
    for (cn in names(sub_subclone_list)) {
        private_desc <- sub_subclone_list[[cn]]
        clone_dic <- clone_dic_list[[cn]]
        passed_mg <- median_genotype_from_list(clone_list=private_desc, mat)
        sub_median_mat_list[[cn]] <- passed_mg

        # Compute intra distance matrix
        passed_mg_dist <- compute_dist_mat(
            mg_mat=passed_mg, use_hamming=use_hamming)

        # Keep track of distance between subclones
        if (length(private_desc) > 1 & sum(passed_mg_dist) > 0) {
            tmp_indx <- which(passed_mg_dist > 0, arr.ind=TRUE)
            tmp_internal_dist <- data.frame(
                src_clone=clone_dic$letters[tmp_indx[, 1]],
                trg_clone=clone_dic$letters[tmp_indx[, 2]],
                sub_dist=passed_mg_dist[tmp_indx],
                original_clone=cn
            )
        } else {
            tmp_internal_dist <- data.frame(
                src_clone=names(private_desc),
                trg_clone=names(private_desc),
                sub_dist=NA,
                original_clone=cn
            )
        }

        if (is.null(internal_dist)) {
            internal_dist <- tmp_internal_dist
        } else {
            internal_dist <- rbind(internal_dist, tmp_internal_dist)
        }
    }

    list(
        internal_dist=internal_dist,
        sub_median_mat_list=sub_median_mat_list,
        sub_subclone_list=sub_subclone_list,
        clone_dic_list=clone_dic_list
    )
}

get_pretty_names_for_loci <- function(str_array) {
    str_array <- gsub('locus_', '', str_array)
    chr <- gsub('([0-9]+|X|Y)_.*', '\\1', str_array)
    str_array <- gsub('.*_([0-9]+)_.*', '\\1', str_array)
    str_array <- substr(str_array, 1, 4)
    str_array <- paste0('chr', chr, '_', str_array, '')
    gsub('chrroot_root', 'root', str_array)
}

get_clone_dic_for_seq <- function(a_seq) {
    data.frame(
        old_K=a_seq,
        is_ref=FALSE, letters=LETTERS[seq(length(a_seq))],
        pretty_names=get_pretty_names_for_loci(a_seq), stringsAsFactors=FALSE
    )
}

median_genotype_from_list <- function(clone_list, mat) {
    clones <- names(clone_list)
    nClones <- length(clones)
    res <- matrix(NA, nrow=nClones, ncol=nrow(mat))
    rownames(res) <- clones

    for (clone_i in seq_along(clones)) {
        clone <- clones[[clone_i]]
        sub_mat <- mat[, colnames(mat) %in% clone_list[[clone]]]
        res[clone_i, ] <- apply(as.matrix(sub_mat), 1, median)
    }
    colnames(res) <- rownames(mat)
    res
}


compute_dist_mat <- function(mg_mat, use_hamming=FALSE) {
    n_clones <- nrow(mg_mat)
    dist_to_median <- matrix(0, nrow=n_clones, ncol=n_clones)
    for (i in seq(n_clones)) {
        for (j in seq(n_clones)) {
            if (i < j) {
                if (use_hamming) {
                    dist_to_median[i, j] <- mean(mg_mat[i, ] != mg_mat[j, ])
                } else {
                    dist_to_median[i, j] <- mean(
                        abs(mg_mat[i, ] - mg_mat[j, ]))
                }
            }
        }
    }
    dist_to_median
}

is_outlier <- function(x, prob, type='t') {
    tt <- outliers::scores(x=x, type=type, prob=prob)
    ff <- outliers::scores(x=x, type=type) < 0
    tt & ff
}


merge_clones <- function(cdat) {
    list_of_clones <- list()
    for (i_cn in seq(nrow(cdat))) {
        xx <- unname(unlist(cdat[i_cn, , drop=TRUE]))
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
                list_of_clones[[p1]] <- unique(c(
                    list_of_clones[[p1]], list_of_clones[[p2]]))
                list_of_clones[[p2]] <- NULL
            }
            list_of_clones[[p1]] <- unique(c(list_of_clones[[p1]], xx))
        }
    }
    list_of_clones
}

who_has_it <- function(it, list_of_clones) {
    for (i in seq_along(list_of_clones)) {
        if (it %in% list_of_clones[[i]])
            return(i)
    }
    return(NA)
}

is_consecuitive <- function(a_list) {
    # TODO: check for B_p...

    if (length(a_list) < 2) return(TRUE)

    a_list <- sort(a_list)
    i1 <- a_list[[1]]
    i2 <- a_list[[length(a_list)]]

    ideal_string <- gsub(sprintf(
        '.*(%s.*%s).*', i1, i2), '\\1', paste0(LETTERS, collapse=''))
    if (nchar(ideal_string) != length(a_list))
        return(FALSE)

    return(TRUE)
}


find_non_consecutive_ones <- function(a_list, orig_name) {
    if (length(a_list) < 2) {
        return(NULL)
    }
    a_list <- sort(a_list)
    i1 <- a_list[[1]] # start
    i2 <- a_list[[length(a_list)]] # end

    ideal_string <- gsub(sprintf(
        '.*(%s.*%s).*', i1, i2), '\\1', paste0(LETTERS, collapse=''))
    ideal_list <- strsplit(ideal_string, '')[[1]]

    # Use the non-existing ones as delimiter to divide
    # -- all after should not appear
    del <- ideal_list[which(!(ideal_list %in% a_list))]
    if (length(del) > 1) {
        print('WARNING! NOT IMPLEMENTED!!!')
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
                new_name <- sprintf(
                    '%s_(%s)', orig_name,
                    paste0(remaining_list[1:(ii_i - 1)], collapse='+')
                )
                new_list[[new_name]] <- remaining_list[1:(ii_i - 1)]
            }
            remaining_list <- remaining_list[-c(1:ii_i)]
        } else {
            print('ERROR!')
        }
    }
    if (length(remaining_list) > 0) {
        new_name <- sprintf(
            '%s_(%s)', orig_name, paste0(remaining_list, collapse='+'))
        new_list[[new_name]] <- remaining_list
    }

    return(new_list)
}


list_to_data_frame <- function(clust_list, clone_dic=NULL) {
    # TODO: incorporate the clone_dic
    dat <- data.frame(
        cell_id=unlist(clust_list, use.names=FALSE),
        clone_id=NA,
        stringsAsFactors=FALSE
    )
    for (nl in names(clust_list)) {
        dat$clone_id[dat$cell_id %in% clust_list[[nl]]] <- nl
    }

    if (!is.null(clone_dic)) {
        dat <- dplyr::left_join(dat, clone_dic, by=c("clone_id"="old_K"))
        dat$clone_id <- dat$letters
        dat <- dat[, c('cell_id', 'clone_id')]
    }

    dat
}

rename_clones <- function(cell_clones) {
    cell_clones$clone_id <- as.integer(factor(cell_clones$clone_id))

    if(length(unique(cell_clones$clone_id)) <= 26) {
        cell_clones$clone_id <- LETTERS[cell_clones$clone_id]
    }

    return(cell_clones)
}



main <- function() {
    argv <- get_args()

    tree <- read.tree(argv$newick)
    g <- read_ltm_tree(tree_2_edge_list(tree))

    copynumber <- read.delim(
        argv$copynumber, check.names=FALSE, stringsAsFactors=FALSE)
    copynumber <- format_copynumber_matrix(copynumber)

    tree_split <- split_by_genotype_driver(
        g=g,
        minimum_fraction=argv$minimum_fraction,
        maximum_fraction=argv$maximum_fraction,
        mat=copynumber
    )

    cell_clones <- list_to_data_frame(tree_split$new_aug)
    cell_clones <- rename_clones(cell_clones)
    write.table(
        cell_clones, argv$output, quote=FALSE, sep="\t", row.names=FALSE)
}

#main()


