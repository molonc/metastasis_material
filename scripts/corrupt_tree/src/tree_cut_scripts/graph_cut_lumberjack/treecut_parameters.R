#!/usr/bin/env Rscript


## Installation
## See notes at run_treecut_parameters.R 
## How to run program, see main() function at run_treecut_parameters.R

suppressPackageStartupMessages({
    require("ape")
    require("gtools")
    require("igraph")
    # require("optparse")
    require("data.table")
})
# initial.options <- commandArgs(trailingOnly = FALSE)
# file.arg.name <- "--file="
# script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
# script.basename <- dirname(script.name)
# 
# source(paste0(script.basename, "/utils.R"))

split_by_genotype_driver <- function(the_graph,
                                     minimum_fraction,
                                     maximum_fraction,
                                     copy_number_matrix,
                                     clone_size_threshold,
                                     edge_diff_threshold,
                                     adjust_ploidy,
                                     prob,
                                     min_merge_back_fraction
                                     ) {

    # TODO remove... but removing changes results
    options(stringsAsFactors = FALSE)

    # 1. Split clones by their fractions
    res <- split_by_genotype(the_graph,
                             minimum_fraction,
                             maximum_fraction,
                             clone_size_threshold = clone_size_threshold,
                             edge_diff_threshold = edge_diff_threshold
                             )

    sub_subclone_list <- res$sub_subclone_list
    clone_dic_list <- res$clone_dic_list

    # 2. Split subclone by their ploidy
    if (adjust_ploidy == 1) {
        res <- split_by_ploidy(sub_subclone_list, clone_dic_list)
        sub_subclone_list <- res$sub_subclone_list
        clone_dic_list <- res$clone_dic_list
    }

    # 3. Compute distance between subclones in each clone
    res <- split_compute_intra_dist(sub_subclone_list = sub_subclone_list,
                                    clone_dic_list = clone_dic_list,
                                    copy_number_matrix = copy_number_matrix
                                    )

    internal_dist <- res$internal_dist
    sub_median_mat_list <- res$sub_median_mat_list
    prob <- 0.975
    # 4. Merge similar clones and drop too small ones
    new_aug <- split_merge_drop_clones(internal_dist = internal_dist,
                                       sub_subclone_list = sub_subclone_list,
                                       prob = prob,
                                       min_merge_back_fraction = min_merge_back_fraction
                                       )

    res <- list(new_aug = new_aug,
                internal_dist = internal_dist,
                clone_dic_list = clone_dic_list,
                sub_subclone_list = sub_subclone_list
                )

    return(res)
}

split_by_genotype <- function(the_graph,
                              minimum_fraction,
                              maximum_fraction,
                              clone_size_threshold,
                              edge_diff_threshold,
                              use_hamming=FALSE
                              ) {
    # get the subgraphs for each clade
    aug_cut <- get_aug_cut(the_graph = the_graph, minimum_fraction, maximum_fraction)
    subgraphs <- get_full_subgraphs_for_aug_cut(the_graph, aug_cut)

    cn <- as.list(names(subgraphs))
    names(cn) <- LETTERS[1:length(cn)]

    # check each clone for split
    internal_dist <- NULL
    sub_subclone_list <- NULL
    clone_dic_list <- NULL
    sub_median_mat_list <- NULL

    for (ii in names(subgraphs)) {
        gg <- subgraphs[[ii]]

        # Edges that are at least 20% of the clade
        ts <- get_edge_data(the_graph = gg)
        candid_edges <- ts[ts$frac > clone_size_threshold, ]

        # Further filter edges to ignore too close parental edges
        # ss is ordered so that closests to the root is first
        ss <- get_desc_names(gg, as.list(candid_edges$clone_id))

        ss_temp <- ss # Keep the updated fractions here
        total_cells <- g.count.cells(gg)
        n_test_edges <- nrow(candid_edges)
        delta_frac <- matrix(0, nrow = n_test_edges, ncol = n_test_edges)
        for (i in seq(n_test_edges)) {
            for (j in seq(n_test_edges)) {
                if (i < j) {
                    delta_frac[i, j] <- length(setdiff(ss[[i]], ss[[j]])) / total_cells
                }
            }
        }

        # Start from all the hanging-leaf-like loci
        test_candids <- candid_edges$clone_id

        # @Debug
        #annotate_nodes_on_tree(p = NULL, edge_list_path = NULL, node_names = test_candids, tree_node_dic = graph_to_tree_dic(gg))

        cc <- candid_edges[, c("clone_id", "frac", "height")]

        # While moving up the tree, remove the cells from the frac of a
        # parent that belong to the already seen ones
        bad_cadidates <- c()
        for (ii in rev(test_candids)) {
            p1 <- find_all_parents(gg, ii, max_height = max(candid_edges$height))
            if (length(p1) > 0) {
                p2 <- test_candids[test_candids %in% p1]
                # Find immediate parent
                immediate_parent <- rev(candid_edges$clone_id[candid_edges$clone_id %in% p2])[1]
                i_index <- which(candid_edges$clone_id == ii)
                j_index <- which(candid_edges$clone_id == immediate_parent)
                if (delta_frac[j_index, i_index] < edge_diff_threshold) {
                    bad_cadidates <- c(bad_cadidates, ii)
                } else {
                    # For all of j's children except i, update the difference
                    # by the negative of size of i
                    ch_in_list <- ss[[j_index]][ss[[j_index]] %in% test_candids[-c(i_index)]]
                    if (length(ch_in_list)) {
                        for (x in ch_in_list) {
                            x_index <- which(candid_edges$clone_id == x)
                            delta_frac[j_index, x_index] <- delta_frac[j_index, x_index] - length(ss[[i_index]])
                        }
                    }
                }
            }
        }

        passed_candids <- setdiff(test_candids, bad_cadidates)
        # Sort by height
        passed_candids <- rev(candid_edges$clone_id[candid_edges$clone_id %in% passed_candids])
        private_desc <- list()
        i_index <- 0
        for (ii in passed_candids) {
            i_index <- i_index + 1
            if (i_index > 1) {
                private_desc[[ii]] <- setdiff(ss[[ii]], unlist(ss[passed_candids[seq(i_index - 1)]]))
            } else {
                private_desc[[ii]] <- ss[[ii]]
            }
        }

        # Drop all loci
        for (ii in names(private_desc)) {
            private_desc[[ii]] <- private_desc[[ii]][!grepl("locus", private_desc[[ii]])]
            if (length(private_desc[[ii]]) == 0) {
                private_desc[[ii]] <- NULL
            }
        }

        # Add subclone letters ...
        clone_dic <- get_clone_dic_for_seq(a_seq = names(private_desc))
        clone_dic_list[[names(cn)[cn == ii]]] <- clone_dic

        # @Debug
        #p <- fast_tree_aug(g = gg, aug_cut = private_desc, outdir = '~/projects/tree_cutting/outputs/debug/July11/', title = names(cn)[cn == ii])

        # Keep track of sub_subclones for merging, dividing later on
        tmp_private_desc <- private_desc
        names(tmp_private_desc) <- clone_dic$letters
        sub_subclone_list[[names(cn)[cn == ii]]] <- tmp_private_desc
    }

    list(sub_subclone_list = sub_subclone_list, clone_dic_list = clone_dic_list)
}

split_by_ploidy <- function(sub_subclone_list, clone_dic_list) {

    # no datatag to generate copy_number_matrix from load_new_cn_data() function?
    #mat <- load_new_cn_data(datatag)

    new_sub_subclone_list <- list()
    new_clone_dic_list <- list()
    for (cn in names(sub_subclone_list)) {
        private_desc <- sub_subclone_list[[cn]]
        clone_dic <- clone_dic_list[[cn]]
        passed_mg <- median_genotype_from_list(private_desc, copy_number_matrix)
        pa <- get_ploidy_for_mat(t(passed_mg))
        ploidies <- unique(pa$ploidy)

        if (length(ploidies) > 1) {
            for (pl in ploidies) {
                new_sub_subclone_list[[paste0(cn, "-", "p", pl)]] <- private_desc[which(pa$ploidy == pl)]
                new_clone_dic_list[[paste0(cn, "-", "p", pl)]] <- clone_dic[which(pa$ploidy == pl), ]
            }
        } else {
            new_sub_subclone_list[[cn]] <- private_desc
            new_clone_dic_list[[cn]] <- clone_dic
        }
    }

    list(sub_subclone_list = new_sub_subclone_list, clone_dic_list = new_clone_dic_list)
}

split_compute_intra_dist <- function(sub_subclone_list,
                                     clone_dic_list,
                                     copy_number_matrix,
                                     use_hamming = FALSE
                                     ) {
    sub_median_mat_list <- list()
    internal_dist <- list()
    # print(length(sub_subclone_list))
    # print(names(sub_subclone_list)[1:3])
    for (cn in names(sub_subclone_list)) {
        private_desc <- sub_subclone_list[[cn]]
        # print(paste0("Debug: private ",length(private_desc)))
        clone_dic <- clone_dic_list[[cn]]
        passed_mg <- median_genotype_from_list(clone_list = private_desc, copy_number_matrix)
        sub_median_mat_list[[cn]] <- passed_mg
        # print(paste0("Debug: pass: ",is.null(passed_mg)))
        # print(paste0("Debug: mtx: ",dim(passed_mg)))
        # print(paste0("Debug: n_clones: ",nrow(passed_mg)))
        
        # Compute intra distance matrix
        
        passed_mg_dist <- compute_dist_mat(mg_mat = passed_mg, use_hamming = use_hamming)
        # print(paste0("Debug: passed mg dist: ",dim(passed_mg_dist)))
        # print(paste0("Debug: passed_mg_dist ",sum(passed_mg_dist)))
        # print(head(passed_mg_dist))
        # Keep track of distance between subclones
        # & !is.na(sum(passed_mg_dist))
        if (length(private_desc) > 1 & sum(passed_mg_dist) > 0) {
            tmp_indx <- which(passed_mg_dist > 0, arr.ind = TRUE)

            tmp_internal_dist <- data.frame(src_clone = clone_dic$letters[tmp_indx[, 1]],
                                            trg_clone = clone_dic$letters[tmp_indx[, 2]],
                                            sub_dist = passed_mg_dist[tmp_indx],
                                            original_clone = cn
                                            )
        } else {
            tmp_internal_dist <- data.frame(src_clone = names(private_desc),
                                            trg_clone = names(private_desc),
                                            sub_dist = NA,
                                            original_clone = cn
                                            )
        }

        if (is.null(internal_dist)) {
            internal_dist <- tmp_internal_dist
        } else {
            internal_dist <- rbind(internal_dist, tmp_internal_dist)
        }
    }

    list(internal_dist = internal_dist,
         sub_median_mat_list = sub_median_mat_list,
         sub_subclone_list = sub_subclone_list,
         clone_dic_list = clone_dic_list
         )
}

split_merge_drop_clones <- function(internal_dist,
                                    sub_subclone_list,
                                    prob=0.75,
                                    min_clone_fraction=0.01,
                                    min_merge_back_fraction=0.05
                                    ) {
    # Check for duplicates
    stopifnot(!any(duplicated(unlist(unname(sub_subclone_list)))))

    # 1.a. Use t-test scores to merge subclades within each clade
    new_clades <- list()
    for (cn in unique(internal_dist$original_clone)) {
        sub_internal_dist <- internal_dist[internal_dist$original_clone == cn, ]
        if (nrow(sub_internal_dist) > 1) {
            tt <- is_outlier(x = sub_internal_dist$sub_dist, prob = prob)
        } else {
            tt <- is_outlier(x = internal_dist$sub_dist, prob = prob)
            if (cn %in% internal_dist$original_clone[tt]) {
                tt <- c(TRUE)
            } else {
                tt <- c(FALSE)
            }
        }
        
        sub_rows <- sub_internal_dist[tt, c(1, 2)]
        all_subclones <- unique(c(sub_internal_dist$src_clone, sub_internal_dist$trg_clone))

        if (nrow(sub_rows) > 0) {
            tmp_clones <- unique(c(sub_rows$src_clone, sub_rows$trg_clone))
            non_merge_subclones <- setdiff(all_subclones, tmp_clones)
            new_clades[[cn]] <- list(merged_clones = list(merge_clones(sub_rows)), indiv_clones = non_merge_subclones)
        } else {
            new_clades[[cn]] <- list(indiv_clones = all_subclones)
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
                if (!is_consecuitive(a_list = jj)) {
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
                            new_new_clades[[cn]]$merged_clones[[1]] <- append(new_new_clades[[cn]]$merged_clones[[1]], ff)
                        } else {
                            # print(sprintf('Adding to indiv %s', qq))
                            # print(the_new_split[[qq]])
                            new_new_clades[[cn]]$indiv_clones <- append(new_new_clades[[cn]]$indiv_clones, the_new_split[[qq]])
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
            new_name <- paste0(cn, "_", ii)
            new_aug[[new_name]] <- sub_subclone_list[[cn]][[ii]]
        }

        # Now add the merged ones
        for (ii in new_clades[[cn]]$merged_clones) {
            # There may be multiple merged subgroups
            for (jj in ii) {
                new_name <- paste0(cn, "_(", paste0(jj, collapse = "+"), ")")
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
        print("Dropping ")
        print(drop_clones)
        new_aug[drop_clones] <- NULL
    }

    # Check for duplicates
    stopifnot(!any(duplicated(unlist(unname(new_aug)))))

    new_new_aug <- new_aug
    # 3.
    # print("Merge the too small ones back to the main clone...")
    # # print("ERROR! NOT MERGING TINY CLONES...")
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

    str(new_new_aug)

    # Check for double assignment
    stopifnot(!any(duplicated(unlist(unname(new_new_aug)))))

    return(new_new_aug)
}

get_edge_data <- function(the_graph) {
    # Number of edges=Number of nodes - number of leaves
    # TODO: generalise this...

    # TODO: @Sohrab: remove this
    # out_ts <- sprintf('%s_%s_ts_dat.rds', './outputs/', edge_density(g))
    # if (file.exists(out_ts)) return(readRDS(out_ts))

    n_nodes <- length(V(the_graph)$name)
    n_leaves <- length(V(the_graph)$name[!grepl("root", V(the_graph)$name) & !grepl("locus", V(the_graph)$name)])
    n_edges <- n_nodes - n_leaves

    # internal_nodes
    i_nodes <- V(the_graph)$name[grepl("root", V(the_graph)$name) | grepl("locus", V(the_graph)$name)]

    # Account for GM cells
    # TODO REMOVE
    # GM_nodes <- grep("SA928", V(the_graph)$name, value = TRUE)

    # Sanity check
    # stopifnot(sum(tcc$Freq) == (n_leaves - length(GM_nodes)))

    # Find all descendents
    desc <- get_decendents(clone_roots = i_nodes, the_graph = the_graph, min_cell_per_clone = 1)
    # Find out which one is a cells and which is a locus (will be NA)
    res <- retrieve_time_point_per_clade(the_graph = the_graph, decendents = desc)

    # Remove all NAs entries in res
    for (cn in names(res)) {
        if (all(is.na(res[[cn]]))) {
            res[[cn]] <- NULL
        } else {
            res[[cn]] <- res[[cn]][!is.na(res[[cn]])]
        }
    }

    # Count the number of cells in each subgroup
    ts <- lapply(names(res), function(x) {length(res[[x]])})
    ts <- data.frame(clone_id = names(res), N = unlist(ts), stringsAsFactors = FALSE)
    stopifnot(ts$N[ts$clone_id == "root"] == n_leaves)

    ts$frac <- ts$N / n_leaves
    stopifnot(max(ts$frac) <= 1 & min(ts$frac) >= 0)

    # Add height that
    height_dat <- get_height_dat(the_graph)
    height_dat <- height_dat[height_dat$id %in% ts$clone_id, ]
    colnames(height_dat) <- c("clone_id", "height")

    ts <- dplyr::right_join(ts, height_dat, by = c("clone_id"))

    ts <- ts[order(ts$height), ]
    ts$id <- paste0("(height=", ts$height, ", frac=", format(ts$frac, digits=2), ")")

    # TODO:@Sohrab remove this
    # saveRDS(ts, out_ts)

    return(ts)
}

collapse_tree_by_edge <- function(the_graph,
                                  nodes,
                                  do_plot = FALSE,
                                  annotate_nodes = NULL
                                  ) {
    h <- get_height_dat(the_graph)
    h <- h[h$id %in% nodes, ]
    h <- h[order(h$dist, decreasing = TRUE), ]

    # The neighbour with the smallest height is the parent
    edge_list <- matrix(NA, nrow = nrow(h), ncol = 2, dimnames = list(c(), c("source", "target")))
    for (i in seq(nrow(h))) {
        p <- find_all_parents(the_graph, h$id[i], h$dist[i])
        parent <- "root"
        if (nrow(h[h$id %in% p, ]) > 0) {
            parent <- h$id[h$id %in% p][1]
        }

        edge_list[i, ] <- c(parent, h$id[i])
    }

    orig_edge_list <- edge_list

    # Convert cell name to clone number
    clone_g <- igraph::graph_from_edgelist(el = edge_list)
    gg <- igraph::graph_from_edgelist(el = orig_edge_list)
    V(gg)$id <- seq(vcount(gg))
    list(edge_list = orig_edge_list, graph = gg)
}

h1_auto_cut <- function(the_graph,
                        candiate_edges,
                        minimum_fraction,
                        maximum_fraction = 0.5
                        ) {

    # Sort edges by height
    candiate_edges <- candiate_edges[order(candiate_edges$height, decreasing = FALSE), ]
    tau <- collapse_tree_by_edge(the_graph = the_graph, nodes = candiate_edges$clone_id)

    qeueu <- get_dfs_queue(tau$graph, "root")
    # Remove the root
    qeueu <- qeueu[-c(1)]
    covered_nodes <- list()
    the_cut <- list()
    for (edge in qeueu) {
        if (edge %in% covered_nodes) {
            next
        }
        tmp <- candiate_edges[candiate_edges$clone_id == edge, ]

        # If passed
        filter <- (tmp$frac < maximum_fraction)

        if (filter) {
            the_cut <- append(the_cut, edge)
            # Remove all children from the queue
            covered_nodes <- append(covered_nodes, find_all_descendents(tau$graph, edge, max(candiate_edges$height)))
        } else {
            # If failed at a leaf, pick the first node that is a child of a
            # parent with 2 or more children
            if (length(find_all_descendents(tau$graph, edge, 1)) == 0) {
                t1 <- find_all_parents(tau$graph, edge, tmp$height)
                t2 <- candiate_edges[candiate_edges$clone_id %in% t1, ]
                degrees <- igraph::degree(graph = tau$graph, v = t1, mode = "out")
                degrees <- data.frame(clone_id = names(degrees), degree = as.vector(degrees))
                t2 <- dplyr::right_join(t2, degrees)
                t2 <- t2[order(t2$height, decreasing = TRUE), ]
                branch_parent_index <- which(t2$degree > 1)[1]
                if (branch_parent_index > 1) {
                    edge <- t2[branch_parent_index - 1, ]$clone_id
                    the_cut <- append(the_cut, edge)
                    # Remove all children from the queue
                    covered_nodes <- append(covered_nodes, find_all_descendents(tau$graph, edge, max(candiate_edges$height)))
                }
                # Otherwise ignore this part of the tree
            }
        }
    }
    return(the_cut)
}

get_the_cut <- function(the_graph, minimum_fraction, maximum_fraction) {
    ts <- get_edge_data(the_graph)

    candiate_edges <- ts[ts$frac > minimum_fraction,]

    the_cut <- h1_auto_cut(the_graph = the_graph,
                           candiate_edges = candiate_edges,
                           minimum_fraction = minimum_fraction,
                           maximum_fraction = maximum_fraction
                           )

    return(the_cut)
}

add_siblings_cut <- function(the_cut,
                             the_graph,
                             minimum_fraction,
                             maximum_fraction,
                             oscilation=0
                             ) {
    # Now use the manually curated clones. Drop the subclades that threaten
    # the purity of the clade
    ts <- get_edge_data(the_graph)
    # heights
    q <- ts %>%
            dplyr::filter(clone_id %in% unlist(the_cut)) %>%
            dplyr::select(clone_id, height) %>%
            dplyr::arrange(desc(height)) %>%
            as.data.frame()

    # Add the descendents
    new_cut <- get_desc_names(the_graph, the_cut)

    # The clones should not include anyone else or any children that are a
    # parent to another existing clone
    cell_in_cut <- unlist(unname(get_desc_names(the_graph = the_graph, the_node = the_cut)))
    cell_in_cut <- cell_in_cut[!is.na(cell_in_cut)]

    # TODO: Check and remove all immediate children that are a
    # parent of someone else
    for (i in seq(nrow(q))) {
        node <- q$clone_id[i]; node_height <- q$height[i]
        parent <- find_parent(the_graph, node = node)

        # Don't overright a parent that has already been processed
        # TODO: should we generlaly break here?
        if (parent %in% names(new_cut)) {
            next
        }

        p_desc <- get_non_overlapping_desc(the_graph = the_graph, parent = parent, node = node, cell_in_cut = cell_in_cut)
        if (length(p_desc) != 0) {
            new_cut[[parent]] <- p_desc
            cell_in_cut <- c(cell_in_cut, p_desc)
        }
    }

    # TODO: replace with a black list
    # Remove GM cells
    # for (nc in names(new_cut)) {
    #     new_cut[[nc]] <- new_cut[[nc]][!grepl("SA928", new_cut[[nc]])]
    # }

    # Remove loci
    for (nc in names(new_cut)) {
        new_cut[[nc]] <- new_cut[[nc]][!grepl("locus_", new_cut[[nc]])]
    }

    # Remove NAs
    for (nc in names(new_cut)) {
        new_cut[[nc]] <- new_cut[[nc]][!is.na(new_cut[[nc]])]
    }

    # Impose min and max fraction
    new_clades <- setdiff(names(new_cut), unlist(the_cut))
    bad_clades <- c()
    n_cells <- ts$N[ts$clone_id == "root"]
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

get_aug_cut <- function(the_graph, minimum_fraction, maximum_fraction) {
    the_cut <- get_the_cut(the_graph = the_graph, minimum_fraction = minimum_fraction, maximum_fraction = maximum_fraction)
    add_siblings_cut(the_cut = the_cut, the_graph = the_graph, minimum_fraction = minimum_fraction, maximum_fraction = maximum_fraction)
}

get_full_subgraphs_for_aug_cut <- function(the_graph, aug_cut) {
    # For each clone, find the edge'
    # Gather all loci below immediate children
    children <- list()
    for (ii in names(aug_cut)) {
        children[[ii]] <- find_immediate_children(the_graph = the_graph, node = ii)
        children[[ii]] <- children[[ii]][!is.na(children[[ii]])]
        children[[ii]] <- children[[ii]][!(children[[ii]] %in% names(aug_cut))]
    }

    # Purge level 2
    bad_children <- list()
    for (ii in names(children)) {
        just_loci <- children[[ii]][grepl("locus_", children[[ii]])]
        temp <- get_desc_names(the_graph = the_graph, the_node = just_loci)
        bad_children[[ii]] <- c()
        for (tmp in names(temp)) {
            if (any(names(aug_cut) %in% temp[[tmp]])) {
                bad_children[[ii]] <- c(bad_children[[ii]], tmp)
            }
        }
    }
    for (ii in names(bad_children)) {
        children[[ii]] <- children[[ii]][!(children[[ii]] %in% bad_children[[ii]])]
    }

    # All loci: children of immediate loci + immediate children that are cells
    desc_loci <- list()
    for (ii in names(children)) {
        loci_children <- grep("locus_", children[[ii]], value = TRUE)
        desc <- unname(unlist(get_desc_names(the_graph = the_graph, the_node = loci_children)))
        desc <- desc[!(desc %in% loci_children)]

        # Add all immediate children
        desc <- c(desc, children[[ii]])
        desc_loci[[ii]] <- desc
    }

    # Remove tip loci
    tip_loci <- get_leaves_names_from_graph(the_graph = the_graph, only_loci = TRUE)

    if (length(tip_loci) > 0) {
        for (ii in names(desc_loci)) {
            desc_loci[[ii]] <- desc_loci[[ii]][!(desc_loci[[ii]] %in% tip_loci)]
        }
    }

    # TODO fix this
    # Remove GM cells
    if (length(tip_loci) > 0) {
        for (ii in names(desc_loci)) {
            desc_loci[[ii]] <- desc_loci[[ii]][!(grepl("SA928", desc_loci[[ii]]))]
        }
    }

    # Create subgraphs
    subgraphs <- list()
    for (ii in names(desc_loci)) {
        vids <- c(desc_loci[[ii]], ii)
        subgraphs[[ii]] <- igraph::induced.subgraph(graph = the_graph, vids = vids)
        subgraphs[[ii]] <- setup_graph(subgraphs[[ii]])
    }

    return(subgraphs)
}

launch_run <- function(desc, tree, copy_number, 
                       grouping_file,save_dir, 
                       minimum_fraction,maximum_fraction,clone_size_threshold, 
                       edge_diff_threshold, adjust_ploidy, 
                       prob, min_merge_back_fraction) {
        
    print("Read data")
    
    the_graph <- read_ltm_tree(tree_2_edge_list(tree))
    copy_number_matrix <- format_copynumber_matrix(copy_number)

    print("Split tree")
    tree_split <- split_by_genotype_driver(the_graph = the_graph,
                                           minimum_fraction = minimum_fraction,
                                           maximum_fraction = maximum_fraction,
                                           copy_number_matrix = copy_number_matrix,
                                           clone_size_threshold = clone_size_threshold,
                                           edge_diff_threshold = edge_diff_threshold,
                                           adjust_ploidy = adjust_ploidy,
                                           prob = prob,
                                           min_merge_back_fraction = min_merge_back_fraction
                                           )

    cell_clones <- list_to_data_frame(tree_split$new_aug)
    cell_clones <- rename_clones(cell_clones)
    print("Write clones info to file")
    data.table::fwrite(cell_clones, paste0(save_dir, 'cell_clones_',desc,'.csv.gz'))
    
    print("Plot heatmap")
    output_hm <- paste0(save_dir, 'heatmap_',desc,'.png')
    # copynumber_df <- read_tsv(copynumber)
    # clones <- read_tsv(output)
    # pdf(output_hm, width=10)
    png(output_hm, height = 2*800, width=2*1400, res = 2*72)
    make_cell_copynumber_tree_heatmap(
        tree, copy_number, cell_clones, NULL, grouping_file
    )
    dev.off()
    print("DONE!!")
}

# option_list <- list(make_option(c("-t", "--newick"),
#                                 type = "character",
#                                 default = NULL,
#                                 help = "output from corrupt_tree",
#                                 metavar = "character"
#                                 ),
#                     make_option(c("-n", "--copynumber"),
#                                 type = "character",
#                                 default = NULL,
#                                 help = "filtered_cell_cn",
#                                 metavar = "character"
#                                 ),
#                     make_option(c("-o", "--output"),
#                                 type = "character",
#                                 default = NULL,
#                                 help = "tree_heatmap.png",
#                                 metavar = "character"
#                                 ),
#                     make_option(c("-m", "--minimum_fraction"),
#                                 type = "double",
#                                 default = 0.02,
#                                 help = "minimum fraction of cells in initial clone selection"
#                                 ),
#                     make_option(c("-x", "--maximum_fraction"),
#                                 type = "double",
#                                 default = 0.38,
#                                 help = "maximum fraction of cells in initial clone selection"
#                                 ),
#                     make_option(c("-c", "--clone_size_threshold"),
#                                 type = "double",
#                                 default = 0.13,
#                                 help = "threshold for clone size"
#                                 ),
#                     make_option(c("-e", "--edge_diff_threshold"),
#                                 type = "double",
#                                 default = 0.04,
#                                 help = "threshold for edge difference"
#                                 ),
#                     make_option(c("-a", "--adjust_ploidy"),
#                                 type = "integer",
#                                 default = 0,
#                                 help = "if 1, adjust ploidy, if 0, do not adjust ploidy"
#                                 ),
#                     make_option(c("-p", "--prob"),
#                                 type = "double",
#                                 default = 0.75,
#                                 help = "probability threshold to call clones"
#                                 ),
#                     make_option(c("-f", "--min_merge_back_fraction"),
#                                 type = "double",
#                                 default = 0.03,
#                                 help = "minimum fraction to merge back"
#                                 )
#                )
# 
# opt_parser <- OptionParser(option_list = option_list)
# opt <- parse_args(opt_parser)
# launch_run(opt$newick,
#     opt$copynumber,
#     opt$output,
#     opt$minimum_fraction,
#     opt$maximum_fraction,
#     opt$clone_size_threshold,
#     opt$edge_diff_threshold,
#     opt$adjust_ploidy,
#     opt$prob,
#     opt$min_merge_back_fraction)





