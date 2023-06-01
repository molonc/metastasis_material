# source("/datadrive/tree_viz/ss_viz/general_utils.R")

split_by_genotype_driver <- function(aug_path=NULL) {
    outpath <- sprintf('%s_split_results.rds', aug_path) 
    if (file.exists(outpath)) {
        return(readRDS(outpath)$new_aug)
    }
    
    use_genotype_cuts <<- FALSE

    # 1. Split clones by their fractions
    res <- split_by_genotype(clone_size_treshold = .15, edge_diff_treshold = .04)
    sub_subclone_list <- res$sub_subclone_list
    clone_dic_list <- res$clone_dic_list; str(clone_dic_list)
    str(clone_dic_list$A)

    # 2. Split subclone by their ploidy -- only for cell lines
    adjust_ploidy <- FALSE
    if (datatag %in% c('SA906a', 'SA906b', 'SA922')) adjust_ploidy <- TRUE
    
    if (adjust_ploidy) {
        res <- split_by_ploidy(sub_subclone_list, clone_dic_list)
        sub_subclone_list <- res$sub_subclone_list
        clone_dic_list <- res$clone_dic_list
    }

    # 3. Compute distance between subclones in each clone
    res <- split_compute_intra_dist(sub_subclone_list, clone_dic_list)
    internal_dist <- res$internal_dist
    sub_median_mat_list <- res$sub_median_mat_list
    
    # 4. Merge similar clones and drop too small ones
    new_aug <- split_merge_drop_clones(internal_dist = internal_dist, sub_subclone_list = sub_subclone_list, prob = .75, min_merge_back_fraction = .03); str(new_aug)
    res <- list(new_aug = new_aug, 
                            internal_dist = internal_dist, 
                            clone_dic_list = clone_dic_list, 
                            sub_subclone_list = sub_subclone_list)
    
    use_genotype_cuts <<- TRUE
    saveRDS(res, outpath)
    return(new_aug)
}

split_by_genotype <- function(clone_size_treshold = .15, edge_diff_treshold = .04) {
    # SA1000; clone_size_treshold <- .15; edge_diff_treshold <- .09
    use_hamming <- FALSE

    # get the subgraphs for each clade
    g <- read_ltm_tree(edge_list_path = edge_list_path)
    subgraphs <- get_full_subgraphs_for_aug_cut(g = g, datatag = datatag, edge_list_path = edge_list_path)
    str(subgraphs, max.level = 1)
    cn <- get_clone_letters(datatag, edge_list_path)
    
    # check each clone for split
    internal_dist <- NULL
    sub_subclone_list <- NULL
    clone_dic_list <- NULL
    sub_median_mat_list <- NULL
    for (ii in names(subgraphs)) {
        gg <- subgraphs[[ii]]
        
        # Edges that are at least 20% of the clade
        ts <- get_edge_data(g = gg, tcc = get_tcc(gg), datatag = NULL, edge_list_path = NULL)
        candid_edges <- ts[ts$frac > clone_size_treshold, ]
        
        # Further filter edges to ignore too close parental edges 
        ss <- get_desc_names(NULL, gg, as.list(candid_edges$clone_id)) # ss is ordered so that closests to the root is first
        
        ss_temp <- ss # Keep the updated fractions here
        total_cells <- g.count.cells(gg)
        n_test_edges <- nrow(candid_edges)
        delta_frac <- matrix(0, nrow = n_test_edges, ncol = n_test_edges)
        for (i in seq(n_test_edges)) {
            for (j in seq(n_test_edges)) {
                if (i < j) {
                    delta_frac[i, j] <- length(setdiff( ss[[i]], ss[[j]] ))/total_cells
                }
            }
        }
        
        # Start from all the hanging-leaf-like loci
        test_candids <- candid_edges$clone_id
        cc <- candid_edges[, c('clone_id', 'frac', 'height')]
        cc$pretty_names <- get_pretty_names_for_loci(cc$clone_id)
        
        # While moving up the tree, remove the cells from the frac of a parent that belong to the already seen ones
        bad_cadidates <- c()
        for (ii in rev(test_candids)) {
            # ii = rev(test_candids)[1]
            p1 <- find_all_parents(gg, ii, max_height = max(candid_edges$height))
            if (length(p1) > 0) {
                p2 <- test_candids[test_candids %in% p1]
                # Find immediate parent
                immediate_parent <- rev(candid_edges$clone_id[candid_edges$clone_id %in% p2])[1]
                i_index <- which(candid_edges$clone_id == ii)
                j_index <- which(candid_edges$clone_id == immediate_parent)
                if (delta_frac[j_index, i_index] < edge_diff_treshold) {
                    bad_cadidates <- c(bad_cadidates, ii)
                } else {
                    # For all of j's children except i, update the difference by the negative of size of i
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
        passed_candids <- rev(candid_edges$clone_id[candid_edges$clone_id %in% passed_candids]) # Sort by height
        private_desc <- list()
        i_index <- 0
        for (ii in passed_candids) {
            i_index <- i_index + 1
            if (i_index > 1) {
                private_desc[[ii]] <- setdiff(ss[[ii]], unlist( ss[passed_candids[seq(i_index-1)]] )  )
            } else {
                private_desc[[ii]] <- ss[[ii]]
            }
        }
        
        # Drop all loci
        for (ii in names(private_desc)) {
            private_desc[[ii]] <- private_desc[[ii]][!grepl('locus', private_desc[[ii]])]
        }
        
        # Add subclone letters ...
        clone_dic <- get_clone_dic_for_seq(a_seq = names(private_desc))
        clone_dic_list[[names(cn)[cn == ii]]] <- clone_dic
        
        # Keep track of sub_subclones for merging, dividing later on 
        tmp_private_desc <- private_desc
        names(tmp_private_desc) <- clone_dic$letters
        sub_subclone_list[[names(cn)[cn == ii]]] <- tmp_private_desc
    }
    
    list(sub_subclone_list = sub_subclone_list, clone_dic_list = clone_dic_list)
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
                new_sub_subclone_list[[paste0(cn, '-', 'p', pl)]] <- private_desc[which(pa$ploidy == pl)]
                new_clone_dic_list[[paste0(cn, '-', 'p', pl)]] <- clone_dic[which(pa$ploidy == pl), ]
            }
        } else {
            new_sub_subclone_list[[cn]] <- private_desc
            new_clone_dic_list[[cn]] <- clone_dic
        }
    }
    
    list(sub_subclone_list = new_sub_subclone_list, clone_dic_list = new_clone_dic_list)
}

split_compute_intra_dist <- function(sub_subclone_list, clone_dic_list) {
    mat <- load_new_cn_data(datatag)
    sub_median_mat_list <- list()
    internal_dist <- list()
    for (cn in names(sub_subclone_list)) {
        private_desc <- sub_subclone_list[[cn]]
        clone_dic <- clone_dic_list[[cn]]
        passed_mg <- median_genotype_from_list(private_desc, mat)
        sub_median_mat_list[[cn]] <- passed_mg  
        
        # Compute intra distance matrix
        passed_mg_dist <- compute_dist_mat(mg_mat = passed_mg, use_hamming = use_hamming)
            
        # Keep track of distance between subclones
        if (length(private_desc) > 1 & sum(passed_mg_dist) > 0) {
            tmp_indx <- which(passed_mg_dist > 0, arr.ind = T)
            tmp_internal_dist <- data.frame(src_clone = clone_dic$letters[tmp_indx[, 1]],
                                                                            trg_clone = clone_dic$letters[tmp_indx[, 2]],
                                                                            sub_dist = passed_mg_dist[tmp_indx], 
                                                                            original_clone = cn)
        } else {
            tmp_internal_dist <- data.frame(src_clone = names(private_desc),
                                                                            trg_clone = names(private_desc),
                                                                            sub_dist = NA, 
                                                                            original_clone = cn)
        }
        if (is.null(internal_dist)) {
            internal_dist <- tmp_internal_dist
        } else {
            internal_dist <- rbind(internal_dist, tmp_internal_dist)
        }
    }
    list(internal_dist = internal_dist, sub_median_mat_list = sub_median_mat_list, sub_subclone_list = sub_subclone_list, clone_dic_list = clone_dic_list)
}

split_merge_drop_clones <- function(internal_dist, sub_subclone_list, prob = .75, min_clone_fraction = .01, min_merge_back_fraction = .05) {
    options(stringsAsFactors = F)
    
    # Check for duplicates
    stopifnot(!any(duplicated(unlist(unname(sub_subclone_list)))))

    # 1.a. Use t-test scores to merge subclades within each clade
    new_clades <- list()
    for (cn in unique(internal_dist$original_clone)) {
        print(cn)
        sub_internal_dist <- internal_dist[internal_dist$original_clone == cn, ]
        if (nrow(sub_internal_dist) > 1) {
            tt <- is_outlier(x=sub_internal_dist$sub_dist, prob=prob)
        } else {
            tt <- is_outlier(x=internal_dist$sub_dist, prob=prob)
            if (cn %in% internal_dist$original_clone[tt]) {
                tt = c(TRUE)
            } else {
                tt = c(FALSE)
            }
        }
        
        sub_rows <- sub_internal_dist[tt, c(1,2)]
        all_subclones <- unique(c(sub_internal_dist$src_clone, sub_internal_dist$trg_clone))
        
        if (nrow(sub_rows) > 0) {
            tmp_clones <- unique(c(sub_rows$src_clone, sub_rows$trg_clone))
            non_merge_subclones <- setdiff(all_subclones, tmp_clones)
            new_clades[[cn]] <- list(merged_clones = list(merge_clones(sub_rows)), indiv_clones = non_merge_subclones)
        } else {
            new_clades[[cn]] <- list(indiv_clones = all_subclones)
        }
    }
    str(new_clades)
    
    new_new_clades <- new_clades
    # 1.a.0 Do not merge subclones that envelope a third one
    for (cn in names(new_clades)) {
        print(cn)
        for (ii in new_clades[[cn]]$merged_clones) {
            print(ii)
            i_index <- 0
            for (jj in ii) {
                i_index <- i_index + 1
                print(jj)
                if (!is_consecuitive(a_list = jj)) {
                    # Remove the offending ones
                    print(jj)
                    the_new_split <- find_non_consecutive_ones(jj, cn)
                    print('the_new_split')
                    print(the_new_split)
                    # remove this group from merged groups
                    new_new_clades[[cn]]$merged_clones[[1]][[i_index]] <- NA
                    # for each one in the new group, either add them as a new merged or a new nonmerged
                    for (qq in names(the_new_split)) {
                        if (nchar(qq) > 5) { # e.g., B_(E)
                            ff <- unname(the_new_split[qq])
                            new_new_clades[[cn]]$merged_clones[[1]] <- append(new_new_clades[[cn]]$merged_clones[[1]], ff)
                        } else {
                            print(sprintf('Adding to indiv %s', qq))
                            print(the_new_split[[qq]])
                            new_new_clades[[cn]]$indiv_clones <- append(new_new_clades[[cn]]$indiv_clones, the_new_split[[qq]])
                        }
                    }
                    print(sprintf('----------------------- begin{State at the end of %s} --------------------', cn))
                    print(new_new_clades[[cn]])
                    print(sprintf('----------------------- end{State at the end of %s} --------------------', cn))
                }
            }
        }
        print(str(new_new_clades[[cn]]))
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
                new_name <- paste0(cn, '_(', paste0(jj, collapse = '+'), ')')
                for (kk in jj) {
                    new_aug[[new_name]] <- c(new_aug[[new_name]], sub_subclone_list[[cn]][[kk]])
                }
            }
        }
    }
    str(new_aug)
    
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
    
    # 3. Merge the too small ones back to the main clone
    tiny_clones <- names(lens)[lens < min_merge_back_fraction]
    new_new_aug <- new_aug
    for (tc in tiny_clones) {
        tokens <- strsplit(tc, '_')[[1]]
        orig_clone_name <- lens[grep(paste0(tokens[[1]], '_'), names(lens))]
        orig_clone_name <- names(orig_clone_name)[which.max(orig_clone_name)]
        if (orig_clone_name != tc) {
            new_new_aug[[orig_clone_name]] <- c(new_new_aug[[orig_clone_name]], new_aug[[tc]])
            new_new_aug[[tc]] <- NULL
        }
    }
    str(new_new_aug)
    # Check for double assignment
    stopifnot(!any(duplicated(unlist(unname(new_new_aug)))))
    new_new_aug
}

get_edge_data <- function(g, tcc, datatag=NULL, edge_list_path=NULL) {    
    if (!is.null(edge_list_path)) {
        outpath <- sprintf('%s_ts_%s.rds', edge_list_path, datatag) 
        if (file.exists(outpath))
            return(readRDS(outpath))
    }
    # Number of edges = Number of nodes - number of leaves
    n_nodes = length(V(g)$name)
    n_leaves = length(V(g)$name[!grepl('root', V(g)$name) & !grepl('locus', V(g)$name)])
    n_edges = n_nodes - n_leaves
    # internal_nodes 
    i_nodes = V(g)$name[grepl('root', V(g)$name) | grepl('locus', V(g)$name)]
    # Account for GM cells
    GM_nodes = grep('SA928', V(g)$name, value=TRUE)
    # Below the threshold of detection
    zero_threshold = .01
    # Sanity check
    stopifnot(sum(tcc$Freq) == (n_leaves - length(GM_nodes)))
    
    # For each edge, 
    # 1. Number of cells (N)
    # 2. Table of trajectory (X)
    # 3. Oscillation
    # 4. Dip to zero?
    # 5. Dip down and up?
    
    # Use functions in cladedivide
    res <- NULL
    
    # A bunch of node ids
    desc <- get_decendents(i_nodes, the_graph = g, min_cell_per_clone = 1)
    res <- retrieve_time_point_per_clade(the_graph = g, decendents = desc, data_tag = datatag)
    print('WARNING - not removing all NA...')
    ts <- get_time_series_for_cluster_cut(res = res, total_cell_counts_break_down = tcc)
    
    # Sort the timepoints (used in oscilation compuation and down and up)
    ts <- order_cols(ts)
    ts_counts <- get_time_series_for_cluster_cut(res, tcc, FALSE)
    ts_counts <- order_cols(ts_counts)
    x_cols <- grep('X[0-9][0-9]*', colnames(ts_counts))
    colnames(ts_counts)[x_cols] <- gsub('X', 'N', colnames(ts_counts)[x_cols])
    ts <- dplyr::right_join(ts, ts_counts, by = c('clone_id'))
    ts$N = unlist(lapply(names(res), function(x) (length(res[[x]]))))
    ts$frac <- ts$N / sum(tcc$Freq)
    
    if (!(max(ts$frac, na.rm = T) <= 1 & min(ts$frac, na.rm = T) >= 0)) {
        print('Warning - fraction boundaries violated -- check for NAs in graph to timepoint res')
    }
    if (is.null(datatag)) {
        ts$oscillation <- 0.0
    } else if (datatag == 'SA666') {
        ts$oscillation <- 0.0
    } else {
        ts$oscillation <- oscillation_for_ts(ts[, x_cols])
    }
    
    ts$dip_to_zero_and_rise = FALSE
    ts$dip_to_zero_and_rise_soft = FALSE
    ts$down_and_up_for_ts_hard = TRUE
    ts$down_and_up_for_ts_soft = TRUE
    
    # Add height that
    height_dat <- get_height_dat(g)
    height_dat <- height_dat[height_dat$id %in% ts$clone_id, ]
    colnames(height_dat) <- c('clone_id', 'height')
    ts <- dplyr::right_join(ts, height_dat, by=c('clone_id'))
    ts <- ts[order(ts$height), ]
    ts$id = paste0('(height=', ts$height, ', frac=', format(ts$frac, digits = 2), ')')
    if (!is.null(edge_list_path)) {
        saveRDS(ts, outpath)
    } 
    ts
}

get_the_cut <- function(datatag, edge_list_path, g = NULL) {
    cut_configs <- load_cut_configs(datatag = datatag)
    outpath <- get_cut_outpath(datatag = datatag, edge_list_path = edge_list_path, cut_configs = cut_configs)
    print(outpath)
    if (file.exists(outpath)) {
        return(readRDS(outpath))
    }

    if (is.null(g)) {
        g <- read_ltm_tree(edge_list_path)
    }
    tcc <- get_tcc(g)
    ts <- get_edge_data(g, tcc, datatag, edge_list_path)
    minimum_fraction <- cut_configs$minimum_fraction
    maximum_frac <- cut_configs$maximum_frac
    oscilation <- cut_configs$oscilation

    # Set oscilation after min has been set
    if (oscilation != 0) {
        oscilation <- summary(ts$oscillation[ts$frac > minimum_fraction])[oscilation]
    } 
    
    # TODO: filter by those trajectories that have hit a minimum at last at some timepoin
    candiate_edges <- ts[ts$frac > minimum_fraction & ts$oscillation > oscilation, ]
    the_cut <- h1_auto_cut(the_g = g, candiate_edges = candiate_edges,
                                                minimum_fraction = minimum_fraction, maximum_frac = maximum_frac, 
                                                use_soft_down_up = T, 
                                                ignore_down_up = T, ignore_dip_to_zero = T)

    saveRDS(the_cut, outpath)
    return(the_cut)
}

h1_auto_cut <- function(the_g, candiate_edges, minimum_fraction, maximum_frac = .5, use_soft_down_up=FALSE, ignore_down_up=TRUE, ignore_dip_to_zero=FALSE) {
    # Sort edges by height
    candiate_edges <- candiate_edges[order(candiate_edges$height, decreasing = F),]
    
    candiate_edges$down_and_up <- candiate_edges$down_and_up_for_ts_hard
    if (use_soft_down_up) {
        candiate_edges$down_and_up = candiate_edges$down_and_up_for_ts_soft
    }
    
    tau <- collapse_tree_by_edge(the_g = the_g, nodes = candiate_edges$clone_id)
    
    qeueu <- get_dfs_queue(tau$graph, 'root')
    # Remove the root
    qeueu <- qeueu[-c(1)]
    covered_nodes <- list()
    the_cut <- list()
    for (edge in qeueu) {
        if(edge %in% covered_nodes) next
        tmp <- candiate_edges[candiate_edges$clone_id == edge, ]
        if (ignore_dip_to_zero)
            tmp$dip_to_zero_and_rise = FALSE
        filter = (tmp$dip_to_zero_and_rise == FALSE & tmp$down_and_up == FALSE & tmp$frac < maximum_frac)
        if (ignore_down_up)
            filter = (tmp$dip_to_zero_and_rise == FALSE & tmp$frac < maximum_frac)
        if (filter) {
            the_cut <- append(the_cut, edge)
            # Remove all children from the queue
            covered_nodes <- append(covered_nodes, find_all_descendents(tau$graph, edge, max(candiate_edges$height)))
        } else {
            # If failed at a leaf, pick the first node that is a child of a parent with 2 or more children
            if (length(find_all_descendents(tau$graph, edge, 1)) == 0) {
                t1 = find_all_parents(tau$graph, edge, tmp$height)
                t2 = candiate_edges[candiate_edges$clone_id %in% t1, ]
                degrees = igraph::degree(graph = tau$graph, v = t1, mode = 'out')
                degrees = data.frame(clone_id=names(degrees), degree=as.vector(degrees))
                t2 <- dplyr::right_join(t2, degrees)
                t2 <- t2[order(t2$height, decreasing = T), ]
                branch_parent_index = which(t2$degree > 1)[1]
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
    the_cut
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

add_siblings_cut <- function(datatag, edge_list_path, the_cut, the_g) {
    cut_path <- get_cut_outpath(datatag = datatag, edge_list_path = edge_list_path)
    outpath <- sprintf('%s_aug_cut.rds', cut_path) 
    
    # Now use the manually curated clones. Drop the subclades that threaten the purity of the clade
    if (file.exists(outpath)) {
        if (use_genotype_cuts == FALSE) {
            aug_cut <- readRDS(outpath)
            return(aug_cut)
        }
        return(split_by_genotype_driver(aug_path=outpath))
    }

    tcc <- get_tcc(the_g)
    ts <- get_edge_data(the_g, tcc, datatag, edge_list_path)
    cut_configs <- load_cut_configs(datatag = datatag)
    minimum_fraction <- cut_configs$minimum_fraction
    maximum_frac <- cut_configs$maximum_frac
    oscilation <- cut_configs$oscilation
    if (oscilation != 0) {
        oscilation <- summary(ts$oscillation[ts$frac > minimum_fraction])[oscilation]
    } 
    
    # heights
    q <- ts %>% dplyr::filter(clone_id %in% unlist(the_cut)) %>% dplyr::select(clone_id, height) %>% dplyr::arrange(desc(height)) %>% as.data.frame()
    
    # Add the descendents
    new_cut <- get_desc_names(datatag, the_g, the_cut)

    # The clones should not include anyone else or any children that are a parent to another existing clone
    cell_in_cut <- unlist(unname(get_desc_names(the_graph = the_g, datatag = datatag, the_node = the_cut)))
    cell_in_cut <- cell_in_cut[!is.na(cell_in_cut)]
    
    # TODO: Check and remove all immediate children that are a parent of someone else
    for (i in seq(nrow(q))) {
        node <- q$clone_id[i]; node_height <- q$height[i]
        parent <- find_parent(the_g, node = node)
        
        # Don't overright a parent that has already been processed
        # TODO: should we generlaly break here?
        if (parent %in% names(new_cut)) next
    
        p_desc <- get_non_overlapping_desc(the_graph = the_g, parent = parent, node = node, cell_in_cut = cell_in_cut)
        if (length(p_desc) != 0) {
            new_cut[[parent]] <- p_desc
            cell_in_cut <- c(cell_in_cut, p_desc)
        }
    }
    
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

    # Impose min and max fraction (no oscilation for now...)
    new_clades <- setdiff(names(new_cut), unlist(the_cut))
    bad_clades <- c()
    n_cells <- sum(tcc$Freq)
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
    
    # TODO: add oscilation...
    res <- retrieve_time_point_per_clade(NULL, new_cut, datatag)
    tser <- order_cols(get_time_series_for_cluster_cut(res))
    os <- oscillation_for_ts(as.matrix(tser[, -c(1)]))
    names(os) <- tser$clone_id
    bad_clades <- names(os)[os < oscilation]
    if (length(bad_clades) > 0) {
        for (nc in bad_clades) {
            new_cut[nc] <- NULL
        }
    }

    saveRDS(new_cut, outpath)
    
    if (use_split_cuts == FALSE)
        return(new_cut)
    
    return(drop_impure_splits(aug_cut=new_cut, aug_path=outpath))
}

get_full_subgraphs_for_aug_cut <- function (g = NULL, datatag=NULL, edge_list_path=NULL) {
    # For each clone, find the edge'
    # Gather all loci below immediate children
    aug_cut <- get_aug_cut(g = g, datatag = datatag, edge_list_path = edge_list_path)
    children <- list()
    for (ii in names(aug_cut)) {
        children[[ii]] <- find_immediate_children(the_g = g, node = ii)
        children[[ii]]<- children[[ii]][!is.na(children[[ii]])]
        children[[ii]] <- children[[ii]][!(children[[ii]] %in% names(aug_cut))]
    }
    
    # Purge level 2
    bad_children <- list()
    for (ii in names(children)) {
        just_loci <- children[[ii]][grepl('locus_', children[[ii]])]
        temp <- get_desc_names(datatag = datatag, the_graph = g, the_node = just_loci)
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
        loci_children <- grep('locus_', children[[ii]], value = T)
        desc <- unname(unlist(get_desc_names(datatag = datatag, the_graph = g, the_node = loci_children)))
        desc <- desc[!(desc %in% loci_children)]
        
        # Add all immediate children
        desc <- c(desc, children[[ii]])
        desc_loci[[ii]] <- desc
    }
    
    # Remove tip loci
    tip_loci <- get_leaves_names_from_graph(the_g = g, only_loci=TRUE)
    
    if (length(tip_loci) > 0) {
        for (ii in names(desc_loci)) {
            desc_loci[[ii]] <- desc_loci[[ii]][!(desc_loci[[ii]] %in% tip_loci)]
        }
    }
    
    # Remove GM cells
    if (length(tip_loci) > 0) {
        for (ii in names(desc_loci)) {
            desc_loci[[ii]] <- desc_loci[[ii]][!(grepl('SA928', desc_loci[[ii]]))]
        }
    }
    
    # Create subgraphs
    subgraphs <- list()
    for (ii in names(desc_loci)) {
        vids <- c(desc_loci[[ii]], ii)
        subgraphs[[ii]] <- igraph::induced.subgraph(graph = g, vids = vids)
        subgraphs[[ii]] <- setup_graph(g = subgraphs[[ii]])
    }
    subgraphs
}

get_clone_letters <- function(datatag, edge_list_path) {
    clone_dic <- cosmic_clone_dic(datatag, edge_list_path)
    cn <- as.list(clone_dic$old_K)
    names(cn) <- clone_dic$letters
    cn
}

drop_impure_splits <- function(aug_cut, aug_path) {
    outpath <- sprintf('%s_drop_split.rds', aug_path) 

    if (file.exists(outpath)) {
        return(readRDS(outpath))
    }
    
    use_split_cuts <<- FALSE
    
    # For the sub-clade original letter
    dummy_clone_dic <- cosmic_clone_dic(datatag, edge_list_path)
    res <- get_split_cut(datatag, edge_list_path, .3)
    splitted_clades <- res$splitted_clades
    
    print(length(unname(unlist(unlist(splitted_clades)))))
    
    # For the sub-subclades to be removed
    clone_dic <- res$clone_dic
    config <- load_split_cut_configs(datatag = datatag)
    clones_to_drop <- config$drop
    tmp <- list()
    for (ii in clones_to_drop) {
        orig_clone_letter <- gsub('[0-9]+', '', ii)
        orig_old_K <- dummy_clone_dic$old_K[dummy_clone_dic$letters == orig_clone_letter]
        old_K <- clone_dic$old_K[clone_dic$letters == ii]
        tmp[[orig_old_K]][[old_K]] <- splitted_clades[[old_K]]
        splitted_clades[[orig_old_K]][[old_K]] <- NULL
    }
    
    str(splitted_clades)
    
    # sub-subclones to combine
    clones_to_combine <- config$combine
    for (ii in names(clones_to_combine)) {
        tmp <- list()
        for (jj in clones_to_combine[[ii]]) {
            orig_clone_letter <- gsub('[0-9]+', '', jj)
            orig_old_K <- dummy_clone_dic$old_K[dummy_clone_dic$letters == orig_clone_letter]
            old_K <- clone_dic$old_K[clone_dic$letters == jj]
            tmp[[jj]] <- splitted_clades[[orig_old_K]][[old_K]]
            splitted_clades[[orig_old_K]][[old_K]] <- NULL
        }
        splitted_clades[[ii]] <- unlist(tmp)
    }

    # Merge the sub-sub-clades back
    for (ii in names(splitted_clades)) {
        splitted_clades[[ii]] <- unname(unlist(splitted_clades[[ii]]))
    }
    
    str(splitted_clades)
    saveRDS(splitted_clades, outpath)
    use_split_cuts <<- TRUE
    splitted_clades
}

cosmic_clone_dic <- function(datatag, edge_list_path, set_colours=FALSE) {
    g <- read_ltm_tree(edge_list_path = edge_list_path)
    tcc <- get_tcc(g)
    ts <- get_edge_data(g, tcc, datatag, edge_list_path)
    the_cut <- get_the_cut(datatag, edge_list_path, g=g)
    aug_cut <- add_siblings_cut(datatag = datatag, edge_list_path = edge_list_path, the_cut = the_cut, the_g = g)
    new_clades <- setdiff(names(aug_cut), unlist(the_cut))
    dummy_clone_dic <- universal_clone_dic(ts = ts, aug_cut = aug_cut, new_clades = new_clades, tcc = tcc)
    
    if (set_colours) {
        # Set the colours in ... with the legacy colours  
        candiate_edges <- get_candi_edge(the_cut = the_cut, aug_cut = aug_cut, ts = ts, tcc = tcc, new_clades)
        myColors <- get_cluster_colours(length(unique(candiate_edges$clone_id)))
        myColors <- myColors[1:length(unique(candiate_edges$clone_id))]
        names(myColors) <- unique(as.character(candiate_edges$clone_id))[order(get_pretty_names_for_loci(unique(candiate_edges$clone_id)))]
        cdf <- data.frame(old_K = names(myColors), color = unname(myColors))
        dummy_clone_dic <- dplyr::right_join(dummy_clone_dic, cdf) 
    }
    dummy_clone_dic
}

universal_clone_dic <- function(ts, aug_cut, new_clades, tcc) {
    config_path <- file.path(exp_path, 'config.yaml')
    if (!is.null(exp_path) & file.exists(config_path)) {
        fitness_exp_path <- file.path(dirname(exp_path), fitclone_exp_dir)
        clone_dic <- get_clone_dic(datatag, fitness_exp_path)
        clone_dic$colour <- get_cluster_colours(nrow(clone_dic))
        clone_dic$K <- as.character(clone_dic$K)
        dummy_clone_dic <- clone_dic
    } else {
        dummy_clone_dic <- sugar_dummy_clone_dic(ts = ts, aug_cut = aug_cut, new_clades = new_clades, tcc = tcc)
    }
    dummy_clone_dic
}

get_aug_cut <- function(g = NULL, datatag, edge_list_path) {
    if (is.null(g)) {
        g <- read_ltm_tree(edge_list_path)
    }
    the_cut <- get_the_cut(datatag, edge_list_path, g=g)
    add_siblings_cut(datatag = datatag, edge_list_path = edge_list_path, the_cut = the_cut, the_g = g)
}

sugar_dummy_clone_dic <- function(ts, aug_cut, new_clades, tcc) {
    candi_edge <- get_candi_edge(the_cut = NULL, aug_cut = aug_cut, ts = ts, tcc = tcc, new_clades = new_clades)
    get_dummy_clone_dic(candi_edge)
}

get_clone_dic <- function(datatag, exp_path, edge_list_path=NULL, use_cashe = TRUE) {

    if (datatag == 'SA1000') {
        return(get_clone_dic('SA609', exp_path = get_fitness_exp_path_for_datatag('SA609'), edge_list_path = get_edge_list_path_for_datatag('SA609')))
    } 
    
    if (datatag == 'SA922') {
        return(cosmic_clone_dic(datatag = datatag, edge_list_path = edge_list_path, set_colours = T))
    }
    
    if (datatag %in% get_treatment_datanames()) {
        outpath <- sprintf('%s/%s_clone_dic.rds', exp_path, 'dummy') 
    } else if (edge_list_path != '') {
            cut_path <- get_cut_outpath(datatag = datatag, edge_list_path = edge_list_path)
            aug_path <- sprintf('%s_aug_cut.rds', cut_path) 
            split_gen_path <- sprintf('%s_split_results.rds', aug_path) 
            outpath <- sprintf('%s_clone_dic.rds', split_gen_path) 
    } else {
        outpath <- sprintf('%s/%s_clone_dic.rds', get_exp_path_for_datatag(datatag), 'dummy') 
    }
    
    if (file.exists(outpath) & use_cashe) {
        return(readRDS(outpath))
    }
    
    config <- read_yaml(file.path(exp_path, 'config.yaml'))
    if (!is.null(config$K_prime)) {
        K = config$K_prime + 1
    } else {
        K = config$K + 1
    }
        
    # Find the original data file
    original_dat <- read_yaml(file.path(exp_path, 'config.yaml'))$original_data
    ss <- read.table(paste0(original_dat, '.gz'), header = T)
    
    # Find original-observed times (remove the last two)
    observed_times <- get_obsered_times(exp_path)
    ss <- ss[ss$time %in% observed_times, ]
    
    # Match by X to find which clone is which (the original name to help with colouring and naming, etc)
    # Use the version here: /Users/sohrabsalehi/projects/fitclone/figures/raw/supp/raw_times/{datatag}_dlp_time_dropped_no_time_autocut_{}_clones.tsv.gz
    raw_dat_path <- file.path('~/projects/fitclone/figures/raw/supp/raw_times', paste0(tolower(datatag), '_dlp_time_dropped_no_time_autocut_', K, '_clones.tsv.gz'))
    if (!file.exists(raw_dat_path)) {
        # Did we add a month tag?
        tokens <- strsplit(config$original_data, '/')[[1]]
        tag <- tokens[length(tokens) - 2]
        raw_dat_path <- file.path('~/projects/fitclone/figures/raw/supp/raw_times', paste0(tolower(datatag), '_', tag, '_dlp_time_dropped_no_time_autocut_', K, '_clones.tsv.gz'))
    }
    raw_dat <- read.table(raw_dat_path, stringsAsFactors = F)
        
    # Restric to one time point, hoping there aren't identical (mostly zero) clones to throw this off

    # A better way, look at the reference directory number
    # But just to make sure, use O(N^2) timeseries comparison
    # Find which of the original coordinates was removed
    
    found_orig_Ks <- matrix(NA, nrow=length(unique(raw_dat$K)), ncol = 1)
    for (j1 in seq(unique(raw_dat$K))) {
        for (i1 in seq(unique(ss$K))) {
            # Get x in observation sorted by time
            t1 <- ss[order(ss$time), ]
            t1 <- t1$X[t1$K == unique(ss$K)[i1]]
            
            # Get x in original_dat sorted by time
            t2 <- raw_dat[order(raw_dat$time), ]
            t2 <- t2$X[t2$K == unique(raw_dat$K)[j1]]
            
            # Check if the two timeseries are the same
            if (all(t1 == t2)) {
                found_orig_Ks[j1, 1] <- unique(ss$K)[i1]
                break
            }
        }
    }
    ref_K_name <- raw_dat$old_K[which(is.na(found_orig_Ks))]
    clone_dic <- data.frame(old_K = unique(raw_dat$old_K), K = as.numeric(found_orig_Ks), stringsAsFactors = F)
    clone_dic$K[is.na(clone_dic$K)] <- K-1
    
    # Join by X
    #@Sohrab, new
    #ss$K <- NULL
    # lookup_time <- 0
    # temp <- dplyr::right_join(ss[ss$time == lookup_time, ], raw_dat[raw_dat$time == lookup_time, c('old_K', 'X')], by=('X'))
    # 
    # for (ti in seq(unique(ss$time))) {
    #   print(ti)
    #   if (any(duplicated(temp$old_K))) {
    #     lookup_time <- unique(ss$time)[ti]
    #     print(lookup_time)
    #     temp <- dplyr::right_join(ss[ss$time == lookup_time, ], raw_dat[raw_dat$time == unique(raw_dat$time)[ti], c('old_K', 'X')], by=('X'))
    #   } else {
    #     break
    #   }
    #   # print(ti)
    # }
    
    # temp1 <- ss[ss$time == lookup_time, ]
    # colnames(temp1) <- c('nottime', 'notK', 'X')
    # raw_dat <- dplyr::right_join(temp1, temp, by=('X'))
    # 
    # # Find the reference K
    # ref_K_name <- unique(as.character(raw_dat$old_K[which(is.na(raw_dat$time))]))
    # 
    # raw_dat$K[is.na(raw_dat$K)] <- K-1
    
    # Generate get_clone_dic_for_datatag
    # clone_dic <- raw_dat[, c('K', 'old_K')]
    
    clone_dic$is_ref <- clone_dic$old_K == ref_K_name
    # New names
    clone_dic <- clone_dic[order(unique(clone_dic$old_K)), ]
    clone_dic$letters = LETTERS[1:nrow(clone_dic)]
    clone_dic$pretty_names <- get_pretty_names_for_loci(clone_dic$old_K)

    # New
    # Reconcile with the clone_dic used to generate tree plots
    get_dummar_clone_dic <- function(datatag, edge_list_path) {
        print(edge_list_path)
        g <- read_ltm_tree(edge_list_path)
        the_cut <- get_the_cut(datatag, edge_list_path, g=g)
        aug_cut <- add_siblings_cut(datatag = datatag, edge_list_path = edge_list_path, the_cut = the_cut, the_g = g)
        new_clades <- setdiff(names(aug_cut), unlist(the_cut))
        tcc <- get_tcc(g)
        ts <- get_edge_data(g, tcc, datatag, edge_list_path)
        sugar_dummy_clone_dic(ts, aug_cut, new_clades, tcc)
    }
    if (datatag != 'SA501') {
        if (datatag %in% get_treatment_datanames()) {
            sugar_dummy_clone_dic <- get_dummar_clone_dic('SA609', get_edge_list_path_for_datatag('SA609'))
        } else {
            sugar_dummy_clone_dic <- get_dummar_clone_dic(datatag, edge_list_path)
        }
    } else {
        sugar_dummy_clone_dic <- get_clone_dic_for_sa501_bulk()
    }
    clone_dic$letters <- NULL
    sugar_dummy_clone_dic$is_ref <- NULL
    sugar_dummy_clone_dic$pretty_names <- NULL
    clone_dic <- dplyr::inner_join(clone_dic, sugar_dummy_clone_dic, by='old_K')
    clone_dic <- clone_dic[order(clone_dic$letters), ]
    saveRDS(clone_dic, outpath)
    return(clone_dic)
}

get_split_cut <- function(datatag, edge_list_path, min_tresh) {
    g <- read_ltm_tree(edge_list_path = edge_list_path)
    subgraphs <- get_full_subgraphs_for_aug_cut(g = g, datatag = datatag, edge_list_path = edge_list_path)
    cn <- get_clone_letters(datatag, edge_list_path)
    splitted_clades <- list()
    clone_dic <- NULL
    for (sg in names(subgraphs)) {
        clade_letter <- names(unlist(cn))[unlist(cn) == sg] 
        print(sprintf('Processing %s = %s', clade_letter, sg))
        splitted_clade <- split_clade(sg = subgraphs[[sg]], min_tresh = min_tresh)
        splitted_clades[[sg]] <- splitted_clade[[sg]]
        dummy_clone_dic <- get_clone_dic_for_spliited_clade(g = g, datatag = datatag, edge_list_path = edge_list_path, splitted_clade = splitted_clade)
        dummy_clone_dic <- dummy_clone_dic[grepl(paste0(clade_letter, '[0-9]+'), dummy_clone_dic$letters), ]  
        if (is.null(clone_dic)) {
            clone_dic <- dummy_clone_dic
        } else {
            clone_dic <- rbind(clone_dic, dummy_clone_dic)  
        }
    }
    str(splitted_clades)
    str(clone_dic)
    list(splitted_clades=splitted_clades, clone_dic=clone_dic)
}

split_clade <- function(sg, min_tresh=.2) {
    tcc <- get_tcc(g = sg)
    count_tresh <- min_tresh * sum(tcc$Freq)
    
    # 1. Find the lowest min edge (LME) -- the major clade in each branch
    new_clades <- find_subclones(some_g = sg, new_clades = list(), count_tresh=count_tresh)
    stopifnot(!any(duplicated(unlist(new_clades))))
    
    # 2. Find the reminder structure
    # Find all parents, except for the root
    if (is.null(V(sg)$height)) {
        h <- get_height_dat(sg)
        stopifnot(all(h$id == V(sg)$name))
        V(sg)$height <- as.integer(h$dist)
    }
    root <- V(sg)$name[V(sg)$height == 0]
    the_parents <- list()
    for (ii in names(new_clades)) {
        the_parents[[ii]] <- setdiff(find_all_parents(the_g = sg, node = ii, max_height = V(sg)$height[V(sg)$name == ii]), root)
        if (length(the_parents[[ii]]) == 0) {
            the_parents[[ii]] <- NULL
        }
    }
    # For the nodes that are parent/child, ditch the child; 
    children <- c()
    if (length(the_parents) > 1 & length(new_clades) > 1) {
        children <- get_children_from_parents(g = sg, parents = the_parents)
    }
    if (length(children) > 0) {
        for (ii in children) {
            the_parents[[ii]] <- NULL
        }
    }
    private_parents <- list()
    for (ii in names(the_parents)) {
        private_parents[[ii]] <- the_parents[[ii]]
        for (jj in names(the_parents)) {
            if (ii != jj) {
                private_parents[[ii]] <- setdiff(private_parents[[ii]], c(the_parents[[jj]], jj))
                if (length(private_parents[[ii]]) == 0) {
                    private_parents[[ii]] <- NULL
                }
            }
        }
    }
    for (ii in names(private_parents)) {
        # Pick the elder parent
        if (length(private_parents[[ii]]) > 0) {
            ep <- pick_elder_parent(g=sg, nodes = private_parents[[ii]])
            new_clades[[ep]] <- setdiff(get_desc_names(datatag = datatag, the_graph = sg, the_node = ep)[[1]], unlist(unname(new_clades)) )
            if (length(new_clades[[ep]]) == 0) {
                new_clades[[ep]] <- NULL
            }
        }
    }
    
    # 3. Everyone else, including the immediate cell-children of the root will form the last clade
    all_else <- setdiff(V(sg)$name, unlist(new_clades))
    if (length(all_else) > 0) {
        if (!(root %in% names(new_clades))) {
            new_clades[[root]] <- all_else
        }
    }
    new_clades <- drop_non_cells(new_clades)
    # Drop empty clades
    for (ii in names(new_clades)) {
        if (length(new_clades[[ii]]) == 0) {
            new_clades[[ii]] <- NULL
        }
    }
    res <- list()
    res[[root]] <- new_clades
    res
}

get_clone_dic_for_spliited_clade <- function(g = NULL, datatag, edge_list_path, splitted_clade) {
    # Setup the colours as shades of the original colour
    if (is.null(g))
        g <- read_ltm_tree(edge_list_path)
    
    tcc <- get_tcc(g)
    ts <- get_edge_data(g, tcc, datatag, edge_list_path)
    the_cut <- get_the_cut(datatag, edge_list_path, g=g)
    aug_cut <- add_siblings_cut(datatag = datatag, edge_list_path = edge_list_path, the_cut = the_cut, the_g = g)
    new_clades <- setdiff(names(aug_cut), unlist(the_cut))
    dummy_clone_dic <- cosmic_clone_dic(datatag = datatag, edge_list_path = edge_list_path)
    candiate_edges <- get_candi_edge(the_cut = the_cut, aug_cut = aug_cut, ts = ts, tcc = tcc, new_clades = new_clades)
    
    # Set the colours in ... with the legacy colours  
    myColors <- get_cluster_colours(length(unique(candiate_edges$clone_id)))
    myColors <- myColors[1:length(unique(candiate_edges$clone_id))]
    names(myColors) <- unique(as.character(candiate_edges$clone_id))[order(get_pretty_names_for_loci(unique(candiate_edges$clone_id)))]
    cdf <- data.frame(old_K = names(myColors), color = unname(myColors))
    dummy_clone_dic <- dplyr::right_join(dummy_clone_dic, cdf) 
    
    ## b.1 remove the old clade, and add the new clades as subsets of that
    sc <- splitted_clade[[1]]
    n_split_clade <- length(names(sc))
    old_clade_name <- names(splitted_clade)
    old_clade_dic <- dummy_clone_dic[dummy_clone_dic$old_K == old_clade_name, ]
    
    # Generate different opacities of the same colour to differentiate
    print('Using distinguishable colours...')
    
    alpha_colours <- viridis_pal(option = "C")(n_split_clade)
    tmp_dic <- data.frame(old_K = names(sc), 
                                                is_ref = old_clade_dic$is_ref, 
                                                letters = paste0(old_clade_dic$letters, seq(length(sc))),
                                                pretty_names = get_pretty_names_for_loci(names(sc)), 
                                                color = alpha_colours)
    dummy_clone_dic <- dummy_clone_dic[dummy_clone_dic$old_K != old_clade_name, ]
    dummy_clone_dic <- dplyr::bind_rows(dummy_clone_dic, tmp_dic)
    dummy_clone_dic
}

find_subclones <- function(some_g, new_clades, count_tresh) {
    print(length(new_clades))
    n_tcc <- get_tcc(some_g)
    new_ts <- get_edge_data(g = some_g, tcc = n_tcc, datatag = NULL, edge_list_path = NULL)
    new_ce <- new_ts %>% dplyr::filter(N > count_tresh) %>% dplyr::arrange(desc(height))
    # see if anyone matches the criterion
    if (nrow(new_ce) > 1) { # There should be one more than the root
        # remove ii from the graph
        node <- new_ce$clone_id[1]
        print(node)
        new_clades[[node]] <- get_desc_names(datatag = datatag, the_graph = some_g, the_node = node)[[1]]
        # TODO: check this out: if the parent has only one locus child, cover all your cell siblings too...
        pr <- find_parent(some_g, node)
        fic <- find_immediate_children(some_g, pr)
        if (length(grep('locus_', fic)) == 1) {
            new_clades[[pr]] <- c(new_clades[[node]], fic[!grepl('locus_', fic)])
            new_clades[[node]] <- NULL
            node <- pr
            print('Moving up to the parent')
        }     
        vids <- setdiff(V(some_g)$name , new_clades[[node]])
        if (length(vids) > 1 & sum(grepl('locus_', vids)) > 1) {
            new_sg <- setup_graph(igraph::induced.subgraph(graph = some_g, vids = vids))
            new_clades <- find_subclones(some_g = new_sg, new_clades = new_clades, count_tresh = count_tresh)
        } else {
            return(new_clades)
        }
    } else {
        return(new_clades)
    }
}