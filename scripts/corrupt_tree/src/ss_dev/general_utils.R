### Parsing inputs
################################################################################################################################################
read_ltm_tree <- function(edge_list_path) {
    ss <- read.csv(edge_list_path, stringsAsFactors = F, header=T)
    # Find the root
    dd <- (as.matrix(ss))
    g <- graph_from_edgelist(dd)
    V(g)$id <- seq(vcount(g))
    return(g)
}

get_tcc <- function(g) {
    tcc <- as.data.frame(table(parse_cell_names(cell_names = V(g)$name[!grepl('root', V(g)$name) & !grepl('locus', V(g)$name)])))
    # Remove GM cells
    tcc <- tcc[grep('X[0-9][0-9]*', tcc$Var1), ]
    tcc[order(as.numeric(gsub('X', '', tcc$Var1))), ]
}

parse_cell_names <- function(cell_names) {
    # Find libids
    libids <- libid_from_cell_id(cell_names)
    cell.dat <- data.frame(library_id = libids, cellids = cell_names, stringsAsFactors = F)
    if (datatag %in% c('SA922', 'SA922n')) {
        core_dat <- load_core_dat(datatag)
    } else {
        core_dat <- load_core_dat()
    }
    
    cell.dat <- dplyr::left_join(cell.dat, core_dat[, c('timepoint', 'library_id')], by=c('library_id'))
    stopifnot(all(cell_names == cell.dat$cellids, na.rm = T))
    
    # Set SA928 to itself
    cell.dat$timepoint[grepl('SA928', cell.dat$cellids)] <- cell.dat$cellids[grepl('SA928', cell.dat$cellids)]
    
    cell.dat$timepoint
}

libid_from_cell_id <- function(cell_id) {
    if (any(grepl('SA609', cell_id)) | any(grepl('SA532', cell_id))) {
        gsub('SA([0-9]|[A-Z]|[a-z])+-(A[0-9|A-Z]+)-.*', '\\2', cell_id)
    } else {
        gsub('SA([0-9|a-z]+)-(A[0-9|A-Z]+)-.*', '\\2', cell_id)
    } 
}

load_core_dat <- function(datatag = NULL) {
    core_dat <- readRDS(get_core_dat_path())
    
    if (!is.null(datatag)) {
        if (datatag == 'SA666') {
            core_dat <- core_dat %>% dplyr::filter(!is.na(drug))
        } else {
            core_dat <- core_dat %>% dplyr::filter(datatag == !!datatag) %>% as.data.frame()
        }
    }
    
    return(core_dat)
}

get_core_dat_path <- function() '~/Desktop/SC-1'

load_new_cn_data <- function(datatag, filter_by_timepoint=NULL) {
    in_dir <- get_in_dir_for_datatag(datatag)
    
    out_path <- file.path(in_dir,  'cnv_data.rds')
    if (file.exists(out_path)) 
        return(readRDS(out_path))
    
    dat <- as.data.frame(fread(file.path(in_dir, 'cnv_data.csv')))

    # Filter for cells in timepoint
    if (!is.null(filter_by_timepoint)) {
        dat <- dat[grepl(paste0(datatag, timepoint), dat$single_cell_id), ]
    }
    
    # tidy to wide
    value.var = 'copy_numer'
    if (is.null(dat[[value_var]])) 
        value.var = 'state'
    
    if (is.null(dat$single_cell_id)) {
        dat$single_cell_id <- dat$cell_id
        dat$cell_id <- NULL
    }
    
    mat <- reshape2::dcast(dat, chr+start+end ~ single_cell_id, value.var = c(value.var))
    rownames(mat) <- paste0(mat$chr, '_', mat$start, '_', mat$end)
    mat <- mat[, -c(1:3)]
    mat <- mat[, grepl(datatag, colnames(mat))]
    saveRDS(mat, out_path)
    mat
}

get_in_dir_for_datatag <- function(datatag) {
    sprintf('~/Desktop/SC-1311/%s/processed_data', datatag)
}

load_main_configs <- function(datatag = NULL, backuppath = NULL) {
    configs <- yaml.load_file(get_main_config_path())
    if (!is.null(backuppath)) {
        configs <- yaml.load_file(backuppath)
        print('WARNING! USING backup path to load the configs...')
    }
    
    if (!is.null(datatag)) {
        configs <- configs[[datatag]]
    }
    configs
}

get_main_config_path <- function() return('~/projects/fitness/fitclone_path_configs.yaml')

get_clone_dic_for_sa501_bulk <- function() {
    dat <- load_sa501_snv_dat()
    K <- nrow(dat)
    clones <- dat$K
    clone_dic <- data.frame(old_K = clones, is_ref = FALSE, letters = clones, pretty_names = clones, stringsAsFactors = F)
    clone_dic
}

load_sa501_snv_dat <- function() {
    in_path <- '~/Desktop/SC-1311/SA501/processed_data/sa501_clonal_prev.rds'
    if (file.exists(in_path)) return(readRDS(in_path))
    
    dat <- read.csv('/Users/sohrabsalehi/projects/fitclone/figures/raw/supp/bulk_estimated_clonal_prevalences.csv', stringsAsFactors = F)
    dat$X <- NULL
    
    # Convert to compatible foramt (Clone by time)
    times <- dat$time
    dat$time <- NULL
    dat
    dat <- as.data.frame(t(dat))
    
    colnames(dat) <- times
    dat$K <- rownames(dat)
    dat <- dat[, c(7, 1:6)]
    
    saveRDS(dat, in_path)
    dat
}

load_split_cut_configs <- function(datatag = NULL) {
    configs <- yaml.load_file('~/projects/fitness/fitclone_split_cut.yaml')
    if (!is.null(datatag)) {
        configs <- configs[[datatag]]
    }
    configs
}

### Graph general
################################################################################################################################################
get_height_dat <- function(g, edge_list_path=NULL) {
    if (!is.null(edge_list_path)) {
        outpath <- sprintf('%s_height.rds', edge_list_path)
        if (file.exists(outpath))
            return(readRDS(outpath))
    }
    the_root <- find_root(g)
    height_search <- bfs(g, root=c(the_root), order=TRUE, neimode = 'out', unreachable=FALSE, dist = TRUE, rank=TRUE, succ=TRUE, father=TRUE, pred=TRUE)
    res <- data.frame(id=names(height_search$dist), dist=as.numeric(height_search$dist), stringsAsFactors = F)
    if (!is.null(edge_list_path)) {
        saveRDS(res, outpath)
    }
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

get_leaves_names_from_graph <- function(the_g, only_loci=FALSE) {
    edge_list <- as.data.frame(igraph::as_edgelist(the_g), stringsAsFactors = F)
    colnames(edge_list) <- c('source', 'target')
    res <- get_leaves_names(edge_list)
    if (only_loci) {
        res <- grep('locus_', res, value = T)
    }
    res
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
retrieve_time_point_per_clade <- function(the_graph, decendents=NULL, data_tag=NULL, trim_cell_name=TRUE, include_loci=FALSE) {
    if (is.numeric(unlist(decendents)))
        decendents <- node_number_to_label(the_graph, decendents)
    res <- list()
    for (clade_root in names(decendents)) {
        cells <- decendents[[clade_root]]
        if (trim_cell_name) {
            cells <- parse_cell_names(cell_names = cells)
        }
        res[[clade_root]] <- cells
    }
    res
}

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

get_time_series_for_cluster_cut <- function(res, total_cell_counts_break_down=NULL, convert_to_freq=TRUE) {
    time_series_dat <- get_time_series(clade_desc_list = res, timepoints = as.character(total_cell_counts_break_down$time))
    if (convert_to_freq)
        time_series_dat = convert_to_freq(time_series_dat = time_series_dat, total_cell_counts_breakdown = total_cell_counts_break_down)
    order_cols(time_series_dat)
}

get_time_series <- function(clade_desc_list, timepoints=NULL) {
    # Also take the clade root into account
    time_series_dat <- NULL
    for (clade_root in names(clade_desc_list)) {
        break_down <- table(clade_desc_list[[clade_root]])
        temp <- data.frame(clone_id = clade_root, as.list(break_down), stringsAsFactors = F)
        
        if (is.null(time_series_dat)) {
            time_series_dat <- temp
        } else
            time_series_dat <- plyr::rbind.fill(time_series_dat, temp)
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
        time_series_dat <- time_series_dat[, c(1, 1+order(colnames(time_series_dat)[2:ncol(time_series_dat)]))]
    }
    
    for (clade_root in names(clade_desc_list)) {
        root_label <- parse_cell_names(clade_root)
        if (!grepl('SA', root_label)) {
            if (root_label %in% colnames(time_series_dat))
                time_series_dat[time_series_dat$clone_id == clade_root, root_label] <- time_series_dat[time_series_dat$clone_id == clade_root, root_label] + 1
            else {
                #print(sprintf('Warning! %s was not added to the count table.', root_label))
            }
        }
    }
    
    order_cols(ts = time_series_dat)
}

oscillation_for_ts <- function(ts) {
    diff = abs(ts[, 2:ncol(ts)] - ts[, 1:(ncol(ts)-1)])
    rowSums(diff)
}

### Naming files
################################################################################################################################################
get_pretty_names_for_loci <- function(str_array) {
    str_array <- gsub('locus_', '', str_array)
    chr = gsub('([0-9]+|X|Y)_.*', '\\1', str_array)
    str_array <- gsub('.*_([0-9]+)_.*', '\\1', str_array)
    str_array <- substr(str_array, 1, 4)
    str_array <- paste0('chr', chr, '_', str_array, '')
    gsub('chrroot_root', 'root', str_array)
}

get_pretty_str_from_funny_str <- function(input_str = NULL, keep_timestamp = TRUE) {
    if (is.null(input_str))
        input_str <- generate_random_funny_string()
    # Remove underscore and add the tag in paranthesis
    tokens <- strsplit(input_str, '__')[[1]]
    pretty_str <- stringr::str_to_title(strsplit(tokens[[2]], '_')[[1]])
    res_str <- paste0(pretty_str, collapse = ' ')
    if (keep_timestamp) {
        timestr <- tokens[length(tokens)]
        timestr <- get_pretty_str_for_time_str(timestr)
        res_str <- paste0(res_str, ' (', timestr,  ')')
    } 
    res_str
}

generate_random_funny_string <- function() {
    r1 <- as.character(base::sample(1:1000, 1))
    r2 <- as.character(base::sample(1:100, 1))
    
    name_str <- gsub(' ', '_', paste0('__', r1, '_', r2, '__'))
    paste0(name_str, generate_time_str())
}

get_pretty_str_for_time_str <- function(timestr) {
    posix_time <- as.POSIXct(x = timestr, format = "%Y%m%d%H-%M-%S")
    format(posix_time, "%Y%m%d%H-%M-%S")
    format(posix_time, "%a %b %X %Y")
}

generate_time_str <- function(n = 1) {
    format(Sys.time(), "%Y%m%d%H-%M-%S")
}

generate_random_str <- function(n = 1) {
    a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
    paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
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

load_cut_configs <- function(datatag = NULL) {
    configs <- yaml.load_file('~/projects/fitness/fitclone_cut_configs.yaml')
    if (!is.null(datatag)) {
        configs <- configs[[datatag]]
    }
    configs
}

get_cut_outpath <- function(datatag, edge_list_path, cut_configs = NULL) {
    if (is.null(cut_configs)) cut_configs <- load_cut_configs(datatag = datatag)
    sprintf('%s_cut_%s_%.2f_%.2f_%s.rds', edge_list_path, datatag, cut_configs$minimum_fraction, cut_configs$maximum_frac, tolower(gsub('( |\\.)', '', cut_configs$oscilation)))
}

order_by_chr_name <- function(the_array) {
    the_array[the_array == 'X'] <- 44
    the_array[the_array == 'Y'] <- 100
    
    order(as.numeric(the_array))
}

strip_chr_names <- function(the_array) {
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

g.count.cells <- function(g=NULL, edge_list_path=NULL) {
    length(g.get.cells(g, edge_list_path))
}

g.get.cells <- function(g=NULL, edge_list_path=NULL) {
    if (is.null(g)) g <- read_ltm_tree(edge_list_path = edge_list_path)
    get_cells(V(g)$name)
}

get_cells <- function(str_arry) str_arry[!grepl('locus_', str_arry)]

get_cluster_colours <- function(nClusters) {
    if (nClusters > 8) {
        clust_colours <- colorRampPalette(brewer.pal(8, "Set2"))(nClusters)
    } else {
        clust_colours <- brewer.pal(nClusters, "Set2")
    }
    clust_colours
}

legacy_colours_for_CN_heatmap <- function() {
    # From CN = 0 to CN = 11
    c('#4880B8', '#A7C9DF','#CCCCCC','#F5CE93','#ED9364','#D2553E','#A42116','#8B1A43','#CB3576','#D06CAD','#C196C4','#D0BAD8')
}

###Utils utils
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

parse_bin_names <- function(bin_names) {
    # Remove corrupt_tree locus tag if it's there
    bin_names <- gsub('locus_', '', bin_names)
    chr <- gsub('([0-9]+|X|Y)_[0-9]+_[0-9]+', '\\1', bin_names)
    start <- as.numeric(gsub('([0-9]+|X|Y)_([0-9]+)_[0-9]+', '\\2', bin_names))
    end <- as.numeric(gsub('([0-9]+|X|Y)_([0-9]+)_([0-9]+)', '\\3', bin_names))
    data.frame(chr = chr, start = start, end = end)
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
                    dist_to_median[i, j] <- mean(abs(mg_mat[i, ] - mg_mat[j, ]))
                }
            }
        }
    }
    dist_to_median
}

is_outlier <- function(x, prob, type='t') {
    tt <- outliers::scores(x = x, type=type, prob=prob)
    ff <- outliers::scores(x = x, type=type) < 0
    tt & ff
}

merge_clones <- function(cdat) {
    print(cdat)
    list_of_clones <- list()
    for (i_cn in seq(nrow(cdat))) {
        xx <- unname(unlist(cdat[i_cn, , drop=TRUE]))
        p1 <- who_has_it(xx[1], list_of_clones)
        p2 <- who_has_it(xx[2], list_of_clones)
        if (is.na(p1) & is.na(p2)) {
            list_of_clones[[length(list_of_clones)+1]] <- xx
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
    
    ideal_string <- gsub(sprintf('.*(%s.*%s).*', i1, i2), '\\1', paste0(LETTERS, collapse = ''))
    if (nchar(ideal_string) != length(a_list))
        return(FALSE)
        
    return(TRUE)
}

find_non_consecutive_ones <- function(a_list, orig_name) {
    if (length(a_list) < 2)  {
        return(NULL)
    }

    a_list <- sort(a_list)
    i1 <- a_list[[1]] # start
    i2 <- a_list[[length(a_list)]] # end
    
    ideal_string <- gsub(sprintf('.*(%s.*%s).*', i1, i2), '\\1', paste0(LETTERS, collapse = ''))
    ideal_list <- strsplit(ideal_string, '')[[1]]
    
    # Use the non-existing ones as delimiter to divide -- all after should not appear
    del <- ideal_list[which(!(ideal_list %in% a_list))]
    if (length(del) > 1) {
        print('WARNING! NOT IMPLEMENTED!!!')
    }
    
    del <- sort(del)
    new_list <- list()
    remaining_list <- ideal_list

    for (ii in del) {
        print(ii)
        ii_i <- which(remaining_list == ii)
        if (length(ii_i) > 0) {
            if (ii_i == 1) {
                # don't add this one
            } else {
                new_name <- sprintf('%s_(%s)', orig_name, paste0(remaining_list[1:(ii_i-1)], collapse = '+'))
                print(new_name)
                new_list[[new_name]] <- remaining_list[1:(ii_i-1)]
            }
            remaining_list <- remaining_list[-c(1:ii_i)]
        } else {
            print('ERROR!')
        }
    }
    if (length(remaining_list) > 0) {
        new_name <- sprintf('%s_(%s)', orig_name, paste0(remaining_list, collapse = '+'))
        new_list[[new_name]] <- remaining_list
    }

    return(new_list)
}

get_dfs_queue <- function(the_g, edge) {
    s1 = dfs(the_g, root=edge, order=TRUE, neimode = 'out', unreachable=FALSE)
    qq = na.omit(s1$order)
    unlist(lapply(seq(length(qq)), function(x) qq[x]$name))
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

get_edge_list_path_for_datatag <- function(dt) {
    sa609_config <- load_main_configs(datatag = dt)
    sa609_config$all_cells$edge_list_path
}

get_fitness_exp_path_for_datatag <- function(dt) {
    sa609_config <- load_main_configs(datatag = dt)
    dt_batch_path <- file.path('~/Desktop/SC-1311', dt, 'batch_runs', sa609_config$all_cells$batch_name)
    dt_exp_path <-  file.path(dt_batch_path, 'outputs', sa609_config$all_cells$exp_dir)
    file.path(dirname(dt_exp_path), sa609_config$all_cells$fitclone_exp_dir) 
}

get_treatment_datanames <- function() {
    c("SA609aRx8p", "SA609aRx4p", "SA609bRx8p")
}

get_obsered_times <- function(exp_path, keep_all_tps = FALSE) {
    original_dat <- read_yaml(file.path(exp_path, 'config.yaml'))$original_data
    ss = read.table(paste0(original_dat, '.gz'), header = T)
    
    # Find original-observed times (remove the last two)
    observed_times <- unique(ss$time)
    observed_times
}

get_dummy_candid_edge <- function(g, aug_cut) {
    ntotal <- sum(get_tcc(g)$Freq)
    frac <- unlist(lapply(names(aug_cut), function(x) length(aug_cut[[x]])/ntotal))
    data.frame(clone_id = names(aug_cut), frac = frac, stringsAsFactors = F)
}

get_exp_path_for_datatag <- function(dt) {
    sa609_config <- load_main_configs(datatag = dt)
    dt_batch_path <- file.path('~/Desktop/SC-1311', dt, 'batch_runs', sa609_config$all_cells$batch_name)
    file.path(dt_batch_path, 'outputs', sa609_config$all_cells$exp_dir)
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

get_dummy_clone_dic <- function(candi_edge) {
    data.frame(old_K = candi_edge$clone_id, 
                         is_ref = FALSE, letters = LETTERS[seq(nrow(candi_edge))],
                         pretty_names = get_pretty_names_for_loci(candi_edge$clone_id), stringsAsFactors = F)
}