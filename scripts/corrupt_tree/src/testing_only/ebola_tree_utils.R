source_dir <- '/home/htran/Projects/farhia_project/rscript/dlp/fitness/'
source(paste0(source_dir,'fitclone_corrupt_utils.R'))
source(paste0(source_dir,'fitclone_autocut.R'))
source(paste0(source_dir,'fitclone_tree.R'))
source(paste0(source_dir,'fitclone_heatmaps.R'))
source(paste0(source_dir,'fitclone_heatmap_utils.R'))
source(paste0(source_dir,'fitclone_s.R'))
source(paste0(source_dir,'fitclone_trajectories.R'))
source(paste0(source_dir,'fitclone_private_cnv.R'))
source(paste0(source_dir,'fitclone_paths.R'))
source(paste0(source_dir,'fitclone_cell_assignment_utils.R'))
source(paste0(source_dir,'fitclone_genotype_cut.R'))

get_aug_clustering <- function(results_dir, datatag, edge_list_path, 
                               the_mat=NULL, remove_gm_cells=TRUE,
                               minimum_fraction=0.02,maximum_fraction=0.55) {
  g <- read_ltm_tree(edge_list_path)
  # the_cut <- get_the_cut(datatag, edge_list_path, g=g)
  the_cut <- get_the_cut(g, minimum_fraction, maximum_fraction)
  # aug_cut <- add_siblings_cut(datatag = datatag, edge_list_path = edge_list_path, the_cut = the_cut, the_g = g)
  aug_cut <- add_siblings_cut(the_cut, g, minimum_fraction,maximum_fraction,oscilation=0)
  
  
  new_clades <- setdiff(names(aug_cut), unlist(the_cut))
  tcc <- get_tcc(g, results_dir)
  # ts <- get_edge_data(g, tcc, datatag, edge_list_path)
  ts <- get_edge_data(g)
  
  clone_dic <- sugar_dummy_clone_dic(ts, aug_cut, new_clades, tcc)
  # clone_dic <- get_clone_dic(datatag = datatag, exp_path = get_fitness_exp_path_for_datatag(datatag), edge_list_path = get_edge_list_path_for_datatag(datatag), use_cashe = TRUE)
  # function() {
  #   stop('sohrab')
  #   get_ref_clone_dic_sa1000(use_cashe = FALSE)
  # }
  corrupt_ca <- create_cell_assignment_dat(results_dir=results_dir, 
                                           datatag = datatag, 
                                           edge_list_path = edge_list_path, 
                                           clone_dic = clone_dic, g = g, 
                                           the_cut = NULL, aug_cut = aug_cut)
  
  # if (datatag == 'SA666') {
  #   corrupt_ca <- corrupt_ca[is_treatment_cell(corrupt_ca$single_cell_id), ]
  # }
  
  # peekm(corrupt_ca)
  # table(corrupt_ca$letters)
  
  # Add GM cells back
  # if (is.null(the_mat)) {
  #   the_mat <- load_new_cn_data(datatag = datatag)
  # }
  # 
  # gm_cells <- grep('SA928', colnames(the_mat), value = T)
  # if (length(gm_cells) > 0) {
  #   corrupt_ca <- bind_rows(corrupt_ca, data.frame(single_cell_id = gm_cells, genotype=gm_cells, time=gm_cells, is_ref=FALSE, letters='Un', pretty_names=gm_cells, K=NA))
  # }
  # 
  # if (remove_gm_cells) {
  #   corrupt_ca <- corrupt_ca[!(corrupt_ca$single_cell_id %in%  gm_cells), ]
  # }
  
  corrupt_ca <- corrupt_ca[!grepl('locus_', corrupt_ca$single_cell_id), ]
  corrupt_ca <- corrupt_ca[!is.na(corrupt_ca$single_cell_id), ]
  corrupt_ca
  
  # dim(corrupt_ca)
  # saveRDS(corrupt_ca, file = paste0(save_dir,'corrupt_ca.RDS'))
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

get_libid_from_cell_names <- function(cell_names) {
  # TODO: add the condition for SA906 and SA666
  gsub('SA([0-9]|[A-Z]|[a-z])+-(A([0-9]|[A-Z])+)-.*', '\\2', cell_names)
}

get_library_labels <- function(cell_ids) {
  labels <- sapply(strsplit(cell_ids, "-"), function(x) {
    return(x[2])
  })
  return(labels)
}

get_timepoint <- function(results_dir){
  grouping_df <- read.csv(paste0(results_dir,'library_groupings.csv'),header=T, check.names=F, 
                          row.names = 1, stringsAsFactors=F)
  colnames(grouping_df)[which(names(grouping_df) == "pdxid")] <- "timepoint"
  colnames(grouping_df)[which(names(grouping_df) == "mainsite")] <- "treatmentSt"
  colnames(grouping_df)[which(names(grouping_df) == "grouping")] <- "library_id"
  return(grouping_df)
}

parse_cell_names <- function(cell_names, results_dir) {
  # Find libids
  #libids <- libid_from_cell_id(cell_names)
  libids <- get_libid_from_cell_names(cell_names)
  cell.dat <- data.frame(library_id = libids, cellids = cell_names, stringsAsFactors = F)
  # if (datatag %in% c('SA922', 'SA922n')) {
  #   core_dat <- load_core_dat(datatag)
  # } else {
  #   core_dat <- load_core_dat()
  # }

  # cell.dat <- dplyr::left_join(cell.dat, core_dat[, c('timepoint', 'library_id')], by=c('library_id'))
  # stopifnot(all(cell_names == cell.dat$cellids, na.rm = T))

  grouping_df <- get_timepoint(results_dir)
  cell.dat <- dplyr::left_join(cell.dat, grouping_df[, c('timepoint', 'library_id')], by=c('library_id'))

  # Set SA928 to itself
  # cell.dat$timepoint[grepl('SA928', cell.dat$cellids)] <- cell.dat$cellids[grepl('SA928', cell.dat$cellids)]

  cell.dat$timepoint
}
# retrieve_time_point_per_clade <- function(results_dir, the_graph, decendents=NULL, data_tag=NULL, 
#                                           trim_cell_name=TRUE, include_loci=FALSE) {
#   if (is.numeric(unlist(decendents)))
#     decendents <- node_number_to_label(the_graph, decendents)
#   res <- list()
#   for (clade_root in names(decendents)) {
#     #print(clade_root)
#     # clade_root = names(decendents)[[1]]
#     cells <- decendents[[clade_root]]
#     # cells <- grep('SA609', cells, value = T)
#     if (trim_cell_name) {
#       # cells <- parse_cell_names(cell_names = cells)
#       cells <- parse_cell_names(cell_names = cells, results_dir)
#     }
#     res[[clade_root]] <- cells
#   }
#   res
# }

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
  res <- retrieve_time_point_per_clade(the_graph = the_graph, 
                                       decendents = desc)
  
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

h1_auto_cut <- function(the_g, candiate_edges, minimum_fraction, maximum_frac = .5, use_soft_down_up=FALSE, ignore_down_up=TRUE, ignore_dip_to_zero=FALSE) {
  # the_g = g
  # Filter by oscilacion & minium fraction
  #candiate_edges <- ts[ts$frac > minimum_fraction & ts$oscillation > summary(ts$oscillation)['3rd Qu.'], ]
  #& ts$dip_to_zero_and_rise == FALSE & ts$down_and_up_for_ts_soft == FALSE
  
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
    # edge = qeueu[[1]]
    # print(edge)
    tmp <- candiate_edges[candiate_edges$clone_id == edge, ]
    
    if (ignore_dip_to_zero)
      tmp$dip_to_zero_and_rise = FALSE
    
    # If passed
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
        #edge = 'locus_12_131000001_131500000'
        t1 = find_all_parents(tau$graph, edge, tmp$height)
        t2 = candiate_edges[candiate_edges$clone_id %in% t1, ]
        degrees = igraph::degree(graph = tau$graph, v = t1, mode = 'out')
        degrees = data.frame(clone_id=names(degrees), degree=as.vector(degrees))
        #degrees <- degrees[degrees$clone_id != 'root', ]
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




get_tcc <- function(g, results_dir) {
  tcc <- as.data.frame(table(parse_cell_names(cell_names = V(g)$name[!grepl('root', V(g)$name) & !grepl('locus', V(g)$name)],results_dir)))
  # Remove GM cells
  tcc <- tcc[grep('X[0-9][0-9]*', tcc$Var1), ]
  tcc[order(as.numeric(gsub('X', '', tcc$Var1))), ]
}



get_clone_dic <- function(datatag, exp_path, edge_list_path=NULL, use_cashe = TRUE) {
  if (!missing(exp_path)) {
    print(sprintf('exp_path is set to %s', exp_path))
  }
  
  # exp_path <- fitness_exp_path
  #if (datatag == 'SA666' | datatag == 'SA1000') {
  if (datatag == 'SA1000' | datatag == 'SA001') {
    #return(cosmic_clone_dic(datatag = 'SA1000', edge_list_path = get_edge_list_path_for_datatag('SA1000'), set_colours = T))
    if (file.exists(get_fitness_exp_path_for_datatag('SA609'))) {
      return(get_clone_dic('SA609', exp_path = get_fitness_exp_path_for_datatag('SA609'), edge_list_path = get_edge_list_path_for_datatag('SA609')))
    } 
    # else {
    #   return(cosmic_clone_dic(datatag = 'SA1000', edge_list_path = get_edge_list_path_for_datatag('SA1000')))
    # }
  } 
  
  if (datatag == 'SA922') {
    return(cosmic_clone_dic(datatag = datatag, edge_list_path = edge_list_path, set_colours = T))
  }
  
  # Check if the config file exists...
  result <- tryCatch({
    stopifnot(length(file.exists(file.path(exp_path, 'config.yaml'))) > 0)
    if (!missing(exp_path)) {
      stopifnot(file.exists(file.path(exp_path, 'config.yaml')))
    }
    NULL
  }, warning = function(w) {
    # warning
    print(w)
  }, error = function(e) {
    print('Cannot find the config.yaml file for the corresponding fitness analysis.')
    print('Returning the usual clone dic.')
    print(e)
    cosmic_clone_dic(datatag = datatag, edge_list_path = edge_list_path, set_colours = T)
  })
  
  if (!is.null(result)) return(result)
  
  # if (datatag == 'SA501') {
  #   clone_dic <- get_clone_dic_for_sa501_bulk()
  #   if (is.numeric(clone_dic$K)) clone_dic$K <- as.character(clone_dic$K)
  #   return(clone_dic)
  # }
  
  
  if (datatag %in% get_treatment_datanames()) {
    outpath <- sprintf('%s/%s_clone_dic.rds', exp_path, 'dummy') 
  } else {
    if (edge_list_path != '') {
      cut_path <- get_cut_outpath(datatag = datatag, edge_list_path = edge_list_path)
      aug_path <- sprintf('%s_aug_cut.rds', cut_path) 
      split_gen_path <- sprintf('%s_split_results.rds', aug_path) 
      outpath <- sprintf('%s_clone_dic.rds', split_gen_path) 
    } else {
      outpath <- sprintf('%s/%s_clone_dic.rds', get_exp_path_for_datatag(datatag), 'dummy') 
    }
  }
  
  print(sprintf('cashe for clone_dic\n%s', outpath))
  if (file.exists(outpath) & use_cashe) {
    return(readRDS(outpath))
  }
  
  # Fix SA609X3X8a
  function() {
    qq <- readRDS('/Users/sohrabsalehi/Desktop/SC-1311/plots/sep_03/redownload/new_trees/SA1000/tree__541_78__2019101017-47-28_SA609X3X8a.csv_cut_SA609X3X8a_0.04_0.30_0.rds_aug_cut.rds_split_results.rds_clone_dic.rds')
    bu <- readRDS('/Users/sohrabsalehi/projects/fitness/backup_sa1000_old_clone_dics.rds')
    ss <- bu$sa1000
    
    qq <- dplyr::bind_rows(ss[!(ss$letters %in% qq$letters), ], qq)
    qq$K[qq$letters == 'A'] <- 7
    saveRDS(qq, '/Users/sohrabsalehi/Desktop/SC-1311/plots/sep_03/redownload/new_trees/SA1000/tree__541_78__2019101017-47-28_SA609X3X8a.csv_cut_SA609X3X8a_0.04_0.30_0.rds_aug_cut.rds_split_results.rds_clone_dic.rds')
  }
  
  print('Reading config yaml file at:')
  print(file.path(exp_path, 'config.yaml'))
  config <- read_yaml(file.path(exp_path, 'config.yaml'))
  #config <- read_yaml(file.path(fitness_exp_path, 'config.yaml'))
  if (!is.null(config$K_prime)) {
    K = config$K_prime + 1
  } else {
    K = config$K + 1
  }
  
  # Find the original data file
  original_dat <- read_yaml(file.path(exp_path, 'config.yaml'))$original_data
  #original_dat <- read_yaml(file.path(fitness_exp_path, 'config.yaml'))$original_data
  ss <- read.table(paste0(original_dat, '.gz'), header = T)
  
  # Find original-observed times (remove the last two)
  #observed_times <- get_obsered_times(fitness_exp_path)
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
  
  #sa609_dec_29_dlp_time_dropped_no_time_autocut_5_clones.tsv.gz
  #sa609_dlp_time_dropped_no_time_autocut_5_clones.tsv.gz
  
  # Restric to one time point, hoping there aren't identical (mostly zero) clones to throw this off
  #raw_dat = raw_dat[raw_dat$time == 0, ]
  #raw_dat <- raw_dat[, c('old_K', 'X')]
  
  # A better way, look at the reference directory number
  # But just to make sure, use O(N^2) timeseries comparison
  # Find which of the original coordinates was removed
  
  found_orig_Ks <- matrix(NA, nrow=length(unique(raw_dat$K)), ncol = 1)
  for (j1 in seq(unique(raw_dat$K))) {
    # j1 = 1
    for (i1 in seq(unique(ss$K))) {
      # i1 = 1
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
    the_cut <- get_the_cut(results_dir, datatag, edge_list_path, g=g)
    aug_cut <- add_siblings_cut(datatag = datatag, edge_list_path = edge_list_path, the_cut = the_cut, the_g = g)
    new_clades <- setdiff(names(aug_cut), unlist(the_cut))
    tcc <- get_tcc(g)
    ts <- get_edge_data(g, tcc, datatag, edge_list_path)
    sugar_dummy_clone_dic(ts, aug_cut, new_clades, tcc)
  }
  
  #sugar_dummy_clone_dic = ss
  if (datatag != 'SA501') {
    #if (datatag %in% c(get_treatment_datanames(), 'SA000', 'SA609X3X8a')) {
    if (datatag %in% get_SA1000_datatags()) {
      sugar_dummy_clone_dic <- get_dummar_clone_dic('SA609', get_edge_list_path_for_datatag('SA609'))
    } else {
      sugar_dummy_clone_dic <- get_dummar_clone_dic(datatag, edge_list_path)
    }
  } else {
    sugar_dummy_clone_dic <- get_clone_dic_for_sa501_bulk()
  }
  # clone_dic = get_clone_dic(datatag, exp_path)
  clone_dic$letters <- NULL
  sugar_dummy_clone_dic$is_ref <- NULL
  sugar_dummy_clone_dic$pretty_names <- NULL
  
  clone_dic <- dplyr::inner_join(clone_dic, sugar_dummy_clone_dic, by='old_K')
  
  # Add the missing letters for SA1000 subgraph
  if ((datatag %in% get_SA1000_datatags()) & datatag %ni% c('SA1000', 'SA609')) {
    print(sprintf('Checking to ammend SA1000 for dt = %s', datatag))
    cd_1000 <- get_ref_clone_dic_sa1000(TRUE)
    ml <- setdiff(cd_1000$letters, clone_dic$letters)
    if (length(ml) > 0) {
      sub_dic <- cd_1000[cd_1000$letters %in% ml, ]
      sub_dic$K <- -1 
      clone_dic <- dplyr::bind_rows(clone_dic, sub_dic)
    }
  }
  
  clone_dic <- clone_dic[order(clone_dic$letters), ]
  
  # FOR SA000 with all the ones that have just zeros for K
  fix_sa000_clone_dic <- function(clone_dic) {
    if (datatag != 'SA000') return(clone_dic)
    K <- clone_dic$K
    k_set <- seq(length(K)) - 1
    rem_set <- setdiff(k_set, K)
    if (length(rem_set) > 0) {
      ii <- which(duplicated(K))
      clone_dic$K[ii] <- rem_set
    }
    clone_dic
  }
  clone_dic <- fix_sa000_clone_dic(clone_dic)
  
  if (!use_cashe) saveRDS(clone_dic, outpath)
  
  return(clone_dic)
}

read_ltm_tree <- function(edge_list_path) {
  if (edge_list_path == '') {
    return(NULL)
  }
  ss <- read.csv(edge_list_path, stringsAsFactors = F, header=T)
  # Find the root
  dd <- (as.matrix(ss))
  g <- graph_from_edgelist(dd)
  V(g)$id <- seq(vcount(g))
  g
}

create_cell_assignment_dat <- function(results_dir, datatag, edge_list_path, clone_dic, g = NULL, the_cut = NULL, aug_cut = NULL) {
  # All assigned cells
  if (is.null(aug_cut)) {
    desc <- get_decendents(the_cut, the_graph = g, min_cell_per_clone = 1)
    # res <- retrieve_time_point_per_clade(g, desc, datatag, FALSE)
    res <- retrieve_time_point_per_clade(g, desc,FALSE)
  } else {
    # res <- retrieve_time_point_per_clade(g, aug_cut, datatag, FALSE)
    res <- retrieve_time_point_per_clade(g, aug_cut,FALSE)
  }
  
  #names(res) <- gsub('locus_', '', names(res))
  
  dat <- ltm_clust_to_df(res, names(res))
  dat$time <- parse_cell_names(dat$single_cell_id, results_dir)
  dat$genotype <- as.character(dat$genotype)
  clone_dic$old_K <- as.character(clone_dic$old_K)
  dat <- dplyr::right_join(dat, clone_dic, by=c('genotype'='old_K'))
  
  # For unassigned cells
  all_cells <- V(g)$name[!grepl('root', V(g)$name) & !grepl('locus_', V(g)$name)]
  # Remove GM cells
  # TODO: centralised GM cells removal
  if (datatag %in% c('SA906a', 'SA906b')) {
    all_cells <- all_cells[grepl(datatag, all_cells)]
  } else if (datatag == 'SA922') {
    # Just do nothing
  } else {
    all_cells <- all_cells[grep('X[0-9]+', all_cells)]
  }
  
  un_assigned_cell <- setdiff(all_cells, dat$single_cell_id)
  # any(un_assigned_cell %in% dat$single_cell_id)
  #dat1 = data.frame(single_cell_id=un_assigned_cell, genotype='unassigned', time=parse_cell_names(un_assigned_cell), K=-1, is_ref=NA, letters = 'Un', pretty_names='')
  if (length(un_assigned_cell) > 0) {
    dat1 <- data.frame(single_cell_id=un_assigned_cell, genotype='unassigned', 
                       time=parse_cell_names(un_assigned_cell, results_dir), 
                       is_ref=NA, letters = 'Un', pretty_names='', K = NA, 
                       stringsAsFactors = F)
    
    dat <- dplyr::bind_rows(dat, dat1)
  }
  
  dat
}

addrow <- function(orig, tmp) {
  if (is.null(orig)) {
    orig <- tmp
  } else {
    orig <- rbind(orig, tmp)
  }
  orig
}


find_parent <- function(the_g, node) {
  igraph::neighborhood(the_g, nodes=node, mode='in')[[1]]$name[[2]]
}

