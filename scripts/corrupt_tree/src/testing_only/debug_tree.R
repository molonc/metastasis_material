project_dir <- '/home/htran/Projects/hakwoo_project/corrupt_tree'
source(paste0(project_dir, "/src/testing_only/blob_v2.R"))
# results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_new_encode/SA1035_whole_local/'
# results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_new_encode/SA1035_whole_local/'
#
# results_dir <- "/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local_tyler/"
# results_dir <- "/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local_tyler_left_padding/"
# results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_filtered/'
results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Sohrab_filtered'



# cnv_path <- paste0(results_dir, 'total_merged_filtered_states.csv')
# cell_clones <- paste0(results_dir, 'tree_cut_out/cell_clones_if_0.1_af_0.65_p0.75_e0.04.csv')

save_dir <- paste0(results_dir,'debug_tree/')
if (!file.exists(save_dir)){
  dir.create(save_dir)
}
cell_clones <- paste0(results_dir, 'cell_clones.csv')

filtered_df <- read.csv(paste0(results_dir, 'filtered.csv'), check.names = F, stringsAsFactors = FALSE)
dim(filtered_df)

chr1_bins <- grep('locus_1_', filtered_df$loci, value=T)
unique(chr1_bins)

copynumber <- paste0(results_dir, 'total_merged_filtered_states.csv')
copy_number <- read.csv(copynumber, header=T, row.names=1, check.names = F,stringsAsFactors = FALSE)
loci_use <- rownames(copy_number)

padding <- paste0(results_dir, 'bin_cnvs_corrupt_double_padding.csv')
padding_df <- read.csv(padding, header=T, check.names = F,stringsAsFactors = FALSE)
loci_chrs <- paste0(padding_df$chr,'_',padding_df$start,'_',padding_df$end)
loci_chrs <- loci_chrs[!loci_chrs %in% loci_use]
length(loci_chrs)
t <- parse_bin_names(loci_chrs)
chr11_bins <- grep('^1_', loci_use, value=T)
unique(chr11_bins)

padding_added <- as.data.frame(matrix(2, nrow = length(loci_chrs),ncol=ncol(copy_number), 
            dimnames = list(loci_chrs, colnames(copy_number))))
copy_number <- rbind(copy_number,padding_added)
write.csv(copy_number, paste0(save_dir,'total_merged_filtered_states_padding.csv'),  row.names = T, quote=F)

View(unique(filtered_df$loci)[1:30])

# get the list of left and right padding, add to copynumber column



# newick <- file.path(results_dir, 'tree.newick')
# tree <- read.tree(newick)
# tree$tip.label[1:3]
# plot(tree, vertex.size=10, vertex.label=NA) 
# ggtree(tree) + geom_point(aes(shape=isTip, color=isTip), size=1)
# length(tree$edge)
# 
# plot(tree)
# # Create a ring by hand
# graph_from_edgelist(cbind(1:10, c(2:10, 1)))
# 

# debug specific chromosomes
# sample_run <- function(newick, copy_number, cell_clones = NULL, output = NULL, use_all_loci = TRUE, filter_cells = NULL, 
#                        corrupt_input = NULL, posterior_mat = NULL, chr_filter = NULL, dev = 'png', g = NULL, drop_loci_names = FALSE) {
#   


sample_run_main_corrupt_tree(cell_clones, ctdir = results_dir, chr_filter = '7')

chr_ls <- c(8, 11, 19)
for(chr in chr_ls){
  sample_run_main_corrupt_tree(cell_clones, ctdir = results_dir, chr_filter = chr)
}

sample_run_main_corrupt_tree(cell_clones, ctdir = results_dir)

sample_run_grow_tree(ctdir = paste0(results_dir,'corrupt_grow/'), 
                     chr_filter = c(7,8), just_normal = TRUE)




sample_run_fast(ctdir = results_dir,
                copy_number = cnv_path,
                cell_clones = cell_clones,
                use_all_loci = TRUE,
                #filter_cells = bad_cells,
                # filter_cells = setdiff(jump_rank$cell_id, bad_cells),
                dev = 'png',
                # chr_filter = c(7),
                chr_filter = NULL,
                just_normal = TRUE,
                drop_loci_names = FALSE)

sample_run_main_corrupt_tree <- function(cell_clones, ctdir = NULL, chr_filter = NULL, just_normal = FALSE, 
                                         use_greedy_tree = FALSE) {
  if (!is.null(ctdir)) {
    newick <- file.path(ctdir, 'tree.newick')
    corrupt_input <- file.path(ctdir, 'filtered.csv')  # corrupt_filter func, input: straighten_output.csv, output: filtered.csv
    posterior_mat <- file.path(ctdir, 'average.csv')
    cell_clones <- file.path(ctdir, 'cell_clones.csv')
    copy_number <- file.path(ctdir, 'total_merged_filtered_states.csv')
    # copy_number <- file.path(save_dir,'total_merged_filtered_states_padding.csv')
    if (use_greedy_tree) {
      newick <- file.path(ctdir, 'consensus.newick')
    }
  }
  print('Generating the ordinary one')
  if (just_normal) {
    sample_run(newick = newick, corrupt_input = NULL, posterior_mat = NULL, ...)
    
  } else {
    # sample_run(newick = newick, corrupt_input = corrupt_input, posterior_mat = posterior_mat, ...)
    print("Crz!!!")
    sample_run(newick, copy_number, cell_clones, 
               output = NULL, use_all_loci = TRUE, filter_cells = NULL, 
               corrupt_input, posterior_mat, chr_filter, 
               dev = 'png', g = NULL, drop_loci_names = FALSE) 
    
  }
}


sample_run_grow_tree <- function(ctdir = NULL, chr_filter = NULL, just_normal = FALSE, 
                                 use_greedy_tree = FALSE, corrupt_input=NULL,
                                 posterior_mat = NULL ) {
  if (!is.null(ctdir)) {
    newick <- file.path(ctdir, 'grown.newick')
    cell_clones <- file.path(ctdir, 'cell_clones.csv')
    copy_number <- file.path(ctdir, 'total_merged_filtered_states_original.csv')
    if (use_greedy_tree) {
      newick <- file.path(ctdir, 'consensus.newick')
    }
  }
  print('Generating the ordinary one')
  if (just_normal) {
    print("Debug without filtered.csv and average.csv")
    sample_run(newick = newick, copy_number=copy_number, 
               cell_clones=cell_clones, chr_filter=chr_filter,
               corrupt_input = NULL, posterior_mat = NULL)
    
  } else {
    # sample_run(newick = newick, corrupt_input = corrupt_input, posterior_mat = posterior_mat, ...)
    print("Full debug")
    sample_run(newick, copy_number, cell_clones, 
               output = NULL, use_all_loci = TRUE, filter_cells = NULL, 
               corrupt_input, posterior_mat, chr_filter, 
               dev = 'png', g = NULL, drop_loci_names = FALSE) 
    
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
    pmat1 <- pmat
    rownames(pmat) <- pmat$loci
    pmat$loci <- NULL
    pmat1$loci <- NULL
    headm(pmat)
    colnames(pmat) <- gsub('cell_', '', colnames(pmat))
    stopifnot(all(colnames(pmat) %in% colnames(copynumber)))
    colnames(pmat1) <- gsub('cell_', '', colnames(pmat1))
    stopifnot(all(colnames(pmat1) %in% colnames(copynumber)))
    
    # Add +4 to the start value of the loci (just used for sorting)
    options(scipen=999) # prevent scientific notation
    bins <- parse_bin_names(rownames(pmat), as_factors = F) 
    bins$start <- bins$start + 3
    # rownames(pmat) <- collapse_bin_names(bins)
    print("*** DEBUG  ****")
    print(dim(bins))
    print(dim(pmat1))
    rownames(pmat1) <- collapse_bin_names(bins)
    
    # Add to the copynumber and resort
    s1 <- as.data.frame(copynumber)
    
    s1$bin_name <- rownames(s1)
    
    pmat1$bin_name <- rownames(pmat1)
    copynumber1 <- dplyr::bind_rows(s1, pmat1)
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
    fdat1 <- fdat
    rownames(fdat) <- fdat$loci
    fdat$loci <- NULL
    fdat1$loci <- NULL
    headm(fdat)
    
    # Add +2 to the start value of the loci
    options(scipen=999) # prevent scientific notation
    bins <- parse_bin_names(rownames(fdat), as_factors = F) 
    bins$start <- bins$start + 2
    # rownames(fdat) <- collapse_bin_names(bins)
    print("*** DEBUG  ****")
    print(dim(bins))
    print(dim(fdat1))
    rownames(fdat1) <- collapse_bin_names(bins)
    
    # Add to the copynumber and resort
    s1 <- as.data.frame(copynumber)
    print("*** DEBUG s1 ****")
    s1$bin_name <- rownames(s1)
    # fdat$bin_name <- rownames(fdat)
    print("*** DEBUG fdat1 ****")
    fdat1$bin_name <- rownames(fdat1)
    # copynumber1 <- dplyr::bind_rows(s1, fdat)
    print("*** DEBUG concatenate ****")
    print(dim(s1))
    print(dim(fdat1))
    # print(s1[1:3,1:3])
    # print(fdat1[1:3,1:3])
    rownames(s1) <- paste0('cp_',s1$bin_name)
    rownames(fdat1) <- paste0('fdat_',rownames(fdat1))
    
    
    copynumber1 <- dplyr::bind_rows(s1, fdat1)
    # rownames(copynumber1) <- copynumber1$bin_name
    rownames(copynumber1) <- c(rownames(s1),rownames(fdat1))
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
