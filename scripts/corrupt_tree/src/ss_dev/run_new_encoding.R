source_dir <- '/home/htran/Projects/hakwoo_project/corrupt_tree/src/ss_dev/'
source(paste0(source_dir, 'blob.R'))
source(paste0(source_dir, 'general_utils.R'))
source(paste0(source_dir, 'tree_and_heatmap_utils_utils.R'))
source(paste0(source_dir, 'tree_and_heatmap_utils.R'))
source(paste0(source_dir, 'tree_and_heatmap.R'))


run_plot_heatmap <- function(){
  # ## Launch
  # option_list <- list(make_option(c("-n", "--newick"), type="character", default=NULL, help="input_newick", metavar="character"),
  #                     make_option(c("-t", "--datatag"), type="character", default=NULL, help="library_id", metavar="character")          
  # )
  # 
  # opt_parser <- OptionParser(option_list=option_list)
  # opt <- parse_args(opt_parser)
  # edge_list_path <- newick_2_edge_list(opt$newick)
  # g <- read_ltm_tree(edge_list_path)
  # fitclone_plot_tree_heatmap(g, opt$datatag, edge_list)
  
  data_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7/X0847_whole_065MF_20percent/'
  newick_path <- paste0(data_dir,'test_new_encoding/tree.newick')
  datatag <- 'SA919X7'
  edge_list_path <- newick_2_edge_list(newick_path)
  g <- read_ltm_tree(edge_list_path)
  exp_path <- paste0(data_dir,'test_new_encoding/')
  fitclone_plot_tree_heatmap(g, datatag, edge_list_path, out_dir = paste0(data_dir,'test_new_encoding/'))
  
}


test_run <- function(data_dir, prob = .80) {
  # TODO: compute average ploidy too for each cell
  
  cnv_path <- paste0(data_dir,'total_merged_filtered_states_original.csv')
  mat <- read.delim(cnv_path, check.names = FALSE, stringsAsFactors = FALSE, sep = ",")
  # mat <- fasttable(cnv_path)
  # rownames(mat) <- mat$V1; mat$V1 <- NULL
  # colnames(mat)[1:2]
  # rownames(mat)[1:2]
  # dim(mat)
  jump_rank <- compute_jump_cells(mat)
  # 
  # # Filter cells by jump/stuff
  
  #bad_cells <- jump_rank %>% dplyr::filter( njumps > quantile(x=njumps, probs = prob) | avgCNA > quantile(x=avgCNA, probs = prob)) %>% dplyr::select(cell_id)
  bad_cells <- jump_rank %>% dplyr::filter( njumps > quantile(x=njumps, probs = prob)) %>% dplyr::select(cell_id)
  bad_cells <- bad_cells$cell_id
  print("Number of bad cells: ")
  print(length(bad_cells))
  
  # uncomment for enhanced viz
  # sample_run_fast(ctdir = data_dir,
  #                 copy_number = cnv_path,
  #                 cell_clones = NULL,
  #                 use_all_loci = TRUE,
  #                 #filter_cells = bad_cells,
  #                 filter_cells = setdiff(jump_rank$cell_id, bad_cells),
  #                 dev = 'png',
  #                 #chr_filter = c(1),
  #                 chr_filter = NULL,
  #                 just_normal = TRUE,
  #                 drop_loci_names = TRUE)

  
  # Filter cells
  cc <- mat[, colnames(mat) %ni% bad_cells]
  print(paste0("Filtered output: ", dim(cc)[1]," ",dim(cc)[2]))
  write.csv(cc, paste0(data_dir,'total_merged_filtered_states.csv'),  row.names = T)
  pad_mod_cnv_mat(cc, data_dir, pad_left = TRUE)
  
  # The results would be in:
  # ./data/results/bin_cnvs_corrupt_double_padding.csv
  
  # Viz
  # sample_run_fast(ctdir = '??',
  #                 copy_number = cnv_path,
  #                 cell_clones = NULL,
  #                 use_all_loci = TRUE,
  #                 dev = 'png',
  #                 #chr_filter = c(1),
  #                 chr_filter = NULL,
  #                 just_normal = TRUE, 
  #                 drop_loci_names = TRUE)
  
}



cn2binary_newencoding <- function(cnv_path, output_file, prob = .80) {
  # TODO: compute average ploidy too for each cell
  
  # cnv_path <- paste0(data_dir,'total_merged_filtered_states_original.csv')
  mat <- read.delim(cnv_path, check.names = FALSE, stringsAsFactors = FALSE, sep = ",")
  jump_rank <- compute_jump_cells(mat)
  
  
  # # Filter cells by jump/stuff
  #bad_cells <- jump_rank %>% dplyr::filter( njumps > quantile(x=njumps, probs = prob) | avgCNA > quantile(x=avgCNA, probs = prob)) %>% dplyr::select(cell_id)
  bad_cells <- jump_rank %>% dplyr::filter( njumps > quantile(x=njumps, probs = prob)) %>% dplyr::select(cell_id)
  bad_cells <- bad_cells$cell_id
  print("Number of bad cells: ")
  print(length(bad_cells))
  cc <- mat[, colnames(mat) %ni% bad_cells]
  print(paste0("Filtered output: ", dim(cc)[1]," ",dim(cc)[2]))
  write.csv(cc, paste0(data_dir,'total_merged_filtered_states.csv'),  row.names = T)
  pad_mod_cnv_mat(cc, output_file, pad_left = TRUE)
  
  # The results would be in:
  # bin_cnvs_corrupt_double_padding.csv
  
}

data_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035/SA1035_20percent_v3/'
cnv_path <- paste0(data_dir,'total_merged_filtered_states.csv')
cnv_path_1 <- paste0(data_dir,'total_merged_filtered_states_original.csv')

mat <- read.delim(cnv_path_1, check.names = FALSE, stringsAsFactors = FALSE, sep = ",")
t <- read.csv(cnv_path, header=T, row.names=1, check.names = F, stringsAsFactors = FALSE)
dim(t)
dim(mat)
mat[1:2,1:2]
t[1:2,1:2]
cn2binary_newencoding(data_dir)
