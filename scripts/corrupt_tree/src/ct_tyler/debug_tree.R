# load debug utils function
project_dir <- '/home/htran/Projects/hakwoo_project/corrupt_tree'
source(paste0(project_dir, "/src/testing_only/blob_v2.R"))

# For Tyler version, there are some left padding bin with negative values
parse_bin_names <- function(bin_names, as_factors = FALSE) {
  # Remove corrupt_tree locus tag if it's there
  bin_names <- gsub('locus_', '', bin_names)
  # ex: bin_names = "7_-1_0"
  # ex: bin_names = 7_150500001_151000000"
  chr <- gsub('([0-9]+|X|Y)_(-1|[0-9]+)_[0-9]+', '\\1', bin_names)
  start <- as.numeric(gsub('([0-9]+|X|Y)_(-1|[0-9]+)_[0-9]+', '\\2', bin_names))
  end <- as.numeric(gsub('([0-9]+|X|Y)_(-1|[0-9]+)_([0-9]+)', '\\3', bin_names))
  # print(paste0(chr,'_',start,'_',end))
  data.frame(chr = chr, start = start, end = end, stringsAsFactors = as_factors)
}

sample_run_main_corrupt_tree_Tyler <- function(ctdir = NULL, chr_filter = NULL, just_normal = FALSE, 
                                         use_greedy_tree = FALSE) {
  if (!is.null(ctdir)) {
    newick <- file.path(ctdir, 'tree.newick')
    corrupt_input <- file.path(ctdir, 'filtered.csv')  # corrupt_filter func, input: straighten_output.csv, output: filtered.csv
    posterior_mat <- file.path(ctdir, 'average.csv')
    cell_clones <- file.path(ctdir, 'cell_clones.csv')
    # copy_number <- file.path(ctdir, 'total_merged_filtered_states.csv')
    copy_number <- file.path(save_dir,'total_merged_filtered_states_padding.csv')
    if (use_greedy_tree) {
      newick <- file.path(ctdir, 'consensus.newick')
    }
  }
  print('Generating the ordinary one')
  if (just_normal) {
    sample_run(newick = newick, corrupt_input = NULL, posterior_mat = NULL, ...)
    
  } else {
    # sample_run(newick = newick, corrupt_input = corrupt_input, posterior_mat = posterior_mat, ...)
    sample_run(newick, copy_number, cell_clones, 
               output = NULL, use_all_loci = TRUE, filter_cells = NULL, 
               corrupt_input, posterior_mat, chr_filter, 
               dev = 'png', g = NULL, drop_loci_names = FALSE) 
    
  }
}



results_dir <- 'your_dir/SA919_Tyler_filtered/'

copynumber <- paste0(results_dir, 'total_merged_filtered_states.csv')
copy_number <- read.csv(copynumber, header=T, row.names=1, check.names = F,stringsAsFactors = FALSE)
loci_use <- rownames(copy_number)

padding <- paste0(results_dir, 'bin_cnvs_corrupt_double_padding.csv')
padding_df <- read.csv(padding, header=T, check.names = F,stringsAsFactors = FALSE)
loci_chrs <- paste0(padding_df$chr,'_',padding_df$start,'_',padding_df$end)
loci_chrs <- loci_chrs[!loci_chrs %in% loci_use]
length(loci_chrs)

# Add new bins to cn state matrix, with value of 2, just to make sure debug function work
padding_added <- as.data.frame(matrix(2, nrow = length(loci_chrs),ncol=ncol(copy_number), 
                                      dimnames = list(loci_chrs, colnames(copy_number))))
copy_number <- rbind(copy_number,padding_added)
write.csv(copy_number, paste0(save_dir,'total_merged_filtered_states_padding.csv'),  row.names = T, quote=F)

sample_run_main_corrupt_tree_Tyler(ctdir = results_dir, chr_filter = '7')


