  project_dir <- '/home/htran/Projects/hakwoo_project/corrupt_tree'
  source(paste0(project_dir, "/src/testing_only/blob_v2.R"))
  # results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_new_encode/SA1035_whole_local/'
  results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local_v5/'  
  #
  cnv_path <- paste0(results_dir, 'total_merged_filtered_states.csv')
  # cell_clones <- paste0(results_dir, 'tree_cut_out/cell_clones_if_0.02_af_0.75_p0.75_e0.02.csv')
  cell_clones <- paste0(results_dir, 'cell_clones.csv')
  
  newick <- file.path(results_dir, 'tree.newick')
  tree <- read.tree(newick)
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
  
  
 
  
  
  sample_run_main_corrupt_tree(ctdir = results_dir, chr_filter = c(7,8))
  sample_run_main_corrupt_tree(ctdir = results_dir)
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
  sample_run_main_corrupt_tree <- function(ctdir = NULL, chr_filter = NULL, just_normal = FALSE, 
                              use_greedy_tree = FALSE) {
    if (!is.null(ctdir)) {
      newick <- file.path(ctdir, 'tree.newick')
      corrupt_input <- file.path(ctdir, 'filtered.csv')  # corrupt_filter func, input: straighten_output.csv, output: filtered.csv
      posterior_mat <- file.path(ctdir, 'average.csv')
      cell_clones <- file.path(ctdir, 'cell_clones.csv')
      copy_number <- file.path(ctdir, 'total_merged_filtered_states.csv')
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
