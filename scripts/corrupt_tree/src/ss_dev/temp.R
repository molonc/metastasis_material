
# For tyler
#
function() {
  qq = fasttable('/Users/sohrabsalehi/projects/fitness/breaktree/data/SA1051_corrupttree_lf0.001/filtered_cell_cn.tsv')
  tt <- qq[, c(1,2,3)]
  qq <- qq[, -c(1,2,3,4)]
  rownames(qq) <- collapse_bin_names(tt)  
  qq$V1 <- rownames(qq)
  #data.table::fwrite(qq, '/Users/sohrabsalehi/projects/fitness/breaktree/data/SA1051_corrupttree_lf0.001/total_merged_filtered_states.csv')
  write.table(x = qq, file = '/Users/sohrabsalehi/projects/fitness/breaktree/data/SA1051_corrupttree_lf0.001/total_merged_filtered_states.csv', sep = ',', row.names = T)
  
  ss = fasttable('/Users/sohrabsalehi/projects/fitness/breaktree/data/SA1051_corrupttree_lf0.001/total_merged_filtered_states.csv')
  headm(ss)
  ss$V1
  
}
#sample_run_fast(ctdir = '/Users/sohrabsalehi/projects/fitness/breaktree/data/SA1051_corrupttree_lf0.001/corrupt_tree_output/',
#sample_run_fast(ctdir = '/Users/sohrabsalehi/projects/corrupt-nextflow/deliverables/no_pt_fast_default_error_global_cell_filter_full_padding_tyler',
#sample_run_fast(ctdir = '/Users/sohrabsalehi/projects/corrupt-nextflow/deliverables/no_pt_fast_default_error_local_cell_filter_full_padding_tyler/',

#sample_run_fast(ctdir = '/Users/sohrabsalehi/projects/corrupt-nextflow/deliverables/no_pt_fast_default_error_global_cell_filter_full_padding/',
function() {
  sample_run_fast(ctdir = '/Users/sohrabsalehi/projects/corrupt-nextflow/deliverables/no_pt_fast_default_error_global_cell_filter_full_padding_tyler_more_loci/',
                  #copy_number = '/Users/sohrabsalehi/projects/fitness/breaktree/data/SA1051_corrupttree_lf0.001/filtered_cell_cn.tsv',
                  copy_number = '/Users/sohrabsalehi/projects/fitness/breaktree/data/SA1051_corrupttree_lf0.001/total_merged_filtered_states.csv',
                  cell_clones = NULL,
                  use_all_loci = TRUE,
                  dev = 'png',
                  #chr_filter = c(20,21,22,"X"))
                  #chr_filter = c(1))
                  chr_filter = NULL,
                  just_normal = TRUE, 
                  drop_loci_names = TRUE)
  
}

function() {
  
  
  sample_run(newick = '/Users/sohrabsalehi/Downloads/corrupttree_v3_Hoa/tree.newick', 
             copy_number = '/Users/sohrabsalehi/Downloads/corrupttree_v3_Hoa/total_merged_filtered_states.csv',
             cell_clones = NULL,
             use_all_loci = TRUE,
             corrupt_input = '/Users/sohrabsalehi/Downloads/corrupttree_v3_Hoa/filtered.csv',
             posterior_mat = '/Users/sohrabsalehi/Downloads/corrupttree_v3_Hoa/average.csv',
             dev = 'png',
             #chr_filter = c(20,21,22,"X"))
             #chr_filter = c(1))
             chr_filter = NULL)
  
  
  
}

## For Hoi
#sample_run_fast(ctdir = '~/projects/corrupt-nextflow/deliverables/no_pt_fast_work_v3_Hoa/', 
#sample_run_fast(ctdir = '~/projects/corrupt-nextflow/deliverables/no_pt_fast_SA919_Hoa/', 
#sample_run_fast(ctdir = '~/projects/corrupt-nextflow/deliverables/no_pt_fast_SA919_Hoa_long/no_pt_fast/', 
#sample_run_fast(ctdir = '~/projects/corrupt-nextflow/deliverables/no_pt_fast_SA919_Hoa_long/no_pt_fast/', 
#sample_run_fast(ctdir = '/Users/sohrabsalehi/projects/corrupt-nextflow/deliverables/no_pt_fast_less_error', 
#sample_run_fast(ctdir = '/Users/sohrabsalehi/projects/corrupt-nextflow/deliverables/no_pt_fast_less_error_local/', 
#sample_run_fast(ctdir = '/Users/sohrabsalehi/projects/corrupt-nextflow/deliverables/no_pt_fast_less_error_local_cell_filter/', 
#sample_run_fast(ctdir = '/Users/sohrabsalehi/projects/corrupt-nextflow/deliverables/no_pt_fast_default_error_global_cell_filter/',
#sample_run_fast(ctdir = '/Users/sohrabsalehi/projects/corrupt-nextflow/deliverables/no_pt_fast_default_error_global_cell_filter_full_padding',
# sample_run_fast(ctdir = '/Users/sohrabsalehi/projects/corrupt-nextflow/deliverables/no_pt_fast_default_error_global_cell_filter_full_padding/',
#              #copy_number = '~/Downloads/corrupttree_v3_Hoa/total_merged_filtered_states.csv',
#              #copy_number = '~/Downloads/output_SA919_Hoa/total_merged_filtered_states.csv',
#              copy_number = '~/Downloads/SA919_without_noise/total_merged_filtered_states.csv',
#              cell_clones = NULL,
#              use_all_loci = TRUE,
#              dev = 'png',
#              #chr_filter = c(20,21,22,"X"))
#              #chr_filter = c(1))
#              chr_filter = NULL,
#              just_normal = TRUE)




function() {
  
  #ss = fasttable('/Users/sohrabsalehi/projects/corrupt-nextflow/deliverables/no_pt_fast_work_v3_Hoa/filtered.csv')
  ss = fasttable('/Users/sohrabsalehi/Downloads/output_SA919_Hoa/results/bin_cnvs_corrupt.csv')
  ss = fasttable('/Users/sohrabsalehi/Downloads/SA919_without_noise/results/bin_cnvs_corrupt_double_padding.csv')
  ss = fasttable('/Users/sohrabsalehi/projects/fitness/breaktree/data/SA1051_corrupttree_lf0.001/results/bin_cnvs_corrupt_double_padding.csv')
  ss = fasttable('/Users/sohrabsalehi/projects/corrupt-nextflow/deliverables/no_pt_fast_default_error_global_cell_filter_full_padding_tyler_more_loci/filtered.csv')
  
  #ss = fasttable('/Users/sohrabsalehi/projects/corrupt-nextflow/deliverables/no_pt_fast/filtered.csv')
  (loci <- unique(ss$loci))
  ll <- loci[get_pad_bin_index(loci)]
  lchr <- gsub('([0-9]+|X|Y)_[0-9]+_[0-9]+', '\\1', gsub('locus_', '', ll)); table(lchr)
  ll[lchr %in% c(10)]
  ll[lchr %in% c(5)]
  
  
  
  tree <- read.tree('/Users/sohrabsalehi/projects/corrupt-nextflow/deliverables/no_pt_fast_work_v3_Hoa/tree.newick')
  g <- read_ltm_tree(tree_2_edge_list(tree))
  g.get.loci(g)
  
  annotate_nodes_on_tree(edge_list_path = NULL, node_names = loci[get_pad_bin_index(loci)], tree_node_dic = graph_to_tree_dic(g), show_label = F)
}