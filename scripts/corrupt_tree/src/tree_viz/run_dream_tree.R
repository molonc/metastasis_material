# Farhia, SA1035
# results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_new_encode/SA1035_whole_local/'
# cellclones <- paste0(results_dir, 'tree_cut_out/cell_clones_if_0.1_af_0.65_p0.75_e0.04.csv')
# cellclones <- paste0(results_dir, 'tree_cut_out/cell_clones_if_0.02_af_0.45_p0.75_e0.04.csv')
# datatag = 'SA1035_13clones'


# Farhia, SA1035
# results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_new_encode/SA1035_filtered/'
# results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_new_encode/SA1035_filtered_v2/'
# results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_new_encode/SA1035_filtered_FP001/'
# cellclones <- paste0(results_dir, 'cell_clones.csv')
# datatag = 'SA1035_Salehi'


# results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_new_encode/SA1035_filtered_FP001/'
# datatag = 'SA1035_Salehi'
# cellclones <- paste0(results_dir, 'graph_cut/cell_clones_0.01.csv')

# results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_new_encode/SA1035_Tyler/'
# cellclones <- paste0(results_dir, 'cell_clones.csv')
# datatag = 'SA1035_Tyler'

# results_dir <- '/home/htran/storage/datasets/metastasis_results/SA535_new_encoding/Tyler_whole/'
# datatag = 'SA535_Tyler_v2'
# cellclones <- paste0(results_dir, 'graph_cut/cell_clones_0.01.csv')

# results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_filtered/'
# datatag = 'SA919_Tyler'
# cellclones <- paste0(results_dir, 'cell_clones.csv')


# base_dir <- '/home/htran/storage/datasets/metastasis_results/SA535_new_encoding/'
# sub_dir <- c("Tyler_whole")
# sub_dir <- c("Tyler_2361","Tyler_2362","Tyler_2363","Tyler_2364","Tyler_whole")

# base_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/'
# sub_dir <- c("Tyler_2251","Tyler_2252","Tyler_2253","Tyler_2254","Tyler_216","Tyler_2164","SA919_Tyler_filtered")
# # s <- "Tyler_2251"
# # source_dir <- '/home/htran/Projects/hakwoo_project/corrupt_tree/src/tree_viz/'
# # source(paste0(source_dir,'trajectory_utils.R'))
# # 
# for(s in sub_dir){
#   s = sub_dir[7]
#   subthres=100
#   results_dir <- paste0(base_dir,s,"/")
# 
#   cellclones <- paste0(results_dir, 'cell_clones.csv')
#   print(s)
#   # datatag <- gsub('Tyler_','SA535_',s)
#   if(s=="SA919_Tyler_filtered"){
#     datatag <- "SA919_total"
#   }
#   # datatag <-'SA919'
#   save_dir <- paste0(results_dir,'dream_tree/')
#   if (!file.exists(save_dir)){
#     dir.create(save_dir)
#   }
#   generate_tree(results_dir, save_dir, datatag, subthres)
# 
# }

generate_tree <- function(results_dir, save_dir, datatag = 'SA1035',subthres=20) {
  
  
  tree_info <- get_tree_info(results_dir)
  plot_trees(cellclones, tree_info, results_dir, save_dir, subthres,
             rm_threshold=-1, datatag, save_data=T)
  
  plot_metadata(save_dir, datatag)
  
}

# generate_tree(results_dir, datatag)
# meta_data_ls <- readRDS(paste0(save_dir,datatag,'_meta_data_ls.rds'))
# plot_metadata(meta_data_ls, save_dir, datatag)






