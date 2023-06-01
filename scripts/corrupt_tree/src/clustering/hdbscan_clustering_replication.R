## Input data
## tree phylogeny tree 
## cell_clones file cell_id, clone_id 
## copy number matrix, chr_start_end in rowname, cell_id in colnames 
suppressPackageStartupMessages({
  library(dplyr)
  library(signals)
})  
## Reference: please cite the signals package if you use these scripts
## https://github.com/shahcompbio/signals
## I add some scripts to run this function but main functions are from signal packgage, Marc William
## Noted: comment main() function in this file
source('/home/htran/Projects/hakwoo_project/corrupt_tree/src/clustering/make_cell_copynumber_tree_heatmap.R')

# "/home/htran/storage/datasets/metastasis_results/dlp_SA919_mixing_exp/graph_cut/A130854B_filtered_cell_clones.csv"
# '/home/htran/storage/raw_DLP/metastasis_DLP/SA919/A130854B/hmmcopy/filtered_reads_repl_timing.csv.gz'
get_filtered_cells <- function(download_dir, save_dir, library_id, reads_df,
                               cellclones_fn=NULL, metric_fn=NULL){
  if(is.null(metric_fn)){
    metric_fn <- paste0(download_dir,library_id,'/annotation/','metrics.csv.gz')
  }
  
  if(is.null(cellclones_fn)){
    cellclones_fn <- paste0(save_dir,library_id,'_cell_clones.csv')
  }  
  cell_clones <- data.table::fread(cellclones_fn) 
  
  metrics_df <- data.table::fread(metric_fn) 
  print(summary(as.factor(cell_clones$clone_id)))  
  noisy_cells <- cell_clones %>%
    dplyr::filter(clone_id=='A') %>%
    dplyr::pull(cell_id)
  
  cell_clones <- cell_clones %>%
    dplyr::filter(!clone_id %in% c('A')) %>%
    dplyr::mutate(clone_id=
      case_when(
        clone_id %in% c('0','B') ~ 'B',
        TRUE ~ 'C')
    )
  dim(cell_clones)    
  print(summary(as.factor(cell_clones$clone_id)))    
  
  # colnames(metrics_df)[grepl('*gc*', colnames(metrics_df))]
  
  metrics_df <- metrics_df %>%
    dplyr::filter(cell_id %in% cell_clones$cell_id) %>%
    dplyr::mutate(cell_type_status=
           case_when(quality>0.75 & is_s_phase==F ~ 'cn_g',
                     TRUE ~ 'cn_s')) %>%
    dplyr::select(cell_id, cell_type_status)
  # t <- metrics_df %>%
  #   dplyr::filter(cell_id %in% cell_clones$cell_id &
  #                 quality>0.75)
  # 
  dim(metrics_df) 
  
  cell_clones <- cell_clones %>%
    inner_join(metrics_df, by='cell_id')
  print(sum(metrics_df$quality>0.75))    
  print(summary(as.factor(metrics_df$is_s_phase)))    
  # print(summary(as.factor(t$is_s_phase)))    
  table(metrics_df$is_s_phase, metrics_df$quality>0.75)
  dim(cell_clones)
  data.table::fwrite(cell_clones, paste0(save_dir,library_id,'_filtered_cell_clones.csv'))
  
  
  # if(is.null(reads_df)){
  #   reads_fn <- paste0(download_dir,library_id,'/hmmcopy/','reads.csv.gz')
  #   # reads_fn <- paste0(download_dir,library_id,'/hmmcopy/',library_id,'_reads.csv.gz')
  # }
  f <- c('cell_id','chr', 'start', 'end', 'gc','state','reads')
  features_use <- colnames(reads_df)[colnames(reads_df) %in% f]
  print(features_use)
  # gc_col : gc content of a bin gc
  # cn_col : copy number state of a bin (via hmmcopy state or some other caller)state
  # clone_col : phylogenetic clone ID for each cell (via sitka or some other copy number clustering)
  # input_col : read depth of each bin (both reads per million or raw read count work)
  reads_df <- reads_df %>% 
    dplyr::select(all_of(features_use)) %>% 
    dplyr::filter(cell_id %in% cell_clones$cell_id)
  dim(reads_df)  
  data.table::fwrite(reads_df, paste0('/home/htran/storage/raw_DLP/metastasis_DLP/SA919/A130854B/hmmcopy/filtered_reads_repl_timing.csv.gz'))
  # temp_cn_s = cn_s[['cell_id', 'chr', 'start', 'end', argv.gc_col, argv.cn_col, 'library_id', argv.input_col]]
  # temp_cn_g = cn_g[['cell_id', 'chr', 'start', 'end', argv.gc_col, argv.cn_col, 'library_id', argv.clone_col, argv.input_col]]
  
  # print('creating scrt object')
  # # create SPF object with input
  # nb_iterations = 100 #1500
  # scrt = scRT(temp_cn_s, temp_cn_g, input_col=argv.input_col, 
  #             rt_prior_col=None, assign_col=argv.cn_col,
  #             cn_state_col=argv.cn_col, gc_col=argv.gc_col, 
  #             cn_prior_method=argv.cn_prior_method, max_iter=nb_iterations)
  # 
  # # run inference
  # print('running inference')
  # cn_s_with_scrt, supp_s_output, cn_g_with_scrt, supp_g_output = scrt.infer_pyro_model()
  
}

## pctcells: percentage of minimum cells to consider as one clone
hdbscran_viz <- function(download_dir, results_dir, library_id, 
                         grouping_file, copynumber_fn=NULL, reads_fn=NULL, pctcells = 0.05){
  if(is.null(copynumber_fn)){
    copynumber_fn <- paste0(results_dir,library_id,'_filtered_states.csv.gz')
  }
  if(is.null(reads_fn)){
    reads_fn <- paste0(download_dir,library_id,'/hmmcopy/','reads.csv.gz')
    # reads_fn <- paste0(download_dir,library_id,'/hmmcopy/',library_id,'_reads.csv.gz')
  }
  # pctcells = 0.05
  copynumber <- data.table::fread(copynumber_fn) %>% as.data.frame()
  rownames(copynumber) <- copynumber$V1
  copynumber$V1 <- NULL
  print(dim(copynumber))
  
  reads_df <- data.table::fread(reads_fn) %>% as.data.frame()
  print(dim(reads_df))
  # reads_df$copy[1:10]
  
  
  reads_df <- reads_df %>%
    dplyr::filter(cell_id %in% colnames(copynumber))
  # length(unique(reads_df$cell_id))
  reads_df$chr_desc <- paste0(reads_df$chr,'_',reads_df$start,'_',reads_df$end)
  reads_df <- reads_df %>%
    dplyr::filter(chr_desc %in% rownames(copynumber))
  
  print(dim(reads_df))
  ncells <- length(unique(reads_df$cell_id))
  ## Using field = "copy", copy number values as the input for cell clustering
  clusters <- signals::umap_clustering(reads_df,
                                        minPts = max(round(pctcells * ncells), 30),
                                        field = "copy")
  
  # clusters <- signals::umap_clustering(reads_df, 
  #                                      minPts = min(round(pctcells * ncells), 10), 
  #                                      field = "copy")
  
  tree <- clusters$tree
  clones <- clusters$clustering
  unique(clones$clone_id)
  
  
  ## Save output of clones here
  
  # ## For A98166A 
  # clones$clone_id <- ifelse(clones$clone_id=='0','C',clones$clone_id)
  # clones$clone_id <- ifelse(clones$clone_id=='C','A',
  #                           ifelse(clones$clone_id %in% c('NA'),'C','NA'))
  # class(clones)
  # data.table::fwrite(clones, paste0(save_dir, library_id,'_cell_clones.csv'))
  # saveRDS(clones, paste0(save_dir, library_id,'_cell_clones.rds'))
  # ## For A98166B
  # clones$clone_id <- ifelse(clones$clone_id %in% c('A','0'),'C',
  #                           ifelse(clones$clone_id=='B','A','D'))
  # summary(as.factor(clones$clone_id))
  save_dir <- paste0(results_dir,'graph_cut/')
  if(!dir.exists(save_dir)){
    dir.create(save_dir)
  }
  data.table::fwrite(clones, paste0(save_dir,library_id,'_cell_clones.csv'))
  output_fn <- paste0(save_dir, library_id,'_hdbscan_cell_cn_tree_heatmap.png')
  png(output_fn, height = 2*700, width=2*1200, res = 2*72)
  make_cell_copynumber_tree_heatmap(
    tree, copynumber, clones, NULL, grouping_file
  )
  dev.off()
  
  
}

# Mixing experiment with primary cells library A98166A
# experiments: 1:0.25, expected: A: 1, B: 0.25
# Output: A: 250 cells, C: 13 cells. 


# results_dir <- '/home/htran/storage/gm_instability_results/HEK293/data_thrs_10000/'
# download_dir <- '/home/htran/storage/datasets/damian_DLP/'
# grouping_file <- '/home/htran/storage/gm_instability_results/HEK293/library_groupings.csv'
# libs_id <- c('A110616A','A110617A','A118346B','A118825A','A118866A')
# for(library_id in libs_id){
#   hdbscran_viz(results_dir, library_id, grouping_file)
# }

# results_dir <- '/home/htran/storage/datasets/metastasis_results/dlp_SA919_mixing_exp/'
# download_dir <- '/home/htran/storage/raw_DLP/metastasis_DLP/SA919/'
# grouping_file <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/dlp_trees/SA919_mixing_experiment/library_groupings.csv'
# libs_id <- c('A98166A','A98166B')
# libs_id <- c('A130841A','A130841B','A130839B')
# library_id <- 'A130839B'
# for(library_id in libs_id){
#   hdbscran_viz(download_dir, results_dir, library_id, grouping_file)
# }

# results_dir <- '/home/htran/storage/datasets/metastasis_results/SA1142/first_results/'
# download_dir <- '/home/htran/storage/raw_DLP/metastasis_DLP/SA1142/'
# grouping_file <- '/home/htran/storage/datasets/metastasis_results/SA1142/first_results/library_groupings.csv.gz'
# lib_df <- data.table::fread(grouping_file)
# 
# libs_id <- lib_df$grouping
# library_id <- lib_df$grouping[1]
# 
# for(library_id in lib_df$grouping[2:dim(lib_df)[1]]){
#   hdbscran_viz(download_dir, results_dir, library_id, grouping_file)
# }


results_dir <- '/home/htran/storage/datasets/metastasis_results/dlp_SA919_mixing_exp/'
download_dir <- '/home/htran/storage/raw_DLP/metastasis_DLP/SA919/'
grouping_file <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/dlp_trees/SA919_mixing_experiment/library_groupings.csv'
library_id <- 'A130854B'
for(library_id in libs_id){
  hdbscran_viz(download_dir, results_dir, library_id, grouping_file)
}

paste0(results_dir,'graph_cut/',library_id,'_filtered_cell_clones.csv')

input_dir <- '/home/htran/storage/raw_DLP/metastasis_DLP/SA919/A130854B/'
reads_fn <- paste0(input_dir,'hmmcopy/filtered_reads_repl_timing.csv.gz')
reads_df <- data.table::fread(reads_fn) %>% as.data.frame()


metrics_df <- data.table::fread(paste0(input_dir,library_id,'_filtered_cell_clones.csv'))
dim(metrics_df)
summary(as.factor(metrics_df$clone_id))
summary(as.factor(metrics_df$cell_type_status))
metrics_df$cell_id
cn_g_cells <- metrics_df %>% 
  dplyr::filter(cell_type_status=='cn_g') %>% 
  dplyr::pull(cell_id)

cn_s_cells <- metrics_df %>% 
  dplyr::filter(cell_type_status=='cn_s') %>% 
  dplyr::pull(cell_id)

reads_df_g <- reads_df %>% 
  dplyr::filter(cell_id %in% cn_g_cells)
reads_df_s <- reads_df %>% 
  dplyr::filter(cell_id %in% cn_s_cells)

length(unique(reads_df_g$cell_id))
length(unique(reads_df_s$cell_id))
length(cn_g_cells)
length(cn_s_cells)
data.table::fwrite(reads_df_g, paste0(input_dir,'hmmcopy/filtered_reads_RT_g_cells.csv.gz'))
data.table::fwrite(reads_df_s, paste0(input_dir,'hmmcopy/filtered_reads_RT_s_cells.csv.gz'))


