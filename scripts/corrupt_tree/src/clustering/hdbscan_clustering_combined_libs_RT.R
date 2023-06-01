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
## Hoa Tran: I add some scripts to run this function but main functions are from signal packgage, Marc William
## Hoa Tran: I add some visualization scripts here
source('/home/htran/Projects/hakwoo_project/corrupt_tree/src/clustering/make_cell_copynumber_tree_heatmap.R')

load_hmmcopy_reads <- function(obs_libs, download_dir, datatag='SA',save_data=F, cells_id=NULL, reads_fn=NULL){
  print('Loading copy number data from library ids: ')
  print(paste(obs_libs, collapse = ' '))
  obs_reads <- list()
  
  for(l in obs_libs){
    lib_fn1 <- paste0(download_dir,l,'/hmmcopy/','reads.csv.gz') # csv.gz or csv
    # lib_fn <- paste0(download_dir,l,'/',l,'/hmmcopy/','reads.csv.gz') # csv.gz or csv
    lib_fn2 <- paste0(download_dir,l,'/hmmcopy/',l,'_reads.csv.gz') # csv.gz or csv
    lib_fn <- ''
    if(file.exists(lib_fn1)){
      lib_fn <- lib_fn1
    } else{
      if(file.exists(lib_fn2)){
        lib_fn <- lib_fn2
      } else{
        print(paste0('***ERROR: library: ', l))
        stop('Do not exist hmmcopy reads for library, double check input data')
      } 
    } 
    if(lib_fn!=''){
      tmp <- data.table::fread(lib_fn)
      # View(head(tmp))
      # dim(tmp)
      if(!is.null(cells_id)){
        tmp <- tmp %>%
          dplyr::filter(cell_id %in% cells_id)  
      }
      
      tmp <- tmp %>%
        dplyr::select(cell_id, chr,start,end,gc,state,reads,copy) # 
      print(l)
      print(dim(tmp))
      if(dim(tmp)[1]>0){
        obs_reads[[l]] <- tmp  
      }else{
        print(paste0('***Warning:  Do not have any output for clone: ',c))
      }
      
    }
    
  }
  reads_df <- as.data.frame(dplyr::bind_rows(obs_reads))
  features_use <- c('cell_id','chr', 'start', 'end', 'gc','state','reads','copy')
  if(sum(colnames(reads_df) %in% features_use) !=length(features_use)){
    stop("Reads mtx do not contain used features, pls double check input reads mtx")
  }
  print(features_use)
  # gc_col : gc content of a bin gc
  # cn_col : copy number state of a bin (via hmmcopy state or some other caller)state
  # clone_col : phylogenetic clone ID for each cell (via sitka or some other copy number clustering)
  # input_col : read depth of each bin (both reads per million or raw read count work)
  reads_df <- reads_df %>% 
    dplyr::select(all_of(features_use)) #%>% 
  # dplyr::filter(cell_id %in% cell_clones$cell_id)
  print(dim(reads_df))  
  
  if(save_data){
    if(is.null(reads_fn)){
      reads_fn <- paste0(download_dir, datatag,'_total_reads.csv.gz')
    }
    data.table::fwrite(reads_df, reads_fn)  
  }
  
  print(dim(reads_df))
  return(reads_df)
}


load_filtered_data <- function(obs_libs, input_dir, datatag, filtered_fn=NULL){
  print('Loading filtered copy number states from library ids: ')
  print(paste(obs_libs, collapse = ' '))
  obs_reads <- list()

  for(l in obs_libs){
    lib_fn <- paste0(input_dir,l,'_filtered_states.csv.gz')
    
    if(file.exists(lib_fn)){
      tmp <- data.table::fread(lib_fn)
      rownames(tmp) <- tmp$V1
      tmp$V1 <- NULL
      
      print(l)
      print(dim(tmp))
      if(dim(tmp)[1]>0){
        obs_reads[[l]] <- tmp  
      }else{
        print(paste0('***Warning:  Do not have any output for clone: ',c))
      }
      
    }else{
      print(paste0('***ERROR: library: ', l))
      # print('Do not exist hmmcopy reads for library, double check input data')
      stop('Do not exist hmmcopy reads for library, double check input data')
    }
    
  }
  filtered_df <- as.data.frame(dplyr::bind_cols(obs_reads))
  print(dim(filtered_df))
  rownames(filtered_df) <- rownames(tmp)
  if(is.null(filtered_fn)){
    filtered_fn <- paste0(input_dir, datatag,'_total_merged_filtered_states.csv.gz')
  }
  data.table::fwrite(filtered_df, filtered_fn, row.names = TRUE)
  
  return(filtered_df)
}

load_metrics_data <- function(obs_libs, download_dir, results_dir, datatag='', cells_id=NULL, metrics_fn=NULL){
  print('Loading copy number hmmcopy metrics from library ids: ')
  print(paste(obs_libs, collapse = ' '))
  obs_reads <- list()
  
  for(l in obs_libs){
    
    lib_fn1 <- paste0(download_dir,l,'/annotation/','metrics.csv.gz') # csv.gz or csv
    # lib_fn <- paste0(download_dir,l,'/',l,'/hmmcopy/','reads.csv.gz') # csv.gz or csv
    lib_fn2 <- paste0(download_dir,l,'/annotation/',l,'_metrics.csv.gz') # csv.gz or csv
    lib_fn <- ''
    if(file.exists(lib_fn1)){
      lib_fn <- lib_fn1
    } else{
      if(file.exists(lib_fn2)){
        lib_fn <- lib_fn2
      } else{
        print(paste0('***ERROR: library: ', l))
        stop('Do not exist annotation metrics for library, double check input data')
      } 
    } 
    if(lib_fn!=''){
      tmp <- data.table::fread(lib_fn)
      print(l)
      print(dim(tmp))
      if(dim(tmp)[1]>0){
        obs_reads[[l]] <- tmp  
      }else{
        print(paste0('***Warning:  Do not have any output for clone: ',c))
      }
      
    }
    
  }
  metrics_df <- as.data.frame(dplyr::bind_rows(obs_reads))
  print(dim(metrics_df))
  obs_metrics <- c('cell_id','experimental_condition','is_s_phase', 'quality')
  obs_metrics <- colnames(metrics_df)[colnames(metrics_df) %in% obs_metrics]
  # gc_col : gc content of a bin gc
  # cn_col : copy number state of a bin (via hmmcopy state or some other caller)state
  # clone_col : phylogenetic clone ID for each cell (via sitka or some other copy number clustering)
  # input_col : read depth of each bin (both reads per million or raw read count work)
  metrics_df <- metrics_df %>% 
    dplyr::select(all_of(obs_metrics))
  
  if(!is.null(cells_id)){
    metrics_df <- metrics_df %>%
      dplyr::filter(cell_id %in% cells_id)  
  }
  if(is.null(metrics_fn)){
    metrics_fn <- paste0(results_dir, datatag,'_total_cells_metrics.csv.gz')
  }
  data.table::fwrite(metrics_df, metrics_fn)
  
  return(metrics_df)
}

## pctcells: percentage of minimum cells to consider as one clone
hdbscran_viz <- function(download_dir, results_dir, obs_libs, grouping_file,
                         reads_df, copynumber=NULL,
                         copynumber_fn=NULL, pctcells = 0.05){
  # pctcells = 0.05
  ## for single library
  # copynumber <- data.table::fread(paste0(results_dir,library_id,'_filtered_states.csv.gz'))%>% as.data.frame()
  
  if(is.null(copynumber)){
    if(is.null(copynumber_fn)){
      # copynumber_fn <- paste0(results_dir,library_id,'_filtered_states.csv.gz')
      copynumber_fn <- paste0(results_dir, 'total_merged_filtered_states.csv.gz')
    }
    ## for combined libraries
    copynumber <- data.table::fread(copynumber_fn) %>% as.data.frame()
    rownames(copynumber) <- copynumber$V1
    copynumber$V1 <- NULL
    print(dim(copynumber))
    print(rownames(copynumber)[1:5])  
  }
  # copynumber <- filtered_df
  ## for single library
  # reads_df <- data.table::fread(paste0(download_dir,library_id,'/hmmcopy/','reads.csv.gz')) %>% as.data.frame()
  # reads_df <- data.table::fread(paste0(download_dir,library_id,'/hmmcopy/',library_id,'_reads.csv.gz')) %>% as.data.frame()
  
  ## for combined libraries
  cells_id <- colnames(copynumber)
  cells_id[1]
  # cells_id <- metrics_df$cell_id
  # length(cells_id)
  # reads_df <- load_hmmcopy_reads(obs_libs, download_dir, cells_id)
  print(dim(reads_df))
  
  # reads_df <- load_hmmcopy_reads(obs_libs, download_dir, cells_id=NULL)
  # print(dim(reads_df))
  # reads_df$copy[1:10]
  
  # reads_df <- reads_df %>%
  #   dplyr::filter(cell_id %in% cells_id)
  
  reads_df$chr_desc <- paste0(reads_df$chr,'_',reads_df$start,'_',reads_df$end)
  reads_df <- reads_df %>%
    dplyr::filter(chr_desc %in% rownames(copynumber)) %>%
    dplyr::select(-chr_desc)
  
  print(dim(reads_df))
  ncells <- length(unique(reads_df$cell_id))
  ncells
  
  ## Using field = "copy", copy number values as the input for cell clustering
  clusters <- signals::umap_clustering(reads_df,
                                       minPts = max(round(pctcells * ncells), 20),
                                       field = "copy")
  
  # clusters <- signals::umap_clustering(reads_df, 
  #                                      minPts = min(round(pctcells * ncells), 10), 
  #                                      field = "copy")
  
  tree <- clusters$tree
  clones <- clusters$clustering
  print(unique(clones$clone_id))
  
  
  ## Save output of clones here
  # data.table::fwrite(clones, paste0(save_dir, 'combined_cell_clones.csv'))
  # write.tree(phy = tree, file = paste0(save_dir, 'combined_tree.newick'), tree.names = F)
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
  tree_fn <- paste0(save_dir,datatag,'_tree.newick')
  ape::write.tree(phy = tree, file = tree_fn, tree.names = F)
  
  data.table::fwrite(clones, paste0(save_dir,datatag,'_cell_clones.csv'))
  
  output_fn <- paste0(save_dir,datatag,'_hdbscan_heatmap.png')
  print(output_fn)
  png(output_fn, height = 2*700, width=2*1200, res = 2*72)
  make_cell_copynumber_tree_heatmap(
    tree, copynumber, clones, NULL, grouping_file
  )
  dev.off()
  
  
  ## NOTED: remove clone A from the list
  # clones <- clones %>%
  #   dplyr::filter(clone_id!='A')
  # dim(clones)
  # reads_df <- reads_df %>%
  #   dplyr::filter(cell_id %in% clones$cell_id)
  # data.table::fwrite(reads_df, paste0(results_dir, 'total_reads.csv.gz'))
  return(clones)
}

# reads_fn <- paste0(results_dir, 'total_reads.csv.gz')
# metrics_fn <- paste0(results_dir, 'SA535X4XB05649_total_cells_clones.csv.gz')
prepare_RT_input_data <- function(cell_clones_fn, results_dir, download_dir, datatag, obs_libs, 
                                  grouping_file, reads_fn=NULL, metrics_fn=NULL, filtered_fn=NULL){
  
  # Depending on the results of hdbscan, you may want to remove the 'very' noisy cells first
  
  if(is.null(reads_fn)){
    reads_fn <- paste0(download_dir, datatag,'_total_reads.csv.gz')
  }
  if(is.null(metrics_fn)){
    metrics_fn <- paste0(results_dir, datatag,'_total_cells_metrics.csv.gz')
  }
  if(is.null(filtered_fn)){
    filtered_fn <- paste0(results_dir, datatag,'_total_merged_filtered_states.csv.gz')
  }
  if(file.exists(filtered_fn)){
    filtered_df <- data.table::fread(filtered_fn) %>% as.data.frame()  
  }else{
    filtered_df <- load_filtered_data(obs_libs, paste0(results_dir,'filtered_data/'), datatag, filtered_fn) 
  }
  print(dim(filtered_df))
  print(colnames(filtered_df)[1])
  print(rownames(filtered_df)[1])
  cells_id <- colnames(filtered_df)
   
  
  if(file.exists(reads_fn)){
    reads_df <- data.table::fread(reads_fn) %>% as.data.frame()  
  }else{
    reads_df <- load_hmmcopy_reads(obs_libs, download_dir, datatag, 
                                   save_data=T, cells_id, reads_fn)
  }
  dim(reads_df)  
  
  # metrics_df <- data.table::fread(paste0(input_dir,library_id,'_filtered_cell_clones.csv'))
  if(file.exists(metrics_fn)){
    metrics_df <- data.table::fread(metrics_fn)
  }else{
    metrics_df <- load_metrics_data(obs_libs, download_dir, results_dir, 
                                    datatag, cells_id, metrics_fn)
  }
  print(dim(metrics_df))
  lids <- sapply(strsplit(metrics_df$cell_id,'-'), function(x){
    return(x[2])
  })
  metrics_df$library_id <- as.character(lids)
  
  
  ## Run hdbscan first
  ## Then removing noisy cells from the list. 
  # download_dir, results_dir, obs_libs, grouping_file,
  # reads_df, 
  
  hdbscan_cell_clones <- hdbscran_viz(download_dir, results_dir, obs_libs, grouping_file,
             reads_df, copynumber=filtered_df,
             copynumber_fn=NULL, pctcells = 0.05)
  ## 
  datatag
  # hdbscan_cell_clones$cell_id[1]
  # excluded_clones <- c('F') #very noisy cells for sample id: SA535X4XB05391
  # excluded_clones <- c('C') #very noisy cells for sample id: SA535X4XB05462
  excluded_clones <- c('E') #very noisy cells for sample id: SA535X4XB02498
  hdbscan_cell_clones <- hdbscan_cell_clones %>%
    dplyr::filter(!clone_id %in% excluded_clones)
  dim(hdbscan_cell_clones)
  reads_df <- reads_df %>%
    dplyr::filter(cell_id %in% hdbscan_cell_clones$cell_id)
  dim(reads_df)
  data.table::fwrite(reads_df, reads_fn)
  
  metrics_df <- metrics_df %>%
    dplyr::filter(cell_id %in% hdbscan_cell_clones$cell_id)
  dim(metrics_df)
  data.table::fwrite(metrics_df, metrics_fn)
  
  filtered_df <- filtered_df[,hdbscan_cell_clones$cell_id]
  data.table::fwrite(filtered_df, filtered_fn, row.names = TRUE)
  
  if(file.exists(cell_clones_fn)){
    cell_clones <- data.table::fread(cell_clones_fn)
  }else{
    stop('Double check input cell_clones_fn!!!')
  }
  cell_clones <- cell_clones %>%
    dplyr::filter(cell_id %in% metrics_df$cell_id & clone_id!='None') %>%
    dplyr::mutate(is_s_phase=FALSE, cell_type_status='cn_g')
  # data.table::fwrite(cell_clones, paste0(results_dir,datatag,'_main_cells_clones.csv.gz'))
  dim(cell_clones)
  metrics_df1 <- metrics_df %>%
    dplyr::filter(!cell_id %in% cell_clones$cell_id) %>%
    dplyr::mutate(clone_id='unassigned', cell_type_status='cn_s') %>%
    dplyr::select(cell_id, clone_id, is_s_phase, cell_type_status)
  dim(metrics_df1)
  
  cell_clones <- dplyr::bind_rows(cell_clones, metrics_df1)
  lids <- sapply(strsplit(cell_clones$cell_id,'-'), function(x){
    return(x[2])
  })
  cell_clones$library_id <- as.character(lids)
  data.table::fwrite(cell_clones, paste0(results_dir,datatag,'_total_cells_clones.csv.gz'))
  dim(cell_clones)
  cn_g_cells <- cell_clones %>% 
    dplyr::filter(cell_type_status=='cn_g') %>% 
    dplyr::pull(cell_id)
  
  cn_s_cells <- cell_clones %>% 
    dplyr::filter(cell_type_status=='cn_s') %>% 
    dplyr::pull(cell_id)
  print(length(cn_g_cells))
  print(length(cn_s_cells))
  reads_df_g <- reads_df %>% 
    dplyr::filter(cell_id %in% cn_g_cells)
  reads_df_s <- reads_df %>% 
    dplyr::filter(cell_id %in% cn_s_cells)
  # 
  # length(unique(reads_df_g$cell_id))
  # length(unique(reads_df_s$cell_id))
  # length(cn_g_cells)
  # length(cn_s_cells)
  data.table::fwrite(reads_df_g, paste0(results_dir,datatag,'_filtered_reads_RT_g_cells.csv.gz'))
  data.table::fwrite(reads_df_s, paste0(results_dir,datatag,'_filtered_reads_RT_s_cells.csv.gz'))
  
}

correct_cell_clones <- function(){
  results_dir <- '/home/htran/storage/datasets/metastasis_results/replication_timing/'
  download_dir <- '/home/htran/storage/raw_DLP/metastasis_DLP/SA535/'
  input_dir <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/dlp_trees/'
  grouping_file <- paste0(input_dir, 'SA535/library_groupings.csv.gz')
  grouping <- data.table::fread(paste0(input_dir, 'SA535/library_groupings.csv.gz'))
  grouping <- grouping %>%
    dplyr::rename(library_id=grouping)
  sample_ids <- c('SA535X4XB02498','SA535X4XB05391','SA535X4XB05462','SA535X4XB05649')
  for(datatag in sample_ids){
    clones_fn <- paste0(results_dir,datatag,'_total_cells_clones.csv.gz')
    cell_clones <- data.table::fread(clones_fn)  
    print(dim(cell_clones))
    lids <- sapply(strsplit(cell_clones$cell_id,'-'), function(x){
      return(x[2])
    })
    cell_clones$library_id <- as.character(lids)
    data.table::fwrite(cell_clones, clones_fn)
  }
  
}
# cell_clones_fn <- paste0(input_dir, 'SA535/cell_clones.csv.gz')
# metaSamples_fn <- paste0(input_dir, 'SA535/library_groupings.csv.gz')
get_library_ids <- function(cell_clones_fn, metaSamples_fn, obs_sample, mouse_id){
  cell_clones <- data.table::fread(cell_clones_fn)
  grouping <- data.table::fread(metaSamples_fn)
  sids <- sapply(strsplit(cell_clones$cell_id,'-'), function(x){
    return(x[1])
  })
  lids <- sapply(strsplit(cell_clones$cell_id,'-'), function(x){
    return(x[2])
  })
  cell_clones$library_id <- as.character(lids)
  cell_clones$sample_id <- as.character(sids)
  unique(cell_clones$sample_id)
  cell_clones <- cell_clones %>%
    left_join(grouping, by=c('library_id'='grouping', 'sample_id'))
  dim(cell_clones)
  meta_df <- cell_clones %>%
    dplyr::filter(clone_id!='None' & sample_id==obs_sample,
                  mainsite=="Primary" & pdxid==mouse_id)%>%
    dplyr::group_by(library_id, sample_id) %>%
    dplyr::summarise(nb_cells=n()) %>%
    dplyr::select(-nb_cells)
  data.table::fwrite(meta_df, paste0(output_dir,'obs_library_ids.csv'))
  return(meta_df)
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


results_dir <- '/home/htran/storage/datasets/metastasis_results/replication_timing/'
download_dir <- '/home/htran/storage/raw_DLP/metastasis_DLP/SA535/'
input_dir <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/dlp_trees/'
grouping_file <- paste0(input_dir, 'SA535/library_groupings.csv.gz')
cell_clones_fn <- paste0(input_dir, 'SA535/cell_clones.csv.gz')
# cell_clones <- data.table::fread(cell_clones_fn)
grouping <- data.table::fread(paste0(input_dir, 'SA535/library_groupings.csv.gz'))
# obs_libs <- c("A98165A","A98197B","A98261A","A98232A","A98277B")
# ticket_ids <- c("SC-3377","SC-3403","SC-3213","SC-3173","SC-2716")
# View(grouping)
# colnames(grouping)
# t$library_id
t <- grouping %>%
  dplyr::rename(library_id=grouping)
# data.table::fwrite(t, paste0(input_dir, 'SA535/library_groupings_test.csv.gz'))
# t <- data.table::fread(paste0(input_dir, 'SA535/jira_tickets_library_groupings.csv.gz'))
# # View(t)
# dim(t)

library(dplyr)
t1 <- t %>%
  dplyr::filter(mainsite=="Primary" & sample_id!="SA535X4XB05649") #%>%
  # dplyr::rename(library_id=grouping)
  # dplyr::group_by(sample_id, pdxid) %>%
  # dplyr::summarise(nb_libraries=n())
t1
# sample_id      pdxid      nb_libraries
# <chr>          <chr>             <int>
# 1 SA535X4XB02498 X0011_2364            3
# 2 SA535X4XB05391 X0011_2361            1
# 3 SA535X4XB05462 X0011_2363            3


datatag <- 'SA535X4XB05391'
datatag <- 'SA535X4XB05462'
datatag <- 'SA535X4XB02498'
obs_libs  <- t %>%
  dplyr::filter(sample_id==datatag) %>%
  dplyr::pull(library_id)
obs_libs
cell_clones_fn
input_dir
download_dir
results_dir

reads_fn=NULL
metrics_fn=NULL
filtered_fn=NULL











# data.table::fwrite(t1, paste0(input_dir, 'SA535/jira_tickets_library_groupings_M1_M3_M4.csv.gz'))
# t <- t %>%
#   dplyr::filter(library_id %in% obs_libs)
# dim(t)
# t$sample_id[1]
# data.table::fwrite(t, paste0(input_dir, 'SA535/jira_tickets_library_groupings_M2.csv.gz'))
# View(t)
# 
# input_dir <- '/home/htran/storage/datasets/metastasis_results/replication_timing/'
# load_filtered_data(obs_libs, input_dir)
# results_dir <- input_dir
# 
# metrics_df <- load_metrics_data(obs_libs, download_dir, results_dir)
# metrics_df <- metrics_df %>%
#   dplyr::filter(cell_id %in% clones$cell_id)
# summary(as.factor(metrics_df$is_s_phase))
# table(metrics_df$is_s_phase, metrics_df$quality>0.75)
# 
# 
# metrics_df <- data.table::fread(paste0(results_dir, 'total_metrics.csv.gz'))
# dim(metrics_df)
# summary(metrics_df$cell_id)
# 
# head(cell_clones)
# dim(cell_clones)
# summary(as.factor(cell_clones$clone_id))
# 
# unique(metrics_df$is_s_phase)
# # cells are not assigned to clonal labels: s-phase cells,..
# metrics_df1 <- metrics_df %>%
#   dplyr::filter(!cell_id %in% cell_clones$cell_id) %>%
#   dplyr::mutate(clone_id='unassigned', cell_type_status='cn_s') %>%
#   dplyr::select(cell_id, clone_id, is_s_phase, cell_type_status)
# head(metrics_df1)
# colnames(cell_clones)
# cell_clones_total <- dplyr::bind_rows(cell_clones, metrics_df1)
# dim(cell_clones_total)
# summary(as.factor(metrics_df1$experimental_condition))
# data.table::fwrite(cell_clones_total, paste0(results_dir, 'SA535X4XB05649_total_cells_clones.csv.gz'))
# 
# clonal_df <- data.table::fread(paste0(results_dir, 'RT_results/SA535X4XB05649_clonal_RT.csv.gz'))
# 
# clonal_df <- clonal_df %>%
#   dplyr::filter(cell_type_status=='cn_s')
# summary(as.factor(clonal_df$clone_id))
# dim(clonal_df)
# View(head(clonal_df))
# ## To Do: 
# ## Cells that are not assigned to clone --> cn_s
# ## Cells that are assigned to clone --> cn_g
# ##  Run functions
# 
# 
# results_dir <- '/home/htran/storage/datasets/metastasis_results/replication_timing/RT_results/'
# out_s_df <- data.table::fread(paste0(results_dir, 'SA535X4XB05649_cn_s_with_scrt.csv.gz'))
# 
# colnames(out_s_df)
# 
# frac_df <- data.table::fread(paste0(results_dir, 'SA535X4XB05649_RT_cell_metrics_RT_frac_clones.csv.gz'))
# dim(frac_df)
# # 53 cells
# outlier_thrs <- 0.001
# frac_df <- frac_df %>%
#   dplyr::filter(cell_type_status=='cn_s' & cell_frac_rep <= (1-outlier_thrs) & cell_frac_rep >= outlier_thrs) #
# dim(frac_df)
# # 26 cells
# frac_df <- frac_df %>%
#   dplyr::filter(cell_type_status=='cn_g' & cell_frac_rep <= (1-outlier_thrs) & cell_frac_rep >= outlier_thrs)
# 
# thrs <- 0.05
# frac_s <- frac_df %>%
#   dplyr::filter(cell_type_status=='cn_s' & cell_frac_rep <=thrs)
# thrs <- 0.95
# frac_s <- frac_df %>%
#   dplyr::filter(cell_type_status=='cn_s' & cell_frac_rep >= thrs)
# 
# 
# thrs <- 0.05
# frac_g <- frac_df %>%
#   dplyr::filter(cell_type_status=='cn_g' & cell_frac_rep <=thrs)
# thrs <- 0.95
# frac_g <- frac_df %>%
#   dplyr::filter(cell_type_status=='cn_g' & cell_frac_rep >= thrs)
# 
# dim(frac_g)
# summary(as.factor(frac_df$cell_type_status))
# summary(frac_df$cell_frac_rep)
# 
# summary(bk_df$rep_bk)
# library(ggplot2)
# frac_df <- frac_df %>%
#   dplyr::filter(cell_type_status=='cn_s')
# frac_df$clone_id <- as.factor(frac_df$clone_id)
# frac_df <- frac_df %>%
#   dplyr::group_by(clone_id) %>%
#   dplyr::summarise(nb_cells=n())
# p <- ggplot(data=frac_df, aes(x=clone_id, y=nb_cells)) +
#     geom_bar(stat="identity", width=0.5, aes(fill=clone_id))
# p
# df <- frac_g
# p <- ggplot(data=df, aes(x=clone_id, y=cell_frac_rep)) +
#     geom_jitter(position=position_jitter(0.2), aes(color=clone_id))  +
#     stat_summary(fun=mean, geom="point", shape=18,
#                size=3, color="black")    
# p
# 
# frac_df <- frac_df %>%
#   dplyr::filter(cell_type_status=='cn_s' & (cell_frac_rep<0.25 | cell_frac_rep>0.75)) 
#   # dplyr::group_by(clone_id) %>%
#   # dplyr::mutate(rep_bk_025= quantile(rep_bk, prop=0.25), rep_bk_075=quantile(rep_bk, prop=0.75))
# 
# bk_df <- data.table::fread(paste0(results_dir, 'SA535X4XB05649_RT_cell_metrics_clones.csv.gz'))
# dim(bk_df)
# bk_df <- bk_df %>%
#   dplyr::filter(cell_type_status=='cn_s')
# bk_df$clone_id <- as.factor(bk_df$clone_id)
# p <- ggplot(data=bk_df, aes(x=clone_id, y=rep_bk, color=clone_id)) +
#   geom_jitter(position=position_jitter(0.2)) +
#   stat_summary(fun=median, geom="point", shape=18,
#                  size=3, color="black")
# # facet_wrap(~ cell_type_status, strip.position = "top")
# p
