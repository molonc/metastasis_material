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

load_hmmcopy_reads <- function(obs_libs, download_dir, cells_id=NULL, datatag='SA', save_data=FALSE){
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
        dplyr::select(chr, start, end, copy, cell_id, state) # 
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
  reads_df <- as.data.frame(dplyr::bind_rows(obs_reads))
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
  
  nb_hq_cells <- c()
  for(l in obs_libs){
    lib_fn <- paste0(input_dir,l,'_filtered_states.csv.gz')
    
    if(file.exists(lib_fn)){
      tmp <- data.table::fread(lib_fn)
      rownames(tmp) <- tmp$V1
      tmp$V1 <- NULL
      
      print(l)
      print(dim(tmp))
      nb_hq_cells <- c(nb_hq_cells, dim(tmp)[2])
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
  metacells <- tibble::tibble(library_id=obs_libs, nb_filtered_cells=nb_hq_cells)
  
  data.table::fwrite(metacells, paste0(input_dir, datatag,'_QC_cells.csv.gz'), row.names = TRUE)
  
  # filtered_df <- as.data.frame(dplyr::bind_cols(obs_reads))
  # print(dim(filtered_df))
  # rownames(filtered_df) <- rownames(tmp)
  # if(is.null(filtered_fn)){
  #   filtered_fn <- paste0(input_dir, datatag,'_total_merged_filtered_states.csv.gz')
  # }
  # data.table::fwrite(filtered_df, filtered_fn, row.names = TRUE)
  # 
  # return(filtered_df)
  return(metacells)
}

## pctcells: percentage of minimum cells to consider as one clone
hdbscran_viz <- function(datatag, download_dir, results_dir, obs_libs, grouping_file, 
                         copynumber_fn=NULL, pctcells = 0.05){
  # pctcells = 0.05
  ## for single library
  # copynumber <- data.table::fread(paste0(results_dir,library_id,'_filtered_states.csv.gz'))%>% as.data.frame()
  # if(is.null(copynumber_fn)){
  #   copynumber_fn <- paste0(results_dir,library_id,'_filtered_states.csv.gz')
  # }
  if(is.null(copynumber_fn)){
    # copynumber_fn <- paste0(results_dir,library_id,'_filtered_states.csv.gz')
    copynumber_fn <- paste0(results_dir, 'total_merged_filtered_states.csv.gz')
    
  }
  metacells <- load_filtered_data(obs_libs, results_dir, datatag, filtered_fn=copynumber_fn)
  
  # copynumber <- load_filtered_data(obs_libs, results_dir, datatag, filtered_fn=copynumber_fn)
  ## for combined libraries
  copynumber <- data.table::fread(copynumber_fn)%>% as.data.frame()
  rownames(copynumber) <- copynumber$V1
  copynumber$V1 <- NULL
  print(dim(copynumber))
  print(rownames(copynumber)[1])  
  print(colnames(copynumber)[1])  
  ## for single library
  # reads_df <- data.table::fread(paste0(download_dir,library_id,'/hmmcopy/','reads.csv.gz')) %>% as.data.frame()
  # reads_df <- data.table::fread(paste0(download_dir,library_id,'/hmmcopy/',library_id,'_reads.csv.gz')) %>% as.data.frame()
  
  ## for combined libraries
  cells_id <- colnames(copynumber)
  reads_df <- load_hmmcopy_reads(obs_libs, download_dir, cells_id, datatag)
  print(dim(reads_df))
  # reads_df$copy[1:10]
  
  reads_df <- reads_df %>%
    dplyr::filter(cell_id %in% cells_id)
  
  reads_df$chr_desc <- paste0(reads_df$chr,'_',reads_df$start,'_',reads_df$end)
  reads_df <- reads_df %>%
    dplyr::filter(chr_desc %in% rownames(copynumber))
  
  print(dim(reads_df))
  ncells <- length(unique(reads_df$cell_id))
  print(ncells)
  ## Using field = "copy", copy number values as the input for cell clustering
  clusters <- signals::umap_clustering(reads_df,
                                       minPts = min(round(pctcells * ncells), 50),
                                       field = "copy")
  
  # clusters <- signals::umap_clustering(reads_df, 
  #                                      minPts = min(round(pctcells * ncells), 10), 
  #                                      field = "copy")
  
  tree <- clusters$tree
  clones <- clusters$clustering
  unique(clones$clone_id)
  
 
  
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
  
  
  output_fn <- paste0(save_dir,datatag,'_combined_hdbscan_heatmap.png')
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
# obs_libs <- unique(lib_df$grouping)
# library_id <- lib_df$grouping[1]
# 
# for(library_id in lib_df$grouping[2:dim(lib_df)[1]]){
#   hdbscran_viz(download_dir, results_dir, library_id, grouping_file)
# }




results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919_X8/'
download_dir <- '/home/htran/storage/raw_DLP/metastasis_DLP/SA919/'
grouping_file <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/dlp_trees/SA919/SA919_passageX8.csv'
lib_df <- data.table::fread(grouping_file)

obs_libs <- unique(lib_df$library_id)
obs_libs
datatag <- 'SA919_X8'
copynumber_fn=NULL
pctcells = 0.05

lib_df <- lib_df %>%
  left_join(metacells, by='library_id')
dim(lib_df)
lib_df$grouping <- lib_df$library_id
data.table::fwrite(lib_df,grouping_file)
lib_df$library_id <- NULL
