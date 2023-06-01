
## library_groupings.csv: cells metadata, sampleid: contain passage info, pdx id: original mouse id
## main site: met or primary
## library labels: ignore this column - unique index for computation 
## library id: id of sequencing input data
# library_labels library_id   mainsite        origin      pdxid      sample_id
# 1        A98254A    A98254A Metastasis Left_Axillary X0011_2362 SA535X4XB09109
# 2       A118855A   A118855A Metastasis      Axillary X0011_2364 SA535X4XB05985


# head(cell_clones): contain clonal infos, output of phylogenetic tree
## cell_id: sampleid-libraryid-rowid-columnid
# cell_id clone_id
# 1 SA535X4XB02498-A98163A-R09-C11        I
# 2 SA535X4XB02498-A98163A-R11-C16        I

# results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA609/'
# datatag <- 'SA609'
# 
# t <- cell_clones %>%
#   dplyr::filter(library_id=='A98289B')
# dim(cell_clones)
# print(summary(as.factor(t$clone_id)))
# t <- cell_clones %>%
#   dplyr::filter(clone_id=='R')

library(dplyr)
# results_dir <- '/home/htran/storage/datasets/metastasis_results/SA535_new_encoding/SA535_wholedata_v2'
# datatag <- 'SA535'
# outputfile <- paste0(results_dir,'/prevalences/clone_samples_Gurdeep/samples.csv')


results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_wholedata/'
datatag <- 'SA919'
outputfile <- paste0(results_dir,'/prevalences/samples_cloneId_mapping.csv')



get_library_id <- function(cell_ids) {
  labels <- sapply(strsplit(cell_ids, "-"), function(x) {
    return(x[2])
  })
  return(as.character(labels))
}
get_sample_id <- function(cell_ids) {
  labels <- sapply(strsplit(cell_ids, "-"), function(x) {
    return(x[1])
  })
  return(as.character(labels))
}


# Combining cells clone info + metadata info
get_meta_data <- function(cell_clones, results_dir, library_grouping_fn=NULL){
  if(is.null(library_grouping_fn)){
    library_grouping_fn <- paste0(results_dir,'library_groupings.csv')
  }
  
  grouping_df <- read.csv(library_grouping_fn, header=T,check.names=F, 
                          stringsAsFactors = F)
  colnames(grouping_df)[which(colnames(grouping_df) == "grouping")] <- "library_id"
  print(colnames(grouping_df))
  print(head(grouping_df))
  cell_clones$library_id <- get_library_id(cell_clones$cell_id)
  cell_clones$sample_id <- get_sample_id(cell_clones$cell_id)
  cell_clones <- cell_clones %>% left_join(grouping_df, by = c("library_id","sample_id"))
  return(cell_clones)
}

get_prevalence <- function(results_dir, outputfile=NULL, cellclones=NULL, datatag = ''){
  
  print("Generating prevalence plots...")
  results_dir <- paste0(results_dir,'/')
  if(is.null(cellclones)){
    cellclones <- paste0(results_dir, 'cell_clones.csv')
  }
  save_dir <- paste0(dirname(outputfile),'/')
  if (!file.exists(save_dir)){
    dir.create(save_dir, recursive=T)
  }
  cell_clones <- read.delim(cellclones, stringsAsFactors = FALSE, sep = ",")
  excluded_clones <- c('Un','unassigned','None')
  cell_clones <- cell_clones %>%
    dplyr::filter(!clone_id %in% excluded_clones)
  print(dim(cell_clones))
  print(summary(as.factor(cell_clones$clone_id)))
  # cell_clones_backup <- cell_clones
  library_grouping_fn <- paste0(results_dir,'library_groupings.csv')
  meta_data <- get_meta_data(cell_clones, results_dir, library_grouping_fn)
  print(dim(meta_data))
  
  # Get samples - clones mapping for DE analysis between clones - bulk RNAseq 
  tmp <- meta_data %>%
    dplyr::filter(clone_id=='C') %>%
    dplyr::group_by(sample_id, clone_id) %>%
    dplyr::summarise(nb_cells=n())
  tmp1 <- tmp %>%
    dplyr::group_by(sample_id) %>%
    dplyr::filter(nb_cells==max(nb_cells))
  dim(tmp1)
  head(tmp1)
  data.table::fwrite(tmp1, outputfile)
  
  tmp1 <- tmp %>%
    dplyr::filter(sample_id=='SA919X7XB05692')
  
  
  # Extract data for each pdx, actually pdx here is mouse id
  pdxids <- unique(meta_data$pdxid)
  for(pd in pdxids){
    tmp <- meta_data %>%
      dplyr::filter(pdxid==pd)
    data.table::fwrite(tmp, paste0(save_dir,'clones_metainfo_mouse_',pd,'.csv'))
  }  
}
base_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/snakemake_10x_v2/'
meta_10x <- data.table::fread(paste0(base_dir,'SA535_10x_metadata.csv')) %>% as.data.frame()
View(meta_10x)

meta_10x <- meta_10x %>%
  dplyr::filter(analysis_status=='included')
dim(meta_10x)
View(head(meta_data))
obs_feature <- 'sample_id'
obs_feature_10x <- 'mouse_id'
pd1 <- get_stat_prevalence_same_clone(meta_data, meta_10x, obs_feature,obs_feature_10x)


obs_feature <- 'origin'
obs_feature_10x <- 'Site_origin'
pd1 <- get_stat_prevalence_same_clone(meta_data, meta_10x, obs_feature,obs_feature_10x)

# meta_data$pdxid <- gsub('X0011_','M',meta_data$pdxid)
obs_feature <- 'pdxid'
obs_feature_10x <- 'pdxid'
pd1 <- get_stat_prevalence_same_clone(meta_data, meta_10x, obs_feature,obs_feature_10x)



get_stat_prevalence_same_clone <- function(meta_data, meta_10x, obs_feature,obs_feature_10x){
  fts <- unique(meta_10x[,obs_feature_10x])
  obs1 <- 'clone_id'
  pd <- as.data.frame(table(meta_data[,obs1], meta_data[,obs_feature]))
  # head(pd)
  colnames(pd) <- c(obs1,obs_feature,'nbcells')
  pd <- pd %>%
    dplyr::filter(nbcells>30 & !!sym(obs_feature) %in% fts)
  clone_used <- pd %>%
    dplyr::group_by(!!sym(obs1))%>%
    dplyr::summarise(nb_occurence=n())%>%
    dplyr::filter(nb_occurence>1)%>%
    dplyr::pull(!!sym(obs1))
  clone_used
  pd <- pd %>%
    dplyr::filter(!!sym(obs1) %in% clone_used)
  
  print(dim(pd))
  # head(pd)
  # summary(as.factor(pd[,obs_feature]))
  pd1 <- pd %>% tidyr::pivot_wider(names_from = !!sym(obs_feature), values_from = 'nbcells')
  
  # obs_clones <- c('C','A','H','D','I','N','R')
  # pd1 <- pd %>%
  #   dplyr::filter(clone %in% obs_clones)  
  data.table::fwrite(pd1, paste0(save_dir, obs_feature,'_same_clone.csv'))
  return(pd1)
}
