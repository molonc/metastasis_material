library(dplyr)


# R24. The manuscript would benefit from the inclusion of a figure showing the raw single cell data, including sequencing read counts and inferred copy number profiles. Such a figure would be important for assessing the noise levels in the data and for evaluating confidence in the identified clones, which underpin many of the key conclusions.
# Hoa: 
# Inferred copy number states based on read counts plot
# Boxplot for total reads
# Heatmap for all cells 
# Copy number profiles dot plots for representative cells - median profile in each clone
# Hamming distance between each cell versus median copy number profile showing small variance, i.e 2% of total mutation events due to noises, small background mutation

# R25. Many of the central results rely on the definition of clones derived from single cell clustering. 
# Given the inherent noise and uncertainty associated with single cell clustering, 
# it would be important to assess the robustness of these findings. 
# For example, the authors could test whether key conclusions like the identification of metastasizing clones 
# are preserved when using alternative clustering approaches, subsampling cells or genomic regions, 
# or applying other methods to evaluate clustering stability.
# Boxplot of Hamming/Mahattan distances between clones showing large change/driver of mutation events. 

# R14- Are there significant differences in copy number signatures between different clones?
  

get_heatmap_CNA_change <- function(){
  # git dir
  # load cna file
  # load cell clones, library grouping 
  input_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/'
  sa919_metasamples <- data.table::fread(paste0(input_dir, 'materials/dlp_trees/SA919/library_groupings.csv.gz'))
  head(sa919_metasamples)
  sa919_metasamples$pdxid <- gsub('X0847','X0847-',sa919_metasamples$pdxid)
  sa919_metasamples
  sa535_metasamples <- data.table::fread(paste0(input_dir, 'materials/dlp_trees/SA535/library_groupings.csv.gz'))
  head(sa535_metasamples)
  sa535_metasamples$grouping
  sa535_metasamples$pdxid <- gsub('X0011_','X0011-',sa535_metasamples$pdxid)
  dplyr::bind_rows()
  dim(sa535_metasamples)
  dim(sa919_metasamples)
  # copynumber_fn <- paste0(input_dir, 'revision/clustering_eval/total_merged_filtered_states.csv')
  
  input_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/'
  copynumber <- data.table::fread(paste0(input_dir, 'revision/clustering_eval/total_merged_filtered_states.csv'))%>% as.data.frame()
  print(dim(copynumber))
  datatag <- 'SA919'
  
  input_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/'
  copynumber <- data.table::fread(paste0(input_dir, 'dlp_trees/SA535/total_merged_filtered_states.csv.gz'))%>% as.data.frame()
  print(dim(copynumber))
  datatag <- 'SA535'
  
  
  
  copynumber <- copynumber %>%
    dplyr::rename(chr_desc=V1) 
  set.seed(42)
  sampling_pct <- 0.2
  cell_ids <- colnames(copynumber)
  cell_ids <- cell_ids[!cell_ids %in% c('chr_desc')]
  sampled <- sample(cell_ids, size = length(cell_ids) * sampling_pct)
  length(sampled)
  sampled[1]
  copynumber <- copynumber[,c('chr_desc',sampled)]
  dim(copynumber)
  # View(copynumber[1:3,1:3])
  copynumber <- get_chr_infos(copynumber)
  copynumber_rand <- copynumber %>%
    dplyr::select(-chr_desc) %>%
    # slice_sample(prop = sampling_prop) %>% 
    tidyr::pivot_longer(cols=-c(chrom, start, end), 
                        names_to = 'CELL', values_to='CN states')%>% 
    dplyr::select(all_of(c('CELL', 'chrom', 'start', 'end', 'CN states'))) %>%
    as.data.frame()
  dim(copynumber_rand)
  head(copynumber_rand)
  # copynumber <- copynumber %>%
  #   dplyr::select(-chr_desc) %>%
  #   tidyr::pivot_longer(cols=-c(chrom, start, end), 
  #                       names_to = 'CELL', values_to='CN states')%>% 
  #   dplyr::select(all_of(c('CELL', 'chrom', 'start', 'end', 'CN states'))) %>%
  #   as.data.frame()
  # 
  # head(copynumber)
  # dim(copynumber)
  
  # data.table::fwrite(copynumber, paste0(input_dir, datatag,'_total_merged_filtered_states_longformat.csv.gz'))
  data.table::fwrite(copynumber_rand, paste0(input_dir, 'revision/clustering_eval/', datatag,'_total_merged_filtered_states_longformat_rand_',sampling_pct,'.tsv'), sep = '\t')
  data.table::fwrite(copynumber, paste0(input_dir, 'revision/clustering_eval/', datatag,'_total_merged_filtered_states_longformat.tsv'), sep = '\t')
  dim(copynumber_rand)
  length(unique(copynumber_rand$CELL))
}

library(dplyr)
get_chr_infos <- function(median_cnv) {
  chr_desc <- as.character(median_cnv$chr_desc)
  obs_chrs = as.character(c(paste0(1:22), "X"))
  
  chrs <- sapply(strsplit(chr_desc, "_"), function(x) {
    return(x[1])
  })
  starts <- sapply(strsplit(chr_desc, "_"), function(x) {
    return(x[2])
  })
  ends <- sapply(strsplit(chr_desc, "_"), function(x) {
    return(x[3])
  })
  median_cnv$chrom <- as.character(chrs)
  median_cnv$start <- as.numeric(starts)
  median_cnv$end <- as.numeric(ends)
  median_cnv <- median_cnv %>% 
    dplyr::filter(chrom %in% obs_chrs)
  median_cnv$chrom <- paste0('chr',median_cnv$chrom)
  # print(unique(median_cnv$chrom))
  # return(list(chr=as.character(chrs),start=as.character(starts),end=as.character(ends)))
  return(median_cnv)
}



sampling_data <- function(datatag, save_dir, copynumber_fn, sampling_pct=0.05, 
                          cellclone_fn=NULL){ #, library_grouping_fn=NULL
  input_dir <- paste0(dirname(copynumber_fn),'/')
  if(is.null(cellclone_fn)){
    cellclone_fn <- paste0(input_dir,'cell_clones.csv.gz')  
  }
  # if(is.null(library_grouping_fn)){
  #   library_grouping_fn <- paste0(input_dir,'library_groupings.csv.gz')  
  # }
  save_dir_data <- paste0(save_dir, datatag, '_',sampling_pct,'/')
  if(!dir.exists(save_dir_data)){
    dir.create(save_dir_data, recursive = T)
  }
  
  # copynumber <- read.csv(copynumber_fn, header=T, row.names = 1, check.names = F,stringsAsFactors = FALSE)
  copynumber <- as.data.frame(data.table::fread(copynumber_fn))
  copynumber <- copynumber %>%
    dplyr::rename(chr_desc=V1) 
  # cell_clones contain 2 columns of cell_id, and clone_id
  # ex:           cell_id                      clone_id
  # 1   SA535X4XB02498-A98163A-R09-C11          C
  cell_clones <- data.table::fread(cellclone_fn) %>% as.data.frame()
  print(dim(cell_clones))
  cell_clones <- cell_clones %>%
    dplyr::filter(!clone_id %in% c('None','unassigned'))
  print(unique(cell_clones$clone_id))
  
  ## Subsampling data per clone
  set.seed(42)
  sampled <- c()
  for(cl in unique(cell_clones$clone_id)){
    cell_ids <- cell_clones %>%
      dplyr::filter(clone_id==cl) %>%
      dplyr::pull(cell_id)
    sampled_per_clone <- sample(cell_ids, size = length(cell_ids) * sampling_pct)
    sampled <- c(sampled, sampled_per_clone)
  }
  print('Number of randomized cells is: ')
  print(length(sampled))
  
  copynumber <- copynumber[,c('chr_desc',sampled)]
  print(dim(copynumber))
  # View(copynumber[1:3,1:3])
  copynumber <- get_chr_infos(copynumber)
  copynumber_rand <- copynumber %>%
    dplyr::select(-chr_desc) %>%
    # slice_sample(prop = sampling_prop) %>% 
    tidyr::pivot_longer(cols=-c(chrom, start, end), 
                        names_to = 'CELL', values_to='CN states')%>% 
    dplyr::select(all_of(c('CELL', 'chrom', 'start', 'end', 'CN states'))) %>%
    as.data.frame()
  # dim(copynumber_rand)
  data.table::fwrite(copynumber_rand, paste0(save_dir_data, datatag,'_copynumber_longformat_',sampling_pct,'.tsv'), sep = '\t')
  
  cell_clones <- cell_clones %>%
    dplyr::filter(cell_id %in% sampled)
  dim(cell_clones)
  data.table::fwrite(cell_clones, paste0(save_dir_data, datatag,'_cell_clones_',sampling_pct,'.csv.gz'))
  
}


input_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/'
save_dir <- paste0(input_dir,'revision/clustering_eval/DICE_SA919/')
datatag <- 'SA919'
copynumber_fn <- paste0(input_dir,'materials/dlp_trees/',datatag,'/total_merged_filtered_states.csv.gz')
# sampling_pct=0.05
sampling_pct=0.5
cellclone_fn=NULL
sampling_data(datatag, save_dir, copynumber_fn, sampling_pct,cellclone_fn)



input_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/'
save_dir <- paste0(input_dir,'revision/clustering_eval/DICE_SA535/')
datatag <- 'SA535'
copynumber_fn <- paste0(input_dir,'materials/dlp_trees/',datatag,'/total_merged_filtered_states.csv.gz')
# sampling_pct=0.05
sampling_pct=0.6
cellclone_fn=NULL
sampling_data(datatag, save_dir, copynumber_fn, sampling_pct,cellclone_fn)

