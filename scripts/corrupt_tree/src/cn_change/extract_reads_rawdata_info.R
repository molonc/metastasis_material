"
  Extract GC corrected reads (the copy column) for a list of libraries
"

library(dplyr)
library(tidyr)
library(data.table)
require("ggtree")

f <- function(...) sprintf(...)
u <- function(...) unique(...)


extract_gc_corrected_reads <- function(libs_listpath, rootpath, outpath) {
  dir_names <- readLines(libs_listpath)
  if(!dir.exists(dirname(outpath))){
    dir.create(dirname(outpath), recursive = T)
  }
  full_mat <- NULL
  # t <- c("A108757B")
  for (dn in dir_names) {
    print(dn)
    lib_id <- dn
    file_path <- file.path(rootpath, dn, 'hmmcopy', f('%s_reads.csv', lib_id))
    if (!file.exists(file_path)) {
      print(sprintf('File %s does not exist. Ignoring...', file_path))
      next
    }
    dat <- data.table::fread(file_path, header = TRUE) %>% as_tibble()
    
    mat <- dat %>% dplyr::mutate(V1 = paste0(chr, '_', start, '_', end)) %>% dplyr::select(V1, copy, cell_id) %>% tidyr::pivot_wider(id_cols = c(V1), names_from = cell_id, values_from = copy)
    
    if (is.null(full_mat)) {
      full_mat <- mat
    } else {
      full_mat <- dplyr::left_join(full_mat, mat, by = c('V1'))
    }
  }
  # full_mat <- as.data.frame(full_mat)
  # print(dim(full_mat))
  data.table::fwrite(x = full_mat, file=outpath)
  print(outpath)
  # write.csv(full_mat, file = outpath)
}

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

get_sampleid_from_cell_names <- function(cell_names) {
  gsub('([A-Z]*|[0-9]*)-.*', '\\1', cell_names)
}


summary_filtered_cells <- function(grouping_fn, input_dir, 
                                   foutpath, datatag = NULL,newick_fn=NULL){
  tree <- read.tree(newick_fn)
  # all_cells <- grep('cell', tree$tip.label, value = T)
  all_cells <- tree$tip.label
  print(length(all_cells))
  all_cells <- gsub('cell_','',all_cells)
  print(length(all_cells))
  
  
  grouping_df <- read.csv(grouping_fn, header=T,check.names=F, 
                          stringsAsFactors = F)#row.names = 1, 
  colnames(grouping_df)[which(names(grouping_df) == "grouping")] <- "library_id"
  colnames(grouping_df)[which(names(grouping_df) == "sample_id")] <- "sample"
  res <- tibble()
  for(i in seq(1, length(grouping_df$library_id))) {
    lib_id <- grouping_df$library_id[i]
    print(lib_id)
    # lib_id <- dn
    # file_path <- file.path(rootpath, lib_id, 'hmmcopy', f('%s_reads.csv', lib_id))
    file_path <- file.path(input_dir, f('%s_filtered_states.csv', grouping_df$library_labels[i]))
    
    if (!file.exists(file_path)) {
      print(sprintf('***** File %s does not exist. Ignoring...', file_path))
      next
    }
    
    dat <- read.csv(file_path, header = TRUE, row.names = 1, check.names = F, stringsAsFactors = F) 
    dat <- dat %>% as_tibble()
    # dat <- data.table::fread(file_path, header = TRUE) %>% as_tibble()
    cells <- colnames(dat)
    print(length(cells))
    print(cells[1])
    cells <- cells[cells %in% all_cells]
    metacell_df <- data.frame(cell_id=cells, stringsAsFactors = F)
    metacell_df$library_id <- get_library_id(cells)
    metacell_df$sample <- get_sample_id(cells)
    
    # cells <- cells[grepl(datatag, sample_ids)]
    # metacell_df <- metacell_df %>% inner_join(grouping_df, by = c("library_id","sample"))
    metacell_df <- metacell_df %>% 
      dplyr::filter(library_id==lib_id, sample==grouping_df$sample[i])
    print(dim(metacell_df))
    
    print("")
    tmp <- tibble(library_id = lib_id, sample=grouping_df$sample[i], #PDX = grouping_df$PDX[i],
                  ncells = length(metacell_df$cell_id))
    dim(tmp)
    # tmp <- tibble(library_id = lib_id, ncells = length(cells))
    res <- dplyr::bind_rows(res, tmp)
  }
  dir.create(dirname(foutpath))
  data.table::fwrite(x = res, foutpath)
  res <- as.data.frame(res)
  print(dim(res))
  return(res)
}


extract_cells_features <- function(cellclones, grouping_fn, 
                                   rootpath, outpath, datatag = NULL) {
  cell_clones <- read.csv(cellclones, header=T,check.names=F, stringsAsFactors = F)
  print(head(cell_clones))
  grouping_df <- read.csv(grouping_fn, header=T,check.names=F, 
                          row.names = 1, stringsAsFactors = F)
  colnames(grouping_df)[which(names(grouping_df) == "grouping")] <- "library_id"
  colnames(grouping_df)[which(names(grouping_df) == "sample_id")] <- "sample"
  print(head(grouping_df))
  res <- tibble()
  for(lib_id in grouping_df$library_id) {
    file_path <- file.path(rootpath, lib_id, 'annotation', f('%s_metrics.csv', lib_id))
    
    if (!file.exists(file_path)) {
      print(sprintf('***** File %s does not exist. Ignoring...', file_path))
      next
    }
    metrics <- data.table::fread(file_path, header = TRUE) %>% as.data.frame()
    metrics <- metrics %>%
      dplyr::filter(cell_id %in% cell_clones$cell_id)
    if(dim(metrics)[1]>0){
      res <- dplyr::bind_rows(res, as_tibble(metrics))
      print(dim(metrics))
    }
  }  
  
  # colnames(metrics)
  
  features_use <- c("cell_id","total_mapped_reads",
                    "coverage_depth","total_reads","coverage_breadth","quality")
  total_metrics <- res %>%
    dplyr::select(all_of(features_use))
  
  total_metrics$library_id <- get_library_id(total_metrics$cell_id)
  total_metrics <- total_metrics %>%
    dplyr::group_by(library_id)%>%
    dplyr::summarise(median_nreads_per_cell=median(total_reads),
                     std_nreads_per_cell=sd(total_reads),
                     median_coverage_depth=median(coverage_depth),
                     std_coverage_depth=sd(coverage_depth),
                     median_coverage_breadth=median(coverage_breadth),
                     std_coverage_breadth=sd(coverage_breadth),
                     median_quality=median(quality),
                     std_quality=sd(quality))%>%
    ungroup()
  
  print(summary(total_metrics$median_nreads_per_cell))
  print(summary(total_metrics$median_coverage_depth))
  print(summary(total_metrics$median_coverage_breadth))
  print(summary(total_metrics$median_quality))
  print(dim(total_metrics))
  dir.create(dirname(outpath), recursive = T)
  data.table::fwrite(x = total_metrics, outpath)
}  
  

extract_cells_features_instability <- function(library_ids, results_dir,
                                   rootpath, outpath, datatag = NULL) {
  
  res <- tibble()
  for(lib_id in library_ids) {
    filtered_fn <- file.path(results_dir,'data_thrs_10000',f('%s_filtered_states.csv', lib_id))
    file_path <- file.path(rootpath, lib_id, 'annotation', f('%s_metrics.csv.gz', lib_id))
    
    if (!file.exists(file_path)) {
      print(sprintf('***** File %s does not exist. Ignoring...', file_path))
      next
    }
    metrics <- data.table::fread(file_path, header = TRUE) %>% as.data.frame()
    print(dim(metrics))
    if(file.exists(filtered_fn)){
      filtered_df <- read.csv(filtered_fn, header = TRUE, row.names=1, check.names=F, stringsAsFactors=F)
      dim(filtered_df)
      # colnames(filtered_df)[1]
      # rownames(filtered_df)[1]
      metrics <- metrics %>%
        dplyr::filter(cell_id %in% colnames(filtered_df))
    }else{
      print(f('Attention: do not exist filtered data for library: %s', lib_id))
    }
    print(dim(metrics))
    if(dim(metrics)[1]>0){
      res <- dplyr::bind_rows(res, as_tibble(metrics))
      print(dim(metrics))
    }
  }  
  
  # colnames(metrics)
  
  features_use <- c("cell_id","total_mapped_reads",
                    "coverage_depth","total_reads","coverage_breadth","quality")
  total_metrics <- res %>%
    dplyr::select(all_of(features_use))
  
  total_metrics$library_id <- get_library_id(total_metrics$cell_id)
  total_metrics <- total_metrics %>%
    dplyr::group_by(library_id)%>%
    dplyr::summarise(median_nmapped_reads_per_cell=median(total_mapped_reads),
                     std_nmapped_reads_per_cell=sd(total_mapped_reads),
                     median_nreads_per_cell=median(total_reads),
                     std_nreads_per_cell=sd(total_reads),
                     median_coverage_depth=median(coverage_depth),
                     std_coverage_depth=sd(coverage_depth),
                     median_coverage_breadth=median(coverage_breadth),
                     std_coverage_breadth=sd(coverage_breadth),
                     median_quality=median(quality),
                     std_quality=sd(quality))%>%
    ungroup()
  
  print(summary(total_metrics$median_nreads_per_cell))
  print(summary(total_metrics$median_coverage_depth))
  print(summary(total_metrics$median_coverage_breadth))
  print(summary(total_metrics$median_quality))
  print(dim(total_metrics))
  dir.create(dirname(outpath), recursive = T)
  data.table::fwrite(x = total_metrics, outpath)
}  


# extract number of cells
extract_n_cells <- function(grouping_fn, rootpath, outpath, datatag = NULL) {
  grouping_df <- read.csv(grouping_fn, header=T,check.names=F, 
                          row.names = 1, stringsAsFactors = F)
  colnames(grouping_df)[which(names(grouping_df) == "grouping")] <- "library_id"
  colnames(grouping_df)[which(names(grouping_df) == "sample_id")] <- "sample"
  # unique(grouping_df$PDX)
  # View(grouping_df)
  # grouping_df <- grouping_df %>% 
  #                dplyr::filter(PDX %in% c("SA535_untreated","SA535_cisplatin"))
  # dim(grouping_df)
  stopifnot(!is.null(datatag))
  # "/home/htran/storage/datasets/hakwoo_metastasis/SA919/A108836B/annotation/A108836B_metrics.csv"
  # dir_names <- readLines(libs_listpath)
  res <- tibble()
  for(i in seq(1, length(grouping_df$library_id))) {
    lib_id <- grouping_df$library_id[i]
    print(lib_id)
    # lib_id <- dn
    # file_path <- file.path(rootpath, lib_id, 'hmmcopy', f('%s_reads.csv', lib_id))
    file_path <- file.path(rootpath, lib_id, 'annotation', f('%s_metrics.csv', lib_id))
    
    if (!file.exists(file_path)) {
      print(sprintf('***** File %s does not exist. Ignoring...', file_path))
      next
    }
    dat <- data.table::fread(file_path, header = TRUE) %>% as_tibble()
    cells <- unique(dat$cell_id)
    print(length(cells))
    metacell_df <- data.frame(cell_id=cells, stringsAsFactors = F)
    metacell_df$library_id <- get_library_id(cells)
    metacell_df$sample <- get_sample_id(cells)
    
    # cells <- cells[grepl(datatag, sample_ids)]
    # metacell_df <- metacell_df %>% inner_join(grouping_df, by = c("library_id","sample"))
    metacell_df <- metacell_df %>% 
      dplyr::filter(library_id==lib_id, sample==grouping_df$sample[i])
    print(dim(metacell_df))
    
    print("")
    tmp <- tibble(library_id = lib_id, sample=grouping_df$sample[i], #PDX = grouping_df$PDX[i],
                  ncells = length(metacell_df$cell_id))
    # tmp <- tibble(library_id = lib_id, ncells = length(cells))
    res <- dplyr::bind_rows(res, tmp)
  }
  dir.create(dirname(outpath))
  data.table::fwrite(x = res, outpath)
  res <- as.data.frame(res)
  stat_res <- res %>%
    group_by(sample) %>%
    dplyr::summarise(total_cells=sum(ncells))
  stat_res <- stat_res %>%
    dplyr::summarise(mean_cells=mean(total_cells),
                     sd_cells=sd(total_cells),
                     max_cells=max(total_cells),
                     min_cells=min(total_cells))
  print(paste0("for a total of ",sum(res$ncells)," single cells (mean = ",round(stat_res$mean_cells,2),
               ", sigma = ",round(stat_res$sd_cells,2),", max = ",stat_res$max_cells,", min = ",stat_res$min_cells," per sample)."))
  # for(pdx in unique(res$PDX)){
  #   pdx <- c("SA1035_CY","SA1035_VC")  #"SA1035_Un"
  #   res_tmp <- res %>%
  #              dplyr::filter(PDX %in% pdx)
  #   
  #   res_tmp1 <- res_tmp %>%
  #              group_by(sample) %>%
  #              dplyr::summarise(total_cells=sum(ncells))
  #   res_tmp1 <- res_tmp1 %>%
  #     dplyr::summarise(total=sum(total_cells),
  #                      mean_cells=mean(total_cells),
  #                      sd_cells=sd(total_cells),
  #                      max_cells=max(total_cells),
  #                      min_cells=min(total_cells))
  #   print(pdx)
  #   print(paste0("for a total of ",sum(res_tmp$ncells)," single cells (mean = ",res_tmp1$mean_cells,
  #                ", sigma = ",res_tmp1$sd_cells,", max = ",res_tmp1$max_cells,", min = ",res_tmp1$min_cells," per sample)."))
  #   
  # }
  colnames(res)
  res <- res %>% left_join(grouping_df, by=c("library_id","sample"))
  dim(res)
  dim(grouping_df)
  for(mt in unique(res$mainsite)){
    res_tmp <- res %>%
      dplyr::filter(mainsite==mt)
    
    res_tmp1 <- res_tmp %>%
      group_by(origin) %>%
      dplyr::summarise(total_cells=sum(ncells))
    res_tmp1 <- res_tmp1 %>%
      dplyr::summarise(total=sum(total_cells),
                       mean_cells=mean(total_cells),
                       sd_cells=sd(total_cells),
                       max_cells=max(total_cells),
                       min_cells=min(total_cells))
    print(paste0("for a total of ",sum(res_tmp$ncells)," ",tolower(mt)," single cells ; (mean = ",round(res_tmp1$mean_cells,2),
                 ", sigma = ",round(res_tmp1$sd_cells,2),", max = ",res_tmp1$max_cells,", min = ",res_tmp1$min_cells," per origin)."))
    
  }
  
  for(mt in unique(res$mainsite)){
    res_tmp <- res %>%
      dplyr::filter(mainsite==mt)
    
    res_tmp1 <- res_tmp %>%
      group_by(sample) %>%
      dplyr::summarise(total_cells=sum(ncells))
    res_tmp1 <- res_tmp1 %>%
      dplyr::summarise(total=sum(total_cells),
                       mean_cells=mean(total_cells),
                       sd_cells=sd(total_cells),
                       max_cells=max(total_cells),
                       min_cells=min(total_cells))
    print(paste0("for a total of ",sum(res_tmp$ncells)," ",tolower(mt)," single cells ; (mean = ",round(res_tmp1$mean_cells,2),
                 ", sigma = ",round(res_tmp1$sd_cells,2),", max = ",res_tmp1$max_cells,", min = ",res_tmp1$min_cells," per sample)."))
    
  }
  
  for(pd in unique(res$pdxid)){
    res_tmp <- res %>%
      dplyr::filter(pdxid==pd)
    
    res_tmp1 <- res_tmp %>%
      group_by(sample) %>%
      dplyr::summarise(total_cells=sum(ncells))
    res_tmp1 <- res_tmp1 %>%
      dplyr::summarise(total=sum(total_cells),
                       mean_cells=mean(total_cells),
                       sd_cells=sd(total_cells),
                       max_cells=max(total_cells),
                       min_cells=min(total_cells))
    print(paste0("For a total of ",sum(res_tmp$ncells)," single cells in ",pd,": (mean = ",round(res_tmp1$mean_cells,2),
                 ", sigma = ",round(res_tmp1$sd_cells,2),", max = ",res_tmp1$max_cells,
                 ", min = ",res_tmp1$min_cells," per sample)."))
    
    
  }
  
  
  return(res)
}


# libs_listpath <- '/home/htran/Projects/farhia_project/fitness_material/cnv_data/SA535_all_libs.txt'
#rootpath <- '/home/htran/storage/datasets/drug_resistance_DLP/SA535'
# rootpath <- '/home/htran/storage/datasets/drug_resistance_DLP/SA535'
# # outpath <- '/home/htran/storage/datasets/drug_resistance_DLP/SA535_Sohrab/%s/'
# # '/home/salehi/drug_data/SA535/all_cells_unfiltered_map_mat.csv.gz'
# results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535/SA535_total_v2/'
# # outpath <- '/home/htran/storage/datasets/drug_resistance_DLP/SA535_Sohrab/SA535_all_cells_unfiltered_map_mat.csv.gz' #gc_corrected_reads_SA535.csv
# grouping_fn <- paste0(results_dir,'library_groupings.csv')
# outpath <- '/home/htran/storage/datasets/drug_resistance_DLP/SA535_Sohrab/SA535_raw_data_info.csv'

# rootpath <- '/home/htran/storage/datasets/drug_resistance_DLP/SA1035/'
# results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_new_encode/SA1035_Tyler_v2/'
# grouping_fn <- paste0(results_dir,'library_groupings.csv')
# outpath <- '/home/htran/storage/datasets/drug_resistance_DLP/SA535_Sohrab/SA1035_raw_data_info.csv'
# extract_n_cells(grouping_fn, rootpath, outpath, datatag = 'SA1035')

# stat <- read.csv(outpath, check.names = F)
# dim(stat)  
# 
# View(stat)
# nb_cells <- sum(stat$ncells)
# mean(stat$ncells)
# sd(stat$ncells)



# SA919 metastasis project
rootpath <- '/home/htran/storage/datasets/hakwoo_metastasis/SA919/'
results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_wholedata/'
grouping_fn <- paste0(results_dir,'library_groupings.csv')
cellclones <- paste0(results_dir,'cell_clones.csv')
newick_fn <- paste0(results_dir,'tree.newick')
outpath <- '/home/htran/storage/datasets/hakwoo_metastasis/SA919/stat/SA919_raw_data_metrics.csv'
input_dir <- paste0(results_dir,'data/')
foutpath <- '/home/htran/storage/datasets/hakwoo_metastasis/SA919/stat/SA919_filtered_data_metrics_info.csv'
datatag = 'SA919'
# res <- extract_n_cells(grouping_fn, rootpath, outpath, datatag)
fres <- summary_filtered_cells(grouping_fn, input_dir, foutpath, datatag, newick_fn)
res <- fres

extract_cells_features(cellclones, grouping_fn, rootpath, outpath, datatag = NULL) 
  
# SA535 metastasis project
rootpath <- '/home/htran/storage/datasets/hakwoo_metastasis/SA535/'
results_dir <- '/home/htran/storage/datasets/metastasis_results/SA535_new_encoding/SA535_wholedata_v2/'
newick_fn <- paste0(results_dir,'tree.newick')
cellclones <- paste0(results_dir,'cell_clones.csv')
grouping_fn <- paste0(results_dir,'library_groupings.csv')
outpath <- '/home/htran/storage/datasets/hakwoo_metastasis/SA535/stat/SA535_filtered_data_metrics_info.csv'
input_dir <- paste0(results_dir,'data/')
foutpath <- '/home/htran/storage/datasets/hakwoo_metastasis/SA535/stat/SA535_filtered_data_info.csv'
datatag = 'SA535'
res <- extract_n_cells(grouping_fn, rootpath, outpath, datatag)
fres <- summary_filtered_cells(grouping_fn, input_dir, foutpath, datatag, newick_fn)

extract_cells_features(cellclones, grouping_fn, rootpath, outpath, datatag = NULL) 



rootpath <- '/home/htran/storage/datasets/damian_DLP'
results_dir <- '/home/htran/storage/gm_instability_results/HEK293'
outpath <- file.path(results_dir,'HEK293_qc.csv')
datatag <- 'HEK293'
library_ids <- c('A110616A','A110617A','A118346B','A118825A','A118866A','A98230A')
extract_cells_features_instability(library_ids, results_dir,
                                  rootpath, outpath, datatag = NULL)


# SA535 combined tree raw data
# We acquired sc-WGS on 7 consecutively transplanted timepoints (UX4, X5, X6, X7, X8, X9, X10) 
# for a total of 44,618 single cells (mean = 1274.8, $\sigma$ = 300.8 per library).

# extract_gc_corrected_reads(libs_listpath = libs_listpath, 
#                            rootpath = rootpath, 
#                            outpath = outpath)