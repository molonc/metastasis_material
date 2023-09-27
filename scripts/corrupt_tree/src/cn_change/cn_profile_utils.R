suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
}) 

compute_dist_mat <- function(datatag, mg_mat, results_dir, use_hamming = FALSE) {
  # print("Testing")
  if(is.null(mg_mat)){
    # print("Testing 1")
    # return(NULL)
    stop('Check input data!!!')
  } else if(nrow(mg_mat)==1){
    # print("Testing 2")
    # return(matrix(0, nrow = 1, ncol = 1))
    stop('Only one row in mtx, check input data!!!')
  } else{
    # print("Testing 3")
    clone_lbs <- colnames(mg_mat)
    out_mtx <- tibble::tibble()
    n_clones <- ncol(mg_mat)
    dist_to_median <- matrix(NA, nrow = n_clones, ncol = n_clones)
    for (j in seq(n_clones)) {
      for (i in seq(n_clones)) {
        if (i<j) {  #(i + j)<=(n_clones+1)
          if (use_hamming) {
            distance_type='Hamming'
            dist_to_median[j,i] <- round(mean(mg_mat[ ,i] != mg_mat[ ,j]),3)
          } else {
            distance_type='Manhattan'
            dist_to_median[j,i] <- round(mean(abs(mg_mat[ ,i] - mg_mat[ ,j])),3) #The Manhattan distance as the sum of absolute differences
          }
          distji <- c(clone_lbs[j],clone_lbs[i], dist_to_median[j,i])
          names(distji) <- c('SourceClone','TargetClone','CNA_Distance')
          out_mtx <- dplyr::bind_rows(out_mtx,distji)
        }
      }
    }
    rownames(dist_to_median) <- colnames(mg_mat)
    colnames(dist_to_median) <- colnames(mg_mat)
    # data.table::fwrite(as.data.frame(dist_to_median), paste0(results_dir,datatag,'_cn_distance_',distance_type,'.csv.gz'), row.names=T)
    data.table::fwrite(out_mtx, paste0(results_dir,datatag,'_cn_distance_',distance_type,'_output.csv.gz'))
    return(list(dist_to_median=dist_to_median,out_mtx=out_mtx))
  }
  
}

get_median_genotype_v3 <- function(copynumber_fn, 
                                   datatag, save_dir,
                                   cell_clones=NULL, library_grouping_fn=NULL){
  if(is.null(cell_clones)){
    cellclone_fn <- paste0(results_dir,'cell_clones.csv')  
    cell_clones <- data.table::fread(cellclone_fn)
  }
  if(is.null(library_grouping_fn)){
    library_grouping_fn <- paste0(results_dir,'library_groupings.csv')  
  }
  if(!dir.exists(save_dir)){
    dir.create(save_dir)
  }
  
  # copynumber <- read.csv(copynumber_fn, header=T, row.names = 1, check.names = F,stringsAsFactors = FALSE)
  copynumber <- as.data.frame(data.table::fread(copynumber_fn))
  rownames(copynumber) <- copynumber$V1
  copynumber$V1 <- NULL
  # cell_clones contain 2 columns of cell_id, and clone_id
  # ex:           cell_id                      clone_id
  # 1   SA535X4XB02498-A98163A-R09-C11          C
  # cell_clones <- data.table::fread(cellclone_fn) %>% as.data.frame()
  # dim(cell_clones)
  # metasample <- data.table::fread(library_grouping_fn) %>% as.data.frame()
  # dim(metasample)
  # head(metasample)
  cell_clones <- cell_clones %>%
    dplyr::filter(!clone_id %in% c('None','unassigned'))
  
  
  # metasample <- metasample %>%
  #   dplyr::rename(library_id=grouping) %>%
  #   dplyr::select(library_id, mainsite, sample_id)

  cell_clones$library_id <- get_library_id(cell_clones$cell_id)
  cell_clones$sample_id <- get_sample_id(cell_clones$cell_id)
  # cell_clones <- cell_clones %>% left_join(metasample, by=c("library_id","sample_id"))
  
  # meta_cells <- cell_clones
  # dim(meta_cells)
  # meta_cells <- meta_cells %>%
  #   dplyr::group_by(clone_id, mainsite) %>%
  #   dplyr::summarise(nb_cells=n()) %>%
  #   dplyr::ungroup()
  # View(meta_cells)
  # meta_cells <- meta_cells %>%
  #   dplyr::filter(nb_cells>=20 & clone_id !='None')
  # dim(meta_cells)
  # colnames(meta_cells)
  # write.csv(meta_cells, paste0(save_dir,'meta_cells.csv'), quote = F, row.names = F)
  
  copynumber <- copynumber[,colnames(copynumber)[colnames(copynumber) %in% cell_clones$cell_id]]
  copynumber$chr_desc <- rownames(copynumber)
  cnv <- copynumber %>%
    pivot_longer(!chr_desc, names_to = "cell_id", values_to = "copy_number")
  
  print(dim(cnv))
  cnv <- cnv %>% inner_join(cell_clones, by=c("cell_id"))
  # length(unique(cnv$cell_id))
  # dim(cnv)
  # colnames(cnv)
  # cnv <- cnv %>% inner_join(meta_cells, by=c("clone_id","mainsite"))
  # dim(cnv)
  # 
  # cnv$clone_label <- paste0(cnv$clone_id,'_',cnv$mainsite)
  # length(unique(cnv$clone_label))
  
  # get median genotype 
  print("Get median genotype")
  # stat_cnv <- cnv %>%
  #   dplyr::group_by(clone_id, chr_desc) %>%
  #   dplyr::summarise(median_cn=median(copy_number),
  #                    mode_cn=calc_mode(copy_number))
  # 
  stat_cnv <- cnv %>%
    dplyr::group_by(clone_id, chr_desc) %>%
    dplyr::summarise(median_cn=median(copy_number)) %>%
    dplyr::ungroup()
  print(dim(stat_cnv))
  # median_cnv <- stat_cnv %>%
  #   dplyr::select(chr_desc, median_cn, clone_id, clone_label) %>%
  #   group_by(clone_label) 
  # head(stat_cnv)
  
  median_cnv_pivot <- stat_cnv %>%
    pivot_wider(names_from = clone_id, values_from = median_cn)
  
  median_cnv_pivot <- as.data.frame(median_cnv_pivot)
  dim(median_cnv_pivot)
  data.table::fwrite(median_cnv_pivot, paste0(save_dir,datatag,'_median_cnv.csv.gz'))
  # data.table::fwrite(stat_cnv, paste0(save_dir,datatag,'_stat_median_cnv.csv.gz'))
  return(median_cnv_pivot)
}  
get_library_id <- function(cell_ids) {
  
  labels <- lapply(strsplit(cell_ids, "-"), function(x) {
    return(x[2])
  })
  return(as.character(labels))
}
get_sample_id <- function(cell_ids) {
  labels <- lapply(strsplit(cell_ids, "-"), function(x) {
    return(x[1])
  })
  return(as.character(labels))
}
