suppressPackageStartupMessages({
  # library("optparse")##BiocManager::install("argparse")
  library("dplyr")##BiocManager::install("dplyr")
  # library("ggplot2")
  library("SingleCellExperiment")##BiocManager::install("SingleCellExperiment")
})


sce_cbind_func_v3 <- function(sce_list, min_cells = 50, 
                              exprs = c("counts", "logcounts","normcounts"), 
                              colData_names = NULL, save_raw=T, save_dir='',tag='SA') { #, meta_data=NULL
  n_batch <- length(sce_list)
  # method = "intersect"
  
  assay_list <- list()
  for (i in seq_len(length(exprs))) {
    assay_list[[i]] <- do.call(cbind, lapply(sce_list, 
                                             function(y) assay(y, exprs[i])))
  }
  names(assay_list) <- exprs
  if(is.null(colData_names)){ 
    # print('Combining all meta cells features in colData(sce)')
    colData_names <- colnames(colData(sce_list[[1]]))
    for(i in rep(2:length(sce_list),1)){
      colData_names <- intersect(colData_names, colnames(colData(sce_list[[i]])))
    }
    # colData_list <- do.call(DelayedArray::rbind, 
    #                         lapply(sce_list, function(y) colData(y)))
  }
  colData_list <- do.call(DelayedArray::rbind, 
                          lapply(sce_list, function(y) colData(y)[, colData_names, drop = FALSE]))
  # genes_df <- as.data.frame(rowData(sce_list[[1]]))
  # genes_df <- genes_df %>%
  #   dplyr::rename(geo_strand=strand,
  #                 geo_chr=chr,
  #                 geo_start=start,
  #                 geo_end=end)%>%
  #   as.matrix()
  sce_combine <- SingleCellExperiment::SingleCellExperiment(assay = assay_list, 
                                                            colData = colData_list,
                                                            rowData=rowData(sce_list[[1]]))
  
  
  
  # sce_combine$mouse_id <- meta_data$mouse_id
  # sce_combine$passage <- meta_data$passage
  # sce_combine$treatmentSt <- meta_data$treatmentSt
  
  print(paste0("Dim sce combine: ",dim(sce_combine)))
  # print(colnames(colData(sce_combine)))
  print(paste0("sce combine assay name: ",assayNames(sce_combine)))
  if(save_raw){
    saveRDS(sce_combine, file = paste0(save_dir,"filtered_combined_",tag,".rds"))
  }
  if(min_cells > 0){
    nonzero_cbind <- DelayedArray::rowSums(assay(sce_combine, exprs[1]) > 0)
    sce_combine <- sce_combine[names(nonzero_cbind[nonzero_cbind >= min_cells]), ]
  }
  
  return(sce_combine)
}

load_data <- function(meta, input_dir, output_dir=NULL, tag=''){
  if (is.null(output_dir)){
    output_dir <- input_dir
  }
  
  if (!file.exists(output_dir)){
    dir.create(output_dir)
  }
  meta <- as.data.frame(meta)
  rownames(meta) <- meta$library_id
  sce_list <- list()
  c <- 0
  for (f in unique(meta$library_id)){
    # norm_fn <- paste0(base_dir,f,".rdata")
    norm_fn <- paste0(input_dir,f,".rds")
    if(file.exists(norm_fn)){
      sce_tmp <- readRDS(norm_fn)
      sce_tmp$library_id <- meta[f,'library_id']
      # print(meta[f,'library_id']==f)
      c <- c + 1
      if(c==1){
        cols_use <- colnames(colData(sce_tmp))
      }else{
        cols_use <- intersect(cols_use, colnames(colData(sce_tmp)))
      }
      print(paste0('DEBUG: ',f,'  ',dim(sce_tmp)[1],' ',dim(sce_tmp)[2]))
      sce_list[[f]] <- sce_tmp
    } else{
      print(paste0('*** Error loading files: ',f))
    }  
    
  }  
  print(length(cols_use))
  
  ## cut_off_overall = 0: do not filter any gene
  ## cut_off_overall = 0.025: filter genes with zero counts values in 97.5% of total cells
  cols_use <- c('Barcode','Sample','library_id')
  # min_cells = 50 at least 50 cells per a given gene --> keep this gene, gene filtering 
  sce_combine <- sce_cbind_func_v3(sce_list, min_cells = 50, exprs = c("counts"), 
                                   colData_names = cols_use, save_raw=F, save_dir=output_dir, tag) 
  print(dim(sce_combine))
  # sce_combine <- readRDS(paste0(output_dir,'combined_total_genes_filtered.rds'))
  # Twice scater normalization
  # output are saved in logcounts exp values
  print(("Normalizing total data..."))
  
  ## cell name = library_id + barcode
  # if('Sample' %in% colnames(colData(sce_combine))){
  #   sce_combine$library_id <- gsub('.cache/','',sce_combine$Sample)
  #   sce_combine$library_id <- gsub('/filtered_feature_bc_matrix','',sce_combine$library_id)
  #   print(unique(sce_combine$library_id))
  #   colnames(sce_combine) <- paste0(sce_combine$library_id,'_',sce_combine$Barcode)
  # }
  # rownames(sce_combine)
  colnames(sce_combine) <- paste0(sce_combine$library_id,'_',sce_combine$Barcode)
  # print(rowData(sce_combine_raw)$Symbol[1:3])
  
  saveRDS(sce_combine, file=paste0(output_dir,"combined_",tag,".rds"))
  
  return(sce_combine)
}

main(){
  meta <- data.table::fread('/home/htran/Projects/hakwoo_project/metastasis_material/materials/10x/SA535_10x_metadata_passage_X4_full_infos.csv')
  
  # meta$nb_cells_introns
  # meta <- meta %>% 
  #   dplyr::filter(mouse_id=='SA535X4XB05462')
  meta <- meta %>% 
    dplyr::filter(pdxid!='M2364')
  colnames(meta)
  View(meta)
  # meta1 <- meta %>% 
  #   dplyr::arrange(pdxid, desc(Grouping)) %>% 
  #   dplyr::select(pdxid, Grouping, Site_origin, nb_cells_introns, Noted)
  # dim(meta1)
  # data.table::fwrite(meta1, paste0(input_dir,'qc_reports/nb_cells_introns.csv'))
  # # meta$library_id
  # View(meta1)
  input_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/filtered_introns/'
  output_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/normalized_introns/'
  tag <- 'SA535' 
}




