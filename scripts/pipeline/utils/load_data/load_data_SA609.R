library(SingleCellExperiment)
# library(DelayedArray) # need to install this package

load_data <- function(meta_cells, base_dir, tag=''){
  sce_list <- list()
  c <- 0
  for (f in unique(meta_cells$mouse_id)){
    # norm_fn <- paste0(base_dir,f,".rdata")
    norm_fn <- paste0(base_dir,f,".rdata")
    if(file.exists(norm_fn)){
      sce_tmp <- readRDS(norm_fn)
      c <- c + 1
      if(c==1){
        cols_use <- colnames(colData(sce_tmp))
      }else{
        cols_use <- intersect(cols_use, colnames(colData(sce_tmp)))
      }
    } else{
      warning('Error loading files')
      print(f)
    }  
    print(paste0('DEBUG: ',f,'  ',dim(sce_tmp)[1],' ',dim(sce_tmp)[2]))
    sce_list[[f]] <- sce_tmp
    
  }  
  length(cols_use)
  
  ## cut_off_overall = 0: do not filter any gene
  ## cut_off_overall = 0.025: filter genes with zero counts values in 97.5% of total cells
  sce_combine <- sce_cbind_func_v2(sce_list, cut_off_overall = 0, exprs = c("counts"), 
                                   colData_names = cols_use, save_raw=F, save_dir=base_dir, tag) 
  print(dim(sce_combine))
  # sce_combine <- readRDS(paste0(output_dir,'combined_total_genes_filtered.rds'))
  # Twice scater normalization
  # output are saved in logcounts exp values
  print(("Normalizing total data..."))
  
  ## cell name = library_id + barcode
  if('Sample' %in% colnames(colData(sce_combine))){
    sce_combine$library_id <- gsub('.cache/','',sce_combine$Sample)
    sce_combine$library_id <- gsub('/filtered_feature_bc_matrix','',sce_combine$library_id)
    print(unique(sce_combine$library_id))
    colnames(sce_combine) <- paste0(sce_combine$library_id,'_',sce_combine$Barcode)
  }
  # rownames(sce_combine)
  # print(rowData(sce_combine_raw)$Symbol[1:3])
  
  saveRDS(sce_combine, file=paste0(base_dir,"combined_",tag,".rds"))
  
  return(sce_combine)
}

sce_cbind_func_v2 <- function(sce_list, cut_off_overall = 0.01, exprs = c("counts", "logcounts","normcounts"), 
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
    saveRDS(sce_combine, file = paste0(save_dir,"raw_combined_",tag,".rds"))
  }
  if(cut_off_overall > 0){
    zero_cbind <- DelayedArray::rowMeans(assay(sce_combine, exprs[1]) == 0)
    sce_combine <- sce_combine[names(zero_cbind[zero_cbind <= (1 - cut_off_overall)]), ]
  }
  
  return(sce_combine)
}

base_dir <- 'yourdir/SA609-v6/'
base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/SA609-v6/'
meta_samples <- data.table::fread(paste0(base_dir,'metadata/SA609_main_line_10x.csv'))
dim(meta_samples)
## Excluding mixture samples first. 
## These samples are sequenced at the revision time, do not totally match with the first sequencing time
excluded_samples <- c('SA609X5XB03235','SA609X4XB03084')
meta_samples <- meta_samples[!meta_samples$mouse_id %in% excluded_samples,]

## Concatenating all samples data to have a total cells file
sce_combine <- load_data(meta_samples, base_dir, tag='SA609')
dim(sce_combine)


## Get counts data
mtx <- as.matrix(counts(sce))
meta_cells <- as.data.frame(colData(sce))
meta_genes <- as.data.frame(rowData(sce))

data.table::fwrite(mtx, paste0(base_dir,'raw_SA609/','SA609_main_line_mtx.csv.gz'))
data.table::fwrite(meta_cells, paste0(base_dir,'raw_SA609/','SA609_main_line_meta_cells.csv.gz'))
data.table::fwrite(meta_genes, paste0(base_dir,'raw_SA609/','SA609_main_line_meta_genes.csv.gz'))

## TO DO: load these matrices into scanpy
# adata = sc.AnnData(mtx, obs=meta_cellls, var=meta_genes)






