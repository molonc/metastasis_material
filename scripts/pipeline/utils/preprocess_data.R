library(tidyverse)
library(stringr)
source(paste0("/home/htran/Projects/farhia_project/rnaseq/pipeline/utils/normalize_utils.R"))


load_sce <- function(metacells_fn, datatag, input_dir=NULL, return_data=T, tag=''){
  
  sample_df <- read.csv(metacells_fn, check.names = F, stringsAsFactors = F)
  # if(tag=='cisplatin' & datatag=='SA535'){
  #   print(tag)
  #   sample_df <- sample_df %>%
  #     dplyr::filter(PDX %in% c('SA535_untreated','SA535_CY'))
  # }
  # if(tag=='CX5461' & datatag=='SA535'){
  #   print(tag)
  #   sample_df <- sample_df %>%
  #     dplyr::filter(PDX %in% c('SA535_untreated','SA535_CX'))
  # }
  
  ## Checking whether all files are included in metasample file
  # fs <- list.files(path=paste0(input_dir,'rnaseq_v6/',datatag,'-v6/'),pattern = '*.rdata',recursive = F)
  # fs <- gsub('.rdata','',fs)
  # fs[!fs %in% sample_df$mouse_id]
  # length(fs)
  # dim(sample_df)
  
  print(dim(sample_df))
  print(sample_df$mouse_id)
  # Just quick check input files
  
  # sample_df <- sample_df %>%
  #   dplyr::filter(treatmentSt %in% grep("*T$",sample_df$treatmentSt, value = T))
  # View(sample_df)
  sce_list <- list()
  c <- 0
  for(s in unique(sample_df$mouse_id)){
    if(datatag=="SA609" & s=="SA609X4XB03083"){
      s <- "SA609X4XB003083"
    }
    sce_fn <- paste0(input_dir,'rnaseq_v6/',datatag,'-v6/', s,'.rdata')
    if(file.exists(sce_fn)){
      sce <- readRDS(sce_fn)
      print(s)
      print(dim(sce))
      if(!'ID' %in% colnames(rowData(sce))){
        stop('Double check input data!!!')
      }
      if(nrow(sce)>0){
        c <- c + 1
        rownames(sce) <- rowData(sce)$ID
        sce_list[c] <- sce
      }
    }else{
      print(paste0('DEBUG: do not exist sample: ',s))
    }
  }
  
  print(length(sce_list))
  # From 33538 to 13214 genes, remove genes with zeroes in 97.5% cells
  sce_combine <- sce_cbind_func(sce_list, rowData(sce), cut_off_overall = 0.025, exprs = c("counts"), 
                                colData_names = colnames(colData(sce)), meta_data=NULL)
  print(dim(sce_combine))
  
  # print("Removing mito, ribo genes: ")
  # mito_genes <- str_detect(rowData(sce_combine)$Symbol, "^MT\\-")
  # sum(mito_genes==TRUE)
  # ribo_genes <- str_detect(rowData(sce_combine)$Symbol, "^RP(L|S)")  # or ^RP[L|S]
  # sum(ribo_genes==TRUE)
  # sce_combine <- sce_combine[(!mito_genes) & (!ribo_genes), ]
  print(dim(sce_combine))
  print(rownames(sce_combine)[1])
  # rownames(sce_combine) <- rowData(sce_combine)$ID  # using ensemble gene id as gene name here
  
  # saveRDS(sce_combine,paste0(input_dir,'rnaseq_v6/',datatag,'-v6/total_sce_treated.rds'))
  clone_df <- load_metaclone(input_dir, sample_df, datatag, tag)
  print(dim(clone_df))
  summary(as.factor(clone_df$clone))
  # clone_df$cell_id <- paste0(clone_df$library_id,'_',clone_df$Barcode)
  
  colnames(sce_combine) <- paste0(sce_combine$sample,'_',sce_combine$Barcode)
  sce_combine$cell_id <- colnames(sce_combine)
  cells_use <- intersect(colnames(sce_combine), clone_df$cell_id)
  print(length(cells_use))
  # sce_combine <- sce_combine[,cells_use]
  sce_combine$clone <- 'unassigned'
  sce_combine$library_id <- get_lib_info(sce_combine)
  rownames(clone_df) <- clone_df$cell_id
  # colData(sce_combine)[cells_use,'library_id'] <- clone_df[cells_use,'library_id']
  colData(sce_combine)[cells_use,'clone'] <- clone_df[cells_use,'clone']
  print(summary(as.factor(sce_combine$clone)))
  print(summary(as.factor(sce_combine$series)))
  sce_combine$series <- gsub(paste0(datatag,'-'),'',sce_combine$series)
  sce_combine$treatmentSt <- get_treatment_status(sce_combine$series)
  print(unique(sce_combine$treatmentSt))
  
  sample_df <- sample_df %>%
    dplyr::select(library_id, mouse_id, batch_info) %>%
    dplyr::rename(batch=batch_info)
  rownames(sample_df) <- sample_df$library_id
  sce_combine$batch <- 'None'
  lids <- unique(sce_combine$library_id)
  print(lids[!lids %in% sample_df$library_id])
  # View(sample_df)
  if(datatag=='SA1035'){
    sce_combine$library_id <- ifelse(sce_combine$library_id=='TENX076','SCRNA10X_SA_CHIP0076_000',sce_combine$library_id)
  }
  
  sce_combine$batch <- sample_df[as.character(sce_combine$library_id),'batch'] 
  print(unique(sce_combine$batch))
  sce_combine$original_clone <- sce_combine$clone
  sce_combine$clone <- NULL
  
  sce_combine$clone <- get_unique_clone_id(sce_combine$original_clone)
  # unique(sce_combine$clone)
  print(summary(as.factor(sce_combine$clone)))
  # Get clone, if cell is assign to 2 clones, get random clone label
  # meta_clones <- data.table::fread('/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/umap_figs/clonealign_labels.csv') %>% as.data.frame()
  # dim(meta_clones)
  # colnames(meta_clones)
  # unique(meta_clones$datatag)
  # meta_clones <- meta_clones %>%
  #   dplyr::filter(datatag==as.name(datatag))
  # unique(meta_clones$clone)
  
  if(tag!=''){
    tag <- paste0(tag,'_')
  }
  saveRDS(sce_combine,paste0(input_dir,'rnaseq_v6/',datatag,'-v6/',datatag, '_',tag, 'total_sce.rds'))
  if(return_data){
    return(sce_combine)
  }
}


load_sce_reduced_version <- function(mouse_ids, datatag, input_dir=NULL, return_data=T, tag=''){
  
  
  sce_list <- list()
  c <- 0
  for(s in unique(mouse_ids)){
    sce_fn <- paste0(input_dir, s,'.rdata')
    if(file.exists(sce_fn)){
      sce <- readRDS(sce_fn)
      print(s)
      print(dim(sce))
      if(nrow(sce)>0){
        c <- c + 1
        sce_list[c] <- sce
      }
    }else{
      print(paste0('DEBUG: do not exist sample: ',s))
    }
  }
  
  print(length(sce_list))
  # From 33538 to 13214 genes, remove genes with zeroes in 97.5% cells
  sce_combine <- sce_cbind_func(sce_list, rowData(sce), cut_off_overall = 0.025, exprs = c("counts"), 
                                colData_names = colnames(colData(sce)), meta_data=NULL)
  print(dim(sce_combine))
  print("Removing mito, ribo genes: ")
  mito_genes <- str_detect(rowData(sce_combine)$Symbol, "^MT\\-")
  sum(mito_genes==TRUE)
  
  ribo_genes <- str_detect(rowData(sce_combine)$Symbol, "^RP(L|S)")  # or ^RP[L|S]
  sum(ribo_genes==TRUE)
  sce_combine <- sce_combine[(!mito_genes) & (!ribo_genes), ]
  
  print(dim(sce_combine))
  rownames(sce_combine) <- rowData(sce_combine)$ID
  print(sce_combine$sample[1])
  print(sce_combine$Barcode[1])
  colnames(sce_combine) <- paste0(sce_combine$sample,'_',sce_combine$Barcode)
  sce_combine$cell_id <- colnames(sce_combine)
  print(summary(as.factor(sce_combine$series)))
  # sce_combine$series <- gsub(paste0(datatag,'-'),'',sce_combine$series)
  # sce_combine$treatmentSt <- get_treatment_status(sce_combine$series)
  
  if(tag!=''){
    tag <- paste0(tag,'_')
  }
  saveRDS(sce_combine,paste0(input_dir, datatag, '_',tag, 'total_sce.rds'))
  if(return_data){
    return(sce_combine)
  }
}


load_metaclone <- function(input_dir, sample_df, datatag, tag){
  clones_list <- list()
  c <- 0
  for(s in unique(sample_df$mouse_id)){
    if(datatag=="SA609" & s=="SA609X4XB03083"){
      s <- "SA609X4XB003083"
    }
    cls_fn <- paste0(input_dir,'rnaseq_v6/',datatag,'-v6/', s,'.csv')
    if(file.exists(cls_fn)){
      sclone_df <- read.csv(cls_fn, check.names = F, stringsAsFactors = F, colClasses=c("clone"="character"))
      # read_csv(input_file,col_types = cols(.default = "?", clone = "c"))
      print(s)
      print(dim(sclone_df))
      if(nrow(sclone_df)>0){
        c <- c + 1
        clones_list[[c]] <- sclone_df
      }
    }else{
      print(paste0('Do not exist sample: ',s))
    }
  }  
  length(clones_list)
  clone_df <- do.call(rbind, clones_list)
  
  clone_df$Sample <- gsub('(.cache/)','',clone_df$Sample)
  clone_df$Sample <- gsub('(/filtered_feature_bc_matrix)','',clone_df$Sample)
  colnames(clone_df)[which(names(clone_df) == "Sample")] <- "library_id"
  colnames(clone_df)[which(names(clone_df) == "id")] <- "mouse_id"
  # head(clone_df)
  print(unique(clone_df$library_id))
  print(dim(clone_df))
  if(tag!=''){
    tag <- paste0(tag,'_')
  }
  write.csv(clone_df, paste0(input_dir,'rnaseq_v6/',datatag,'-v6/',datatag,'_',tag,'total_clones.csv'), row.names = F, quote = F)
  
  # clone_df <- read.csv(paste0(input_dir,'rnaseq_v6/',datatag,'-v6/total_clones.csv'), check.names=F, stringsAsFactors=F)
  # clone_df$clone
  print(dim(clone_df))
  print(head(clone_df))
  return(clone_df)
  
}

load_metaclone_from_directory <- function(input_dir, sample_df, datatag, tag){
  input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/clonealign_plot/clonealign/'
  datatag <- 'SA535'
  fns <- list.files(input_dir, pattern='*.csv')
  fns <- fns[grepl(datatag, fns)]

  clone_df <- tibble::tibble()
  for(s in fns){
    if(datatag=="SA609" & s=="SA609X4XB03083"){
      s <- "SA609X4XB003083"
    }
    cls_fn <- paste0(input_dir, s)
    if(file.exists(cls_fn)){
      sclone_df <- read.csv(cls_fn, check.names = F, stringsAsFactors = F, colClasses=c("clone"="character"))
      # read_csv(input_file,col_types = cols(.default = "?", clone = "c"))
      print(s)
      print(dim(sclone_df))
      if(nrow(sclone_df)>0){
        clone_df <- dplyr::bind_rows(clone_df, sclone_df)
      }
    }else{
      print(paste0('Do not exist sample: ',s))
    }
  }  
  dim(clone_df)
  # clone_df <- do.call(rbind, clones_list)
  
  clone_df$Sample <- gsub('(.cache/)','',clone_df$Sample)
  clone_df$Sample <- gsub('(/filtered_feature_bc_matrix)','',clone_df$Sample)
  colnames(clone_df)[which(names(clone_df) == "Sample")] <- "library_id"
  colnames(clone_df)[which(names(clone_df) == "id")] <- "mouse_id"
  # head(clone_df)
  print(unique(clone_df$library_id))
  print(dim(clone_df))
  if(tag!=''){
    tag <- paste0(tag,'_')
  }
  tag <- ''
  base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
  write.csv(clone_df, paste0(base_dir,'rnaseq_v6/',datatag,'-v6/',datatag,'_',tag,'total_clones.csv'), row.names = F, quote = F)
  
  # clone_df <- read.csv(paste0(input_dir,'rnaseq_v6/',datatag,'-v6/total_clones.csv'), check.names=F, stringsAsFactors=F)
  # clone_df$clone
  print(dim(clone_df))
  print(head(clone_df))
  unique(clone_df$clone)
  return(clone_df)
  
}


# input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
# datatag <- 'SA609'
# metacells_fn <-  paste0(input_dir,'SA609_rna/snakemake_10x/SA609_10x.csv')
# sce <- load_sce(metacells_fn, datatag, input_dir, return_data=T)
# fs <- list.files(path=paste0(input_dir,'rnaseq_v6/',datatag,'-v6/'),pattern = '*.rdata',recursive = F)
# fs <- gsub('.rdata','',fs)
# fs[!fs %in% sample_df$mouse_id]
# output_dir <- paste0(input_dir,'rnaseq_v6/',datatag, '-v6/')
# # dir.create(output_dir)
# print("SCTransform normalization")
# normalize_SCTransform(sce, output_dir, datatag, return_data=F)
# unique(sce_raw$sample)
# dim(sce_raw)
# filtered_sce_raw <- sce_raw[,!grepl('*TU$',sce_raw$treatmentSt)]
# filtered_sce_raw <- filtered_sce_raw[,!filtered_sce_raw$treatmentSt %in% c("UTTTT","UUUUU","U")]
# dim(filtered_sce_raw)
# saveRDS(filtered_sce_raw, paste0(input_dir,'rnaseq_v6/',datatag,'-v6/','subset_total_sce_v3.rds'))

input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
datatag <- 'SA1035'
metacells_fn <-  paste0(input_dir,'SA1035_rna/snakemake_10x_SA1035/metasample_SA1035.csv')
dim(sce)
sce <- load_sce(metacells_fn, datatag, input_dir, return_data=T)
output_dir <- paste0(input_dir,'rnaseq_v6/',datatag, '-v6/')
dir.create(output_dir)
print("SCTransform normalization")
normalize_SCTransform(sce, output_dir, datatag, return_data=F)
tag <- ''
sce <- readRDS(paste0(input_dir,'rnaseq_v6/',datatag,'-v6/',datatag, '_',tag, 'total_sce_v3.rds'))

saveRDS(sce, paste0(input_dir,'rnaseq_v6/',datatag,'-v6/',datatag, '_',tag, 'total_sce_v3.rds'))
length(unique(sce$batch))
dim(sce)
sample_df[sample_df$library_id=='SCRNA10X_SA_CHIP0076_000','library_id'] <- 'TENX076'
rownames(sample_df) <- sample_df$library_id      

sce$batch <- sample_df[as.character(sce$library_id),'batch']
print(unique(sce$batch))
View(sample_df)
print(unique(sce$library_id))

table(sce$batch, sce$library_id)
unique(colData(sce[is.na(sce$batch),])$library_id)
dim(sce)

print("Twice scran normalization")
source(paste0("/home/htran/Projects/farhia_project/rscript/pipeline/utils/twice_scran_normalization_corrected.R"))
sce_scran <- twice_scran_normalize_v2(sce, input_dir, output_dir, paste0(datatag,'_',tag), return_data=F)

print("Seurat normalization")
sce_seurat <- normalize_Seurat(sce, input_dir, output_dir, paste0(datatag,'_',tag), return_data=F)


# input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
# # input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/others_dataset/SA535_cisplatin/'
# datatag <- 'SA535'
# tag <- 'cisplatin'
# metacells_fn <-  paste0(input_dir,'SA535_total_rna_v2/snakemake/metasample_SA535.csv')
# sce <- load_sce(metacells_fn, datatag, input_dir, return_data=T, tag)

# sce <- load_sce_reduced_version(mouse_ids, datatag, input_dir, return_data=T, tag)


# sample_df$mouse_id[!sample_df$mouse_id %in% fs]
# mouse_ids <- fs
# output_dir <- paste0(input_dir,'rnaseq_v6/',datatag, '-v6/')
# dir.create(output_dir)
# output_dir <- input_dir
# print("SCTransform normalization")
# normalize_SCTransform(sce, output_dir, paste0(datatag,'_',tag), return_data=F)



# input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
# datatag <- 'SA535'
# # tag <- 'CX5461'
# tag <- 'total'
# metacells_fn <-  paste0(input_dir,'SA535_total_rna_v2/snakemake/metasample_SA535.csv')
# # sce <- load_sce(metacells_fn, datatag, input_dir, return_data=T, tag)
# sce <- readRDS(paste0(input_dir,'rnaseq_v6/',datatag,'-v6/',datatag, '_',tag, '_total_sce_v3.rds'))
# output_dir <- paste0(input_dir,'rnaseq_v6/',datatag, '-v6/')
# dir.create(output_dir)
# 
# print("Twice scran normalization")
# source(paste0("/home/htran/Projects/farhia_project/rscript/pipeline/utils/twice_scran_normalization_corrected.R"))
# sce_scran <- twice_scran_normalize_v2(sce, input_dir, output_dir, paste0(datatag,'_',tag), return_data=F)
# 
# print("Seurat normalization")
# sce_seurat <- normalize_Seurat(sce, input_dir, output_dir, paste0(datatag,'_',tag), return_data=F)
# 
# print("SCTransform normalization")
# normalize_SCTransform(sce, output_dir, paste0(datatag,'_',tag), return_data=F)




# input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/others_dataset/SA906a/'
# datatag <- 'SA906a'
# tag <- ''
# mouse_ids <- c('SA906_p11a','SA906_p57a')
# sce <- load_sce_reduced_version(mouse_ids, datatag, input_dir, return_data=T, tag)
# dim(sce)
# output_dir <- input_dir
# # dir.create(output_dir)
# print("SCTransform normalization")
# normalize_SCTransform(sce, output_dir, paste0(datatag,'_',tag), return_data=F)


# input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/others_dataset/SA906b/'
# datatag <- 'SA906b'
# tag <- ''
# mouse_ids <- c('SA906_p15b','SA906_p50b')
# sce <- load_sce_reduced_version(mouse_ids, datatag, input_dir, return_data=T, tag)
# dim(sce)
# output_dir <- input_dir
# # dir.create(output_dir)
# print("SCTransform normalization")
# normalize_SCTransform(sce, output_dir, datatag, return_data=F)




