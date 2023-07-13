suppressPackageStartupMessages({
  # require("optparse")
  # require("scater")
  # require("argparse")
  require("SingleCellExperiment")
  require("stringr")
  require("dplyr")
  require("data.table")
  # require("tidyverse")
})


## Run TreeAlign joint tumour samples: a given sample rnaseq align to given sample in dlp+ clonal labels
## scrnaseq: get libraries from same sample --> combine them into one mtx 
## dlp+: get clone with the minimum of 15 cells per clone for the given sample --> clonal csv profile
## Run TreeAlign
## Get normalized mtx using sctransform (normalized counts, and log normalized)
## Get umap results
## Plot umap for each mouse
## Check if cells from same clone are located at the same neighbourhood areas. 
##
## 
##
##
##
## 
##
##
##

get_rnaExp <- function(){
  input_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/filtered_introns/'
  meta <- data.table::fread('/home/htran/Projects/hakwoo_project/metastasis_material/materials/10x/SA535_10x_metadata_passage_X4_full_infos.csv')
  sids <- unique(meta$mouse_id)
  ## Get cells ids 
  cells_df <- tibble::tibble()
  for(sid in sids){
    print(sid)
    meta_sid <- meta %>% 
      dplyr::filter(mouse_id==sid)
    
    print(meta_sid$library_id)
    
    sce_list <- list()
    for(lid in unique(meta_sid$library_id)){
      tmp_fn <- paste0(input_dir, lid, '.rds')
      if(file.exists(tmp_fn)){
        sce <- readRDS(tmp_fn)
        tmp <- colData(sce) %>%
          as.data.frame() %>%
          dplyr::select(Barcode, Sample) #%>%
          # dplyr::rename(library_id=Sample)
        print(dim(tmp))
        tmp$mouse_id <- sid
        tmp$library_id <- lid
        cells_df <- dplyr::bind_rows(cells_df, tmp)  
      }else{
        print('Do not exist sce file, check this library!!!')
        print(tmp_fn)
        print(sid)
      }
      
    }
  }  
  dim(cells_df)
  sum(cells_df$Sample==cells_df$library_id)
  cells_df$Sample <- NULL
  data.table::fwrite(cells_df, paste0(save_dir,'meta_cells.csv.gz'))
  
  input_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/filtered_introns/'
  save_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/clonealign_introns/'
  
  meta <- data.table::fread('/home/htran/Projects/hakwoo_project/metastasis_material/materials/10x/SA535_10x_metadata_passage_X4_full_infos.csv')
  
  sids <- unique(meta$mouse_id)
  # paste(sids, collapse = ',')
  View(meta)
  cnv <- data.table::fread(paste0(save_dir, 'mapped_wholedata_SA535_v2.csv.gz'))
  dim(cnv)
  head(cnv)
  for(sid in sids){
    print(sid)
    meta_sid <- meta %>% 
      dplyr::filter(mouse_id==sid)
    
    print(meta_sid$library_id)
    
    sce_list <- list()
    for(lid in unique(meta_sid$library_id)){
      tmp_fn <- paste0(input_dir, lid, '.rds')
      if(file.exists(tmp_fn)){
        sce <- readRDS(tmp_fn)
        sce_list[[lid]] <- sce  
      }else{
        print('Do not exist sce file, check this library!!!')
        print(tmp_fn)
        print(sid)
      }
      
    }
    assay_counts <- do.call(cbind, lapply(sce_list, function(y) assay(y, 'counts')))
    colData_names <- c('Barcode','Sample')
    
    colData_list <- do.call(DelayedArray::rbind, 
                            lapply(sce_list, function(y) colData(y)[, colData_names, drop = FALSE]))
    sce_combine <- SingleCellExperiment::SingleCellExperiment(assays = list(counts=assay_counts), 
                                                              colData = colData_list)
    print(dim(sce_combine))
    
    colnames(sce_combine) <- sce_combine$Barcode
    if(is.null(rownames(sce_combine))){
      rownames(sce_combine) <- rowData(sce_combine)$ID
    }
    # rownames(sce_combine)[1]
    # saveRDS(sce_combine, paste0(save_dir, sid, '_introns_sce.rds'))
    # print(assayNames(sce_combine))
    
    mtx <- as.data.frame(counts(sce_combine))
    # dim(mtx)
    print(rownames(mtx)[1])
    print(colnames(mtx)[1])
    
    
    genes_used <- intersect(cnv$ensembl_gene_id, rownames(mtx))
    print(length(genes_used))
    mtx <- mtx[genes_used,]
    
    data.table::fwrite(mtx, paste0(save_dir, sid, '_expr_introns.csv.gz'), row.names = T)
    
    
  }
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

get_cnv_sid <- function(min_cells_per_clone=10){
  save_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/clonealign_introns/'
  
  cnv <- data.table::fread('/home/htran/storage/raw_DLP/metastasis_DLP/SA535/reads_clones_v3/mapped_wholedata_SA535_v2.csv.gz')
  dim(cnv)
  colnames(cnv)
  cnv1 <- cnv %>%
    dplyr::select(-ensembl_gene_id)
  
  
  var_genes <- sapply(seq(1:dim(cnv1)[1]), function(x){
    var(as.numeric(cnv1[x,]))
  })
  
  sum(var_genes>0)
  length(var_genes)
  meta_genes <- tibble::tibble(ensembl_gene_id=cnv$ensembl_gene_id, variance=var_genes)
  meta_genes <- meta_genes %>%
    dplyr::filter(variance>0)
  dim(meta_genes)
  
  cnv <- cnv %>%
    dplyr::filter(ensembl_gene_id %in% meta_genes$ensembl_gene_id)
  dim(cnv)
  data.table::fwrite(cnv, paste0(save_dir, 'mapped_wholedata_SA535_v2.csv.gz'))
  
  cells_clone_fn <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/dlp_trees/SA535/cell_clones.csv.gz'
  cells_clone <- data.table::fread(cells_clone_fn)
  cells_clone$library_id <- get_library_id(cells_clone$cell_id)
  cells_clone$sample_id <- get_sample_id(cells_clone$cell_id)
  unique(cells_clone$sample_id)
  
  meta <- data.table::fread('/home/htran/Projects/hakwoo_project/metastasis_material/materials/10x/SA535_10x_metadata_passage_X4_full_infos.csv')
  
  sids <- unique(meta$mouse_id)
  # 
  dim(cnv)
  ## for shell script
  writeLines(sids, paste0(save_dir, "samples_run_list.txt"))  
  save_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/clonealign_introns/input_data/'
  unassigned <- c('None','Un','unassigned','Unassigned')
  for(sid in sids){
    print(sid)
    
    clone_used <- cells_clone %>% 
      dplyr::filter(sample_id==sid & !clone_id %in% unassigned) %>%
      dplyr::group_by(clone_id) %>%
      dplyr::summarise(nb_cells=n()) %>%
      dplyr::filter(nb_cells>=min_cells_per_clone)
    
    if(dim(clone_used)[1]>0){
      cnv1 <- cnv %>%
        dplyr::select(all_of(c('ensembl_gene_id', clone_used$clone_id)))
      colnames(cnv1) <- c('ensembl_gene_id', paste0('cell_',clone_used$clone_id))
      data.table::fwrite(cnv1, paste0(save_dir, sid, '_clones_cnv.csv.gz'))  
      
      clones_df <- data.frame(cell_id=paste0('cell_',clone_used$clone_id),
                              clone_id=clone_used$clone_id)
      
      data.table::fwrite(clones_df, paste0(save_dir,sid,'_cell_clones.csv.gz'))
      
    }
    
  }  
}


get_clone_labels <- function(){
  meta <- data.table::fread('/home/htran/Projects/hakwoo_project/metastasis_material/materials/10x/SA535_10x_metadata_passage_X4_full_infos.csv')
  # View(meta)
  sids <- unique(meta$mouse_id)
  sids
  base_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/clonealign_introns/'
  save_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/clonealign_introns/TreeAlign_result_v1/'
  fns <- list.files(pattern = '*_clone_assign.csv.gz', path = save_dir)
  fns_sid <- gsub('_clone_assign.csv.gz','',fns)
  # sids[!sids %in% fns_sid]  
  
  ## Note: remove SA535X4XB05989 from the list, small number of cells
  fns_sid <- fns_sid[fns_sid != 'SA535X4XB05989']
  fns_sid
  clonal_df <- tibble::tibble()
  for(s in fns_sid){
    assigned <- data.table::fread(paste0(save_dir,s,'_clone_assign.csv.gz'))
    assigned$sample_id <- s
    print("----------------------")
    print(s)
    print(summary(as.factor(assigned$clone_id)))
    clonal_df <- dplyr::bind_rows(clonal_df, assigned)
  }
  dim(clonal_df)
  data.table::fwrite(clonal_df, paste0(save_dir,'total_clonal_assignment.csv.gz'))
  colnames(meta)
  head(meta)
  cells_df <- data.table::fread(paste0(save_dir,'meta_cells.csv.gz'))
  dim(cells_df)
  colnames(cells_df)
  
  meta1 <- meta %>%
    dplyr::select(mouse_id, pdxid, Site_origin, Grouping, library_id)
  cells_df <- cells_df %>%
    dplyr::left_join(meta1, by = c('mouse_id','library_id'))
  
  clonal_df <- clonal_df %>%
    dplyr::left_join(cells_df, by = c('sample_id'='mouse_id','cell_id'='Barcode'))
  dim(clonal_df)
  summary(as.factor(clonal_df$Grouping))
  
  colnames(clonal_df)
  ## Same clone across different sites: primary, met
  stat1 <- clonal_df %>%
    dplyr::filter(clone_id!='' & !is.na(Grouping)) %>%
    dplyr::group_by(clone_id, Grouping) %>%
    dplyr::summarise(nb_cells=n())%>%
    dplyr::filter(nb_cells>=10)
  dim(stat1)    
  # View(stat1)
  
  ## Same clone across different mice : M1, M2, M3, M6, M64,..
  stat2 <- clonal_df %>%
    dplyr::filter(clone_id!='' & !is.na(Grouping)) %>%
    dplyr::group_by(clone_id, pdxid) %>%
    dplyr::summarise(nb_cells=n()) %>%
    dplyr::filter(nb_cells>=10)
  dim(stat2)    
  # View(stat2)
  
  ## Same clone across different mice : M1, M2, M3, M6, M64,..
  stat3 <- clonal_df %>%
    dplyr::filter(clone_id!='' & !is.na(Grouping)) %>%
    dplyr::group_by(clone_id, Site_origin) %>%
    dplyr::summarise(nb_cells=n()) %>%
    dplyr::filter(nb_cells>=10)
  dim(stat3)    
  View(stat3)
  stat_dir <- paste0(base_dir,'summary_results/')
  # dir.create(stat_dir)
  data.table::fwrite(stat1, paste0(stat_dir,'clone_assign_main_site.csv.gz'))
  data.table::fwrite(stat2, paste0(stat_dir,'clone_assign_pdxids.csv.gz'))
  data.table::fwrite(stat3, paste0(stat_dir,'clone_assign_Site_origin.csv.gz'))
  
  stat1 <- data.table::fread(paste0(stat_dir,'clone_assign_main_site.csv.gz'))
  stat2 <- data.table::fread(paste0(stat_dir,'clone_assign_pdxids.csv.gz'))
  stat3 <- data.table::fread(paste0(stat_dir,'clone_assign_Site_origin.csv.gz'))
  
  table(stat1$Grouping, stat1$clone_id)
  s1 <- stat1 %>%
    tidyr::pivot_wider(names_from = 'clone_id', values_from = 'nb_cells')
  View(s1)
  data.table::fwrite(s1, paste0(stat_dir,'clone_assign_main_site_wideformat.csv'))
  
  
  s2 <- stat2 %>%
    tidyr::pivot_wider(names_from = 'clone_id', values_from = 'nb_cells')
  View(s2)
  data.table::fwrite(s2, paste0(stat_dir,'clone_assign_pdxids_wideformat.csv'))
  
  
  s3 <- stat3 %>%
    tidyr::pivot_wider(names_from = 'clone_id', values_from = 'nb_cells')
  View(s3)
  data.table::fwrite(s3, paste0(stat_dir,'clone_assign_Site_origin_wideformat.csv'))
  
  
  
  # assigned$cell_id
  s <- 'SA535X4XB05662'
  df <- data.table::fread(paste0(base_dir,'input_data/',s,'_expr_introns.csv.gz'))
  dim(df)
  df$V1 <- NULL
  df1 <- tibble::tibble(cell_id=colnames(df), clone_id='C')  
  data.table::fwrite(df1, paste0(save_dir,s,'_clone_assign.csv.gz'))
  
}





input_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/clonealign_v3/'
sid <- 'SA535X4XB05462'

sce <- readRDS(paste0(input_dir, '/',sid, '/',sid,'_sce.rds'))
dim(sce)
sce[1:3,1:3]
dim(counts(sce))
data.table::fwrite(as.matrix(counts(sce)), paste0(input_dir, '/',sid, '/',sid,'_expr.csv.gz'))

sce <- readRDS(paste0(input_dir, '/',sid, '/',sid,'_processed_data.rds'))
dim(sce$gene_expression_data)
expr <- t(sce$gene_expression_data)
class(expr)
View(expr[1:3,1:3])
dim(expr)
expr <- as.data.frame(expr)
data.table::fwrite(expr, paste0(input_dir, '/',sid, '/',sid,'_expr.csv.gz'), row.names=T)
t <- data.table::fread(paste0(input_dir, '/',sid, '/',sid,'_expr.csv.gz'))
t[1:3,1:3]
cnv <- sce$copy_number_data
class(cnv)
dim(cnv)
cnv[1:5,1:5]

clones_df <- data.frame(cell_id=paste0('cell_',colnames(cnv)),clone_id=colnames(cnv))
colnames(cnv) <- clones_df$cell_id
data.table::fwrite(clones_df, paste0(input_dir, '/',sid, '/',sid,'_cell_clones.csv.gz'))
t <- data.table::fread(paste0(input_dir, '/',sid, '/',sid,'_cell_clones.csv.gz'))
head(t)

cnv <- as.data.frame(cnv)
data.table::fwrite(cnv, paste0(input_dir, '/',sid, '/',sid,'_clones_cnv.csv.gz'),row.names = T)
t <- data.table::fread(paste0(input_dir, '/',sid, '/',sid,'_clones_cnv.csv.gz'))
t[1:3,1:3]
input_dir <- '/home/htran/storage/datasets/metastasis_results/SA535_new_encoding/SA535_wholedata_v2/'
cnv <- data.table::fread(paste0(input_dir, 'total_merged_filtered_states.csv'))
rownames(cnv) <- cnv$V1
cnv$V1 <- NULL
dim(cnv)
clones <- data.table::fread(paste0(input_dir, 'cell_clones.csv'))
dim(clones)
sum(colnames(cnv) %in% clones$cell_id)
cells_use <- colnames(cnv)[colnames(cnv) %in% clones$cell_id]
length(cells_use)
cnv <- as.data.frame(cnv)
cnv <- cnv[,cells_use]
dim(cnv)
