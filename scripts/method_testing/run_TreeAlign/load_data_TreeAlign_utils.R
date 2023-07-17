suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(RColorBrewer)
  # library(Seurat)
  library(dplyr)
  library(ggplot2)
})


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

load_10x_TreeAlign <- function(sce_fn, ens_obs_genes, sid, save_dir){
  sce <- readRDS(sce_fn)
  print(dim(sce))
  
  colnames(sce) <- sce$Barcode
  if(is.null(rownames(sce)) || !grepl('EN',rownames(sce)[1])){
    rownames(sce) <- rowData(sce)$ID
  }
  # rownames(sce)[1]
  # saveRDS(sce, paste0(save_dir, sid, '_introns_sce.rds'))
  # print(assayNames(sce))
  
  mtx <- as.data.frame(as.matrix(counts(sce)))
  # dim(mtx)
  print(rownames(mtx)[1])
  print(colnames(mtx)[1])
  
  
  genes_used <- intersect(ens_obs_genes, rownames(mtx))
  print(length(genes_used))
  mtx <- mtx[genes_used,]
  
  data.table::fwrite(mtx, paste0(save_dir, sid, '_expr.csv.gz'), row.names = T)
  print(dim(mtx))
}

get_cnv_sid <- function(cnv, save_dir, sid, obs_clones,
                        purity_thrs=0.4, avg_mappability_thrs=0.7){
  # save_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/clonealign_introns/'
  # cnv_fn <- '/home/htran/storage/raw_DLP/metastasis_DLP/SA535/reads_clones_v3/mapped_wholedata_SA535_v2.csv.gz'
  # cnv <- data.table::fread(cnv_fn) %>% as.data.frame()
  print(dim(cnv))
  colnames(cnv)
  
  cols_use <- c("ensembl_gene_id", "gene_symbol", 
                paste0('pct_pure_',obs_clones),
                paste0('avg_mappability_',obs_clones),
                paste0('cn_median_',obs_clones))
  cnv <- cnv %>% 
    dplyr::select(all_of(cols_use))
  print(dim(cnv))
  if(avg_mappability_thrs>0){
    map_cols <- colnames(cnv)[grepl('avg_mappability',colnames(cnv))]
    for(j in seq(1:length(map_cols))){
      if(j==1){
        cond1 <- cnv[,map_cols[j]] > purity_thrs   
      }else{
        cond1 <- cond1  & (cnv[,map_cols[j]] > avg_mappability_thrs)
      }
      
    }
    
    cnv <- cnv[cond1,]
    print(dim(cnv))
  }
  ## Filtering genomic regions by purity of copy number
  if(purity_thrs>0){
    purity_cols <- colnames(cnv)[grepl('pct_pure',colnames(cnv))]
    for(i in seq(1:length(purity_cols))){
      if(i==1){
        cond2 <- cnv[,purity_cols[i]] > purity_thrs   
      }else{
        cond2 <- cond2  & (cnv[,purity_cols[i]] > purity_thrs)
      }
      
    }
    cnv <- cnv[cond2,]
    print(dim(cnv))
  }
  
  
  cnv1 <- cnv %>%
    dplyr::select(all_of(paste0('cn_median_',obs_clones)))
  colnames(cnv1)
  var_genes <- sapply(seq(1:dim(cnv1)[1]), function(x){
    var(as.numeric(cnv1[x,]))
  })
  names(var_genes) <- cnv$ensembl_gene_id
  var_genes <- var_genes[var_genes>0]
  print(sum(var_genes>0))
  
  cnv <- cnv %>%
    dplyr::filter(ensembl_gene_id %in% names(var_genes))
  print(dim(cnv))
  # head(cnv)
  data.table::fwrite(cnv, paste0(save_dir, sid, '_cnv_details.csv.gz'))
  cols_use <- c("ensembl_gene_id", paste0('cn_median_',obs_clones))
  cnv <- cnv %>% 
    dplyr::select(all_of(cols_use))
  cols_use <- colnames(cnv)
  cols_use <- gsub('cn_median_','',cols_use)
  colnames(cnv) <- cols_use
  return(cnv)
  
}

cells_clone_fn <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cell_clones/SA609_cell_clones_metadata.csv'
cnv_fn <- '/home/htran/storage/raw_DLP/drug_resistance_DLP/SA609/reads_clones_v3/mapped_wholedata_SA609_grch38_v2.csv.gz'
save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/TreeAlign_testing/'
input_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/SA609-v6/"
# sids <- c('SA609X5XB03223','SA609X7XB03554')
sids <- c('SA609X3XB01584','SA609X4XB003083','SA609X4XB03080',
          'SA609X5XB03230','SA609X5XB03231','SA609X6XB03447')


# cnv contains ensembl_gene_id and list of clones, each clone is one column
get_clones_per_sample <- function(cnv_fn, cells_clone_fn, sids, 
                                  save_dir, min_cells_per_clone=10){
  # cells_clone_fn <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/dlp_trees/SA535/cell_clones.csv.gz'
  # save_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/clonealign_introns/input_data/'
  min_cells_per_clone=10
  if(!dir.exists(save_dir)){
    dir.create(save_dir)  
  }
  cnv <- data.table::fread(cnv_fn) %>% as.data.frame()
  cols <- colnames(cnv)
  cols <- gsub('_R','_A',cols)
  colnames(cnv) <- cols
  cells_clone <- data.table::fread(cells_clone_fn)
  colnames(cells_clone)
  dim(cells_clone)
  cells_clone$library_id <- get_library_id(cells_clone$cell_id)
  cells_clone$sample_id <- get_sample_id(cells_clone$cell_id)
  unique(cells_clone$sample_id)
  unique(cells_clone$library_id)
  
  # meta <- data.table::fread('/home/htran/Projects/hakwoo_project/metastasis_material/materials/10x/SA535_10x_metadata_passage_X4_full_infos.csv')
  # sids <- unique(meta$mouse_id)
  # # 
  
  ## for TreeAlign shell script
  writeLines(sids, paste0(save_dir, "samples_run_list.txt"))  
  
  unassigned <- c('None','Un','unassigned','Unassigned')
  for(sid in sids){
    print(sid)
    
    clone_used <- cells_clone %>% 
      dplyr::filter(sample_id==sid & !clone_id %in% unassigned) %>%
      dplyr::group_by(clone_id) %>%
      dplyr::summarise(nb_cells=n()) %>%
      dplyr::filter(nb_cells>=min_cells_per_clone)
    
    if(dim(clone_used)[1]>0){
      obs_clones <- unique(clone_used$clone_id)
      print(obs_clones)
      cnv_sid <- get_cnv_sid(cnv, save_dir, sid, obs_clones,
                  purity_thrs=0.4, avg_mappability_thrs=0.7)
      colnames(cnv_sid)
      cnv1 <- cnv_sid %>%
        dplyr::select(all_of(c('ensembl_gene_id', obs_clones)))
      colnames(cnv1) <- c('ensembl_gene_id', paste0('cell_',obs_clones))
      data.table::fwrite(cnv1, paste0(save_dir, sid, '_clones_cnv.csv.gz'))  
      
      clones_df <- data.frame(cell_id=paste0('cell_',obs_clones),
                              clone_id=obs_clones)
      
      data.table::fwrite(clones_df, paste0(save_dir,sid,'_cell_clones.csv.gz'))
      sce_fn <- paste0(input_dir, sid,'.rdata')
      if(file.exists(sce_fn)){
        ens_obs_genes <- unique(cnv1$ensembl_gene_id)
        load_10x_TreeAlign(sce_fn, ens_obs_genes, sid, save_dir)  
      }
      
    }
  }  
    
}


