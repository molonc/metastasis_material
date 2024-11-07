suppressPackageStartupMessages({
  library(tidyverse)
  library(annotables)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  # require(scales)
  # require(ggrepel)
  require(stringr)
  library(scran)
  library(SingleCellExperiment)
  library(tximport)
})

get_meta_samples_SA919_mixing_exp <- function(meta_samples_fn, datatag, save_dir){
  df <- data.table::fread(meta_samples_fn) 
  # dim(df)
  # sum(df$`Bulk RNA ATID`!='')
  df <- df %>% 
    dplyr::filter(`Bulk RNA ATID`!='')
  colnames(df) <- gsub(' ','_', colnames(df))
  df <- df %>%
    dplyr::mutate(mouse_id=stringr::str_sub(Sample_Name,1, 5)) %>%
    dplyr::select(mouse_id, Bulk_RNA_ATID, Anatomical_Site, everything())
  head(df)
  data.table::fwrite(df, paste0(save_dir, datatag, '_meta_samples.csv.gz'))
  return(df)
}
## Tutorial here: https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html
## Reading txt files
# files <- file.path(dir, "kallisto", samples$run, "abundance.tsv.gz")

load_raw_counts_kallisto <- function(input_dir, datatag, 
                                     save_dir, sample_ids=NULL){
  
  if(is.null(sample_ids)){
    sample_ids <- list.files(input_dir)
  }
  
  input_fns <- c()
  sids <- c()
  for(s in sample_ids){
    bulk_fn <- paste0(input_dir, s,'/abundance.h5')
    if(file.exists(bulk_fn)){
      print(paste0('Sample: ',s,', read file: ', bulk_fn))
      input_fns <- c(input_fns,bulk_fn)
      sids <- c(sids, s)
    }
  }
  sample_ids <- sids
  ## input files
  # input_fns <- paste0(input_dir, sample_ids,'/abundance.h5') #.h5 or .tsv both are fine
  # input_fns <- paste0(input_dir, sample_ids,'/abundance.h5') #.h5 or .tsv both are fine
  
  ## Best version of mapping reference between transcript and genes from annotables package
  tx2gene <- annotables::grch38_tx2gene
  colnames(tx2gene) <- c('TXNAME','GENEID')
  #meta genes, the mapping from transcript to genes
  
  # The most updated version: tx2gene_updated_06_Dec_2021.csv.gz
  # ref_tx2gene_fn <- 'mydir/convert_transcript_gene/tx2gene.csv'
  # ref_tx2gene_fn <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/tx2gene_updated_06_Dec_2021.csv.gz'
  # tx2gene <- data.table::fread(ref_tx2gene_fn, header=T) %>% as.data.frame()
  # tx2gene$V1 <- NULL
  # # head(tx2gene)
  # dim(tx2gene)
  # for(f in input_fns){
  #   if(!file.exists(f)){
  #     stop(paste0('Do not exist ',f))
  #   }
  # }
  txi <- tximport(input_fns, type = "kallisto", tx2gene = tx2gene, 
                  ignoreAfterBar = TRUE, ignoreTxVersion=TRUE)
  # class(txi$counts)
  # head(txi$abundance)
  # class(txi$length)
  colnames(txi$counts) <- sample_ids
  colnames(txi$length) <- sample_ids
  colnames(txi$abundance) <- sample_ids
  
  print(paste0('# of genes: ',dim(txi$counts)[1], ', #samples: ',dim(txi$counts)[2]))
  txi$counts <- round(txi$counts,2) # cut off numbers after ,
  
  df_counts <- as.data.frame(txi$counts)
  df_counts$ens_gene_id <- rownames(df_counts)
  # dim(df_counts)
  
  gene_length_df <- as.data.frame(txi$length)
  gene_length_df$ens_gene_id <- rownames(gene_length_df)
  abundance_df <- as.data.frame(txi$abundance)
  abundance_df$ens_gene_id <- rownames(abundance_df)
  
  data.table::fwrite(df_counts, paste0(save_dir, datatag,'_total_raw_counts.csv.gz'))
  data.table::fwrite(gene_length_df, paste0(save_dir, datatag,'_total_raw_length.csv.gz'))
  data.table::fwrite(abundance_df, paste0(save_dir, datatag,'_total_raw_abundance.csv.gz'))
  
  ## Normalizing TPM: 
  ## https://www.biostars.org/p/335187/
  ## Or https://rdrr.io/github/skimlab/CCSBUtils/man/counts_to_tpm.html
  return(df_counts)
}

normalize_by_size_factor <- function(df_counts_fn, datatag, save_dir){
  # input_dir <- '/Users/htran/Documents/storage_tmp/metastasis_trees/SA919_10x/'
  
  df <- data.table::fread(df_counts_fn) %>% as.data.frame()
  print(dim(df))
  # head(df)
  if('ens_gene_id' %in% colnames(df)){
    meta_genes <- data.frame(ens_gene_id=df$ens_gene_id)
    rownames(df) <- df$ens_gene_id
    df$ens_gene_id <- NULL  
  }
  if(is.null(rownames(df))){
    stop('Need gene id in row names')
  }
  print(df[1:2,])
  print('Sequencing depth (total UMI counts) for each sample is: ')
  print(colSums(df))
  print(colnames(df))
  print(rownames(df)[1:3])
  ref <- annotables::grch38 %>%
    dplyr::select(ensgene,symbol,chr) %>%
    dplyr::rename(ens_gene_id=ensgene) %>%
    dplyr::filter(ens_gene_id %in% meta_genes$ens_gene_id)
  ref <- ref[!duplicated(ref$ens_gene_id),]
  dim(ref)
  meta_genes <- meta_genes %>%
    dplyr::left_join(ref, by='ens_gene_id')
  dim(meta_genes)
  sum(rownames(df)==meta_genes$ens_gene_id)==dim(df)[1]
  
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts=as.matrix(df)))
  rowData(sce) <- meta_genes
  # sce
  if(is.null(rownames(sce))){
    rownames(sce) <- rownames(df)
  }
  
  sce <- scran::computeSumFactors(sce)
  print(sce$sizeFactor)
  sce$size_factor <- sizeFactors(sce)
  sce <- scuttle::logNormCounts(sce, log=FALSE, exprs_values="counts", 
                                           assay.type="counts",
                                           name='normcounts', size_factors=sce$size_factor)
  
  saveRDS(sce, paste0(save_dir, datatag, '_sizefactor_normalized.rds'))
  ## Manual normalization instead of using logNormCounts function - same output
  # norm_counts_mtx <- sweep(raw_counts,2,sce$sizeFactor,FUN="/")
  # s2 <- colSums(norm_counts_mtx)
  
  raw_counts <- as.data.frame(counts(sce))
  rownames(raw_counts)[1:3]
  print('Sequencing depth (total UMI counts) for each sample at raw data is: ')
  s <- colSums(raw_counts)
  print(s)
  norm_counts_mtx <- normcounts(sce)
  s2 <- colSums(norm_counts_mtx)
  
  print('Sequencing depth (total UMI counts) for each sample after size factor normalization is: ')
  print(s2)
  dim(norm_counts_mtx)
  norm_counts_mtx <- as.data.frame(norm_counts_mtx)
  norm_counts_mtx$ens_gene_id <- rownames(norm_counts_mtx)
  norm_counts_mtx$ens_gene_id[1:2]
  
  data.table::fwrite(norm_counts_mtx, paste0(save_dir, datatag, '_total_sizefactor_normalized.csv.gz'))
  return(norm_counts_mtx)
}


main <- function(){
  input_dir <- '/home/htran/storage/rnaseq_datasets/bulk_metastasis/mixing_exp_SA919_bulkRNAseq/'  
  datatag <- 'SA919_mixing'
  save_dir <- '/home/htran/storage/datasets/metastasis_results/bulk_SA919/mixing_SA919/'
  df_counts <- load_raw_counts_kallisto(input_dir, datatag, save_dir, sample_ids=NULL)
  
  ## Loading counts - output of above function, and normalize data
  df_counts_fn <- paste0(save_dir, datatag,'_total_raw_counts.csv.gz')
  df_normalized <- normalize_by_size_factor(df_counts_fn, datatag, save_dir)
  
  ## Meta data
  save_dir <- '/home/htran/storage/datasets/metastasis_results/bulk_SA919/mixing_SA919/'
  meta_samples_fn <- paste0(save_dir,'Samples-Metastasis_Hakwoo_bulkRNA_mixing_exp.csv')
  get_meta_samples_SA919_mixing_exp(meta_samples_fn, datatag, save_dir)
  
}

main()