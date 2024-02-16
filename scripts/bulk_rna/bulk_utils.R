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
get_normalized_TPM <- function(counts, lengths) {
  print(sum(lengths==0))
  lengths <- ifelse(lengths==0, 0.01,lengths)
  rate = counts / lengths
  norm_vals <- sapply(rate,function(x){
    return(1e6*x/sum(x))
  }) 
  return(norm_vals)  
}  
get_normalized_FPKM <- function(counts, effLen)
{
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

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
  # head(df)
  data.table::fwrite(df, paste0(save_dir, datatag, '_meta_samples.csv.gz'))
  return(df)
}
## Tutorial here: https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html
## Reading txt files
# files <- file.path(dir, "kallisto", samples$run, "abundance.tsv.gz")
get_geneId_without_version <- function(gene_ids) {  
  labels <- sapply(strsplit(gene_ids, "\\."), function(x) {  
    return(x[1])  
  })  
  return(as.character(labels))  
}

build_mapping_reference <- function(){
  ## Best version of mapping reference between transcript and genes from annotables package
  tx2gene1 <- annotables::grch38_tx2gene
  colnames(tx2gene1) <- c('TXNAME','GENEID')
  # head(tx2gene1)
  # dim(tx2gene1)
  
  # The most updated version: tx2gene_updated_06_Dec_2021.csv.gz
  # ref_tx2gene_fn <- 'mydir/convert_transcript_gene/tx2gene.csv'
  # ref_tx2gene_fn <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/tx2gene_sets/tx2gene_updated_06_Dec_2021.csv.gz'
  
  ref_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/tx2gene_sets/'
  ref_tx2gene_fn <- paste0(ref_dir, 'tx2gene.csv.gz')
  tx2gene <- data.table::fread(ref_tx2gene_fn, header=T) %>% as.data.frame()
  tx2gene$V1 <- NULL
  tx2gene$TXNAME <- get_geneId_without_version(tx2gene$TXNAME)
  tx2gene$GENEID <- get_geneId_without_version(tx2gene$GENEID)
  extra_mapping_genes <- tx2gene %>%
    dplyr::filter(!GENEID %in% tx2gene1$GENEID)
  # head(extra_mapping_genes)
  # dim(extra_mapping_genes)
  tx2gene <- dplyr::bind_rows(tx2gene1, extra_mapping_genes)
  dim(tx2gene)
  head(tx2gene)
  length(unique(tx2gene$GENEID))
  data.table::fwrite(tx2gene, paste0(ref_dir, 'tx2gene_updated_23_Oct_2023.csv.gz'), 
                     row.names = F)
  # sum(!extra_mapping_genes$GENEID %in% tx2gene1$GENEID)
  # sum(!extra_mapping_genes$TXNAME %in% tx2gene1$TXNAME)
  # # head(tx2gene)
  # dim(tx2gene)
  # for(f in input_fns){
  #   if(!file.exists(f)){
  #     stop(paste0('Do not exist ',f))
  #   }
  # }
  
  return(tx2gene)
  
}
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
  
  ## Just for testing
  # input_fns <- input_fns[1]
  # sample_ids <- sids[1]
  ref_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/tx2gene_sets/'
  ref_tx2gene_fn <- paste0(ref_dir, 'tx2gene_updated_23_Oct_2023.csv.gz')
  tx2gene <- data.table::fread(ref_tx2gene_fn, header=T) %>% as.data.frame()
  head(tx2gene)
  
  txi <- tximport(input_fns, type = "kallisto", tx2gene = tx2gene, 
                  ignoreAfterBar = T, ignoreTxVersion=T)
  dim(txi$counts)
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
  
  # gene_length_df <- as.data.frame(txi$length)
  # gene_length_df$ens_gene_id <- rownames(gene_length_df)
  # abundance_df <- as.data.frame(txi$abundance)
  # abundance_df$ens_gene_id <- rownames(abundance_df)
  
  # data.table::fwrite(df_counts, paste0(save_dir, datatag,'_total_raw_counts.csv.gz'))
  # data.table::fwrite(gene_length_df, paste0(save_dir, datatag,'_total_raw_length.csv.gz'))
  # data.table::fwrite(abundance_df, paste0(save_dir, datatag,'_total_raw_abundance.csv.gz'))
  
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
  print(colSums(df, na.rm = T))
  print(colnames(df))
  print(rownames(df)[1:3])
  # df_backup <- df
  # df <- df[,]
  
  ## Removing NA rows before normalization
  na_rows <- apply(is.na(df),1,sum)
  names(na_rows) <- rep(1:dim(df)[1])
  na_rows <- na_rows[na_rows>0]
  length(na_rows)
  selected_rows <- rep(1:dim(df)[1])
  selected_rows <- selected_rows[!selected_rows %in% names(na_rows)]
  length(selected_rows)
  df <- df[selected_rows,]
  dim(df)
  
  # ref <- annotables::grch38 %>%
  #   dplyr::select(ensgene,symbol,chr) %>%
  #   dplyr::rename(ens_gene_id=ensgene) %>%
  #   dplyr::filter(ens_gene_id %in% meta_genes$ens_gene_id)
  # ref <- ref[!duplicated(ref$ens_gene_id),]
  # dim(ref)
  # meta_genes <- meta_genes %>%
  #   dplyr::left_join(ref, by='ens_gene_id')
  # dim(meta_genes)
  # sum(rownames(df)==meta_genes$ens_gene_id)==dim(df)[1]
  
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts=as.matrix(df)))
  # rowData(sce) <- meta_genes
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
  
  # saveRDS(sce, paste0(save_dir, datatag, '_sizefactor_normalized.rds'))
  
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
  
  data.table::fwrite(norm_counts_mtx, paste0(save_dir, datatag, '_sizefactor_normalized.csv.gz'))
  return(norm_counts_mtx)
}

get_mito_genes <- function(meta_genes){
  mito_gene_regex = '^MT-' # for human
  # mito_gene_regex = '^Mt-' # for mouse
  mt_genes_symbols <- meta_genes$symbol[grepl(mito_gene_regex,meta_genes$symbol)]
  # [1] "MT-ND6"  "MT-CO2" 
  # [3] "MT-CYB"  "MT-ND2" 
  # [5] "MT-ND5"  "MT-CO1" 
  # [7] "MT-ND3"  "MT-ND4" 
  # [9] "MT-ND1"  "MT-ATP6"
  # [11] "MT-CO3"  "MT-ND4L"
  # [13] "MT-ATP8"
  
  mt_genes_ens <-meta_genes$ens_gene_id[grepl(mito_gene_regex,meta_genes$symbol)]
  # [1] "ENSG00000198695"
  # [2] "ENSG00000198712"
  # [3] "ENSG00000198727"
  # [4] "ENSG00000198763"
  # [5] "ENSG00000198786"
  # [6] "ENSG00000198804"
  # [7] "ENSG00000198840"
  # [8] "ENSG00000198886"
  # [9] "ENSG00000198888"
  # [10] "ENSG00000198899"
  # [11] "ENSG00000198938"
  # [12] "ENSG00000212907"
  # [13] "ENSG00000228253"
  return(mt_genes_ens)
  
}
