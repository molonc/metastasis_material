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


get_geneId_without_version <- function(gene_ids) {  
  labels <- sapply(strsplit(gene_ids, "\\."), function(x) {  
    return(x[1])  
  })  
  return(as.character(labels))  
}
build_mapping_reference_GRCh38_96 <- function(){
  
  ref_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/tx2gene_sets/'
  # library(rtracklayer)
  
  # gtf <- rtracklayer::import("~/Downloads/Homo_sapiens.GRCh38.96.gtf")
  # 
  # # View(head(gtf))
  # t2g <- as.data.frame(gtf) %>%
  #   filter(type =="transcript") %>%
  #   mutate(transcript_id_with_version = paste0(transcript_id, ".", transcript_version),
  #          ensembl_gene_with_version = paste0(gene_id, ".", gene_version))
  # 
  # t2g <- t2g[, c("transcript_id_with_version", "ensembl_gene_with_version", "gene_name")]
  # 
  # dim(t2g)
  # head(t2g)
  # data.table::fwrite(t2g, paste0(ref_dir, 'Homo_sapiens_GRCh38_96.csv.gz'))
  t2g <- data.table::fread(paste0(ref_dir, 'Homo_sapiens_GRCh38_96.csv.gz'))
  t2g$TXNAME <- get_geneId_without_version(t2g$transcript_id_with_version)
  t2g$GENEID <- get_geneId_without_version(t2g$ensembl_gene_with_version)
  # head(t2g)
  data.table::fwrite(t2g, paste0(ref_dir, 'Homo_sapiens_GRCh38_96.csv.gz'))
}  

load_raw_counts_kallisto <- function(input_dir, datatag, 
                                     save_dir, sample_ids=NULL){
  
  if(is.null(sample_ids)){
    sample_ids <- list.files(input_dir)
  }
  print(sample_ids)
  
  ## Only taking into account samples with abundance.h5 or abundance.tsv files in the folder
  ## .h5 file is faster
  input_fns <- c()
  sids <- c()
  for(s in sample_ids){
    bulk_fn <- paste0(input_dir, s,'/abundance.h5') ## you can use .tsv file as well
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
  # ref_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/tx2gene_sets/'
  ref_dir <- save_dir
  ref_tx2gene_fn <- paste0(ref_dir, 'tx2gene_updated_23_Oct_2023.csv.gz')
  tx2gene <- data.table::fread(ref_tx2gene_fn, header=T) %>% as.data.frame()
  head(tx2gene)
  dim(tx2gene)
  
  txi <- tximport(input_fns, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = T, ignoreTxVersion=T)
  dim(txi$counts)
  # class(txi$counts)
  # head(txi$abundance)
  # class(txi$length)
  
  ## Tested with Homo_sapiens_GRCh38_96 reference but obtained less number of mapping transcripts to genes
  # ref_tx2gene_fn <- paste0(ref_dir, 'Homo_sapiens_GRCh38_96.csv.gz')
  # tx2gene <- data.table::fread(ref_tx2gene_fn, header=T) %>% as.data.frame()
  # head(tx2gene)
  # tx2gene <- tx2gene %>%
  #   dplyr::select(TXNAME, GENEID)
  # tx2gene <- tx2gene %>%
  #   dplyr::select(transcript_id_with_version, ensembl_gene_with_version) %>%
  #   dplyr::rename(TXNAME=transcript_id_with_version, GENEID=ensembl_gene_with_version)
  # txi2 <- tximport(input_fns, type = "kallisto", tx2gene = tx2gene) #ignoreAfterBar = T, ignoreTxVersion=T)
  # dim(txi2$counts)
  
  colnames(txi$counts) <- sample_ids
  colnames(txi$length) <- sample_ids
  colnames(txi$abundance) <- sample_ids
  
  print(paste0('# of genes: ',dim(txi$counts)[1], ', #samples: ',dim(txi$counts)[2]))
  head(txi$counts)
  txi$counts <- round(txi$counts,2) # cut off numbers after , only keep 2 numbers after ,
  head(txi$counts)
  
  df_counts <- as.data.frame(txi$counts)
  # head(df_counts)
  # dim(df_counts)
  
  ## Removing NA rows 
  na_rows <- apply(is.na(df_counts),1,sum)
  names(na_rows) <- rep(1:dim(df_counts)[1])
  na_rows <- na_rows[na_rows>0]
  length(na_rows)
  selected_rows <- rep(1:dim(df_counts)[1])
  selected_rows <- selected_rows[!selected_rows %in% names(na_rows)]
  length(selected_rows)
  df_counts <- df_counts[selected_rows,]
  df_counts$ens_gene_id <- rownames(df_counts)
  dim(df_counts)
  data.table::fwrite(df_counts, paste0(save_dir, datatag,'_total_raw_counts.csv.gz'))
  
  gene_length_df <- as.data.frame(txi$length)
  gene_length_df$ens_gene_id <- rownames(gene_length_df)
  # head(gene_length_df)
  gene_length_df <- gene_length_df %>%
    dplyr::filter(ens_gene_id %in% df_counts$ens_gene_id)
  data.table::fwrite(gene_length_df, paste0(save_dir, datatag,'_total_raw_length.csv.gz'))
  
  # abundance_df <- as.data.frame(txi$abundance)
  # abundance_df$ens_gene_id <- rownames(abundance_df)
  # data.table::fwrite(abundance_df, paste0(save_dir, datatag,'_total_raw_abundance.csv.gz'))
  
  
  return(df_counts)
}

## gene_filter_thresh=30: Sleuth use the thres=5, edgeR use 10 or 15, depending number of samples, you can use 
## different threshold, here with 20 samples, you can use 30
normalize_by_size_factor <- function(df_counts_fn, datatag, save_dir, gene_filter_thresh=30){
  # input_dir <- '/Users/htran/Documents/storage_tmp/metastasis_trees/SA919_10x/'
  
  df <- data.table::fread(df_counts_fn) %>% as.data.frame()
  print(dim(df))
  # head(df)
  if('ens_gene_id' %in% colnames(df)){
    # meta_genes <- data.frame(ens_gene_id=df$ens_gene_id)
    rownames(df) <- df$ens_gene_id
    df$ens_gene_id <- NULL  
  }
  if(is.null(rownames(df))){
    stop('Need gene id in row names')
  }
  print(df[1:2,])
  print('Sequencing depth (total UMI counts) for each sample in raw data is: ')
  print(colSums(df, na.rm = T))
  print(colnames(df))
  # print(rownames(df)[1:3])
  
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
  print('Sequencing depth (total UMI counts) for each sample at raw data is: ')
  s <- colSums(raw_counts)
  print(s)
  
  norm_counts_mtx <- normcounts(sce)
  s2 <- colSums(norm_counts_mtx)
  print('Sequencing depth (total UMI counts) for each sample after size factor normalization is: ')
  print(s2)
  dim(norm_counts_mtx)
  norm_counts_mtx <- as.data.frame(norm_counts_mtx)
  
  
  ## Do gene filtering here
  ## gene_filter_thresh=30: Sleuth use the thres=5, edgeR use 10 or 15, depending number of samples, you can use 
  ## different threshold, here with 20 samples, you can use 30
  rs <- rowSums(norm_counts_mtx)
  raw_lg <- length(rs)
  rs[1:10]
  rs <- rs[rs>gene_filter_thresh]
  filtered_nb_genes <- length(rs)
  print(paste0('Removing #',(raw_lg-filtered_nb_genes),' genes with low expression from normalized data'))
  print(paste0('# filtered genes: ',filtered_nb_genes))
  norm_counts_mtx$ens_gene_id <- rownames(norm_counts_mtx)
  norm_counts_mtx$ens_gene_id[1:2]
  norm_counts_mtx <- norm_counts_mtx %>%
    dplyr::filter(ens_gene_id %in% names(rs))
  print(dim(norm_counts_mtx))
  
  ## Removing mitochondrial genes from the list
  meta_genes <- data.frame(ens_gene_id=norm_counts_mtx$ens_gene_id)
  ref <- annotables::grch38 %>%
    dplyr::select(ensgene, symbol) %>%
    dplyr::rename(ens_gene_id=ensgene) %>%
    dplyr::filter(ens_gene_id %in% meta_genes$ens_gene_id)
  ref <- ref[!duplicated(ref$ens_gene_id),]
  dim(meta_genes)
  meta_genes <- meta_genes %>%
    dplyr::left_join(ref, by='ens_gene_id')
  dim(meta_genes)
  data.table::fwrite(meta_genes, paste0(save_dir, datatag, '_filtered_genes_list.csv'))
  
  mito_genes <- get_mito_genes(ref)
  print('Removing mito genes from the list')
  norm_counts_mtx <- norm_counts_mtx %>%
    dplyr::filter(!ens_gene_id %in% mito_genes)
  print(dim(norm_counts_mtx))
  
  data.table::fwrite(norm_counts_mtx, paste0(save_dir, datatag, '_sizefactor_normalized.csv.gz'))
  
  # dim(raw_counts)
  # raw_counts$ens_gene_id <- rownames(raw_counts)
  # raw_counts <- raw_counts %>%
  #   dplyr::filter(ens_gene_id %in% norm_counts_mtx$ens_gene_id)
  # dim(raw_counts)
  # ## Raw data but with only filtered genes
  # data.table::fwrite(raw_counts, paste0(save_dir, datatag,'_total_raw_counts_filtered_genes.csv.gz'))
  
  return(norm_counts_mtx)
}

# meta_genes has columns: symbol, ens_gene_id
get_mito_genes <- function(meta_genes, is_human=T){
  if(is_human){
    mito_gene_regex = '^MT-' # for human  
    print('human genes option')
  }else{
    mito_gene_regex = '^Mt-' # for mouse  
    print('mouse genes option')
  }
  
  
  # mt_genes_symbols <- meta_genes$symbol[grepl(mito_gene_regex,meta_genes$symbol)]
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

## Normalizing TPM: 
## https://www.biostars.org/p/335187/
## Or https://rdrr.io/github/skimlab/CCSBUtils/man/counts_to_tpm.html
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

## To Do: 
## Get normalized data, convert it to cpm/tpm/fpkm format if needed
## Selecting only clone C
## Do clustering using ComplexHeatmap with Pearson correlation as distance
## Visualize output
get_clustering <- function(normalized_fn, save_dir, datatag){
  
  normalized_df <- data.table::fread(normalized_fn)
  dim(normalized_df)  
  colnames(normalized_df)
  norm_df <- normalized_df
  rownames(norm_df) <- norm_df$ens_gene_id
  norm_df$ens_gene_id <- NULL
  colnames(norm_df)
  metasamples <- tibble::tibble(sample_id=c("AT22908","AT22909","AT22910","AT22911"),
                                condition=c('UnRx','CX','CX','CX'),
                                time_point=c('0h','6h','16h','30h'))
  metasamples <- as.data.frame(metasamples)
  rownames(metasamples) <- metasamples$sample_id
  
  predefined_cols <- c('#EE220C','#004D80')
  names(predefined_cols) <- c('CX','UnRx')
  
  time_point_cols <- brewer.pal(8, "Dark2")[1:length(unique(metasamples$time_point))]
  names(time_point_cols) <- unique(metasamples$time_point)
  # # colors_use <- predefined_cols[metasamples[colnames(norm_df),'mainsite']]
  # # cols_use=c("#66A61E", "#66A61E", "#666666") # green and grey
  # # names(colors_use) <- colnames(norm_df)
  # # unique(metasamples[colnames(norm_df),'main_clone'])
  # clones_color <- c('A'='#66C2A5','B'='#FC8D62','C'='#8DA0CB')
  # # colors_clone_use <- clones_color[metasamples[colnames(norm_df),'main_clone']]
  # # cols_use=c("#66A61E", "#66A61E", "#666666") # green and grey
  # # names(colors_clone_use) <- colnames(norm_df)
  # mouse_cols <- brewer.pal(8, "Accent")[1:length(unique(metasamples$transplanted_mouse_id))]
  # names(mouse_cols) <- unique(metasamples$transplanted_mouse_id)
  # 
  
  # 
  # 
  top_anno = ComplexHeatmap::HeatmapAnnotation(Condition = factor(metasamples[colnames(norm_df),'condition']),
                                               TimePoint=factor(metasamples[colnames(norm_df),'time_point']),
                                               # MouseId=factor(metasamples[colnames(norm_df),'transplanted_mouse_id']),
                                               # Origin=factor(metasamples[colnames(norm_df),'origin']),
                                               col = list(Condition=predefined_cols,
                                                          TimePoint=time_point_cols#,
                                                          # MouseId=mouse_cols#,
                                                          # Origin=origin_cols
                                               ))
  # library(ComplexHeatmap)
  
  ## NOTE: here I used normalized counts as input data, you can use FPKM, or TPM values as input, see normalized functions above
  # fpkm_df <- get_normalized_data(normalized_df)
  # 
  p <- ComplexHeatmap::Heatmap(as.matrix(norm_df), na_col = "white",
                               # col = col_fun,
                               show_column_names=T,
                               show_row_names = F,
                               cluster_rows=F,
                               cluster_columns=T,
                               clustering_distance_columns = "pearson",
                               name = "Test", 
                               # row_order = sort(rownames(test)),
                               # row_split= samples_use,
                               row_title_rot = 0,
                               row_gap = unit(2, "mm"),
                               # column_split = genes_type$gt,
                               # column_title = paste0("Filtered Normalized Data Clustering ",datatag),
                               column_gap = unit(2, "mm"),
                               column_names_gp = grid::gpar(fontsize = 10),
                               row_names_gp = grid::gpar(fontsize = 10),
                               show_heatmap_legend = T,
                               top_annotation=top_anno,
                               # left_annotation = left_anno,
                               # cell_fun = cell_func,
                               row_dend_reorder=T
  )
  # p
  png(paste0(save_dir,datatag,"_hierarchial_clusters.png"), height = 2*800, width=2*1000, res = 2*72)
  print(p)
  dev.off()  
}  
  
main <- function(){
  base_dir <- '/home/htran/storage/rnaseq_datasets/bulk_metastasis/example_bulkRNAseq_analysis/'
  input_dir <- paste0(base_dir, 'input_Kallisto_files/')
  datatag <- 'example'
  save_dir <- paste0(base_dir, 'preprocessed/')
  if(!dir.exists(save_dir)){
    dir.create(save_dir)
  }
  sample_ids=NULL
  df_counts <- load_raw_counts_kallisto(input_dir, datatag, save_dir, sample_ids)
  
  ## Loading counts - output of above function, and normalize data
  df_counts_fn <- paste0(save_dir, datatag,'_total_raw_counts.csv.gz')
  df_normalized <- normalize_by_size_factor(df_counts_fn, datatag, save_dir)
  
  ## Then do hierarchical clustering using normalized data
  normalized_fn <- paste0(save_dir, datatag, '_sizefactor_normalized.csv.gz')
  get_clustering(normalized_fn, save_dir, datatag)
  
  
  ## DE analysis 
  ## By DRUG TREATMENT CONDITION
  
  ## Get raw counts data, selecting only 6h Rx vs 6h UnRx samples
  ## Do DE analysis using DESeq2
  ## Get genes with abs(logFC) >=0.5 ## too small values may be noise, not real biology effect
  ## Intersecting obtained genes with filtered genes list from normalized step above. 
  ## Filtered genes list is at paste0(save_dir, datatag, '_filtered_genes_list.csv')
  
  ## Get raw counts data, selecting only 16h Rx vs 16h UnRx samples
  ## Do DE analysis using DESeq2
  ## Get genes with abs(logFC) >=0.5 ## too small values may be noise, not real biology effect
  ## Intersecting obtained genes with filtered genes list from normalized step above. 
  ## Filtered genes list is at paste0(save_dir, datatag, '_filtered_genes_list.csv')
  
  
  
  ## By TIME POINT
  ## Get raw counts data, selecting only 168h Rx vs 6h UnRx samples
  ## Do DE analysis using DESeq2
  ## Get genes with abs(logFC) >=0.5 ## too small values may be noise, not real biology effect
  ## Intersecting obtained genes with filtered genes list from normalized step above. 
  ## Filtered genes list is at paste0(save_dir, datatag, '_filtered_genes_list.csv')
  
  ## DE analysis OTHER COMPARISONS
  ## Looking at hierarchical clustering output, pick untreated, and treated from same cluster
  ## Form DE analysis comparison, and selecting genes. 
  
  ## Have fun!!!
}
