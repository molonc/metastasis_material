suppressPackageStartupMessages({
  library(tidyverse)
  library(SingleCellExperiment)
  # library(scater)
  library(data.table)
  # library(methods)
  # library(scran)
  library(parallel)
  # library(feather)
  # library(annotables)
  library(dplyr)
  library(ggplot2)
  # library(snpEnrichment)
  # library(scMerge)
  options(dplyr.summarise.inform = FALSE)
  options(tidyverse.quiet = TRUE)
  # library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(IRanges)
})




get_gene_coordinates <- function(segments){ #, min_pc_overlap=0.1
  # segments$chr <- paste0('chr',chrs_df$chr)
  # library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  # library(org.Hs.eg.db)
  # library(IRanges)
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  g <- genes(txdb, single.strand.genes.only=FALSE)
  
  segments_gr <- makeGRangesFromDataFrame(segments, keep.extra.columns = TRUE)
  overlaps <- findOverlaps(g, segments_gr, ignore.strand = TRUE)
  overlapping_ranges <- pintersect(g[queryHits(overlaps)], segments_gr[subjectHits(overlaps)], ignore.strand = TRUE, drop.nohit.ranges = TRUE)
  
  percentage_overlap <- width(overlapping_ranges) / width(segments_gr[subjectHits(overlaps)])
  # wd_overlap <- width(overlapping_ranges)
  # percentage_overlap <- width(overlapping_ranges) / width(g[queryHits(overlaps)])
  
  percentage_overlap_sum <- sapply(percentage_overlap, sum)
  # wd_overlap_sum <- sapply(wd_overlap, sum)
  gene_cn <- tibble(entrezgene = names(g)[queryHits(overlaps)], percent_overlap=percentage_overlap_sum) #wd_overlap=wd_overlap_sum
  # print(dim(gene_cn))
  # length(unique(gene_cn$entrezgene))
  # t <- gene_cn[duplicated(gene_cn$entrezgene),]
  # dim(t)
  # head(t)
  # t1 <- gene_cn[gene_cn$entrezgene=='10000',]
  gene_cn <- gene_cn %>%
    cbind(mcols(segments_gr[subjectHits(overlaps)]))
  # head(gene_cn)
  # dim(gene_cn)
  # length(unique(gene_cn$entrezgene))
  ext_rows <- gene_cn %>% # Remove the cases where a given gene that map to several regions.
    dplyr::group_by(entrezgene, cluster) %>%
    dplyr::summarise(percent_overlap=max(percent_overlap))%>% 
    dplyr::mutate(desc=paste0(entrezgene, cluster, percent_overlap))%>%
    dplyr::ungroup()%>%
    dplyr::pull(desc)
  # length(ext_rows)
  gene_cn <- gene_cn %>%
    dplyr::mutate(desc=paste0(entrezgene, cluster, percent_overlap))%>%
    dplyr::filter(desc %in% ext_rows)%>%
    dplyr::select(-desc)%>%
    dplyr::mutate(percent_overlap=round(percent_overlap, 4))
  
  
  # dim(gene_cn)
  gene_cn$ensembl_gene_id <- mapIds(org.Hs.eg.db,
                                    keys = gene_cn$entrezgene,
                                    column="ENSEMBL",
                                    keytype="ENTREZID",
                                    multiVals="first") #"CharacterList"
  # ens_mapped <- mapIds(org.Hs.eg.db,
  #             keys = gene_cn$entrezgene,
  #             column="ENSEMBL",
  #             keytype="ENTREZID",
  #             multiVals="list")
  
  ## Other way to do it
  # chrs <- c(as.character(1:22), "X")
  # t <- annotables::grch38 %>%
  #   dplyr::select(ensembl_gene_id = ensgene, entrezgene=entrez)  %>%
  #   dplyr::filter(entrezgene %in% gene_cn$entrezgene) #  & chr %in% chrs
  # #   inner_join(deg_df) %>%
  
  
  # ens_mapped <- unlist(ens_mapped)
  # ens_mapped_out <- data.frame(entrezgene=names(ens_mapped), ensembl_gene_id=ens_mapped)
  # ens_mapped_out <- ens_mapped_out %>%
  #   dplyr::filter(!is.na(ensembl_gene_id))
  # gene_cn <-  gene_cn %>% inner_join(ens_mapped_out, by='entrezgene')
  # print(dim(gene_cn))
  
  # Filter for complete entries
  gene_cn <- gene_cn %>%
    as.data.frame %>%
    drop_na()
  
  # gene_cn <- gene_cn %>%
  #   dplyr::filter(percent_overlap>min_pc_overlap)
  # print(dim(gene_cn))
  # gene_cn <- annotables::grch37 %>% 
  #   dplyr::select(ensembl_gene_id = ensgene, symbol) %>% 
  #   inner_join(gene_cn)
  # print(dim(gene_cn))
  # gene_cn$percent_overlap <- round(gene_cn$percent_overlap, 4)
  return(gene_cn)
}


get_mapping_genes <- function(){
  input_dir <- '/home/htran/storage/datasets/metastasis_results/SA535_new_encoding/SA535_wholedata_v2/CN_profile/'
  cnv <- data.table::fread(paste0(input_dir, 'median_cnv.csv'))
  selected_cols <- c('chr_desc',colnames(cnv)[grepl('_Primary',colnames(cnv))])
  cnv <- cnv %>%
    dplyr::select(all_of(selected_cols))
  dim(cnv)
  colnames(cnv) <- gsub('_Primary','',selected_cols)
  cols_use <- colnames(cnv)[colnames(cnv) != 'chr_desc']
  cnv_df <- cnv %>%
    tidyr::pivot_longer(cols = all_of(cols_use), names_to = 'cluster', values_to = 'copy_number')
  dim(cnv_df)
  head(cnv_df)
  chrs <- unlist(lapply(strsplit(cnv_df$chr_desc,'_'), function(x){
    return(x[1])
  }))
  starts <- unlist(lapply(strsplit(cnv_df$chr_desc,'_'), function(x){
    return(x[2])
  }))
  ends <- unlist(lapply(strsplit(cnv_df$chr_desc,'_'), function(x){
    return(x[3])
  }))
  cnv_df$chr <- as.character(chrs)
  cnv_df$start <- as.character(starts)
  cnv_df$end <- as.character(ends)
  if(!grepl('chr',cnv_df$chr[1])){
    cnv_df <- cnv_df %>%
      dplyr::mutate(chr=paste0("chr", chr))
  }
  
  gene_cn <- get_gene_coordinates(cnv_df)
  dim(gene_cn)
  gene_cn$entrezgene <- NULL
  gene_cn$ensembl_gene_id[1]
  chrs <- c(as.character(1:22), "X")
  meta_genes <- annotables::grch38 %>%
    dplyr::select(ensembl_gene_id = ensgene, gene_symbol=symbol, chr)  %>%
    dplyr::filter(ensembl_gene_id %in% unique(gene_cn$ensembl_gene_id) & chr %in% chrs) #  
  #   inner_join(deg_df) %>%
  gene_cn <- gene_cn %>%
    dplyr::left_join(meta_genes, by='ensembl_gene_id')
  sum(gene_cn$gene_symbol %in% c('MIB-3','MIB-1','MKI67'))
  dim(gene_cn)
  View(head(gene_cn))
  data.table::fwrite(gene_cn, paste0(input_dir, 'median_cnv_genes.csv'))
  
  colnames(gene_cn)
  gene_cn1 <- gene_cn
  gene_cn1$desc <- paste0(gene_cn1$ensembl_gene_id, gene_cn1$cluster)
  gene_cn1 <- gene_cn1[!duplicated(gene_cn1$desc),]
  gene_cn1 <- gene_cn1 %>%
    dplyr::select(cluster, copy_number, ensembl_gene_id, gene_symbol) %>%
    pivot_wider(names_from = 'cluster', values_from = 'copy_number')  
  
  dim(gene_cn1)
  data.table::fwrite(gene_cn1, paste0(input_dir, 'median_cnv_genes_wider.csv'))
  
  
  gene_cn1 <- gene_cn %>%
    dplyr::filter(gene_symbol %in% c('MIB-3','MIB-1','MKI67'))
  View(gene_cn1)
  return(gene_cn)
}

get_mapping_genes_SA919 <- function(){
  input_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_wholedata/CN_profile/'
  cnv <- data.table::fread(paste0(input_dir, 'median_cnv.csv'))
  selected_cols <- c('chr_desc',colnames(cnv)[grepl('_Primary',colnames(cnv))])
  
  selected_cols <- c('B_Primary','C_Metastasis','chr_desc')
  cnv <- cnv %>%
    dplyr::select(all_of(selected_cols))
  dim(cnv)
  colnames(cnv) <- gsub('_Primary','',colnames(cnv))
  colnames(cnv) <- gsub('_Metastasis','',colnames(cnv))
  cols_use <- colnames(cnv)[colnames(cnv) != 'chr_desc']
  cnv_df <- cnv %>%
    tidyr::pivot_longer(cols = all_of(cols_use), names_to = 'cluster', values_to = 'copy_number')
  dim(cnv_df)
  head(cnv_df)
  
  chrs <- unlist(lapply(strsplit(cnv_df$chr_desc,'_'), function(x){
    return(x[1])
  }))
  starts <- unlist(lapply(strsplit(cnv_df$chr_desc,'_'), function(x){
    return(x[2])
  }))
  ends <- unlist(lapply(strsplit(cnv_df$chr_desc,'_'), function(x){
    return(x[3])
  }))
  cnv_df$chr <- as.character(chrs)
  cnv_df$start <- as.character(starts)
  cnv_df$end <- as.character(ends)
  if(!grepl('chr',cnv_df$chr[1])){
    cnv_df <- cnv_df %>%
      dplyr::mutate(chr=paste0("chr", chr))
  }
  
  gene_cn <- get_gene_coordinates(cnv_df)
  dim(gene_cn)
  gene_cn$entrezgene <- NULL
  gene_cn$ensembl_gene_id[1]
  chrs <- c(as.character(1:22), "X")
  meta_genes <- annotables::grch38 %>%
    dplyr::select(ensembl_gene_id = ensgene, gene_symbol=symbol, chr)  %>%
    dplyr::filter(ensembl_gene_id %in% unique(gene_cn$ensembl_gene_id) & chr %in% chrs) #  
  #   inner_join(deg_df) %>%
  gene_cn <- gene_cn %>%
    dplyr::left_join(meta_genes, by='ensembl_gene_id')
  sum(gene_cn$gene_symbol %in% c('MIB-3','MIB-1','MKI67'))
  dim(gene_cn)
  View(head(gene_cn))
  data.table::fwrite(gene_cn, paste0(input_dir, 'median_cnv_genes.csv'))
  
  colnames(gene_cn)
  gene_cn1 <- gene_cn
  gene_cn1$desc <- paste0(gene_cn1$ensembl_gene_id, gene_cn1$cluster)
  gene_cn1 <- gene_cn1[!duplicated(gene_cn1$desc),]
  
  gene_cn2 <- gene_cn1 %>%
    dplyr::filter(gene_symbol %in% de_genes$gene_symb)
  
  # t <- gene_cn2 %>%
  #   dplyr::select(cluster, copy_number, ensembl_gene_id, gene_symbol) %>%
  #   dplyr::group_by(cluster, gene_symbol) %>%
  #   dplyr::summarise(varr=var(copy_number))
  # dim(t)  
  # t
  gene_cn2 <- gene_cn2 %>%
    dplyr::select(cluster, copy_number, ensembl_gene_id, gene_symbol) %>%
    pivot_wider(names_from = 'cluster', values_from = 'copy_number')  
  dim(gene_cn2)
  sum(gene_cn2$B!=gene_cn2$C)
  View(gene_cn2)
  deg_fn <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA919/deg_analysis_X7/SA919_cls1_SupraSpinal_C_vs_cls_0_Primary_B/de_significant_genes.csv.gz'
  de_genes <- data.table::fread(deg_fn)
  dim(de_genes)  
  de_genes$avg_log2FC
  de_genes <- de_genes %>%
    dplyr::filter(avg_log2FC>0.25)
  
  
  
}
