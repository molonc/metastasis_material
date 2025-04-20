suppressPackageStartupMessages({
  library(tidyverse)
  library(annotables)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  # require(scales)
  # require(ggrepel)
  # library("edgeR")
  # library("DESeq2")
  require(stringr)
  library(scran)
  library(SingleCellExperiment)
  library(tximport)
})
options(dplyr.summarise.inform = FALSE)
options(tidyverse.quiet = TRUE)


## To Do: get 80% of samples from cis genes 
## and comparing with the same number of sampling trans genes. 

get_bootstrap_stat_sampling <- function(cis_genes, trans_genes, 
                                        sampling_fraction=0.7, nsamples=1000){
  set.seed(42)
  # if(is.null(genome_genes)){
  #   # ref_dif <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/'
  #   # genome_genes_df <- read.csv(paste0(ref_dif, 'Symbol_ensembl.csv'), check.names = F, stringsAsFactors = F)  
  #   # # dim(genome_genes_df)
  #   # genome_genes <- unique(genome_genes_df$Symbol) # entire genes set
  #   
  #   ref <- annotables::grch38 %>%
  #     dplyr::select(ensgene,symbol) %>%
  #     dplyr::rename(gene_id=ensgene)
  #   ref <- ref[!duplicated(ref$gene_id),]
  #   genome_genes <- ref$gene_id
  # }
  nb_sampled_cis <- round(length(cis_genes) * sampling_fraction)
  nb_sampled_cis <- length(cis_genes)
  # cis_samples <- lapply(1:nsamples, function(i) sample(cis_genes, size=nb_sampled_cis, replace = T))
  # trans_samples <- lapply(1:nsamples, function(i) sample(trans_genes, size=nb_sampled_cis, replace = T))
  cis_samples <- list()
  trans_samples <- list()
  for(i in seq(nsamples)){
    cis_samples[[i]] <- sample(cis_genes, size=nb_sampled_cis, replace = T)
    trans_samples[[i]] <- sample(trans_genes, size=nb_sampled_cis, replace = T)
  }
  
  stat_df <- tibble::tibble()
  
  for(k in seq(nsamples)){
    out_stat <- ks.test(cis_samples[[k]], trans_samples[[k]], alternative='less')
    
    out_vals <- tibble::tibble('idx'=k, 'stat'=out_stat$statistic, 'p_val'=out_stat$p.value)
    stat_df <- dplyr::bind_rows(stat_df, out_vals)
  }
  
  stat_df <- stat_df %>%
    dplyr::mutate(is_signf=
      case_when(
        p_val < 0.05 ~ 'T',
        TRUE ~ 'F'
      )
    )
  # summary(as.factor(stat_df$is_signf))
  summary_df <- stat_df %>%
    dplyr::group_by(is_signf) %>%
    dplyr::summarise(nb_val=n()) %>%
    dplyr::mutate(pct_signf=round(100*nb_val/dim(stat_df)[1]))
  summary_df
  # return(list(CI=as.numeric(r['our']),pval=as.numeric(pval)))
  return(summary_df)
}

get_DE_genes_DESeq2 <- function(dds, DE_comp=c("Metastasis","Primary"),
                                filter_genes=F, min_total_exp_samples=10){
  
  print(dim(dds))
  print(sizeFactors(dds))
  print(colData(dds))
  ## ----prefilter----------------------------------------------------------------
  ## filter genes function pose an issue for next step, need to debug it later
  ## set filter_genes=F first
  if(filter_genes){
    # smallestGroupSize <- 3 # group untreated with 3 samples is the smallest group
    # keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
    keep <- rowSums(counts(dds) >= min_total_exp_samples)
    dds <- dds[keep,]
    print(dim(dds))  
  }
  print(colnames(dds))
  ## ----factorlvl----------------------------------------------------------------
  dds$condition <- factor(dds$condition, levels = DE_comp)
  
  ## ----relevel------------------------------------------------------------------
  # dds$condition <- relevel(dds$condition, ref = "untreated")
  
  ## ----droplevels---------------------------------------------------------------
  # dds$condition <- droplevels(dds$condition)
  
  ## ----deseq--------------------------------------------------------------------
  ## size factors are noted in colData column
  dds <- DESeq(dds)
  contrast_DE_comp <- c("condition")
  contrast_DE_comp <- c(contrast_DE_comp, DE_comp)
  res <- results(dds, contrast=contrast_DE_comp) #c("condition","treated","untreated")
  res$ensembl_gene_id <- rownames(res)
  return(res)
}

## Size factor is kept in coldata$size_factor
create_DESeq2_obj <- function(coldata, cts, use_existing_size_factor=T){
  # library("DESeq2")
  cts <- cts[, rownames(coldata)]
  print(all(rownames(coldata) == colnames(cts)))
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~ condition)
  print(dds)
  
  ## Supplying existing size factor for each sample to deseq2 object
  if(use_existing_size_factor){
    # is.null(sizeFactors(object)
    sizeFactors(dds) <- coldata$size_factor  
  }
  print(sizeFactors(dds))
  return(dds)
  
}
get_vst_clustering_results <- function(dds){
  library("pheatmap")
  library("RColorBrewer")
  vsd <- vst(dds, blind=FALSE)
  # rld <- rlog(dds, blind=FALSE)
  print(head(assay(vsd), 3))
  # this gives log2(n + 1)
  # ntd <- normTransform(dds)
  # # BiocManager::install("vsn")
  # library("vsn")
  # meanSdPlot(assay(ntd))
  # meanSdPlot(assay(vsd))
  # meanSdPlot(assay(rld))
  
  ## ----heatmap------------------------------------------------------------------
  # df <- as.data.frame(colData(dds)[,c("mainsite","main_clone")])
  # p_clust <- pheatmap(assay(vsd), cluster_rows=FALSE, show_rownames=FALSE,
  #          cluster_cols=TRUE, annotation_col=df)
  # p_clust
    
  mainsite_cols <- c('red','blue')
  names(mainsite_cols) <- c('Metastasis','Primary')
  clone_cols <- c('#FC8D62','#8DA0CB')
  names(clone_cols) <- c('B','C')
  
  exp_cols <- c('#808000','#C1E1C1')
  names(exp_cols) <- c('main_exp','mixing_exp')
  
  top_anno = ComplexHeatmap::HeatmapAnnotation(Main_Site = factor(colData(dds)[,"mainsite"]),
                                               Main_Clone=factor(colData(dds)[,"main_clone"]),
                                               Experiment=factor(colData(dds)[,"experiment"]),
                                               col = list(Main_Site=mainsite_cols,
                                                          Main_Clone=clone_cols, 
                                                          Experiment=exp_cols))
  
  p <- ComplexHeatmap::Heatmap(as.matrix(assay(vsd)), na_col = "white",
                               # col = col_fun,
                               show_column_names=T,
                               show_row_names = F,
                               cluster_rows=F,
                               cluster_columns=T,
                               clustering_distance_columns = "pearson",
                               name = "Exp", 
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
  saveRDS(p, paste0(save_dir,'complexheatmap_clustering_allsamples.rds'))
  png(paste0(save_dir,datatag,"_hierarchial_clusters.png"), height = 2*1000, width=2*1000, res = 2*72)
  print(p)
  dev.off()  
  p1 <- grid::grid.grabExpr(ComplexHeatmap::draw(p, padding = unit(c(2, 2, 2, 2), "mm")))
  ggsave(  
    filename = paste0(save_dir,datatag,"_hierarchial_clusters.svg"),  
    plot = p1,  
    height = 8, width = 11, dpi = 150)
  # select <- order(rowMeans(counts(dds, normalized=TRUE)),
  #                 decreasing=TRUE)[1:20]
  # df <- as.data.frame(colData(dds)[,c("condition","type")])
  # pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
  #          cluster_cols=FALSE, annotation_col=df)
  # pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
  #          cluster_cols=FALSE, annotation_col=df)
  # pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
  #          cluster_cols=FALSE, annotation_col=df)
  
  ## ----sampleClust--------------------------------------------------------------
  sampleDists <- dist(t(assay(vsd))) ## transform to the format in ML, samples in rows, features in cols
  
  ## ----figHeatmapSamples, fig.height=4, fig.width=6-----------------------------
  
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsd$mainsite, vsd$main_clone, vsd$experiment, sep="-")
  colnames(sampleDistMatrix) <- NULL
  # colnames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  p <- pheatmap(sampleDistMatrix,
                clustering_distance_rows=sampleDists,
                clustering_distance_cols=sampleDists,
                col=colors)
  return(p)
  # sampleDistMatrix: full mtx
  # sampleDists: half mtx
  
}

# deg_df: contains 2 columns of gene_symbol, and logFC
# gmt_fn <- paste0(ref_dif,'pathway_set/h.all.v7.0.symbols.gmt')

get_fgsea_pathways <- function(deg_df, save_dir, base_name, gmt_fn){
  save_dir <- paste0(save_dir, base_name, '_pathways/')
  if(!dir.exists(save_dir)){
    dir.create(save_dir)
  }
  print(save_dir)
  deg_df <- deg_df %>%
    dplyr::select(gene_symbol, logFC) %>%
    na.omit() %>%
    distinct() %>%
    group_by(gene_symbol) %>%
    summarize(logFC=mean(logFC))
  
  deg_stat <- deg_df$logFC
  names(deg_stat) <- deg_df$gene_symbol
  
  ref_set <- fgsea::gmtPathways(gmt_fn)
  gsea_out <- fgsea(pathways=ref_set, stats=deg_stat)  #, scoreType = "pos"  
  gsea_out$datatag <- base_name
  gsea_out$signf_genes <- ''
  for(i in rep(1:length(gsea_out$pathway), 1)){
    signf_genes <- unlist(gsea_out$leadingEdge[[i]])
    gsea_out$nb_signf_genes[i] <- length(unlist(signf_genes))
    gsea_out$signf_genes[i] <- paste(signf_genes, collapse=',')  
  }
  gsea_out$leadingEdge <- NULL
  gsea_out <- as.data.frame(gsea_out)
  gsea_out <- gsea_out %>%
    dplyr::filter(pval<0.05)
  data.table::fwrite(gsea_out,paste0(save_dir, "signf_pathways_",base_name,".csv"))
  return(gsea_out)
}
get_gprofiler_pathways_obsgenes_ <- function(obs_genes_symb, save_dir, datatag, 
                                            custom_id=NULL, pathway_fn=NULL, 
                                            save_data=F, correction_func='gSCS'){
  library(gprofiler2)
  if(is.null(pathway_fn)){
    # pathway_fn = '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/pathway_set/c2.cp.kegg.v7.1.symbols.gmt'  
    pathway_fn = '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/pathway_set/h.all.v7.0.symbols.gmt'  
  }
  
  ref_set <- fgsea::gmtPathways(pathway_fn)
  # for(s in names(ref_set)){ ## just quick check the size of each set
  #   print(s)
  #   print(length(ref_set[[s]]))
  # }
  if(is.null(custom_id)){
    custom_id <- gprofiler2::upload_GMT_file(pathway_fn)  
    print('Custom id from gprofiler is: ')
    print(custom_id)
  }
  stat <- NULL
  ## correction_method: one of 'fdr', 'gSCS', 'bonferroni' #gSCS is the most popular one
  gostres <- gprofiler2::gost(list(obs_genes_symb), organism = custom_id, 
                              correction_method=correction_func) #, significant = F
  dim(stat)
  # stat <- stat %>%
  #   arrange(p_value)
  # View(stat)
  if(!is.null(gostres$result)){
    stat <- gostres$result
    cols_use <- c('p_value','intersection_size','precision','recall','term_id')
    stat <- stat %>%
      dplyr::select(all_of(cols_use)) %>%
      dplyr::rename(reference_set=term_id, nb_signif_genes=intersection_size) %>%
      dplyr::filter(p_value<0.05) # just to be sure
    
    # Get pathway genes 
    for(i in seq(nrow(stat))){
      pw_set <- stat$reference_set[i]
      ref_genes <- ref_set[[pw_set]]
      # obs_genes <- deg_df$gene_symbol
      intersect_genes <- intersect(obs_genes_symb, ref_genes)
      stat$signif_genes[i] <- paste(intersect_genes, collapse=',')
    }
    if(save_data){
      # added_time <- gsub(':','',format(Sys.time(), "%Y%b%d_%X"))
      # data.table::fwrite(stat, paste0(save_dir, 'pathways_',datatag,'_',added_time,'.csv.gz'))  
      data.table::fwrite(stat, paste0(save_dir, 'pathways_',datatag,'.csv.gz'))  
    }
  }  
  return(stat)
  
}  
get_gprofiler_pathways <- function(genes_df, save_dir, datatag, 
                                   custom_id=NULL, pathway_fn=NULL, save_data=F){
  library(gprofiler2)
  if(is.null(pathway_fn)){
    # pathway_fn = '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/pathway_set/c2.cp.kegg.v7.1.symbols.gmt'  
    pathway_fn = '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/pathway_set/h.all.v7.0.symbols.gmt'  
  }
  
  ref_set <- fgsea::gmtPathways(pathway_fn)
  if(is.null(custom_id)){
    custom_id <- gprofiler2::upload_GMT_file(pathway_fn)  
    print('Custom id from gprofiler is: ')
    print(custom_id)
  }
  
  pathway_stat <- tibble::tibble()
  for(gm in unique(genes_df$gene_type_module)){
    genes_use <- genes_df %>%
      dplyr::filter(gene_type_module==gm) %>%
      dplyr::pull(gene_symbol)
    ## correction_method: one of 'fdr', 'gSCS', 'bonferroni'
    gostres <- gprofiler2::gost(list(genes_use), organism = custom_id, correction_method='gSCS')
    if(!is.null(gostres$result)){
      stat <- gostres$result
      cols_use <- c('p_value','intersection_size','precision','recall','term_id')
      stat <- stat %>%
        dplyr::select(all_of(cols_use)) %>%
        dplyr::rename(reference_set=term_id, nb_signif_genes=intersection_size)
      stat$gene_type_module <- gm
      
      # Get pathway genes 
      for(i in seq(nrow(stat))){
        pw_set <- stat$reference_set[i]
        ref_genes <- ref_set[[pw_set]]
        # obs_genes <- deg_df$gene_symbol
        intersect_genes <- intersect(genes_use, ref_genes)
        stat$signif_genes[i] <- paste(intersect_genes, collapse=',')
      }
      if(save_data){
        data.table::fwrite(stat, paste0(save_dir, 'gene_module_',gm,'_pathways.csv'))  
      }
      pathway_stat <- dplyr::bind_rows(pathway_stat, stat)
    }
  }
  pathway_stat$datatag <- datatag
  return(pathway_stat)
}

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
    bulk_fn1 <- paste0(input_dir, s,'/abundance.h5')
    bulk_fn2 <- paste0(input_dir, s,'_trim/abundance.h5')
    bulk_fn3 <- paste0(input_dir, s,'/abundance.tsv')
    bulk_fn <- ifelse(file.exists(bulk_fn1),bulk_fn1,
                      ifelse(file.exists(bulk_fn2),bulk_fn2, 
                             ifelse(file.exists(bulk_fn3),bulk_fn3,NULL)))
   
   if(!is.null(bulk_fn)){
      print(paste0('Sample: ',s,', read file: ', bulk_fn))
      input_fns <- c(input_fns,bulk_fn)
      sids <- c(sids, s)
   }else{
      stop(paste0('!!!Do not exist file, check input abundance.h5 file for sample',s))
    }
  }
  sample_ids <- sids
  print(sample_ids)
  ## input files
  # input_fns <- paste0(input_dir, sample_ids,'/abundance.h5') #.h5 or .tsv both are fine
  # input_fns <- paste0(input_dir, sample_ids,'/abundance.h5') #.h5 or .tsv both are fine
  
  ## Just for testing
  # input_fns <- input_fns[1]
  # sample_ids <- sids[1]
  ref_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/tx2gene_sets/'
  ref_dir <- '/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/metastasis_material/materials/bulkRNAseq/mapping_transcripts2gene/'
  ref_tx2gene_fn <- paste0(ref_dir, 'tx2gene_updated_23_Oct_2023.csv.gz')
  tx2gene <- data.table::fread(ref_tx2gene_fn, header=T) %>% as.data.frame()
  # print(head(tx2gene))
  # dim(tx2gene)
  tx2gene <- tx2gene[!duplicated(tx2gene$TXNAME),]
  
  ## Just for testing the txtImport read file function
  # for(i in seq(1:length(input_fns))){
  #   print(input_fns[i])
  #   txi <- tximport(input_fns[i], type = "kallisto", tx2gene = tx2gene, 
  #                   ignoreAfterBar = T, ignoreTxVersion=T)
  #   print(dim(txi$counts))
  #   
  # }
  txi <- tximport(input_fns, type = "kallisto", tx2gene = tx2gene, 
                  ignoreAfterBar = T, ignoreTxVersion=T)
    # class(txi$counts)
  # head(txi$abundance)
  # class(txi$length)
  colnames(txi$counts) <- sample_ids
  # colnames(txi$length) <- sample_ids
  # colnames(txi$abundance) <- sample_ids
  
  print(paste0('# of genes: ',dim(txi$counts)[1], ', #samples: ',dim(txi$counts)[2]))
  txi$counts <- round(txi$counts,0) # cut off numbers after ,
  
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

# Returned normalized matrix and size factor
normalize_by_size_factor_v2 <- function(df_counts_fn, datatag, save_dir){
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
  sf <- sizeFactors(sce)
  names(sf) <- colnames(sce)
  print('Size factor for each sample applied Scran size factor function: ')
  print(sf)
  sf_df <- tibble::tibble(sample=names(sf), size_factor=sf)
  data.table::fwrite(sf_df, paste0(save_dir, datatag, '_sizefactors_Scran.csv.gz'))
  
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
  # print(all(colnames(raw_counts)==colnames(norm_counts_mtx)))
  
  stat <- tibble::tibble(library_size=c(s, s2), sample_type=c(names(s), names(s2)),
                         processed=c(rep('counts',length(s)), rep('normalized',length(s2))))
  data.table::fwrite(stat, paste0(save_dir, datatag, '_sizefactor_eval_raw_norm.csv.gz'))
  library("ggplot2")
  p <- ggplot(stat, aes(sample_type, y=library_size, fill=sample_type)) + #
    geom_col(width = 0.3) + 
    facet_grid(~processed) + 
    theme_bw() + 
    theme(axis.text.x = element_text(size=10, angle = 90),
          legend.position = 'none')
  
  png(paste0(save_dir,datatag,"_sizeFactor_raw_normalized.png"), 
      height = 550, width=2*750, res = 2*72)
  print(p)
  dev.off()
  
  # saveRDS(p, paste0(save_dir,datatag,"_sizeFactor_raw_normalized.rds"))
  
  return(list(size_factor=sf, norm_mtx = norm_counts_mtx))
}

# Returned only normalized matrix
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
  
  saveRDS(sce, paste0(save_dir, datatag, '_sizefactor_normalized.rds'))
  meta_samples_df <- as.data.frame(colData(sce))
  data.table::fwrite(meta_samples_df, paste0(save_dir, datatag, '_sizefactors.csv.gz'))
  
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
