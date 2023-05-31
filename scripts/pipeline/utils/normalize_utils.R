# BiocManager::install('extrafont')
suppressPackageStartupMessages({
  # require("optparse")
  require("scater")
  # require("argparse")
  require("SingleCellExperiment")
  require("stringr")
  require("tidyverse")
  # require("scran")
  require("Seurat")
  # require("sctransform")
  require("dplyr")
  
})

# library(extrafont)
# font_import(prompt=F) # import all your fonts
# fonts()
# library(extrafont)
# font_import(prompt=F, paths ='/usr/share/fonts/truetype/myfonts/') # import Helvetica font
# fonts()

# options(scipen = 999)  # disable scientific notation in R
# options(scipen = 0) # numeric option in R
gene_var_stat <- function(sce){
  dec <- scran::modelGeneVar(sce)
  plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
  curve(metadata(dec)$trend(x), col="blue", add=TRUE)
  
  dec3 <- modelGeneVar(sce, block=sce$library_id)
  per.block <- dec3$per.block
  par(mfrow=c(3, 2))
  for (i in seq_along(per.block)) {
    decX <- per.block[[i]]
    plot(decX$mean, decX$total, xlab="Mean log-expression", 
         ylab="Variance", main=names(per.block)[i])
    curve(metadata(decX)$trend(x), col="blue", add=TRUE)
  }
  
  par(mfrow=c(4, 2))
  for (i in seq_along(per.block)) {
    decX <- per.block[[i]]
    plot(decX$mean, decX$total, xlab="Mean log-expression", 
         ylab="Variance", main=names(per.block)[i])
    curve(metadata(decX)$trend(x), col="blue", add=TRUE)
  }
}

# input data for ridge regression - batch effect correction: 
# + normalized mtx csv.gz of genes x cells
# meta data: clone, and batch columns
# output: a mtx of genes x cells correction 
# create a sce file of corrected output
generate_input_data_ridge_regression <- function(norm_sce, save_dir, datatag='SA'){
  # require("SingleCellExperiment")
  
  meta_info <- as.data.frame(colData(norm_sce))
  if(!'cell_id' %in% colnames(meta_info)){
    meta_info$cell_id <- rownames(meta_info)
  }
  # dim(meta_info)
  ts <- colnames(meta_info)[grepl('*treat*',colnames(meta_info))]
  if(ts !='treatmentSt'){
    meta_info$treatmentSt <- meta_info[,ts]
  }
  # meta_info$batch <- sample_df[meta_info$sample,'batch']
  meta_info <- meta_info %>%
    dplyr::select(batch, cell_id, clone, treatmentSt, sample, library_id, timepoint)
  dim(meta_info)
  # SA609
  # meta_info1 <- meta_info %>%
  #   dplyr::filter(clone!='R')
  # dim(meta_info1)
  # meta_info2 <- meta_info %>%
  #   dplyr::filter(clone=='R')
  # dim(meta_info2)
  # meta_info <- dplyr::bind_rows(meta_info1, meta_info2)
  meta_fn <- paste0(save_dir, datatag, "_meta_info.csv")
  write.csv(meta_info, meta_fn, row.names = F, quote = F)
  print(paste0('Save meta data as: ',meta_fn))
  
  # save_dir <- paste0(save_dir,'BE_mtx_v2/')
  norm_sce <- norm_sce[,meta_info$cell_id]
  dim(meta_info)
  dim(norm_sce)
  norm_data <- as.data.frame(as.matrix(logcounts(norm_sce)))
  mtx_fn <- paste0(save_dir, datatag, "_norm.csv.gz")
  data.table::fwrite(norm_data, file=mtx_fn, row.names=T, quote=F)
  print(dim(norm_data))
  print(paste0('Save normalized mtx output as: ',mtx_fn))
}


##corrected_mtx_fn: output of batch correction 
## norm_sce: normalized exp sce of SCTransform, scran,....

convert_corrected_mtx_to_sce <- function(norm_sce, save_dir, 
                             corrected_mtx_fn=NULL, datatag='SA', 
                             return_sce=F, save_data=T){
  ## require("SingleCellExperiment")
  if(is.null(corrected_mtx_fn)){
    # corrected_mtx_fn <- paste0(save_dir, datatag, "_norm.csv.gz")
    corrected_mtx_fn <- paste0(save_dir, datatag, "_corrected.csv.gz")
  }
  corrected_data <- data.table::fread(corrected_mtx_fn) %>% as.data.frame()
  rownames(corrected_data) <- corrected_data$V1
  corrected_data$V1 <- NULL
  norm_sce <- norm_sce[rownames(corrected_data),colnames(corrected_data)]
  print(dim(norm_sce))
  print(dim(corrected_data))
  corrected_sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=counts(norm_sce),
                                                              logcounts=as.matrix(corrected_data)),
                                                              colData = colData(norm_sce), 
                                                              rowData(norm_sce))
  print(corrected_sce)
  if(save_data){
    saveRDS(corrected_sce, paste0(save_dir, datatag, "_norm_BE_sce.rds"))
  }
  if(return_sce){
    return(corrected_sce)  
  }
  
}



sce_normalize_size_factors <- function(sce, min_size=100, rlog=FALSE, exprs="counts"){
  # library(scater)
  ## BiocManager::install('scran')
  # library("scran")
  
  print("Quick clustering")
  if(min_size < dim(sce)[2]){
    qclust <- scran::quickCluster(sce, min.size = min_size, assay.type=exprs)
    print("Compute sum factors")
    sce <- scran::computeSumFactors(sce, clusters = qclust, assay.type=exprs)
    sce$size_factor <- sizeFactors(sce)
    # print(sce$size_factor[1:5])
    print("Normalize data")
    # scater version 1.14.6
    
    # count --> normcounts
    # normcounts --> logcounts 
    # String containing an assay name for storing the output normalized values. 
    # Defaults to "logcounts" when log=TRUE and "normcounts" otherwise.
    sce_normalized <- scater::logNormCounts(sce, log=rlog, exprs_values=exprs, size_factors=sce$size_factor)
    
  } else{
    print("Small nb cells in this library")
    sce_normalized <- logNormCounts(sce, log=rlog, exprs_values=exprs, size_factors=NULL)
  }
  # ?logNormCounts
  # scater version < 1.14.6
  # sce_normalized <- normalize(sce, return_log=rlog, exprs_values=exprs)
  print(assayNames(sce_normalized))
  return(sce_normalized)
}


#' @title Combind several \code{SingleCellExperiment} objects from different batches/experiments
#'
#' @param cut_off_overall A numeric vector  indicating the cut-off for the proportion of 
#' a gene is expressed overall data
# Based on scMerge sce_cbind func, https://sydneybiox.github.io/scMerge/reference/sce_cbind.html
sce_cbind_func <- function(sce_list, genes_df, cut_off_overall = 0.05, exprs = c("counts", "logcounts","normcounts"), 
                      colData_names = NULL, meta_data=NULL) {
  # n_batch <- length(sce_list)
  # method = "intersect"
  genes_df <- rowData(sce)
  if(sum(rownames(sce_list[[1]])==rownames(sce_list[[2]]))!=33538){
    stop('Double check input data')
  }
  
  
  assay_list <- list()
  for (i in seq_len(length(exprs))) {
    assay_list[[i]] <- do.call(cbind, lapply(sce_list, 
                                             function(y) assay(y, exprs[i])))
  }
  names(assay_list) <- exprs
  colData_list <- do.call(DelayedArray::rbind, 
                          lapply(sce_list, function(y) colData(y)[, colData_names, drop = FALSE]))
  sce_combine <- SingleCellExperiment::SingleCellExperiment(assays = assay_list, 
                                                            colData = colData_list)
  
  # genes_df <- as.data.frame(genes_df)
  # dim(genes_df)
  # genes_df$Symbol[genes_df$Symbol!=rownames(sce_combine)]
  # rownames(sce_combine)[rownames(sce_combine)!=genes_df$Symbol]
  # genes_df <- genes_df %>%
  #   dplyr::rename(geo_strand=strand,
  #                 geo_chr=chr,
  #                 geo_start=start,
  #                 geo_end=end)
  if(sum(genes_df$ID==rownames(sce_combine))==33538){
    rowData(sce_combine) <- as.matrix(genes_df)  
  }
  
  # for(c in colnames(genes_df)){
  #   print(class(genes_df[,c]))
  #   print(c)
  #   # if(class(genes_df[,c])=='CompressedGRangesList'){
  #   #   print(c)
  #   # }
  # }
  
  # if (is.null(batches)) {
  #   batches <- paste("batch", seq_len(n_batch))
  # }
  # meta_data$batch_info <- gsub('^CHIP','C_',meta_data$batch_info)
  # print(meta_data$batch_info[1:4])
  
  # meta_data$pdxid <- gsub('^X0847_','',meta_data$pdxid)
  # print(meta_data$pdxid[1:4])
  
  # sce_combine$mouse_id <- meta_data$mouse_id
  # sce_combine$pdxid <- meta_data$pdxid
  # sce_combine$passage <- meta_data$passage
  # sce_combine$treatmentSt <- meta_data$treatmentSt
  # sce_combine$cell_cycle_phases <- meta_data$cphases
  # sce_combine$library_label <- meta_data$library_label
  # sce_combine$batch_info <- meta_data$batch_info
  # sce_combine$batch <- rep(batches, unlist(lapply(sce_list, ncol)))  #id
  # sce_combine$libid <- rep(libid, unlist(lapply(sce_list, ncol)))
  # sce_combine$treatment_status <- rep(treatment_status, unlist(lapply(sce_list, ncol)))
  
  if(cut_off_overall > 0){
    zero_cbind <- DelayedArray::rowMeans(assay(sce_combine, exprs[1]) == 0)
    sce_combine <- sce_combine[names(zero_cbind[zero_cbind <= (1 - cut_off_overall)]), ]
  }
  # rownames(sce_combine)[1]
  # dim(sce_combine)
  print(paste0("Dim sce combine: ",dim(sce_combine)[1],' ',dim(sce_combine)[2]))
  print(paste0("sce combine assay name: ",assayNames(sce_combine)))
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
  sce_combine <- SingleCellExperiment::SingleCellExperiment(assays = assay_list, 
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


# UMI counts in each cell were normalized via division by the cell's total UMI count, scaled by the median UMI count across all 
# cells and log2-transformed after the addition of 1.

median_normalization <- function(counts_df){
  cs <- colSums(t)
  mval <- median(cs) 
  # zseq <- apply(t, 2, function(x) x /colSums(x) * median()
  for(c in rep(1:ncol(t),1)){
    print(c)
    t[,c] <- t[,c] / cs[c]
    t[,c] <- t[,c] * mval
  }
  norm_df <- log2(t + 1)
  return(norm_df)
  
} 

# based on SCTransform source
# See more at: https://github.com/ChristophH/sctransform
# umi: matrix 
inspect_rawdata <- function(umi){
  # dups <- which(duplicated(rownames(umi)))
  gene_attr <- data.frame(mean = rowMeans(umi), 
                          detection_rate = rowMeans(umi > 0),
                          var = apply(umi, 1, var))
  gene_attr$log_mean <- log10(gene_attr$mean)
  gene_attr$log_var <- log10(gene_attr$var)
  rownames(gene_attr) <- rownames(umi)
  # cell_attr <- data.frame(n_umi = colSums(umi),
  #                         n_gene = colSums(umi > 0))
  # rownames(cell_attr) <- colnames(umi)
  
  p <- ggplot(gene_attr, aes(log_mean, log_var)) + 
    geom_point(alpha=0.3, shape=16) + 
    geom_density_2d(size = 0.3) +
    geom_abline(intercept = 0, slope = 1, color='red')
  p
  return(p)
  
}

plot_raw_variance <- function(sce_raw, save_dir, datatag){
  output_fn <- paste0(save_dir,'ref_genes_total.rds')
  if(file.exists(output_fn)){
    print("Output exist in folder, reloading...")
    ref_genes_total <- readRDS(output_fn)
    return(ref_genes_total)
  }
  ref_stable_genes <- load_stable_genes_ls(sce_raw, save_dir, 0.8)  #meta_genes_HK=meta_genes_HK, meta_genes_scSEG
  
  meta_cells <- as.data.frame(colData(sce_raw))
  meta_cells$cell_id <- colnames(sce_raw)
  meta_cells$series[1]
  meta_cells$treatmentSt <- get_treatment_status(meta_cells$series)
  meta_cells <- meta_cells %>%
    dplyr::select(cell_id, treatmentSt)
  
  gene_attr_HK <- inspect_rawdata_std_mean(as.matrix(assay(sce_raw, 'counts')), ref_stable_genes$meta_genes_HK, 
                                           save_dir, 'HK genes variance', T)
  
  gene_attr_scSEG <- inspect_rawdata_std_mean(as.matrix(assay(sce_raw, 'counts')), ref_stable_genes$meta_genes_scSEG, 
                                              save_dir, 'Ref scMerge scSEG genes variance', T)
  
  gene_attr_scSEG_our <- inspect_rawdata_std_mean(as.matrix(assay(sce_raw, 'counts')), ref_stable_genes$meta_genes_scSEG_our, 
                                              save_dir, 'Custom scMerge scSEG genes variance', T)
  
  filtered_gene_attr_HK <- gene_attr_HK$gene_attr %>%
    dplyr::filter(abs(log_var)<=0.5 & classified_gene!='2_other')
  
  filtered_gene_attr_scSEG <- gene_attr_scSEG$gene_attr %>%
    dplyr::filter(abs(log_var)<=0.5 & classified_gene!='2_other')
  
  filtered_our_gene_attr_scSEG <- gene_attr_scSEG_our$gene_attr %>%
    dplyr::filter(abs(log_var)<=0.5 & classified_gene!='2_other')
  
  write.csv(filtered_gene_attr_HK, paste0(save_dir,'filtered_HK_genes.csv'), quote=F, row.names = F)
  write.csv(filtered_gene_attr_scSEG, paste0(save_dir,'filtered_scSEG_genes.csv'), quote=F, row.names = F)
  write.csv(filtered_our_gene_attr_scSEG, paste0(save_dir,'filtered_ours_scSEG_genes.csv'), quote=F, row.names = F)
  
  ptotal <- cowplot::plot_grid(gene_attr_HK$p, gene_attr_scSEG$p, gene_attr_scSEG_our$p, align = 'h', ncol = 3)
  saveRDS(ptotal, paste0(save_dir,datatag,"_HK_genes_versus_ref_scSEG_var_means.rds"))
  
  png(paste0(save_dir,datatag,"_HK_genes_versus_ref_scSEG_var_means.png"), height = 2*300, width=2*1000,res = 2*72)
  print(ptotal)
  dev.off()
  ref_genes_total <- list(ref_stable_genes=ref_stable_genes,
       filtered_gene_attr_HK=filtered_gene_attr_HK,
       filtered_gene_attr_scSEG=filtered_gene_attr_scSEG,
       filtered_our_gene_attr_scSEG=filtered_our_gene_attr_scSEG)
  saveRDS(ref_genes_total, output_fn)
  return(ref_genes_total)
  
}
# umi: matrix with row names is gene names
# classified_gene is a dataframe with genes as rownames and gene_type column
# plottitle='Stably expressed genes verification'
inspect_rawdata_std_mean <- function(umi, classified_gene, save_dir, plottitle='', save_output=F){
  # dups <- which(duplicated(rownames(umi)))
  rowvar <- apply(umi, 1, var)
  rowsd <- apply(umi, 1, sd)
  gene_attr <- data.frame(mean = DelayedArray::rowMeans(umi), 
                          detection_rate = DelayedArray::rowMeans(umi > 0),
                          var = rowvar,
                          std = rowsd)
  gene_attr$log_mean <- log2(gene_attr$mean)
  gene_attr$log_std <- log2(gene_attr$std)
  gene_attr$log_var <- log10(gene_attr$var)
  # gene_attr$log_var <- log10(gene_attr$var)
  # rownames(gene_attr) <- rownames(umi)
  gene_attr$gene_id <- rownames(umi)
  gene_attr$classified_gene <- classified_gene[gene_attr$gene_id,'gene_type']
  plottitle <- paste0(nrow(gene_attr[gene_attr$classified_gene!='2_other',]),' ',plottitle)
  # cell_attr <- data.frame(n_umi = colSums(umi),
  #                         n_gene = colSums(umi > 0))
  # rownames(cell_attr) <- colnames(umi)
  if(save_output){
    write.csv(gene_attr, paste0(save_dir,'gene_attr_',gsub(' ','_',plottitle),'.csv'), quote=F, row.names = F)
  }
  
  p <- ggplot(gene_attr, aes(x=log_mean, y=log_var,color=classified_gene)) + 
    geom_point(alpha=1, shape=1, size=2) +
    geom_hline(yintercept=0.2, linetype="dashed", color = "red", size=0.4) + 
    geom_hline(yintercept=-0.2, linetype="dashed", color = "red", size=0.4) +
    theme(plot.title = element_text(color="black", size=12, hjust = 0.5),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            axis.text.y = element_text(color="black", size=11),
            axis.text.x = element_text(color="black", size=11),
            axis.title = element_text(color="black", size=11),
            legend.title = element_text(color="black", size=8), 
            legend.text = element_text(color="black", size=8))
    # geom_density_2d(size = 0.3) #+
    # geom_abline(intercept = 0, slope = 1, color='red')
  p <- p + labs(x="Mean exp across cells (log2)",y="Variance across cells (log10)",
                title=plottitle)
  
  png(paste0(save_dir, gsub(' ','_',plottitle), "_std_mean.png"), height = 2*330, width=2*550,res = 2*72)
  print(p)
  dev.off()
  
  return(list(gene_attr=gene_attr, p=p))
}


format_theme_plot <- function(p, title='Seurat Clustering'){
  p <- p + labs(title)
  p <- p + theme(
    # legend.title = element_text(size=10), 
    legend.text = element_text(size=10), 
    plot.title = element_text(color="black", size=17, hjust = 0.5)
  )
  return(p)
}



# This function calls sctransform::vst. The sctransform package is available at https://github.com/ChristophH/sctransform. 
# Use this function as an alternative to the NormalizeData, FindVariableFeatures, ScaleData workflow. Results are saved 
# in a new assay (named SCT by default) with counts being (corrected) counts - ATTENTION, data being log1p(counts), 
# scale.data being pearson residuals; sctransform::vst intermediate results are saved in misc slot of new assay. 
# Here using sctransform function in Seurat3 package
normalize_SCTransform <- function(sce, output_dir, datatag, return_data=F, output_fn=NULL){
  if(is.null(output_fn)){
    output_fn <- paste0(output_dir, datatag,'_sctransform_normalized_21_mar_2023.rds')
  }
    # if(file.exists(output_fn)){
  #   print("Output exist in folder, reloading...")
  #   sce2 <- readRDS(output_fn)
  #   # return(sce2)
  # }
  # else{
  #   
  # }
  print("Removing mito genes prior to normalization...")
  print(rowData(sce)$Symbol[1])
  mito_genes <- str_detect(rowData(sce)$Symbol, "^MT\\-")
  print(sum(mito_genes==TRUE))
  # ribo_genes <- str_detect(rowData(sce)$Symbol, "^RP(L|S)")  # or ^RP[L|S]?
  # print(sum(ribo_genes==TRUE))
  # sce <- sce[(!mito_genes) & (!ribo_genes),]
  meta_genes <- rowData(sce)
  irange_cols <- c("seqnames","ranges", "strand", "start", "end", "width", "element")
  cols_use <- colnames(meta_genes)[!colnames(meta_genes) %in% irange_cols]
  meta_genes <- meta_genes[,cols_use]
  rowData(sce) <- as.data.frame(meta_genes)
  
  sce <- sce[!mito_genes,]
  print(dim(sce))
  
    # meta_genes <- meta_genes %>%
  #   dplyr::select(-start, -strand, -end)
  rownames(meta_genes) <- meta_genes$ID
  print(meta_genes$ID[1])
  # sce_counts <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=counts(sce)), 
  #                                                          colData = colData(sce), 
  #                                                          rowData = rowData(sce))
  # print(dim(sce_counts))
  srt <- Seurat::as.Seurat(sce, counts = "counts", data=NULL) #set to NULL if only normalized data are present #, assay = "RNA", project = "SingleCellExperiment"
  
  ncells_use <- 10000 # nb cells use to find residual
  if(dim(sce)[2]<ncells_use){
    ncells_use <- dim(sce)[2]
  }
  # ?Seurat::SCTransform
  srt <- Seurat::SCTransform(object = srt, verbose = T, 
                             do.scale=F, 
                             return.only.var.genes = T, 
                             assay = "originalexp", new.assay.name = "SCT", ncells = ncells_use)
  # print(summary(srt[["RNA"]]@data))
  print(dim(srt[["SCT"]]@data))
  print(dim(srt[["SCT"]]@counts))
  # t <- srt[["SCT"]]@counts
  # norm_mtx <- srt[["SCT"]]@data
  # # sc <- srt[["RNA"]]@scale.data
  # max(norm_mtx)
  # min(norm_mtx)
  # print('Scale sctransform normalized data...')
  # srt <- ScaleData(srt, features = rownames(srt), assay = 'SCT')
  norm_sce <- Seurat::as.SingleCellExperiment(srt, assay = 'SCT')
  print('Debug...')
  print(rownames(norm_sce)[1:10])
  # rownames(sce2) <- rownames(sce)
  # rownames(sce2) <- rownames(srt)
  # rowData(sce2) <- rowData(sce)
  print(dim(norm_sce))
  if(sum(rownames(sce)==rownames(norm_sce))==dim(sce)[1]){
    rowData(norm_sce) <- dplyr::bind_cols(as.data.frame(rowData(sce)),
                                      as.data.frame(rowData(norm_sce)))
    
  }else{
    stop('Double check genes name!')
  }
    # normcounts(sce2) <- srt[["SCT"]]@scale.data  # save scale data to normcounts assay
  print(assayNames(norm_sce))
  # print(dim(srt[["SCT"]]@scale.data))
  # print(max(normcounts(sce2)))
  # # min(normcounts(sce2))
  # print('Save normalized data...')
  saveRDS(norm_sce, output_fn)
  
  # res = 0.3
  # dims = 1:15
  # srt <- RunPCA(object = srt, verbose = FALSE)
  # srt <- FindNeighbors(srt, dims = dims)
  # srt <- FindClusters(srt, resolution = res)
  # srt <- RunUMAP(object = srt, dims = dims, verbose = FALSE) #umap.method = "umap-learn", metric='correlation' 
  # 
  # p21 <- DimPlot(srt, reduction = "umap")
  # p21 <- format_theme_plot(p=p21, title='Seurat Clustering')
  # 
  # p22 <- DimPlot(srt, reduction = "umap", group.by = 'library_id')
  # p22 <- format_theme_plot(p=p22, title='library_id')
  # 
  # # p23 <- DimPlot(srt, reduction = "umap", group.by = 'Grouping')
  # # p23 <- format_theme_plot(p=p23, title='Grouping')
  # # 
  # # p24 <- DimPlot(srt, reduction = "umap", group.by = 'Site_origin')
  # # p24 <- format_theme_plot(p=p24, title='Site_origin')
  # # 
  # # pumap <- p21 + p22 + patchwork::plot_layout(ncol=2)
  # png(paste0(output_dir,datatag,"_clusters_umap.png"), height = 2*350, width=2*450,res = 2*72)
  # print(p21)
  # dev.off()
  # 
  # pumap2 <- p23 + p24 + patchwork::plot_layout(ncol=2)
  # png(paste0(output_dir,datatag,"_clusters_umap_2.png"), height = 2*350, width=2*800,res = 2*72)
  # print(pumap2)
  # dev.off()
  if(return_data){
    # return(list(sce=sce, sctransform=sce2))
    return(norm_sce)
  }  
}  
viz_umap <- function(sce, output_dir, datatag='SA', res = 0.3, dims = 1:25, nb_hvg=3000,
                     return_seurat_obj=F){
  # library(Seurat)
  # sce$library_id <- gsub('SCRNA10X_SA_CHIP','',sce$library_id)
  # sce$Grouping <- ifelse(sce$Site_origin=="Tumor_Recur","Primary",sce$Grouping)
  reducedDim(sce) <- NULL
  srt <- Seurat::as.Seurat(sce, counts = "counts", 
                           data="logcounts") #set to NULL if only normalized data are present
  
  # srt <- Seurat::as.Seurat(sce, counts = "counts", 
  #                          data="logcounts", assay = "RNA", 
  #                          project = "SingleCellExperiment") #set to NULL if only normalized data are present
  print(dim(srt))
  # rownames(srt)[1]
  # class(srt[["RNA"]]@data)
  
  srt <- Seurat::FindVariableFeatures(srt, selection.method = "vst", nfeatures = nb_hvg)
  
  # all.genes <- rownames(srt)
  # srt <- ScaleData(srt, features = all.genes) # scale for all genes, can scale by hvg genes
  srt[["SCT"]]@scale.data <- as.matrix(logcounts(sce))
  srt <- RunPCA(object = srt, verbose = FALSE)
  # srt <- JackStraw(srt, num.replicate = 100)
  # srt <- ScoreJackStraw(srt, dims = 1:20)
  # JackStrawPlot(srt, dims = 1:15)
  # ElbowPlot(srt)
  srt <- FindNeighbors(srt, dims = dims)
  srt <- FindClusters(srt, resolution = res)
  srt <- RunUMAP(object = srt, dims = dims, verbose = FALSE) #umap.method = "umap-learn", metric='correlation'
  
  p21 <- DimPlot(srt, reduction = "umap", group.by ='seurat_clusters')
  # p21
  meta_df <- srt@meta.data
  stat <- meta_df %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::summarise(nb_cells=n())
  meta_df <- meta_df %>% left_join(stat, by='seurat_clusters')
  meta_df$seurat_clusters <- paste0('Cluster ',meta_df$seurat_clusters,' (N=',meta_df$nb_cells,')')
  # srt@meta.data <- meta_df
  summary(srt@meta.data$seurat_clusters)
  # srt@meta.data$library_id[1]
  # cap <- paste(unique(meta_df$seurat_clusters), collapse='; ')
  # p21 <- p21 + labs(caption = cap, title = 'SCRNA10X_SA_CHIP0229_001')
  colnames(meta_df)
  
  # meta_df$batch_info <- NULL
  p21 <- DimPlot(srt, reduction = "umap", group.by = 'seurat_clusters')
  # p22 <- DimPlot(srt, reduction = "umap", group.by = 'lid')
  # p23 <- DimPlot(srt, reduction = "umap", group.by = 'Site_origin')
  # p24 <- DimPlot(srt, reduction = "umap", group.by = 'pdxid')
  # p22 <- format_theme_plot(p=p22, title='library_id')
  # srt$treatment <- srt$treat
  # srt$clone <- 'unassigned'
  # srt$passage <- NULL
  
  # cells_use <- intersect(rownames(clonealign_df),colnames(srt))
  # srt@meta.data[cells_use,'clone'] <- clonealign_df[cells_use,'clone']
  # srt@meta.data[cells_use,'passage'] <- clonealign_df[cells_use,'passage']
  
  # srt@meta.data <- srt@meta.data %>% left_join(meta_cells, by=c('sample'='mouse_id'))
  # srt@meta.data <- meta_df
  # rownames(srt@meta.data) <- colnames(srt)
  # p22 <- DimPlot(srt, reduction = "umap", group.by = 'treatment')
  # p22 <- format_theme_plot(p=p22, title='Treatment')
  # 
  # p22 <- DimPlot(srt, reduction = "pca", group.by = 'treatment')
  # p22 <- format_theme_plot(p=p22, title='Treatment')
  # p22 <- DimPlot(srt, reduction = "umap", group.by = 'clone_desc')
  # p22 <- format_theme_plot(p=p22, title='Subclone')

  # p23 <- DimPlot(srt, reduction = "umap", group.by = 'Grouping')
  # p23 <- format_theme_plot(p=p23, title='Grouping')
  
  # p23 <- DimPlot(srt, reduction = "umap", group.by = 'clone')
  # p23 <- format_theme_plot(p=p23, title='Clone')

  # p24 <- DimPlot(srt, reduction = "umap", group.by = 'passage')
  # p24 <- format_theme_plot(p=p23, title='passage')
  # p24 <- DimPlot(srt, reduction = "umap", group.by = 'Site_origin')
  # p24 <- format_theme_plot(p=p24, title='Site_origin')

  # pumap <- p21 + p22 + patchwork::plot_layout(ncol=2)
  # png(paste0(output_dir,datatag,"_clusters_umap_1.png"), height = 2*350, width=2*800,res = 2*72)
  # print(pumap)
  # dev.off()
  # 
  # pumap2 <- p23 + p24 + patchwork::plot_layout(ncol=2)
  # png(paste0(output_dir,datatag,"_clusters_umap_2.png"), height = 2*350, width=2*800,res = 2*72)
  # print(pumap2)
  # dev.off()

  # pumap_total <- p21 + p24 + p22 + p23  + patchwork::plot_layout(ncol=2)
  pumap_total <- p21 + p22  + patchwork::plot_layout(ncol=2)
  out_fig <- paste0(output_dir,datatag,"_clusters_umap1.png")
  png(out_fig, height = 2*300, width=2*800,res = 2*72)
  print(pumap_total)
  dev.off()
  
  # pumap_total <- p23 + p24  + patchwork::plot_layout(ncol=2)
  # out_fig <- paste0(output_dir,datatag,"_clusters_umap2.png")
  # png(out_fig, height = 2*300, width=2*800,res = 2*72)
  # print(pumap_total)
  # dev.off()
  # print('Save output as: ')
  # print(out_fig)
  if(return_seurat_obj){
    return(srt) #srt@meta.data
  }
}

normalize_Seurat <- function(sce, input_dir, output_dir, datatag, return_data=F){
  sce_counts <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=counts(sce)), 
                                                            colData = colData(sce), 
                                                           rowData = rowData(sce))
  print(dim(sce_counts))
  print('Normalizing data...')
  srt <- Seurat::as.Seurat(sce_counts, counts = "counts", data=NULL) #set to NULL if only normalized data are present
  srt <- NormalizeData(srt, normalization.method = "LogNormalize", scale.factor = 100000)
  srt <- ScaleData(srt, features = rownames(srt))
  # print(summary(srt[["RNA"]]@scale.data))
  # norm_mtx <- srt[["RNA"]]@data
  # sc <- srt[["RNA"]]@scale.data
  
  sce2 <- as.SingleCellExperiment(srt)
  normcounts(sce2) <- srt[["RNA"]]@scale.data  # save scale data to normcounts assay
  print(assayNames(sce2))
  print(dim(sce2))
  print('Save normalized data...')
  saveRDS(sce2, paste0(output_dir,datatag,'_seurat_normalized_v3.rds'))
  if(return_data){
    return(sce2)
  }  
  
  
  # max(normcounts(sce2))
  # min(normcounts(sce2))
  # max(logcounts(sce2))
  # srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = 2000)
  # srt <- RunPCA(srt, features = VariableFeatures(object = srt), verbose = F)
  # srt <- FindNeighbors(srt, dims = 1:15)
  # srt <- FindClusters(srt, resolution = 0.5)
  # srt <- RunUMAP(srt, dims = 1:15)
  
  # p <- DimPlot(srt, reduction = "umap")
  # p1 <- DimPlot(srt, reduction = "umap", group.by = 'series')
  # pumap <- p + p1 + patchwork::plot_layout(ncol=2)
  # png(paste0(output_dir, datatag, "_clusters_umap_Seurat_normalization.png"), height = 2*400, width=2*720,res = 2*72)
  # print(pumap)
  # dev.off()
}
# Basic seurat clustering, based on scrnafitness source code
# 
get_seurat_clusters <- function(base_name, sce, output_file, save_dir, 
                                res = 0.3, dims = 1:15, saveSrt = TRUE, plotRD=TRUE) {

  # High proportion of mito, can affect output of clustering 
  # Remove mito, ribo genes first
  mito_genes <- str_detect(rowData(sce)$Symbol, "^MT\\-")
  sum(mito_genes==TRUE)
  
  ribo_genes <- str_detect(rowData(sce)$Symbol, "^RP(L|S)")  # or ^RP[L|S]?
  sum(ribo_genes==TRUE)
  sce1 <- sce[(!mito_genes) & (!ribo_genes), ]
  print(paste0('Removing mito and ribo genes: ',dim(sce1)[1],'_',dim(sce1)[2]))
  
  
  # reduction = "mnn",
  # algorithm = 1,
  # save_dir <- paste0(dirname(output_file),'/')
  if (!file.exists(save_dir)){
    dir.create(save_dir)
  }
  srt <- as.Seurat(sce1, counts = "counts", data = "logcounts")  
  
  
  if(sum(isUnique(colnames(sce))==FALSE) > 0){
    print("Cells name are not unique, make cells name unique first")
    return(NULL)
  }
  
  srt <- FindVariableFeatures(object = srt, verbose = FALSE)
  var_genes <- VariableFeatures(object = srt)
  print(paste0('Nb var genes: ',length(var_genes)))
  saveRDS(var_genes, file = paste0(save_dir, base_name,"_var_genes.rds")) 
  
  srt <- ScaleData(object = srt, verbose = FALSE)
  srt <- RunPCA(object = srt, verbose = FALSE)
  srt <- FindNeighbors(srt, dims = dims)
  srt <- FindClusters(srt, resolution = res)
  print("Clustering results: ")
  print(summary(srt@meta.data$seurat_clusters))
  
  
  
  if(plotRD){
    srt <- RunTSNE(object = srt, verbose = FALSE)
    srt <- RunUMAP(object = srt, dims = dims, verbose = FALSE) #umap.method = "umap-learn", metric='correlation' 
    p11 <- DimPlot(srt, reduction = "tsne")
    p11 <- format_theme_plot(p=p11, title='seurat_clusters')
    
    p12 <- DimPlot(srt, reduction = "tsne", group.by = 'treatmentSt')
    p12 <- format_theme_plot(p=p12, title='treatmentSt')
    
    # p13 <- DimPlot(srt, reduction = "tsne", group.by = 'Site_origin')
    # p13 <- p13 + labs(title='Site Origin')
    p14 <- DimPlot(srt, reduction = "tsne", group.by = 'batch_info')
    p14 <- format_theme_plot(p=p14, title='Batch Info')
    
    p15 <- DimPlot(srt, reduction = "tsne", group.by = 'passage')
    p15 <- format_theme_plot(p=p15, title='passage')
    
    p21 <- DimPlot(srt, reduction = "umap")
    p21 <- format_theme_plot(p=p21, title='Seurat Clustering')
    
    p22 <- DimPlot(srt, reduction = "umap", group.by = 'treatmentSt')
    p22 <- format_theme_plot(p=p22, title='treatmentSt')
    
    
    # p23 <- DimPlot(srt, reduction = "umap", group.by = 'Site_origin')
    # p23 <- p23 + labs(title='Site Origin')
    p24 <- DimPlot(srt, reduction = "umap", group.by = 'batch_info')
    p24 <- format_theme_plot(p=p24, title='Batch Info')
    
    p25 <- DimPlot(srt, reduction = "umap", group.by = 'passage')
    p25 <- format_theme_plot(p=p25, title='passage')
    
    ls_tsne <- list(p11, p12, p14, p15)
    ls_umap <- list(p21, p22, p24, p25)
    ptsne <- CombinePlots(plots = ls_tsne, ncol=2)
    pumap <- CombinePlots(plots = ls_umap, ncol=2)
    # saveRDS(ls_tsne, file = paste0(save_dir, "plots_tsne.rds")) 
    # saveRDS(ls_umap, file = paste0(save_dir, "plots_umap.rds")) 
    # # pdf(paste0(save_dir,"clusters_treatment_tsne.pdf"), height = 4, width=10)
    png(paste0(save_dir,base_name,"_clusters_tnse.png"), height = 2*500, width=2*620,res = 2*72)
    print(ptsne)
    dev.off()
    png(paste0(save_dir,base_name,"_clusters_umap.png"), height = 2*500, width=2*620,res = 2*72)
    print(pumap)
    dev.off()
  }
  
  if(saveSrt){
    saveRDS(srt, file = paste0(save_dir, base_name,"_clustering_srt.rds")) 
  }
  
  sce2 <- as.SingleCellExperiment(srt)
  sce$cluster_label <- sce2$seurat_clusters
  sce1$cluster_label <- sce2$seurat_clusters
  s = sum(colnames(sce)==colnames(sce2))
  s1 = sum(colnames(sce1)==colnames(sce2))
  print(paste0("Verification: sce: ", s, ' sce1:',s1))
  saveRDS(sce, file = output_file)
  saveRDS(sce1, file = paste0(save_dir, base_name,"_sans_mito_ribo_genes.rds"))
  return(sce)
  
  # if(return_SCE) {
  #   sce2 <- as.SingleCellExperiment(srt)
  #   sce$cluster_label <- sce2$seurat_clusters
  #   if(saveSCE){
  #     saveRDS(sce, file = output_file)
  #   }
  #   # return(sce)
  # } else {
  #   # return(sce$seurat_clusters)
  #   return(srt)
  # }
  
}


run_QC_Cluster <- function(base_name, sce, output_dir=NULL){
  
  mito_genes <- str_detect(rowData(sce)$Symbol, "^MT-")  #"^MT\\-"
  sum(mito_genes==TRUE)
  
  ribo_genes <- str_detect(rowData(sce)$Symbol, "^RP(L|S)")  # or ^RP[L|S]?
  sum(ribo_genes==TRUE)
  
  dim(sce)
  
  # For scater version < 1.14.6
  # sce <- calculateQCMetrics(sce, exprs_values = "logcounts")
  # sce$total_features_by_logcounts
  # For scater version >= 1.14.6 
  per.cell <- perCellQCMetrics(sce, subsets=list(Mito=mito_genes, Ribo=ribo_genes))
  colData(sce) <- cbind(colData(sce), per.cell)
  # print(summary(sce$total_counts))  # sum
  # print(summary(summary(sce$pct_counts_Mito)))
  # print(summary(summary(as.factor(sce$pct_counts_Ribo))))
  colnames(colData(sce))
  
  meta_data <- as.data.frame(colData(sce))
  
  ls_features <- c("detected","total",
                   "subsets_Mito_percent","subsets_Ribo_percent")
  plots <- list()
  idx = 0 
  obs_feature <- "cluster_label"
  for(f in ls_features){
    idx <- idx+1
    p <- ggplot(meta_data, aes_string(x=obs_feature, y=ls_features[idx], color=obs_feature)) + 
      geom_violin()
    p <- p + geom_boxplot(width=0.1)
    plots[[idx]] <- p 
  }
  length(plots) 
  # saveRDS(plots, file=paste0(output_dir,'plots.rds'))
  # pqc <- CombinePlots(plots, ncol=length(plots))
  pqc <- cowplot::plot_grid(plotlist=plots,  align = "h", ncol = 2)
  png(paste0(output_dir,base_name,"_qc_clusters.png"), height = 2*450, width=2*600, res = 2*72)
  print(pqc)
  dev.off()
} 

# based on SCTransform source
# See more at: https://github.com/ChristophH/sctransform
check_normalization_output <- function(umi, normalized_df){
  min_cells = 5
  gmean_eps = 1
  genes_cell_count <- rowSums(umi > 0)
  genes <- rownames(umi)[genes_cell_count >= min_cells]
  genes_log_gmean <- log10(row_gmean(umi, eps = gmean_eps))
  gene_attr <- data.frame(
    detection_rate = genes_cell_count[genes] / ncol(umi),
    gmean = 10 ^ genes_log_gmean,
    variance = row_var(umi))
  if (ncol(normalized_df) > 0) {
    gene_attr$residual_mean = rowMeans(normalized_df)
    gene_attr$residual_variance = sctransform:::row_var(normalized_df)
  }
  p1 <- ggplot(gene_attr, aes(residual_mean)) + geom_histogram(binwidth=0.01)
  p2 <- ggplot(gene_attr, aes(residual_variance)) + geom_histogram(binwidth=0.1) + 
    geom_vline(xintercept=1, color='red') + xlim(0, 10)
  p3 <- ggplot(gene_attr, aes(log10(gmean), residual_variance)) + geom_point(alpha=0.3, shape=16) +
    geom_density_2d(size = 0.3)
  plots <- list(p1, p2, p3)
  return(plots)
}



plot_stable_genes_exp_v2 <- function(sce, stable_genes, use_raw=F, exprs='logcounts', plottitle='Batch effect estimation', 
                                     subtitle='', xlabel='', ylabel="Mean exp of stably expressed genes", yl=c(0,1.6), legend_visible=F){
  print(length(stable_genes))
  stable_genes <- intersect(rownames(sce), as.character(stable_genes))
  print(length(stable_genes))
  sce_SEG <- sce[stable_genes,]
  print(dim(sce_SEG))
  if(use_raw){
    logcounts(sce_SEG) <- log2(counts(sce_SEG)+1)
  }
  
  # summary_SEG <- data.frame(mean_seg = colMeans(logcounts(sce_SEG)), row.names=colnames(sce_SEG))
  # # head(summary_SEG)
  # sce_SEG$mean_seg <- summary_SEG[colnames(sce_SEG),"mean_seg"]
  meta_data <- as.data.frame(colData(sce_SEG))
  meta_data$mean_exp <- DelayedArray::colMeans(assay(sce_SEG, exprs))
  # meta_data$batch_info <- as.factor(meta_data$batch_info)
  meta_data$library_id <- gsub("SCRNA10X_SA_CHIP","",meta_data$library_id)
  if(!is.null(yl)){
    meta_data$mean_exp <- ifelse(meta_data$mean_exp>yl[2],yl[2],
                                 ifelse(meta_data$mean_exp<yl[1],yl[1],meta_data$mean_exp))
  }  
  
  meta_data$series <- factor(meta_data$series, levels = gtools::mixedsort(unique(meta_data$series)))
  # if(legend_visible){
  #   p <- plot_variation_function_with_legend(meta_data, xstring="series", ystring="mean_exp", 
  #                                plottype="batch", plottitle,
  #                                xlabel, ylabel, subtitle)
  # }else{
  #   p <- plot_variation_function(meta_data, xstring="series", ystring="mean_exp", 
  #                                plottype="batch", plottitle,
  #                                xlabel, ylabel, subtitle)
  # }
  # # p <- p + scale_y_continuous(breaks = round(seq(yl[1], yl[2], by = 0.2),1))
  # if(!is.null(yl)){
  #   p <- p + ylim(yl[1],yl[2])
  # }
  # list(p=p,meta_data=meta_data)
  return(meta_data)
  
}
plot_stable_genes_exp_metastasis <- function(sce, stable_genes, 
                                             use_raw=F, exprs='logcounts', yl=NULL){
  if(length(stable_genes)!=sum(stable_genes %in% rownames(sce))){
    stop('Double check input data')
  }
  # stable_genes <- intersect(rownames(sce), as.character(stable_genes))
  sce <- sce[stable_genes,]
  print(dim(sce))
  if(use_raw){
    logcounts(sce) <- log2(counts(sce)+1)
  }
  
  # summary_SEG <- data.frame(mean_seg = colMeans(logcounts(sce_SEG)), row.names=colnames(sce_SEG))
  # # head(summary_SEG)
  # sce_SEG$mean_seg <- summary_SEG[colnames(sce_SEG),"mean_seg"]
  meta_data <- as.data.frame(colData(sce))
  meta_data$mean_exp <- DelayedArray::colMeans(assay(sce, exprs))
  # meta_data$batch_info <- as.factor(meta_data$batch_info)
  meta_data$library_id <- gsub("SCRNA10X_SA_CHIP","C_",meta_data$library_id)
  if(!is.null(yl)){
    meta_data$mean_exp <- ifelse(meta_data$mean_exp>yl[2],yl[2],
                                 ifelse(meta_data$mean_exp<yl[1],yl[1],meta_data$mean_exp))
  }  
  
  meta_data$library_id <- factor(meta_data$library_id, levels = gtools::mixedsort(unique(meta_data$library_id)))
  return(meta_data)
}
plot_stable_genes_exp <- function(sce, stable_genes, use_raw=F, exprs='logcounts', plottitle='Batch effect estimation', 
                                  xlabel='', ylabel="Mean exp of stably expressed genes", yl=c(0,1.6), scale=F){
  print(length(stable_genes))
  stable_genes <- intersect(rownames(sce), as.character(stable_genes))
  sce_SEG <- sce[stable_genes,]
  print(dim(sce_SEG))
  if(use_raw){
    logcounts(sce_SEG) <- log2(counts(sce_SEG)+1)
  }
  
  # summary_SEG <- data.frame(mean_seg = colMeans(logcounts(sce_SEG)), row.names=colnames(sce_SEG))
  # # head(summary_SEG)
  # sce_SEG$mean_seg <- summary_SEG[colnames(sce_SEG),"mean_seg"]
  meta_data <- as.data.frame(colData(sce_SEG))
  corrected <- assay(sce_SEG, exprs)
  if(scale){
    corrected <- corrected %>% 
      # transpose the matrix so genes are as columns
      t() %>% 
      # apply scalling to each column of the matrix (genes)
      scale() %>% 
      # transpose back so genes are as rows again
      t()
  }
  
  
  meta_data$mean_exp <- DelayedArray::colMeans(corrected)
  # meta_data$batch_info <- as.factor(meta_data$batch_info)
  meta_data$library_id <- gsub("SCRNA10X_SA_CHIP","",meta_data$library_id)
  p <- plot_variation_function(meta_data, xstring="series", ystring="mean_exp", 
                               plottype="series", plottitle,
                               xlabel, ylabel)
  if(!is.null(yl)){
    p <- p + ylim(yl[1],yl[2])
  }
  
  return(p)
  
}

plot_scatter_variation <- function(meta_data, xstring="treatment_status", ystring="gene_exp", plottype="batches_info", plottitle="Raw data",
                                    xlabel='', ylabel="Stably express genes") {
  
  p <- ggplot(meta_data, aes_string(x=xstring, y=ystring, fill=plottype)) +
       geom_violin() + 
       geom_jitter(position=position_jitter(0.5), size=0.3) 
  
  p <- p + labs(x=xlabel,y=ylabel,title=plottitle)
  p <- p + theme(legend.title = element_text(color="black", size=12, hjust = 0.5),
                 legend.text = element_text(color="black", size=12, hjust = 0.5),
                 plot.title = element_text(color="black", size=15, hjust = 0.5),
                 # legend.position = "none",
                 # axis.line = element_blank(),
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 # panel.border = element_blank(),
                 axis.text.x = element_text(color = "black", angle = 90, hjust = 1, size=12),
                 axis.text.y = element_text(color = "black", angle = 90, hjust = 1, size=12),
                 #axis.title = element_blank(),
                 # axis.ticks = element_blank()
  )
  # p
  return(p)
  
}


plot_variation_function_with_legend <- function(meta_data, xstring="treatment_status", 
                                                ystring="gene_exp", plottype="batches_info", plottitle="Raw data",
                                                xlabel='', ylabel="Stably express genes", 
                                                subtbl=NULL) {
  # subtbl='Using scTransform log normalized mtx'
  p <- ggplot(meta_data, aes_string(x=xstring, y=ystring, fill=plottype)) +
    stat_boxplot(aes_string(xstring, ystring), geom='errorbar', linetype=1, width=0.4)+ #whiskers
    geom_boxplot(aes_string(xstring, ystring),outlier.shape=1) #+
    # stat_summary(fun=mean, geom="point", size=2) #+
    # stat_summary(fun.data = mean_se, geom = "errorbar")
    
   
  p <- p + labs(x=xlabel,y=ylabel,title=plottitle, subtitle = subtbl) 
           
  p <- p + theme(plot.title = element_text(color="black", size=12, hjust = 0.5, face = "bold"), #
                 plot.subtitle = element_text(color="black", size=10, hjust = 0.5),
                 legend.title = element_text(color="black", size=10, hjust = 0.5),
                 legend.text = element_text(color="black", size=10, hjust = 0.5),
                 # legend.position = "none",
                 # axis.line = element_blank(),
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 # panel.border = element_blank(),
                 axis.text.x = element_text(color = "black", angle = 90, hjust = 1, size=9),
                 axis.text.y = element_text(color = "black", angle = 90, hjust = 1, size=10),
                 axis.title.y = element_text(color = "black", hjust=0.5, size=10),
                 axis.title.x = element_text(color = "black", hjust=0.5, size=8)
                 # axis.ticks = element_blank()
  ) 
  # guides(fill = guide_legend(nrow = 1, override.aes = list(size=0.1)))
  p <- p + guides(fill = guide_legend(override.aes = list(shape = 0, size=0.5)))
  # p
  return(p)
  
}

plot_variation_function <- function(meta_data, xstring="treatment_status", ystring="gene_exp", 
                                    plottype="batches_info", plottitle="Raw data",
                                    xlabel='', ylabel="Stably express genes",
                                    subtbl=NULL) {
  # meta_data <- meta_data[gtools::mixedsort(meta_data[,xstring]),]
  p <- ggplot(meta_data, aes_string(x=xstring, y=ystring, fill=plottype)) +
    stat_boxplot(aes_string(xstring, ystring), geom='errorbar', linetype=1, width=0.4)+ #whiskers
    geom_boxplot(aes_string(xstring, ystring),outlier.shape=1) #+
    #scale_y_continuous(labels = scientific, breaks=c(100000,200000, 30000))
    # stat_summary(fun=mean, geom="point", size=2) +
    # stat_summary(fun.data = mean_se, geom = "errorbar")
  
  p <- p + labs(x=xlabel,y=ylabel,title=plottitle, subtitle = subtbl) 
  # p <- p + theme(plot.title = element_text(color="black", size=10, hjust = 0.5, face="bold"),
  #                # legend.title = element_text(color="black", size=8, hjust = 0.5),
  #                # legend.text = element_text(color="black", size=8, hjust = 0.5),
  #                legend.position = "none",
  #                # axis.line = element_blank(),
  #                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #                # panel.border = element_blank(),
  #                axis.text.x = element_text(color = "black", angle = 90, hjust = 1, size=7),
  #                axis.text.y = element_text(color = "black", angle = 90, hjust = 1, size=7),
  #                axis.title = element_text(color = "black", hjust=0.5, size=9),
  #                # axis.ticks = element_blank()
  # )
  p <- p + theme(plot.title = element_text(color="black", size=12, hjust = 0.5, face = "bold"),
                 plot.subtitle = element_text(color="black", size=10, hjust = 0.5),
                 legend.title = element_text(color="black", size=10, hjust = 0.5),
                 legend.text = element_text(color="black", size=10, hjust = 0.5),
                 legend.position = "none",
                 # axis.line = element_blank(),
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 # panel.border = element_blank(),
                 axis.text.x = element_text(color = "black", angle = 90, hjust = 1, size=9),
                 axis.text.y = element_text(color = "black", angle = 90, hjust = 1, size=10),
                 axis.title.y = element_text(color = "black", hjust=0.5, size=10),
                 axis.title.x = element_text(color = "black", hjust=0.5, size=8)
                 # axis.ticks = element_blank()
  ) 
  
  # p
  return(p)
  
}


plot_variation_function_v3 <- function(meta_data, xstring="treatment_status", ystring="gene_exp", 
                                    plottype="batches_info", plottitle="Raw data",
                                    xlabel='', ylabel="Stably express genes",lg_pos="none") {
  if(ystring=="total_counts"){
    meta_data <- meta_data %>% 
      dplyr::group_by(series) %>%
      dplyr::mutate(total_counts=replace(total_counts, total_counts>quantile(total_counts,probs = 0.95), quantile(total_counts,probs = 0.95)))
  }
  
  # meta_data <- meta_data[gtools::mixedsort(meta_data[,xstring]),]
  p <- ggplot(meta_data, aes_string(x=xstring, y=ystring, fill=plottype)) +
    stat_boxplot(aes_string(xstring, ystring), geom='errorbar', linetype=1, width=0.4)+ #whiskers
    geom_boxplot(outlier.shape=1) +#+   #aes_string(xstring, ystring),
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))
  #scale_y_continuous(labels = scientific, breaks=c(100000,200000, 30000))
  # stat_summary(fun=mean, geom="point", size=2) +
  # stat_summary(fun.data = mean_se, geom = "errorbar")
  
  p <- p + labs(x=xlabel,y=ylabel,title=plottitle,fill="Sample Info") 
  
  # p <- p + theme(plot.title = element_text(color="black", size=12, hjust = 0.5, face = "bold"),
  #                plot.subtitle = element_text(color="black", size=10, hjust = 0.5),
  #                legend.title = element_text(color="black", size=10, hjust = 0.5),
  #                legend.text = element_text(color="black", size=10, hjust = 0.5),
  #                legend.position = "none",
  #                # axis.line = element_blank(),
  #                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #                # panel.border = element_blank(),
  #                axis.text.x = element_text(color = "black", angle = 90, hjust = 1, size=9),
  #                axis.text.y = element_text(color = "black", angle = 90, hjust = 1, size=10),
  #                axis.title.y = element_text(color = "black", hjust=0.5, size=10),
  #                axis.title.x = element_text(color = "black", hjust=0.5, size=8)
  #                # axis.ticks = element_blank()
  # ) 
  my_font <- "Helvetica"
  thesis_theme <- ggplot2::theme(
    text = element_text(size = 8, hjust = 0.5, family=my_font),
    axis.title.x = element_text(size=8, hjust = 0.5, family=my_font),
    axis.title.y = element_text(size=8, hjust = 0.5, family=my_font),
    axis.text.x = element_text(size=7, hjust = 0, family=my_font, angle = 90),
    # axis.text.x = element_blank(),
    axis.text.y = element_text(size=5, hjust = 0.5, family=my_font),
    plot.title = element_text(color = "black", size=9, face="bold", hjust=0.5, family=my_font),
    legend.title=element_text(size=7, hjust = 0.5, family=my_font),
    legend.text=element_text(size=7, hjust = 0, family=my_font),
    strip.text.x = element_text(size=9, family=my_font),
    strip.text.y = element_text(size=9, family=my_font),
    legend.spacing.x = unit(0.1, 'mm'),
    legend.spacing.y = unit(0.06, 'mm'),
    legend.key.height=unit(0.6,"line"),
    legend.position = lg_pos,
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  )
  p <- p + thesis_theme
  return(p)
  
}

get_metadata <- function(meta_data){
  meta_data$tm_st <- ifelse(grepl('TU$',meta_data$treatmentSt),'RxH',
                            ifelse(grepl('T$',meta_data$treatmentSt),'Rx','UnRx'))
  
  b <- unique(meta_data$batch)
  bl <- data.frame(batch=b, batch_label=paste0('B',seq(1:length(b))))
  meta_data <- meta_data %>% inner_join(bl, by='batch')
  # metacells_SA1035$series <- paste0(metacells_SA1035$series,'_',gsub("CHIP0","C",metacells_SA1035$batch))
  meta_data$series <- paste0(meta_data$timepoint,' ',
                             meta_data$tm_st, ' ',meta_data$batch_label)
  
  # meta_data$tag <- factor(meta_data$tag, levels = c("Raw Data","Seurat", "Scran","SCTransform"))
  meta_data$tag <- factor(meta_data$tag, levels = c("Raw Data","Scran","SCTransform","Batch Corrected"))
  meta_data$series <- factor(meta_data$series, levels = gtools::mixedsort(unique(meta_data$series)))
  return(meta_data)
}
plot_batch_effects <- function(meta_data, xstring="treatment_status", ystring="gene_exp", 
                                       plottype="batches_info", plottitle="Raw data",
                                       xlabel='', ylabel="Stably express genes",
                                       lg_pos="none",save_dir) {
  
  # meta_data <- meta_data[gtools::mixedsort(meta_data[,xstring]),]
  # meta_data$tm_st <- ifelse(grepl('TU$',meta_data$treatmentSt),'RxH',
  #                            ifelse(grepl('T$',meta_data$treatmentSt),'Rx','UnRx'))
  # meta_data$tm_st <- ifelse(grepl('XU$',meta_data$treatmentSt),'RxH',
  #                           ifelse(grepl('X$',meta_data$treatmentSt),'Rx','UnRx'))
  # b <- unique(meta_data$batch)
  # bl <- data.frame(batch=b, batch_label=paste0('B',seq(1:length(b))))
  # meta_data <- meta_data %>% inner_join(bl, by='batch')
  # metacells_SA1035$series <- paste0(metacells_SA1035$series,'_',gsub("CHIP0","C",metacells_SA1035$batch))
  # meta_data$series <- paste0(meta_data$timepoint,' ',
  #                            meta_data$tm_st, ' ',meta_data$batch_label)
  
  # meta_data$tag <- factor(meta_data$tag, levels = c("Raw Data","Seurat", "Scran","SCTransform"))
  # meta_data$series <- factor(meta_data$series, levels = gtools::mixedsort(unique(meta_data$series)))
  # 
  p <- ggplot(meta_data, aes_string(x=xstring, y=ystring, fill=plottype)) +
    stat_boxplot(aes_string(xstring, ystring), geom='errorbar', linetype=1, width=0.4)+ #whiskers
    geom_boxplot(outlier.shape=1) +#+   #aes_string(xstring, ystring),
    facet_grid(rows = vars(tag)) + #, space='free', , scales="free"
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) 
  #scale_y_continuous(labels = scientific, breaks=c(100000,200000, 30000))
  # stat_summary(fun=mean, geom="point", size=2) +
  # stat_summary(fun.data = mean_se, geom = "errorbar")
  
  p <- p + labs(x=xlabel,y=NULL,title=plottitle,subtitle = ylabel, fill="Batch") 
  
  # p <- p + theme(plot.title = element_text(color="black", size=12, hjust = 0.5, face = "bold"),
  #                plot.subtitle = element_text(color="black", size=10, hjust = 0.5),
  #                legend.title = element_text(color="black", size=10, hjust = 0.5),
  #                legend.text = element_text(color="black", size=10, hjust = 0.5),
  #                legend.position = "none",
  #                # axis.line = element_blank(),
  #                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #                # panel.border = element_blank(),
  #                axis.text.x = element_text(color = "black", angle = 90, hjust = 1, size=9),
  #                axis.text.y = element_text(color = "black", angle = 90, hjust = 1, size=10),
  #                axis.title.y = element_text(color = "black", hjust=0.5, size=10),
  #                axis.title.x = element_text(color = "black", hjust=0.5, size=8)
  #                # axis.ticks = element_blank()
  # ) 
  my_font <- "Helvetica"
  thesis_theme <- ggplot2::theme(
    text = element_text(size = 8, hjust = 0.5, family=my_font),
    axis.title.x = element_text(size=8, hjust = 0.5, family=my_font),
    axis.title.y = element_text(size=8, hjust = 0.5, family=my_font),
    axis.text.x = element_text(size=8, hjust = 0, family=my_font, angle = 90),
    # axis.text.x = element_blank(),
    axis.text.y = element_text(size=5, hjust = 0.5, family=my_font),
    plot.title = element_text(color="black", size=9, face="bold", hjust=0.5, family=my_font),
    plot.subtitle = element_text(size=7, hjust=0.5, family=my_font),
    legend.title=element_text(size=7, hjust = 0.5, family=my_font, angle = 90),
    legend.text=element_text(size=7, hjust = 0.5, family=my_font),
    strip.text.x = element_text(size=9, family=my_font),
    strip.text.y = element_text(size=9, family=my_font),
    legend.spacing.x = unit(0.05, 'mm'),
    legend.spacing.y = unit(0.05, 'mm'),
    legend.key.height=unit(0.3,"line"),
    legend.position = lg_pos,
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  )
  p <- p + thesis_theme
  # if(!is.null(yl)){
  #   p <- p + ylim(yl[1],yl[2])
  # }
  # saveRDS(p, paste0(save_dir,"batch_effect_HK_plots.rds"))
  # ggsave(paste0(save_dir,"batch_effect_HK_plots.pdf"),
  #        plot = p,
  #        height = 8.5,
  #        width = 2.5,
  #        useDingbats=F
  # )#dpi = 150
  
  
  
  return(p)
  
}


get_mean_value_by_treatment_condition <- function(sce_scMerge, treated_cond, exprs){
  norm_data_treated <- assay(sce_scMerge[,treated_cond], exprs)
  norm_data_treated <- as.data.frame(as.matrix(norm_data_treated))
  norm_data_treated[norm_data_treated==0] <- NA
  norm_data_treated$gene <- rownames(sce_scMerge)
  # norm_data$gene <- rownames(observed_sce)
  dim(norm_data_treated)
  meta_cells_treated <- as.data.frame(colData(sce_scMerge[,treated_cond]))
  meta_cells_treated$cell_id <- colnames(sce_scMerge[,treated_cond])
  meta_cells_treated <- meta_cells_treated %>%
    dplyr::select(cell_id,treatmentSt)
  
  exp_mean_treated <- norm_data_treated %>% 
    # convert to long format
    pivot_longer(!gene, names_to = "cell_id", values_to = "exp")  %>% 
    # join with sample info table
    left_join(meta_cells_treated, by = ("cell_id")) %>%
    # filter to retain only genes of interest
    # filter(gene %in% observed_genes) %>% 
    # for each gene
    group_by(gene) %>% 
    # scale the cts column
    mutate(mean_exp = (exp - mean(exp, na.rm=T))/sd(exp))
  
  
  gexp_treated <- exp_mean_treated %>% 
    # for each gene, summary by treatment_status
    group_by(gene, treatmentSt) %>%
    # calculate the mean (scaled) cts
    summarise(mean_exps= mean(mean_exp, na.rm=T),
              nrep = n()) %>% 
    ungroup()
  
  print(dim(gexp_treated))
  print(head(gexp_treated))
  # dim(sce_scMerge)
  return(gexp_treated)
}

load_stable_genes_ls <- function(sce_raw, save_dir, cutoff_hk=0.8){
  
  input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
  output_dir <- paste0(input_dir,'rnaseq_v6/normalization_evaluation/')
  scHK_df <- readxl::read_excel(paste0(output_dir,'HKGenes.xlsx'))
  print(dim(scHK_df))
  # head(scHK_df)
  scHK_df$stability_index <- scHK_df$`Stability index`
  scHK_df <- scHK_df %>%
    dplyr::filter(stability_index >= quantile(x=stability_index, probs = cutoff_hk))
  
  # scHK_df <- scHK_df %>%
  #   dplyr::filter(stability_index <= quantile(x=stability_index, probs = cutoff_hk))
  
  print(dim(scHK_df))
  print(summary(scHK_df$stability_index))
  
  scHK_df <- scHK_df[scHK_df$GeneSymbol %in% rowData(sce_raw)$Symbol,]
  scHK_df <- scHK_df[order(scHK_df$stability_index,decreasing = T),] 
  # scHK_df <- scHK_df[1:500,]
  # stable_genes <- intersect(scHK_df$GeneSymbol, rowData(sce)$Symbol)
  meta_genes_HK <- data.frame(ens_gene=rowData(sce_raw)$ID, gene_symbol=rowData(sce_raw)$Symbol, row.names = rowData(sce_raw)$ID, stringsAsFactors = F)
  
  meta_genes_HK$gene_type <- '2_other'
  meta_genes_HK[meta_genes_HK$gene_symbol %in% scHK_df$GeneSymbol, 'gene_type'] <- '1_HK_gene'
  print(summary(as.factor(meta_genes_HK$gene_type)))
  
  # Load stably expressed genes
  # meta_genes_scSEG <- data.frame(ens_gene=rowData(sce_raw)$ID, gene_symbol=rowData(sce_raw)$Symbol, row.names = rowData(sce_raw)$ID, stringsAsFactors = F)
  # data("segList_ensemblGeneID", package = "scMerge") 
  # ref_scSEG <- segList_ensemblGeneID$human$human_scSEG
  # meta_genes_scSEG$gene_type <- '2_other'
  # meta_genes_scSEG[intersect(meta_genes_scSEG$ens_gene, ref_scSEG), 'gene_type'] <- '1_ref_scSEG_gene'
  # print(summary(as.factor(meta_genes_scSEG$gene_type)))
  # 
  # scSEG_df <- read.csv(paste0(save_dir,'segIndx_filtered_80.csv'),check.names = F, stringsAsFactors = F)
  # dim(scSEG_df)
  # scSEG_df <- scSEG_df[order(scSEG_df$segIdx,decreasing = T),]
  # meta_genes_scSEG_our <- data.frame(ens_gene=rowData(sce_raw)$ID, gene_symbol=rowData(sce_raw)$Symbol, row.names = rowData(sce_raw)$ID, stringsAsFactors = F)
  # meta_genes_scSEG_our$gene_type <- '2_other'
  # meta_genes_scSEG_our[intersect(meta_genes_scSEG_our$ens_gene, scSEG_df$gene_ens), 'gene_type'] <- '1_ours_scSEG_gene'
  # print(summary(as.factor(meta_genes_scSEG_our$gene_type)))
  
  
  # class(segList_ensemblGeneID$human$human_scSEG)
  # colnames(scSEG_df)
  # View(head(scSEG_df))
  # summary(scSEG_df$segIdx)
  # scSEG_df <- scSEG_df %>%
  #   dplyr::filter(segIdx >= quantile(segIdx, probs = 0.75))
  # print(dim(scSEG_df))
  # return(list(meta_genes_HK=meta_genes_HK, meta_genes_scSEG=meta_genes_scSEG,
  #             meta_genes_scSEG_our=meta_genes_scSEG_our))
  return(list(meta_genes_HK=meta_genes_HK))
  
}

# meta_genes
plot_pca <- function(sce_raw, meta_cells, save_dir, feature_use='treatmentSt',
                     use_raw=F, assay_use='logcounts', plttitle='', npcs=20){
  
  # genes_use <- intersect(meta_genes$ens_gene,rownames(sce_raw))
  # meta_genes <- meta_genes %>%
  #   dplyr::filter(gene_type!='2_other')
  # plttitle <- paste0(nrow(meta_genes),' ',plttitle)
  # 
  # sce_raw <- sce_raw[genes_use,]
  print(dim(sce_raw))
  if(use_raw){
    logcounts(sce_raw) <- log2(counts(sce_raw)+1)
  }
  normalized_data <- as.data.frame(assay(sce_raw, assay_use)) 
  normalized_data_t <- t(normalized_data)
  print("Computing pca")
  # pca_mat <- stats::prcomp(normalized_data_t, rank = npcs, retx=TRUE, center = TRUE, scale. = FALSE)
  pca_mat <- gmodels::fast.prcomp(normalized_data_t, retx=TRUE, center = TRUE, scale. = FALSE)
  pca_mat <- as.data.frame(pca_mat$x)
  pca_mat$cell_id <- rownames(pca_mat)
  pca_mat <- pca_mat %>% left_join(meta_cells, by=c("cell_id"))
  write.csv(pca_mat, paste0(save_dir,'gene_attr_',gsub(' ','_',plttitle),'.csv'), quote=F, row.names = F)
  
  p1 <- ggplot(pca_mat, aes_string(x = 'PC1', y = 'PC2', color = feature_use)) + 
    geom_point(alpha = 0.6, size=2) + theme_bw() + ggtitle(plttitle) + 
    theme(legend.title = element_text(size=17), 
          legend.key.size = unit(1.1, "cm"),
          legend.key.width = unit(0.5,"cm"), 
          legend.text = element_text(size=14), 
          plot.title = element_text(color="black", size=20, hjust = 0.5))
  png(paste0(save_dir,gsub(' ','_',plttitle),"_PCA.png"), height = 2*350, width=2*650,res = 2*72)
  print(p1)
  dev.off()
  return(p1)
}  


plot_tsne <- function(sce, exprs, metadata, output_dir, npcs=20, tsne_perplex=90){
  normalized_data <- as.data.frame(assay(sce, exprs)) 
  normalized_data_t <- t(normalized_data)
  print("Computing pca")
  pca_mat <- stats::prcomp(normalized_data_t, rank = npcs, retx=TRUE, center = TRUE, scale. = FALSE)
  pca_mat <- as.data.frame(pca_mat$x)
  pca_mat$cell_id <- rownames(pca_mat)
  pca_mat <- pca_mat %>% left_join(metadata, by=c("cell_id"))
  write.csv(pca_mat, file=paste0(output_dir, "pca.csv"), quote=F, row.names = F)
  print("Computing tsne")
  out_tsne <- Rtsne::Rtsne(pca_mat, perplexity = tsne_perplex, pca = TRUE, check_duplicates=FALSE)
  tsne_df <- as.data.frame(out_tsne$Y)
  
  ##########################################################
  #tSNE plot
  
  rownames(tsne_df) <- rownames(pca_mat)
  colnames(tsne_df) <- c('tSNE_1', 'tSNE_2')
  tsne_df$cell_id <- rownames(tsne_df)
  tsne_df <- tsne_df %>% left_join(metadata, by=c("cell_id"))
  write.csv(tsne_df, file=paste0(output_dir, "tsne.csv"), quote=F, row.names = F)
  
  p1 <- ggplot(tsne_df, aes_string(x = 'tSNE_1', y = 'tSNE_2', color = 'treatmentSt')) + 
    geom_point(alpha = 0.6) + theme_bw() + ggtitle('scMerge') + 
    theme(legend.title = element_text(size=17), 
          legend.key.size = unit(1.1, "cm"),
          legend.key.width = unit(0.5,"cm"), 
          legend.text = element_text(size=14), 
          plot.title = element_text(color="black", size=20, hjust = 0.5))
  
  p2 <- ggplot(tsne_df, aes_string(x = 'tSNE_1', y = 'tSNE_2', color = 'clone')) + 
    geom_point(alpha = 0.6) + theme_bw() + ggtitle('scMerge') + 
    theme(legend.title = element_text(size=17), 
          legend.key.size = unit(1.1, "cm"),
          legend.key.width = unit(0.5,"cm"), 
          legend.text = element_text(size=14), 
          plot.title = element_text(color="black", size=20, hjust = 0.5))
  
  p <- cowplot::plot_grid(p1, p2, nrow=1)
  png(paste0(output_dir,"tsne_wholedata.png"),width = 2*1000, height = 800, res = 2*72)
  print(p)
  dev.off()
  
}
plot_umap <- function(sce, exprs, metadata, output_dir, datatag, dims = 1:15){
  # assay(sce, exprs)
  # exprs <- 'logcounts'
  srt <- Seurat::as.Seurat(sce, counts = "counts", data = exprs)  
  # srt <- Seurat::FindVariableFeatures(object = srt, verbose = T)
  # srt <- Seurat::ScaleData(object = srt, verbose = FALSE) # problematics
  srt[["RNA"]]@scale.data <- as.matrix(logcounts(sce))
  # class(srt[["RNA"]]@data)
  srt <- Seurat::FindVariableFeatures(object = srt, verbose = FALSE, selection.method = "vst", nfeatures = 3000)
  var_genes <- Seurat::VariableFeatures(object = srt)
  print(paste0('Nb var genes: ',length(var_genes)))
  var_genes_df <- data.frame(var_gene=var_genes)
  # var_genes_df$var_gene[1]
  write.csv(var_genes_df, file=paste0(output_dir,'var_genes.csv'), quote=F, row.names = F)
  srt <- Seurat::RunPCA(object = srt, verbose = FALSE, npcs = 20) #, features = var_genes
  srt <- Seurat::FindNeighbors(srt, dims = dims)
  print("Run UMAP")
  # # srt <- RunTSNE(object = srt, verbose = FALSE)
  srt <- Seurat::RunUMAP(object = srt, dims = dims, verbose = FALSE) #umap.method = "umap-learn", metric='correlation'
  saveRDS(srt, paste0(output_dir,'srt.rds'))
  # # p11 <- DimPlot(srt, reduction = "tsne")
  # # p11 <- format_theme_plot(p=p11, title='Seurat Clustering')
  # # p12 <- DimPlot(srt, reduction = "tsne", group.by = 'treatmentSt')
  # # p12 <- format_theme_plot(p=p12, title='treatmentSt')
  # 
  # # p13 <- DimPlot(srt, reduction = "tsne", group.by = 'Site_origin')
  # # p13 <- p13 + labs(title='Site Origin')
  # # p14 <- DimPlot(srt, reduction = "tsne", group.by = 'batch_info')
  # # p14 <- format_theme_plot(p=p14, title='Batch Info')
  # # p15 <- DimPlot(srt, reduction = "tsne", group.by = 'passage')
  # # p15 <- format_theme_plot(p=p15, title='passage')
  # 
  
  p21 <- Seurat::DimPlot(srt, reduction = "umap", group.by = 'treat')
  # p21 <- Seurat::DimPlot(srt, reduction = "umap", group.by = 'treatmentSt')
  # p21 <- format_theme_plot(p=p21, title='Treatment Status')
  # 
  # p22 <- Seurat::DimPlot(srt, reduction = "umap", group.by = 'clone')
  # p22 <- format_theme_plot(p=p22, title='Clone Id')
  # 
  # 
  # # p23 <- DimPlot(srt, reduction = "umap", group.by = 'Site_origin')
  # # p23 <- p23 + labs(title='Site Origin')
  # # p24 <- DimPlot(srt, reduction = "umap", group.by = 'batch_info')
  # # p24 <- format_theme_plot(p=p24, title='Batch Info')
  # # p25 <- DimPlot(srt, reduction = "umap", group.by = 'passage')
  # # p25 <- format_theme_plot(p=p25, title='passage')
  # 
  # # ls_tsne <- list(p11, p12, p14, p15)
  # # ls_umap <- list(p21, p22)
  # # ptsne <- CombinePlots(plots = ls_tsne, ncol=2)
  # pumap <- cowplot::plot_grid(p21, p22, ncol=2)
  # # saveRDS(ls_tsne, file = paste0(save_dir, "plots_tsne.rds"))
  # # saveRDS(ls_umap, file = paste0(save_dir, "plots_umap.rds"))
  # # # pdf(paste0(save_dir,"clusters_treatment_tsne.pdf"), height = 4, width=10)
  # # png(paste0(save_dir,base_name,"_clusters_tnse.png"), height = 2*500, width=2*620,res = 2*72)
  # # print(ptsne)
  # # dev.off()
  # png(paste0(output_dir,datatag,"_umap_U_T.png"), height = 2*300, width=2*700,res = 2*72)
  # print(pumap)
  # dev.off()
  
  pca_df <- as.data.frame(Seurat::Embeddings(object = srt, reduction = "pca"))
  umap_df <- as.data.frame(Seurat::Embeddings(object = srt, reduction = "umap"))
  class(umap_df)
  class(pca_df)
  head(pca_df)
  dim(umap_df)
  dim(pca_df)
  # head(umap_df)
  # p <- Seurat::DimPlot(srt, reduction = "pca", group.by = 'treat')
  # p
  # srt@meta.data$treatmentSt <- srt@meta.data$treat
  library(dplyr)
  meta_info <- srt@meta.data
  meta_info$cell_id <- rownames(meta_info)
  meta_info <- meta_info %>%
    dplyr::select(cell_id, Barcode, library_id, sample, treatmentSt, timepoint)
  
  umap_df$cell_id <- rownames(umap_df)
  umap_df <- umap_df %>% inner_join(meta_info, by=c("cell_id"))
  umap_df$cell_id <- paste0(umap_df$library_id,'_',umap_df$Barcode)
  
  pca_df$cell_id <- rownames(pca_df)
  pca_df <- pca_df %>% inner_join(meta_info, by=c("cell_id"))
  pca_df$cell_id <- paste0(pca_df$library_id,'_',pca_df$Barcode)
  dim(umap_df)
  dim(pca_df)
  
  data.table::fwrite(pca_df, paste0(save_dir, datatag, "_norm_pca.csv"))
  data.table::fwrite(umap_df, paste0(save_dir, datatag, "_norm_umap.csv"))
  
}
get_treatment_status <- function(status) {
  labels <- sapply(strsplit(status, "-"), function(x) {
    return(x[2])
  })
  return(as.character(labels))
}

get_batch_infos <- function(status) {
  labels <- sapply(strsplit(status, "_"), function(x) {
    return(x[3])
  })
  return(as.character(labels))
}
# col_use <- 'treatmentSt'
# col_use <- 'clone_id'
# cell_clones should contain cell_id 
downsampling_data <- function(cell_clones, col_use, downsample_ratio=0.4, thres_small_clone=600){
  print("Sampling cells, for each clone")
  col_vals <- cell_clones[,col_use]
  clones <- unique(col_vals)
  thres_minority <- 20
  ext_cells <- c()
  for(c in clones){
    
    tmp <- cell_clones[col_vals==c,]
    if(!is.null(tmp) & nrow(tmp)<thres_small_clone & nrow(tmp)>thres_minority){   #if clone contain less than 15 cells, remove this clone from consideration
      cells_to_sample <- tmp$cell_id
    } else if(!is.null(tmp) & nrow(tmp) >= thres_small_clone){
      cells_to_sample <- sample(tmp$cell_id, floor(nrow(tmp)*downsample_ratio),replace = FALSE)
    } else{
      cells_to_sample <- NULL
      
      if(!is.null(tmp)){
        print(paste0('DEBUG: double check this step, clone: ',c, ' nb cells: ',nrow(tmp)))
      } else{
        print(paste0('DEBUG: double check this step, clone: ',c))
      }
    }
    print(paste0('Extract ',length(cells_to_sample),' from clone: ', c))
    
    if(length(cells_to_sample)>0){
      ext_cells <- c(ext_cells, cells_to_sample)
    }
    
  }
  print(paste0('Number of extracted cells: ', length(ext_cells)))
  return(ext_cells)
}

plot_rawdata_qc <- function(sce, output_dir, datatag){
  # options("scipen"=-100, "digits"=4)
  # library(scales)
  meta_cells <- as.data.frame(colData(sce))
  # meta_cells$total_features_by_counts
  # meta_cells$total_counts
  
  meta_cells <- meta_cells %>%
    dplyr::group_by(series) %>%
    dplyr::mutate(nb_cells = n())
  
  meta_cells$sample_feature <- paste0(meta_cells$series,'_',meta_cells$nb_cells,' cells')
  if(!is.null(meta_cells$batch)){
    meta_cells$series <- paste0(meta_cells$series,'_',gsub("CHIP0","C",meta_cells$batch))
    print(meta_cells$series[1])
  }
  meta_cells$series <- factor(meta_cells$series, levels = gtools::mixedsort(unique(meta_cells$series)))
  meta_cells$sample_feature <- factor(meta_cells$sample_feature, levels = gtools::mixedsort(unique(meta_cells$sample_feature)))
  p <- plot_variation_function_with_legend(meta_cells, xstring="series", ystring="total_features_by_counts",
                               plottype="sample_feature", paste0(datatag, ": raw data"),
                               '', '# genes with non-zero counts',subtbl='Nb genes with non-zero counts')
  
  p_c <- plot_variation_function(meta_cells, xstring="series", ystring="total_counts",
                                 plottype="sample_feature", paste0(datatag, ": raw data - total counts"),
                                 '', 'total_counts')
  write.csv(meta_cells, paste0(output_dir,datatag,"_raw_data_features_v2.csv"), row.names = F, quote = F)
  # p_c <- p_c + theme(legend.position = "none")
  p <- cowplot::plot_grid(p, p_c, ncol=2, align = 'h', rel_widths = c(4,3))
  datatag <- gsub(' ','_',datatag)
  saveRDS(p, paste0(output_dir,datatag,"_raw_data_features_v2.rds"))
  png(paste0(output_dir,datatag,"_raw_data_features.png"), height = 2*420, width=2*860,res = 2*72)
  print(p)
  dev.off()
  return(p)  
}

plot_rawdata_qc_metastasis <- function(sce, output_dir, datatag){
  # options("scipen"=-100, "digits"=4)
  # library(scales)
  meta_cells <- as.data.frame(colData(sce))
  # meta_cells$total_features_by_counts
  # meta_cells$total_counts
  
  meta_cells <- meta_cells %>%
    dplyr::group_by(library_id) %>%
    dplyr::mutate(nb_cells = n())
  
  meta_cells$library_id <- gsub('SCRNA10X_SA_','',meta_cells$library_id)
  meta_cells$sample_feature <- paste0(meta_cells$library_id,'_',meta_cells$nb_cells,' cells')
  # if(!is.null(meta_cells$batch)){
  #   meta_cells$series <- paste0(meta_cells$series,'_',gsub("CHIP0","C",meta_cells$batch))
  #   print(meta_cells$series[1])
  # }
  # meta_cells$series <- factor(meta_cells$series, levels = gtools::mixedsort(unique(meta_cells$series)))
  
  meta_cells$sample_feature <- factor(meta_cells$sample_feature, levels = gtools::mixedsort(unique(meta_cells$sample_feature)))
  p <- plot_variation_function_with_legend(meta_cells, xstring="library_id", ystring="total_features_by_counts",
                                           plottype="sample_feature", paste0(datatag, ": raw data"),
                                           '', '# genes with non-zero counts',subtbl='Nb genes with non-zero counts')
  
  p_c <- plot_variation_function(meta_cells, xstring="library_id", ystring="total_counts",
                                 plottype="sample_feature", paste0(datatag, ": raw data - total counts"),
                                 '', 'total_counts')
  write.csv(meta_cells, paste0(output_dir,datatag,"_raw_data_features_v2.csv"), row.names = F, quote = F)
  # p_c <- p_c + theme(legend.position = "none")
  p <- cowplot::plot_grid(p, p_c, ncol=2, align = 'h', rel_widths = c(4,3))
  datatag <- gsub(' ','_',datatag)
  saveRDS(p, paste0(output_dir,datatag,"_raw_data_features_v2.rds"))
  png(paste0(output_dir,datatag,"_raw_data_features.png"), height = 2*420, width=2*860,res = 2*72)
  print(p)
  dev.off()
  return(p)  
}

get_unique_clone_id <- function(clone_labels){
  set.seed(42) 
  cls <- sapply(strsplit(clone_labels,'_'), function(x){
    if(length(x)==1){
      return(x[1])
    }else if(length(x)==2){
      idx <- sample(c(1,2),1)
      return(x[idx])
    }else{
      idx <- sample(c(1,2,3),1)
      return(x[idx])
    }
    
  })
  return(as.character(cls))
}
get_lib_info <- function(sce){
  # sce$library_id <- 'None'
  sce$library_id <- gsub('(.cache/)','',sce$Sample)
  sce$library_id <- gsub('(/filtered_feature_bc_matrix)','',sce$library_id)
  print(unique(sce$library_id))
  return(sce$library_id)
}

# DE analysis: clone in sce1 vs clone in sce2
# desc <- '(Rx - UnRx)'

fix_batch <- function(sce){
  batch_ls <- as.character(sce$batch)
  batch_ls <- ifelse(batch_ls=='CHIP001','CHIP0047',
                     ifelse(batch_ls=='CHIP003' | batch_ls=='CHIP004','CHIP0192',batch_ls))
  return(batch_ls)
}


compute_umap <- function(srt){
  srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = 2000)
  srt <- RunPCA(srt, verbose = F)
  srt <- RunUMAP(srt, dims = 1:30, verbose = FALSE)
  srt <- FindNeighbors(srt, dims = 1:30, verbose = FALSE)
  srt <- FindClusters(srt, verbose = FALSE, resolution = 0.5)
  
  Idents(object = srt) <- "seurat_clusters"
  dim(srt)
  srt1 <- subset(x = srt, idents = c(0,1,2,3))
  dim(srt1)
  srt <- srt1
  p_umap <- DimPlot(srt, label = TRUE, group.by = 'treatmentSt') #+ NoLegend()
  
  p_umap2 <- DimPlot(srt, label = TRUE) #+ NoLegend()
  p_umap3 <- DimPlot(srt, label = TRUE, group.by = 'clone')
  
  pm <- p_umap + p_umap3 + p_umap2 + patchwork::plot_layout(ncol=3)
  png(paste0(save_dir, datatag, "_clusters_umap_sctransform_treated_cells_excluded_cluster4.png"), height = 2*300, width=2*1000,res = 2*72)
  print(pm)
  dev.off()
  meta_info <- srt@meta.data
  class(meta_info)
  summary(as.factor(meta_info$seurat_clusters))
  sum(meta_info$cell_id==pca_df$cell_id)
  dim(srt@reductions$umap)
  # umap_df <- srt[["umap"]]
  pca_df <- as.data.frame(Embeddings(object = srt, reduction = "pca"))
  umap_df <- as.data.frame(Embeddings(object = srt, reduction = "umap"))
  class(umap_df)
  class(pca_df)
  head(pca_df)
  dim(umap_df)
  head(umap_df)
  p <- DimPlot(srt, reduction = "pca")
  p
  
  meta_info <- meta_info %>%
    dplyr::select(cell_id, Barcode, library_id, sample, clone, treatmentSt, timepoint)
  
  umap_df$cell_id <- rownames(umap_df)
  umap_df <- umap_df %>% inner_join(meta_info, by=c("cell_id"))
  umap_df$cell_id <- paste0(umap_df$library_id,'_',umap_df$Barcode)
  
  pca_df$cell_id <- rownames(pca_df)
  pca_df <- pca_df %>% inner_join(meta_info, by=c("cell_id"))
  pca_df$cell_id <- paste0(pca_df$library_id,'_',pca_df$Barcode)
  dim(umap_df)
  dim(pca_df)
  # test_sce <- readRDS(paste0(input_dir,'rnaseq_v6/SA535-v6/SA535X10XB03696.rdata'))
  # t <- colnames(colData(sce))
  # write.table(t, file = paste0(save_dir,'colnames_sce.txt'), sep='\n')
  save_dir <- '/home/htran/storage/python_workspace/velocity/result_SA535/cx5461/'
  
  write.csv(pca_df, paste0(save_dir, datatag, "_norm_pca.csv"), row.names = F, quote = F)
  write.csv(umap_df, paste0(save_dir, datatag, "_norm_umap.csv"), row.names = F, quote = F)
  
  
}

scTransform_v3 <- function(){
  # BiocManager::install("glmGamPoi")
}
load_data <- function(meta_cells, base_dir, tag='CIS_Rx_RxH'){
  sce_list <- list()
  c <- 0
  for (f in unique(meta_cells$mouse_id)){
    # norm_fn <- paste0(base_dir,f,".rdata")
    norm_fn <- paste0(base_dir,f,".rds")
    if(file.exists(norm_fn)){
      sce_normalized <- readRDS(norm_fn)
      c <- c + 1
      if(c==1){
        cols_use <- colnames(colData(sce_normalized))
      }else{
        cols_use <- intersect(cols_use, colnames(colData(sce_normalized)))
      }
    } else{
      warning('Error loading files')
      print(f)
    }  
    print(paste0('DEBUG: ',f,'  ',dim(sce_normalized)[1],' ',dim(sce_normalized)[2]))
    sce_list[[f]] <- sce_normalized
    
  }  
  length(cols_use)
  
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


load_data <- function(samples, cells_use, input_dir, tag='CIS_Rx_RxH'){
  sce_list <- list()
  c <- 0
  for (f in unique(samples)){
    norm_fn <- paste0(input_dir,f,".rdata")
    if(file.exists(norm_fn)){
      sce_normalized <- readRDS(norm_fn)
      print(dim(sce_normalized))
      sce_normalized <- sce_normalized[,sce_normalized$Barcode %in% cells_use]
      print(dim(sce_normalized))
      c <- c + 1
      if(c==1){
        cols_use <- colnames(colData(sce_normalized))
      }else{
        cols_use <- intersect(cols_use, colnames(colData(sce_normalized)))
      }
    } else{
      warning('Error loading files')
      print(f)
    }  
    print(paste0('DEBUG: ',f,'  ',dim(sce_normalized)[1],' ',dim(sce_normalized)[2]))
    sce_list[[f]] <- sce_normalized
    
  }  
  length(cols_use)
  
  sce_combine <- sce_cbind_func_v2(sce_list, cut_off_overall = 0, exprs = c("counts"), 
                                   colData_names = cols_use, save_raw=F, save_dir=input_dir, tag) 
  print(dim(sce_combine))
  # sce_combine <- readRDS(paste0(output_dir,'combined_total_genes_filtered.rds'))
  # Twice scater normalization
  # output are saved in logcounts exp values
  print(("Normalizing total data..."))
  
  ## cell name = library_id + barcode
  sce_combine$library_id <- gsub('.cache/','',sce_combine$Sample)
  sce_combine$library_id <- gsub('/filtered_feature_bc_matrix','',sce_combine$library_id)
  print(unique(sce_combine$library_id))
  # sum(unique(sce_combine$library_id) %in% meta_data$library_id)
  colnames(sce_combine) <- paste0(sce_combine$sample,'_',sce_combine$Barcode)
  # sce_combine <- sce_combine[,colnames(sce_combine) %in% meta_data$cell_id]
  # rownames(sce_combine)
  # print(rowData(sce_combine_raw)$Symbol[1:3])
  
  saveRDS(sce_combine, file=paste0(input_dir,"combined_",tag,".rds"))
  
  return(sce_combine)
}


plot_batch_estimation <- function(sce_raw, sce_scran, sce_seurat, sce_transform, stable_genes, save_dir, 
                                  datatag='', genes_type='', save_data=T, yl=NULL){
  p_raw <- plot_stable_genes_exp_v2(sce_raw, stable_genes, use_raw=T, exprs='logcounts', 
                                    plottitle=paste0(datatag,': Raw Data '),
                                    subtitle = paste0(length(stable_genes),' ',genes_type,' genes avg exp'),
                                    xlabel='', ylabel="Avg. genes exp", yl=yl, legend_visible=T)
  p_scran <- plot_stable_genes_exp_v2(sce_scran, stable_genes, use_raw=F, exprs='logcounts', 
                                      plottitle=paste0(datatag,': Scran Normalization '),
                                      subtitle = paste0(length(stable_genes),' ',genes_type,' genes avg exp'),
                                      xlabel='', ylabel="Avg. genes exp", yl=yl, legend_visible=F)
  
  p_seurat <- plot_stable_genes_exp_v2(sce_seurat, stable_genes, use_raw=F, exprs='logcounts', 
                                       plottitle=paste0(datatag,': Seurat Normalization '),
                                       subtitle = paste0(length(stable_genes),' ',genes_type,' genes avg exp'), 
                                       xlabel='', ylabel="Avg. genes exp", yl=yl, legend_visible=F)
  
  p_sctransform <- plot_stable_genes_exp_v2(sce_transform, stable_genes, use_raw=F, exprs='logcounts', 
                                            plottitle=paste0(datatag,': SCTransform Normalization '),
                                            subtitle = paste0(length(stable_genes),' ',genes_type,' genes avg exp'), 
                                            xlabel='', ylabel="Avg. genes exp", yl=yl, legend_visible=F)
  
  
  
  # p_total <- cowplot::plot_grid(p_raw$p, p_scran$p, p_seurat$p, p_sctransform$p, nrow = 1, align='h', rel_widths = c(1.35,1,1,1))
  # png(paste0(save_dir,datatag,'_',genes_type,"_eval.png"), height = 2*500, width=2*1200,res = 2*72)
  # print(p_total)
  # dev.off()
  
  # if(save_data){
  #   # saveRDS(p_total, paste0(save_dir,gsub(' ','_',datatag),'_',genes_type,"_HK_eval_v2.rds"))
  #   meta_ls <- p_raw$meta_data + p_scran$meta_data + p_seurat$meta_data + p_sctransform$meta_data
  #   saveRDS(meta_ls, paste0(save_dir,gsub(' ','_',datatag),'_',genes_type,"_HK_eval_metadata_v2.rds"))
  # }
  col_use <- c("tag","sample","batch","cell_id","clone","library_id","treatmentSt","mean_exp","timepoint","series")
  p_raw$tag <- 'Raw Data'
  p_scran$tag <- 'Scran'
  p_seurat$tag <- 'Seurat'
  p_sctransform$tag <- 'SCTransform'
  meta_ls <- dplyr::bind_rows(p_raw, p_scran, p_seurat, p_sctransform)
  meta_ls <- meta_ls[,colnames(meta_ls) %in% col_use]
  print(dim(meta_ls))
  saveRDS(meta_ls, paste0(save_dir,gsub(' ','_',datatag),'_HK_eval_metadata_v2.rds'))
  return(meta_ls)
}

#include scMerge
plot_batch_estimation_v2 <- function(sce_scMerge, sce_raw, sce_scran, sce_seurat, sce_transform, stable_genes, save_dir, datatag='', genes_type='', save_data=F){
  p_raw <- plot_stable_genes_exp_v2(sce_raw, stable_genes, use_raw=T, exprs='logcounts', plottitle=paste0(datatag,': Raw - ',length(stable_genes),' ',genes_type,' genes mean exp'),
                                    xlabel='', ylabel="Mean genes exp", yl=c(0.2,1.7), legend_visible=T)
  
  p_scran <- plot_stable_genes_exp_v2(sce_scran, stable_genes, use_raw=F, exprs='logcounts', plottitle=paste0(datatag,': Twice Scran ',length(stable_genes),' ',genes_type,' genes mean exp'),
                                      xlabel='', ylabel="Mean genes exp", yl=NULL, legend_visible=F)
  
  p_seurat <- plot_stable_genes_exp_v2(sce_seurat, stable_genes, use_raw=F, exprs='logcounts', plottitle=paste0(datatag,': Seurat ',length(stable_genes),' ',genes_type,' genes mean exp'),
                                       xlabel='', ylabel="Mean genes exp", yl=NULL, legend_visible=F)
  
  p_sctransform <- plot_stable_genes_exp_v2(sce_transform, stable_genes, use_raw=F, exprs='logcounts', plottitle=paste0(datatag,': SCTransform ',length(stable_genes),' ',genes_type,' genes mean exp'),
                                            xlabel='', ylabel="Mean genes exp", yl=NULL, legend_visible=F)
  
  p_scmerge <- plot_stable_genes_exp_v2(sce_scMerge, stable_genes, use_raw=F, exprs='scMerge_fast', plottitle=paste0(datatag,': scMerge ',length(intersect(stable_genes, rownames(sce_scMerge))),' ',genes_type,' genes mean exp'),
                                        xlabel='', ylabel="Mean genes exp", yl=NULL, legend_visible=F)
  
  p_total <- cowplot::plot_grid(p_raw, p_scran, p_seurat, p_sctransform, p_scmerge, ncol = 1, align='v')
  png(paste0(save_dir,datatag,'_',genes_type,"_eval.png"), height = 2*1000, width=2*400,res = 2*72)
  print(p_total)
  dev.off()
  
  if(save_data){
    saveRDS(p_total, paste0(save_dir,datatag,'_',genes_type,"_HK_eval.rds"))
  }
  return(p_total)
}


#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
# sd = sd(x[[col]], na.rm=TRUE)

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      selower = mean(x[[col]]) - sd(x[[col]])/sqrt(length(x[[col]])),
      seupper = mean(x[[col]]) + sd(x[[col]])/sqrt(length(x[[col]])))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


eval_edgeR_v2 <- function(sce, markers_ls, treatmentSts, clones, de_desc, desc, datatag, output_dir){
  
  sce <- sce[,sce$treatmentSt %in% treatmentSts & sce$clone %in% clones]
  sce <- sce[intersect(rownames(sce), markers_ls$ensembl_gene_id),]
  print(dim(sce))
  cond1 <- sce$treatmentSt==treatmentSts[1] & sce$clone==clones[1]
  cond2 <- sce$treatmentSt==treatmentSts[2] & sce$clone==clones[2]
  sce1 <- sce[,cond1]
  sce2 <- sce[,cond2]
  print(dim(sce1))
  print(dim(sce2))
  norm_data1 <- as.data.frame(as.matrix(logcounts(sce1)))
  norm_data1[norm_data1==0] <- NA
  norm_data1$gene <- rownames(sce1)
  
  exp_mean1 <- norm_data1 %>% 
    # convert to long format
    pivot_longer(!gene, names_to = "cell_id", values_to = "exp")  %>% 
    group_by(gene) %>% 
    summarise(mean_exp = mean(exp, na.rm=T))
  
  norm_data2 <- as.data.frame(as.matrix(logcounts(sce2)))
  norm_data2[norm_data2==0] <- NA
  norm_data2$gene <- rownames(sce2)
  
  exp_mean2 <- norm_data2 %>% 
    # convert to long format
    pivot_longer(!gene, names_to = "cell_id", values_to = "exp")  %>% 
    group_by(gene) %>% 
    summarise(mean_exp = mean(exp, na.rm=T))
  
  dim(exp_mean1)
  dim(exp_mean2)
  sum(exp_mean1$gene==exp_mean2$gene)
  mean_df <- data.frame(gene=exp_mean1$gene, exp_mean=exp_mean1$mean_exp - exp_mean2$mean_exp, stringsAsFactors = F)
  dim(mean_df)
  
  # 
  rownames(markers_ls) <- markers_ls$ensembl_gene_id
  markers_ls <- markers_ls[mean_df$gene,]
  # if(sum(mean_df$gene==markers_ls$ensembl_gene_id)!=nrow(markers_ls)){
  #   stop('Double check analysis')
  # }
  cr <- cor(mean_df$exp_mean, markers_ls$logFC, method = "spearman")
  pltttitle <- paste0(datatag,': edgeR - SCTransform \n Correlation: ',round(cr,3))
  markers_ls$DE_gene <- ifelse(markers_ls$logFC>0,'Up-regulated','Down-regulated')
  # stat <- data.frame(logFC_edgeR=markers_ls$logFC,SCTransform_norm_exp=mean_df$exp_mean, 
  #                    gene_type=markers_ls$DE_gene,
  #                    ensembl_gene_id = markers_ls$ensembl_gene_id,
  #                    gene_symbol = markers_ls$gene_symbol,
  #                    stringsAsFactors = F)
  stat <- markers_ls
  # stat <- stat %>%
  #   dplyr::rename(logFC_edgeR=logFC)
  stat <- stat %>% inner_join(mean_df, by=c('ensembl_gene_id'='gene'))
  # head(stat)
  de_desc <- gsub(' ','_',de_desc)
  write.csv(stat, paste0(output_dir,"edgeR_corr_eval_",de_desc,".csv"), quote = F, row.names = F)
  p <- ggplot(stat, aes(logFC, exp_mean, color=DE_gene, shape=DE_gene)) + 
    geom_point(size=2.3, alpha=1) +
    scale_shape_manual(values=c('Up-regulated'=1, 'Down-regulated'=2)) +
    scale_color_manual(values=c('Up-regulated'='#E69F00', 'Down-regulated'='#56B4E9'))
  p <- p + labs(x='edgeR log2 FC',y=paste0('SCTransform ',desc), title = pltttitle) + 
    theme(plot.title = element_text(color="black", size=13, hjust = 0.5, face = "bold"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text.y = element_text(color="black", size=12),
          axis.text.x = element_text(color="black", size=12),
          axis.title = element_text(color="black", size=12),
          legend.title = element_text(color="black", size=12), 
          legend.text = element_text(color="black", size=12))
  # p
  saveRDS(p, paste0(output_dir,"edgeR_corr_eval_",de_desc,".rds"))
  png(paste0(output_dir,"edgeR_corr_eval_",de_desc,".png"), height = 2*350, width=2*500,res = 2*72)
  print(p)
  dev.off()
  return(p)
}


eval_edgeR_counts <- function(sce, markers_ls, treatmentSts, clones, de_desc, desc, datatag, output_dir){
  
  sce <- sce[,sce$treatmentSt %in% treatmentSts & sce$clone %in% clones]
  sce <- sce[intersect(rownames(sce), markers_ls$ensembl_gene_id),]
  print(dim(sce))
  cond1 <- sce$treatmentSt==treatmentSts[1] & sce$clone==clones[1]
  cond2 <- sce$treatmentSt==treatmentSts[2] & sce$clone==clones[2]
  sce1 <- sce[,cond1]
  sce2 <- sce[,cond2]
  
  logcounts(sce1) <- log2(counts(sce1)+1)
  logcounts(sce2) <- log2(counts(sce2)+1)
  print(dim(sce1))
  print(dim(sce2))
  norm_data1 <- as.data.frame(as.matrix(logcounts(sce1)))
  norm_data1[norm_data1==0] <- NA
  norm_data1$gene <- rownames(sce1)
  
  exp_mean1 <- norm_data1 %>% 
    # convert to long format
    pivot_longer(!gene, names_to = "cell_id", values_to = "exp")  %>% 
    group_by(gene) %>% 
    summarise(mean_exp = mean(exp, na.rm=T))
  
  norm_data2 <- as.data.frame(as.matrix(logcounts(sce2)))
  norm_data2[norm_data2==0] <- NA
  norm_data2$gene <- rownames(sce2)
  
  exp_mean2 <- norm_data2 %>% 
    # convert to long format
    pivot_longer(!gene, names_to = "cell_id", values_to = "exp")  %>% 
    group_by(gene) %>% 
    summarise(mean_exp = mean(exp, na.rm=T))
  
  dim(exp_mean1)
  dim(exp_mean2)
  sum(exp_mean1$gene==exp_mean2$gene)
  mean_df <- data.frame(gene=exp_mean1$gene, exp_mean=exp_mean1$mean_exp - exp_mean2$mean_exp, stringsAsFactors = F)
  dim(mean_df)
  
  # 
  rownames(markers_ls) <- markers_ls$ensembl_gene_id
  markers_ls <- markers_ls[mean_df$gene,]
  # if(sum(mean_df$gene==markers_ls$ensembl_gene_id)!=nrow(markers_ls)){
  #   stop('Double check analysis')
  # }
  cr <- cor(mean_df$exp_mean, markers_ls$logFC, method = "spearman")
  pltttitle <- paste0(datatag,': edgeR - Raw log counts \n Correlation: ',round(cr,3))
  markers_ls$DE_gene <- ifelse(markers_ls$logFC>0,'Up-regulated','Down-regulated')
  # stat <- data.frame(logFC_edgeR=markers_ls$logFC,SCTransform_norm_exp=mean_df$exp_mean, 
  #                    gene_type=markers_ls$DE_gene,
  #                    ensembl_gene_id = markers_ls$ensembl_gene_id,
  #                    gene_symbol = markers_ls$gene_symbol,
  #                    stringsAsFactors = F)
  stat <- markers_ls
  # stat <- stat %>%
  #   dplyr::rename(logFC_edgeR=logFC)
  stat <- stat %>% inner_join(mean_df, by=c('ensembl_gene_id'='gene'))
  # head(stat)
  de_desc <- gsub(' ','_',de_desc)
  write.csv(stat, paste0(output_dir,"edgeR_corr_eval_",de_desc,".csv"), quote = F, row.names = F)
  p <- ggplot(stat, aes(logFC, exp_mean, color=DE_gene, shape=DE_gene)) + 
    geom_point(size=2.3, alpha=1) +
    scale_shape_manual(values=c('Up-regulated'=1, 'Down-regulated'=2)) +
    scale_color_manual(values=c('Up-regulated'='#E69F00', 'Down-regulated'='#56B4E9'))
  p <- p + labs(x='edgeR log2 FC',y=paste0('Raw logcounts ',desc), title = pltttitle) + 
    theme(plot.title = element_text(color="black", size=13, hjust = 0.5, face = "bold"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text.y = element_text(color="black", size=12),
          axis.text.x = element_text(color="black", size=12),
          axis.title = element_text(color="black", size=11),
          legend.title = element_text(color="black", size=12), 
          legend.text = element_text(color="black", size=12))
  # p
  saveRDS(p, paste0(output_dir,"edgeR_corr_eval_",de_desc,".rds"))
  png(paste0(output_dir,"edgeR_corr_eval_",de_desc,".png"), height = 2*350, width=2*500,res = 2*72)
  print(p)
  dev.off()
  return(p)
}

eval_edgeR <- function(sce, markers_ls, de_desc, datatag, output_dir){  #vol_plt, 
  obs_genes <- markers_ls[markers_ls$logFC<0,'ensembl_gene_id']
  print(length(obs_genes))
  exprs <- 'logcounts'
  # de_desc <- paste0(obs_clones,'-',obs_treatment_st,' vs ', obs_clones_untreated,'-', obs_untreated_st)
  # desc <- 'Using scTransform log normalized mtx'
  p_down <- plot_stable_genes_exp_v3(sce, obs_genes, use_raw=F, exprs=exprs, 
                                     plottitle=paste0('edgeR - down-regulated genes \n'), de_desc,
                                     xlabel='', ylabel="Avg exp of down-regulated genes", yl=NULL, scale=F)
  
  obs_genes <- markers_ls[markers_ls$logFC>0,'ensembl_gene_id']
  print(length(obs_genes))
  p_up <- plot_stable_genes_exp_v3(sce, obs_genes, use_raw=F, exprs=exprs, 
                                   plottitle=paste0('edgeR - up-regulated genes \n'), de_desc,
                                   xlabel='', ylabel="Avg exp of up-regulated genes", yl=NULL)
  # p <- cowplot::plot_grid(vol_plt, p_up, p_down, nrow=1, rel_widths = c(1.5, 1, 1))
  p <- cowplot::plot_grid(p_up, p_down, nrow=1, rel_widths = c(1, 1))
  
  # output_dir <- paste0(input_dir,'rnaseq_v6/normalization_evaluation/SA535_cisplatin/')
  saveRDS(p, paste0(output_dir,"edgeR_DE_analysis_",gsub(' ','_',de_desc),".rds"))
  png(paste0(output_dir,"edgeR_DE_analysis_",de_desc,".png"), height = 2*400, width=2*950,res = 2*72)
  print(p)
  dev.off()
  return(p)
}
plot_stable_genes_exp_v3 <- function(sce, stable_genes, use_raw=F, exprs='logcounts', plottitle='Batch effect estimation'
                                     ,desc, 
                                     xlabel='', ylabel="Mean exp of stably expressed genes", yl=c(0,1.6), scale=F){
  print(length(stable_genes))
  stable_genes <- intersect(rownames(sce), as.character(stable_genes))
  plottitle <- paste0(plottitle, desc) #length(stable_genes), 
  sce_SEG <- sce[stable_genes,]
  print(dim(sce_SEG))
  if(use_raw){
    logcounts(sce_SEG) <- log2(counts(sce_SEG)+1)
  }
  
  # summary_SEG <- data.frame(mean_seg = colMeans(logcounts(sce_SEG)), row.names=colnames(sce_SEG))
  # # head(summary_SEG)
  # sce_SEG$mean_seg <- summary_SEG[colnames(sce_SEG),"mean_seg"]
  meta_data <- as.data.frame(colData(sce_SEG))
  corrected <- assay(sce_SEG, exprs)
  if(scale){
    corrected <- corrected %>% 
      # transpose the matrix so genes are as columns
      t() %>% 
      # apply scalling to each column of the matrix (genes)
      scale() %>% 
      # transpose back so genes are as rows again
      t()
  }
  
  
  meta_data$mean_exp <- DelayedArray::colMeans(corrected)
  # meta_data$batch_info <- as.factor(meta_data$batch_info)
  meta_data$library_id <- gsub("SCRNA10X_SA_CHIP","",meta_data$library_id)
  p <- plot_variation_function_with_legend(meta_data, xstring="series", ystring="mean_exp", 
                                           plottype="series", plottitle,
                                           xlabel, ylabel)
  if(!is.null(yl)){
    p <- p + ylim(yl[1],yl[2])
  }
  
  return(p)
  
}
