## edgeR script
## Testing whether edgeR or DESeq2 provide better outputs. 


## To Do: loading data from main exp, clone B, C, and 
## Comparing output of edgeR versus DESeq2


load_input_data <- function(){
  input_dir <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/bulkRNAseq/'
  meta_samples <- data.table::fread(paste0(input_dir, 'metadata/metadata_Hakwoo_bulkRNA_mixing_main_exp.csv'))
  dim(meta_samples)
  
  meta_samples <- meta_samples %>%
    dplyr::filter(main_clone %in% c('B','C'))
  dim(meta_samples)
  
  input_dir <- '/home/htran/storage/rnaseq_datasets/bulk_metastasis/bulk_SA919/'
  datatag <- 'SA919Fig6'
  save_dir <- paste0(input_dir, 'preprocessed/')
  
  normalized_fn <- paste0(save_dir, datatag, '_sizefactor_normalized.csv.gz')
  df_normalized <- data.table::fread(normalized_fn)
  # head(df_normalized)
  sf_df <- data.table::fread(paste0(save_dir, datatag, '_sizefactors_Scran.csv.gz'))
  
  df_counts_fn <- paste0(save_dir, datatag,'_total_raw_counts.csv.gz')
  cts <- data.table::fread(df_counts_fn) %>% as.data.frame()
  print(dim(cts))
  head(cts)
  cts <- cts %>%
    tibble::column_to_rownames('ens_gene_id')
  coldata <- meta_samples %>%
    dplyr::filter(bulk_sid %in% colnames(cts)) %>%
    dplyr::select(mainsite, main_clone, experiment, bulk_sid) %>%
    as.data.frame() %>%
    tibble::column_to_rownames('bulk_sid')
  coldata <- coldata[colnames(cts),]
    
  
}

# meta_df <- coldata
# counts_df <- cts

run_edgeR <- function(datatag, meta_df, input_dir, base_dir){
  library("edgeR")
  groups_use <- c('Rx','UnRx')
  colnames(coldata)
  meta_df$condition <- meta_df$mainsite
  groups_use <- c('Metastasis','Primary')
  #edgeR package use inverted order compared to general convention in other DE tools
  # i.e: Rx vs. UnRx
  # edgeR: 1_UnRx vs. 2_Rx
  meta_df$condition <- ifelse(grepl(groups_use[2], meta_df$condition),paste0("1_",groups_use[2]),
                              paste0("2_",groups_use[1]))
  print(meta_df)
  meta_df$condition <- factor(meta_df$condition,levels = c("1_Primary","2_Metastasis"))
  # rownames(meta_df) <- meta_df$sample
  # counts_df <- load_counts_data(meta_df$sample, input_dir, datatag)
  print(colnames(counts_df))
  # rm(dge)
  dge <- load_dge_counts(counts_df, meta_df$condition)
  # save_dir_de <- paste0(base_dir,paste(meta_df$sample,collapse = '_'),'/')
  save_dir_de <- paste0(save_dir,'edgeR_testing/')
  dir.create(save_dir_de, recursive = T)
  dge_backup <- dge
  de_genes <- edgeR_de_v2(dge, meta_df, save_dir_de, tag=NULL)
  res <- list(counts_df=counts_df,de_genes=de_genes, 
              save_dir=save_dir_de, datatag=datatag,
              meta_df=meta_df)
  return(res)
  
}
edgeR_de_v2 <- function(dge, meta_df, save_dir, save_dge=T, tag=NULL){
  
  # print("Filtering data...")
  # sce_de <- sce_de[, sce_de$total_features_by_counts > 1500]
  # # rs <- rowSums(as.matrix(counts(sce_de)))
  # # qplot(rs, log='x') + geom_vline(xintercept = 100)
  # 
  # sce_de <- sce_de[rowSums(as.matrix(counts(sce_de))) > 100, ]
  # print(dim(sce_de))
  # mycounts <- as.matrix(counts(sce_de))    # zeros are good
  # print("Create DGE edgeR object...")
  
  # dge <- DGEList(counts=mycounts, group=sce_de$clone)
  # dge <- DGEList(counts=mycounts, group=sce_de$treatmentSt)
  
  print("DE Analysis...")
  # which vs which ???
  # design <- model.matrix(~ clone, data = colData(sce_de))
  design <- model.matrix(~condition, data = meta_df)
  if(is.null(rownames(design))){
    rownames(design) <- meta_df$sample
  }
  # 
  # design <- model.matrix(~ treatmentSt, data = colData(sce_de))
  # This describes the edgeR user manual
  # http://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
  # dge <- calcNormFactors(dge, method = "TMM") ## check the update one
  dge <- calcNormFactors(dge, method = "TMM")
  sf_df <- data.table::fread(paste0(save_dir, datatag, '_sizefactors_Scran.csv.gz'))
  sf_df <- sf_df %>%
    tibble::column_to_rownames('sample')
  
  dge$samples$norm.factors <- sf_df[rownames(dge$samples),'size_factor']
  dge <- edgeR::estimateDisp(dge, design = design)
  # dge$common.dispersion
  # dge <- edgeR::estimateCommonDisp(dge, design = design)
  # dge <- edgeR::estimateTagwiseDisp(dge)
  # dge <- edgeR::estimateGLMCommonDisp(dge, design = design)
  # de.com <- exactTest(dge, dispersion=predefined_disp)
  # View(meta_df)
  fit <- edgeR::glmQLFit(dge, design = design)
  qlf <- edgeR::glmQLFTest(fit)
  tt <- edgeR::topTags(qlf, n = Inf)
  dim(tt)
  
  print(summary(decideTests(qlf)))
  png(paste0(save_dir,tag,"_DE_genes.png"), height = 2*350, width=2*500, res = 2*72)
  plotMD(qlf)
  abline(h=c(-1, 1), col="blue")
  dev.off()
  
  
  print("Generate output...")
  tt <- as.data.frame(tt) %>% 
    tibble::rownames_to_column("gene_symbol")
  sum(tt$FDR<0.05 & tt$PValue<0.05 & abs(tt$logFC)>0.25)
  # tt$gene_symbol <- rowData(sce_de[tt$gene_id,])$ID
  # tt$gene_symbol <- rowData(sce_de[tt$gene_id,])$Symbol
  # View(head(tt))
  print(dim(tt))
  # tt$ens_gene_id <- get_ens_gene(tt$gene_id)
  tt$gene_symbol <- get_gene_symbol(tt$gene_id)
  # print(length(unique(tt$gene_symbol)))
  print(summary(tt$logFC))
  print(summary(tt$PValue))
  # Saving the logFC file
  # tt <- tt[abs(tt$logFC)>1,] # tt$FDR<0.01 & tt$PValue<0.05 
  print(head(tt))
  if(is.null(tag)){
    tag <- paste(unique(meta_df$condition),collapse='_')
  }
  write.csv(tt, file=paste0(save_dir, 'edgeR_significant_genes_',tag,'.csv'), row.names=FALSE, quote=FALSE)
  # print("With threshold logFC>0.25, number of significant genes is: ")
  # print(nrow(tt[abs(tt$logFC)>0.25,]))
  if(save_dge){
    saveRDS(dge,paste0(save_dir, 'edgeR_dge',tag,'.rds'))  
  }
  
  return(tt)
}

load_dge_counts <- function(cts, groups_use){
  # Creating a DGEList object for use in edgeR.
  # class(cts)
  # dim(cts)
  # head(cts)
  dge <- edgeR::DGEList(cts, group = groups_use)
  
  # dge <- edgeR::calcNormFactors(cts)
  # dge$genes <- meta_genes
  
  # normMat <- edgeR::cpm(cts, log=TRUE) 
  # dge <- edgeR::scaleOffset(dge, normMat)
  # filtering
  print(dim(dge))
  dge$samples$lib.size <- colSums(dge$counts)
  keep <- edgeR::filterByExpr(dge)
  dge <- dge[keep, ]
  print(dim(dge))
  print(class(dge))
  
  # dge <- calcNormFactors(dge, method = "TMM")
  
  # lcpm_dge <- cpm(dge, log=TRUE)
  # dim(lcpm_dge)
  # View(head(lcpm_dge))
  # y is now ready for estimate dispersion functions see edgeR User's Guide
  return(dge)
}
