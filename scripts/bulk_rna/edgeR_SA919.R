## edgeR script
## Testing whether edgeR or DESeq2 provide better outputs. 

# Figure 6
suppressPackageStartupMessages({
  library("dplyr")
  library("ggplot2")
  library("SingleCellExperiment")
  # library("DESeq2")
})

## Loading utility functions
script_dir <- '/home/htran/Projects/hakwoo_project/metastasis_material/'
script_dir <- '/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/metastasis_material/'
source(paste0(script_dir, 'scripts/bulk_rna/bulk_utils.R'))


## To Do: loading data from main exp, clone B, C, and 
## Get FPKM norm data, and comparing mean exp of C vs. B
## Boxplot of mean exp
## Bootstrap testing

## Get dispersion values for cis genes, for all genome genes 
## Bootstrap testing


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


get_sce_object <- function(dge, meta_df, save_dir, save_dge=T, tag=NULL){
  input_dir <- "/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/metastasis_material/materials/bulkRNAseq/preprocessed_09April2024/"
  input_dir <- "/Users/hoatran/Documents/BCCRC_projects/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/preprocessed_09April2024/"
  input_data_dir <- "/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/bulk_SA919/"
  datatag <- 'SA919Fig6'
  save_dir <- input_dir
  
  meta_samples <- data.table::fread(paste0(input_dir, '../metadata/metadata_Hakwoo_bulkRNA_mixing_main_exp.csv'))
  dim(meta_samples)
  
  meta_samples <- meta_samples %>%
    dplyr::filter(main_clone %in% c('B','C'))
  dim(meta_samples)
  
  df_counts_fn <- paste0(save_dir, datatag,'_total_raw_counts.csv.gz')
  # counts_total <- data.table::fread(df_counts_fn)
  # dim(counts_total)
  
  normalized_fn <- paste0(save_dir, datatag, '_sizefactor_normalized.csv.gz')
  df_normalized <- data.table::fread(normalized_fn)
  # head(df_normalized)
  dim(df_normalized)
  sf_df <- data.table::fread(paste0(save_dir, datatag, '_sizefactors_Scran.csv.gz'))
  
  
  cts <- data.table::fread(df_counts_fn) %>% as.data.frame()
  print(dim(cts))
  # head(cts)
  colnames(meta_samples)
  
  cts <- cts %>%
    tibble::column_to_rownames('ens_gene_id')
  coldata <- meta_samples %>%
    dplyr::filter(bulk_sid %in% colnames(cts)) %>%
    dplyr::select(mainsite, main_clone, experiment, bulk_sid) %>%
    as.data.frame() %>%
    tibble::column_to_rownames('bulk_sid')
  coldata <- coldata[colnames(cts),]
  sfs <- sf_df$size_factor
  names(sfs) <- sf_df$sample
  coldata$size_factor <- sfs[rownames(coldata)]
  
  
  sids <- colnames(cts)
  for(s in sids){
    cts[,s] <- round(cts[,s],0)
  }
  # head(cts)
  # coldata <- meta_samples
  coldata$condition <- coldata$mainsite
  
  head(df_normalized)
  df_normalized <- df_normalized %>%
    tibble::column_to_rownames('ens_gene_id') %>%
    dplyr::select(all_of(colnames(cts)))
  sum(colnames(df_normalized)==colnames(cts))
  dim(df_normalized)
  df_normalized <- df_normalized[rownames(cts),]
  sids <- colnames(df_normalized)
  for(s in sids){
    df_normalized[,s] <- round(df_normalized[,s],1)
  }
  
  
  meta_genes <- tibble::tibble(ens_gene_id=rownames(cts))
  ref <- annotables::grch38 %>%
    dplyr::select(ensgene,symbol,chr) %>%
    dplyr::rename(ens_gene_id=ensgene) %>%
    dplyr::filter(ens_gene_id %in% meta_genes$ens_gene_id)
  ref <- ref[!duplicated(ref$ens_gene_id),]
  dim(ref)
  meta_genes <- meta_genes %>%
    dplyr::left_join(ref, by='ens_gene_id')
  dim(meta_genes) ## some genes do not have gene symbol 
  sum(meta_genes$ens_gene_id==rownames(cts))
  
  
  ## Save data as SingleCellExperiment object for future uses
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts=as.matrix(cts), 
                                                         normcounts=as.matrix(df_normalized)),
                                                    colData=coldata,
                                                    rowData=meta_genes)
  
  print(sce)
  saveRDS(sce, paste0(save_dir, datatag, '_cloneBC_sce.rds'))
  
  
}


create_edgeR_object <- function(sce){
  # sce <- readRDS(paste0(save_dir, datatag, '_cloneBC_sce.rds'))
  meta_df <- as.data.frame(colData(sce))
  counts_df <- as.data.frame(counts(sce))
  library("edgeR")
  meta_df$condition <- meta_df$mainsite
  groups_use <- c('Metastasis','Primary')
  #edgeR package use inverted order compared to general convention in other DE tools
  # i.e: Rx vs. UnRx
  # edgeR: 1_UnRx vs. 2_Rx
  meta_df$condition <- ifelse(grepl('Primary', meta_df$condition),paste0("1_",'Primary'),
                              paste0("2_",'Metastasis'))
  print(meta_df$condition)
  meta_df$condition <- factor(meta_df$condition,levels = c("1_Primary","2_Metastasis"))
  # rownames(meta_df) <- meta_df$sample
  # counts_df <- load_counts_data(meta_df$sample, input_dir, datatag)
  print(colnames(counts_df))
  
  
  # Creating a DGEList object for use in edgeR.
  groups_use <- meta_df$condition
  dge <- edgeR::DGEList(counts_df, group = groups_use)
  
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
  
  
  design <- model.matrix(~condition, data = meta_df)
  if(is.null(rownames(design))){
    rownames(design) <- meta_df$sample
  }
  # dge <- calcNormFactors(dge, method = "TMM")
  # sf_df <- data.table::fread(paste0(save_dir, datatag, '_sizefactors_Scran.csv.gz'))
  # sf_df <- sf_df %>%
  #   tibble::column_to_rownames('sample')
  # dge$samples$norm.factors <- sf_df[rownames(dge$samples),'size_factor']
  
  dge$samples$norm.factors <- colData(sce)[rownames(dge$samples),'size_factor']
  dge <- edgeR::estimateDisp(dge, design = design)
  # dge <- edgeR::estimateCommonDisp(dge, design = design)
  # dge <- edgeR::estimateTagwiseDisp(dge)
  # dge <- edgeR::estimateGLMCommonDisp(dge, design = design)
  
  print(length(dge$tagwise.dispersion))
  return(dge)
  
}

get_bootstrap_stat <- function(){
  input_dir <- "/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/metastasis_material/materials/bulkRNAseq/preprocessed_09April2024/"
  input_dir <- "/Users/hoatran/Documents/BCCRC_projects/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/preprocessed_09April2024/"
  datatag <- 'SA919Fig6'
  save_dir <- input_dir
  meta_genes <- data.table::fread(paste0(save_dir, 'dispersion_cis_trans_genes_cloneCB.csv'))
  dim(meta_genes)
  summary(as.factor(meta_genes$gene_type))
  cis_disp <- meta_genes %>%
    dplyr::filter(gene_type=='cis') %>%
    dplyr::pull(tagwise_dispersion)
  trans_disp <- meta_genes %>%
    dplyr::filter(gene_type=='trans') %>%
    dplyr::pull(tagwise_dispersion)
  cis_genes <- cis_disp
  trans_genes <- trans_disp
  
  cis_disp <- cis_logCPM
  trans_disp <- trans_logCPM
  # F-test
  res.ftest <- var.test(tagwise_dispersion ~ gene_type, data = meta_genes)
  res.ftest
  summary(cis_genes)
  
}

## To Do: report for clone B, and clone C separately
get_avg_exp <- function(sce, cis_genes, meta_samples){
  
  ## Check if we have clone B, C in sce meta samples? 
  
  ## To Do: report for clone B, and clone C separately
  # sce$main_clone
  # genes_used <- rownames(sce)
  # genes_used <- intersect(genes_used, meta_genes$ens_gene_id)
  # cis_genes_used <- genes_used[genes_used %in% cnv$ensembl_gene_id]
  # sce_cis <- sce[cis_genes_used,]
  # trans_genes_used <- genes_used[!genes_used %in% cnv$ensembl_gene_id]
  # length(trans_genes_used)
  # length(cis_genes_used)
  # sce_trans <- sce[trans_genes_used,]
  # 
  # cis_logCPM <- aveLogCPM(as.matrix(normcounts(sce_cis)))
  # length(cis_logCPM)
  # summary(cis_logCPM)
  # 
  # trans_logCPM <- aveLogCPM(as.matrix(normcounts(sce_trans)))
  # length(trans_logCPM)
  # summary(trans_logCPM)
}
get_gene_wise_dispersion_edgeR_cloneBC <- function(sce, datatag, save_dir){
  input_dir <- "/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/metastasis_material/materials/bulkRNAseq/preprocessed_09April2024/"
  datatag <- 'SA919Fig6'
  save_dir <- input_dir
  sce <- readRDS(paste0(save_dir, datatag, '_cloneBC_sce.rds'))
  dim(sce)
  
  dge <- create_edgeR_object(sce)
  saveRDS(dge, paste0(save_dir,'figs/',datatag,"_dge.rds"))
  #Visualize the dispersion estimates with a BCV plot
  png(paste0(save_dir,'figs/',datatag,"_glmQLF_plotBCV_check.png"), height = 2*350, width=2*400, res = 2*72)
  # jpeg("glmQLF_plotBCV.jpg")
  plotBCV(dge)
  dev.off()
  
  # A[1:10]
  # length(A)
  # summary(A)
  # A <- dge$AveLogCPM
  # if (is.null(A)) 
  #   A <- aveLogCPM(dge$counts, offset = getOffset(y))
  if (!is.null(dge$tagwise.dispersion)) {
    print('Get gene wise dispersion')
    print(length(dge$tagwise.dispersion))
    df <- data.frame(ens_gene_id=rownames(dge), 
                     tagwise_dispersion=dge$tagwise.dispersion)
    
    meta_genes <- rowData(sce) %>% as.data.frame() 
    meta_genes <- meta_genes %>%
      dplyr::inner_join(df, by='ens_gene_id')
    dim(meta_genes)
    cnv <- data.table::fread(paste0(save_dir, 'cnv_variance_CB.csv.gz'))
    bc <- data.table::fread(paste0(input_dir, 'SA919_CMetMainMix_Bpri_DESeq2/DE_signif_genes.csv.gz'))
    dim(bc)
    # meta_genes <- meta_genes %>%
    #   dplyr::filter(ens_gene_id %in% bc$ens_gene_id)%>%
    #   dplyr::mutate(gene_type=
    #                   case_when(ens_gene_id %in% cnv$ensembl_gene_id ~ 'cis',
    #                             TRUE ~ 'trans'))
    meta_genes <- meta_genes %>%
      # dplyr::filter(ens_gene_id %in% bc$ens_gene_id)%>%
      dplyr::mutate(gene_type=
                      case_when(ens_gene_id %in% cnv$ensembl_gene_id ~ 'cis',  # not correct, need to intersect with DE genes
                                TRUE ~ 'trans'))
    meta_genes$tagwise_dispersion <- round(meta_genes$tagwise_dispersion, 3)
    data.table::fwrite(meta_genes, paste0(save_dir, 'dispersion_cis_trans_genes_cloneCB.csv'))
    
    
    ## Testing
    st <- meta_genes %>%
      dplyr::group_by(gene_type) %>%
      dplyr::summarise(median_dispersion=median(tagwise_dispersion),
                       # mean_dispersion=mean(tagwise_dispersion),
                       nb_genes=n())
    # View(st)
    # st
      
    dim(meta_genes)
    # View(head(meta_genes))
    cis_df <- meta_genes %>%
      dplyr::filter(gene_type=='cis')
    data.table::fwrite(cis_df, paste0(save_dir, 'dispersion_cis_genes_cloneCB.csv'))
    
    
    library("ggplot2")
    p <- ggplot(meta_genes, aes(x=gene_type, y=tagwise_dispersion)) + #
      geom_boxplot(outlier.shape = NA)+
      geom_jitter(size=0.1, aes(colour = gene_type), position=position_jitter(0.2))  + 
      theme_bw() + 
      labs(x='Gene type', y='Gene wise dispersion')
    
    png(paste0(save_dir,'figs/',datatag,"_cis_trans_disperion_cloneBC.png"), height = 2*300, width=2*300, res = 2*72)
    print(p)
    dev.off()
    
    
    # cis_disp <- meta_genes %>%
    #   dplyr::filter(gene_type=='cis') %>%
    #   dplyr::pull(tagwise_dispersion)
    # trans_disp <- meta_genes %>%
    #   dplyr::filter(gene_type=='trans') %>%
    #   dplyr::pull(tagwise_dispersion)
    # # ks.test(cis_disp, trans_disp, alternative='less')
    # ks.test(cis_disp, trans_disp)
    # ?ks.test
    
  }  
    
  
  #Estimate and view the QL dispersions
  fit <- glmQLFit(list, design, robust=TRUE)
  head(fit$coefficients)
  #Plot to the QL dispersions and write to file
  jpeg("glmQLF_plotQLDisp.jpg")
  plotQLDisp(fit)
  dev.off()
  
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
