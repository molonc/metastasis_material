suppressPackageStartupMessages({
  library("dplyr")
  library("ggplot2")
  library("DESeq2")
  library("SingleCellExperiment")
  
})
# BiocManager::install("fgsea")
source('/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/metastasis_material/scripts/bulk_rna/bulk_utils.R')


## TO DO 
## Get list of pathways and leading genes 
## Check mitosis pw and genes are belong to cis category? 
## get bootstrap stat test 

# 1. mitotic spindle pathway genes - clone C, 
# proportion of cis genes in this pathway compared to trans, 
# split into chr 5, 10. chi-square test for other chrs except 5, and 10.
# 2. overall, cis connected pathway- belong to trans category.



get_cis_genes <- function(save_dir){
  # script_dir <- "/home/htran/Projects/hakwoo_project/metastasis_material/materials/bulkRNAseq/"
  script_dir <- "/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/metastasis_material/materials/bulkRNAseq/SA919_DE_analysis_DESeq2_Hoa_09April2024/"
  save_dir <- paste0(script_dir, 'test_09Nov2024/')
  cnv <- data.table::fread(paste0(script_dir, 'mapping_gene_cnv_SA919.csv.gz'))
  # cnv <- data.table::fread(paste0(script_dir, 'mapped_wholedata_SA919_grch38.csv.gz'))  
  dim(cnv)
  cnv <- cnv[!duplicated(cnv$ensembl_gene_id),]
  # View(head(cnv))
  colnames(cnv)
  cnv <- cnv %>%
    dplyr::filter(pct_pure_A>0.5 & pct_pure_B > 0.5 & pct_pure_C > 0.5)
  dim(cnv) 
  
  rownames(cnv) <- cnv$ensembl_gene_id
  rv <- rowVars(as.matrix(cnv[,c(5, 6)]))  # B, C median copy number profile
  # rv <- rowVars(as.matrix(cnv[,c(7, 10)]))  # B, C median copy number profile
  # length(rv)
  cnv <- as.data.frame(cnv)
  genes_used <- cnv$ensembl_gene_id[rv>0]
  # length(genes_used)
  cnv <- cnv %>%
    dplyr::filter(ensembl_gene_id %in% genes_used)
  dim(cnv)
  # 447 genes with purity and variance with B, C
  # View(cnv)
  data.table::fwrite(cnv, paste0(save_dir, 'cnv_variance_CB.csv.gz'))
  # cnv <- data.table::fread(paste0(save_dir, 'cnv_variance_CB.csv.gz'))
  # sum(cnv$B < cnv$C)
  return(cnv)
}

get_data_cloneBC <- function(){
  input_dir <- "/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/metastasis_material/materials/bulkRNAseq/preprocessed_09April2024/"
  sce <- readRDS(paste0(input_dir, 'sce_cloneBC.rds'))
  dim(sce)
  dim(colData(sce))
  head(colData(sce))
  assayNames(sce)
}
get_cis_genes_DE <- function(){
  input_dir <- "/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/metastasis_material/materials/bulkRNAseq/SA919_DE_analysis_DESeq2_Hoa_09April2024/"
  save_dir <- paste0(input_dir, 'test_09Nov2024/')
  
  cnv <- data.table::fread(paste0(save_dir, 'cnv_variance_CB.csv.gz'))
  dim(cnv)
  bb <- data.table::fread(paste0(input_dir, 'SA919_Bmet_Bpri_DESeq2/DE_signif_genes.csv.gz'))
  bc <- data.table::fread(paste0(input_dir, 'SA919_CMetMainMix_Bpri_DESeq2/DE_signif_genes.csv.gz'))
  # bc <- data.table::fread(paste0(input_dir, 'SA919_CMetMain_Bpri_DESeq2/DE_signif_genes.csv.gz'))
  dim(bb)
  dim(bc)
  
  # genes_used <- intersect(bb$symbol, bc$symbol)
  # genes_used <- intersect(bc$ens_gene_id, cnv$ensembl_gene_id)
  # length(genes_used)
  genes_used <- intersect(bc$symbol, cnv$gene_symbol)
  # genes_used <- intersect(genes_used, cnv$ensembl_gene_id)
  length(genes_used)
  
  
  cis_df <- cnv %>%
    dplyr::filter(gene_symbol %in% genes_used)
  # View(cis_df)
  
  # used_cols <- c('ens_gene_id','logFC') #'symbol',,'chr'
  used_cols <- c('symbol','logFC')
  # bb <- bb %>%
  #   dplyr::select(all_of(used_cols))
  bc <- bc %>%
    dplyr::select(all_of(used_cols))
  # colnames(bb) <- paste0(colnames(bb), '_bb')
  colnames(bc) <- paste0(colnames(bc), '_bc')
  # cis_df <- cis_df %>%
  #   dplyr::inner_join(bb, by=c('gene_symbol'='symbol_bb'))
  cis_df <- cis_df %>%
    dplyr::inner_join(bc, by=c('gene_symbol'='symbol_bc'))
  # cis_df <- cnv %>%
  #   dplyr::inner_join(bc, by=c('ensembl_gene_id'='ens_gene_id'))
  dim(cis_df)
  View(cis_df)
  cis_df$cna_tendency <- ifelse(cis_df$C - cis_df$B>0, 'Gain', 'Loss')
  cis_df$expr_tendency <- ifelse(cis_df$logFC_bc>0, 'Up', 'Down')
  cis_df$desc <- paste0(cis_df$cna_tendency, '-',cis_df$expr_tendency)
  summary(as.factor(cis_df$desc))
  # 100*(66/73) 90.4% positive tendency
  # 100*(7/73) 9.6% negative tendency
  
  data.table::fwrite(cis_df, paste0(save_dir, 'cis_genes_SA919_CB.csv'))
  ## To Do: need to work at gene expression level
}

get_pathways_fgsea_bc <- function(){
  
  input_dir <- "/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/metastasis_material/materials/bulkRNAseq/SA919_DE_analysis_DESeq2_Hoa_09April2024/"
  save_dir <- paste0(input_dir, 'test_09Nov2024/')
  
  # cnv <- data.table::fread(paste0(save_dir, 'cnv_variance_CB.csv.gz'))
  # dim(cnv)
  bb <- data.table::fread(paste0(input_dir, 'SA919_Bmet_Bpri_DESeq2/DE_signif_genes.csv.gz'))
  bc <- data.table::fread(paste0(input_dir, 'SA919_CMetMainMix_Bpri_DESeq2/DE_signif_genes.csv.gz'))
  # bc <- data.table::fread(paste0(input_dir, 'SA919_CMetMain_Bpri_DESeq2/DE_signif_genes.csv.gz'))
  dim(bb)
  dim(bc)
  head(bc)
  bc <- bc %>%
    dplyr::filter(logFC>=1) %>%
    dplyr::rename(gene_symbol=symbol) #logFC=log2FoldChange, 
  bb <- bb %>%
    dplyr::rename(gene_symbol=symbol) #logFC=log2FoldChange, 
  gmt_fn <- '/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/metastasis_material/materials/bulkRNAseq/pathway_set/h.all.v7.0.symbols.gmt'
  base_name <- 'SA919_bc'
  # library("fgsea")
  
  gsea_out <- get_fgsea_pathways(bc, save_dir, base_name, gmt_fn)
  
  base_name <- 'SA919_bb'
  gsea_out <- get_fgsea_pathways(bb, save_dir, base_name, gmt_fn)
  
}
get_genes_pathways <- function(){
  input_dir <- "/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/metastasis_material/materials/bulkRNAseq/SA919_DE_analysis_DESeq2_Hoa_09April2024/"
  save_dir <- paste0(input_dir, 'test_09Nov2024/')
  
  cnv <- data.table::fread(paste0(save_dir, 'cnv_variance_CB.csv.gz'))
  dim(cnv)
  bb <- data.table::fread(paste0(input_dir, 'SA919_Bmet_Bpri_DESeq2/DE_signif_genes.csv.gz'))
  bc <- data.table::fread(paste0(input_dir, 'SA919_CMetMainMix_Bpri_DESeq2/DE_signif_genes.csv.gz'))
  # bc <- data.table::fread(paste0(input_dir, 'SA919_CMetMain_Bpri_DESeq2/DE_signif_genes.csv.gz'))
  dim(bb)
  dim(bc)
  sum(bb$symbol %in% cnv$gene_symbol)
  sum(bc$symbol %in% cnv$gene_symbol)
  obs_genes <- bc$symbol[bc$symbol %in% cnv$gene_symbol]
  obs_genes <- unique(obs_genes)
  length(obs_genes)
  datatag <- 'SA919_cis'
  pw_df <- get_gprofiler_pathways_obsgenes(obs_genes, save_dir, datatag, 
                                           custom_id="gp__GN18_LkQt_dzk", pathway_fn=pathway_fn, 
                                           save_data=T, correction_func='gSCS')
  
  sum(bc$log2FoldChange > 0)
  datatag <- 'bb'
  bb <- bb %>%
    dplyr::filter(logFC>=0.5) 
  obs_genes_symb <- bb$symbol
  pathway_fn <- '/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/metastasis_material/materials/bulkRNAseq/pathway_set/c5.all.v2024.1.Hs.symbols.gmt'
  pathway_fn <- '/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/metastasis_material/materials/bulkRNAseq/pathway_set/h.all.v7.0.symbols.gmt'
  pw_df <- get_gprofiler_pathways_obsgenes(obs_genes_symb, save_dir, datatag, 
                                           custom_id=NULL, pathway_fn=pathway_fn, 
                                           save_data=T, correction_func='gSCS')
  
  ## correction_method: one of 'fdr', 'gSCS', 'bonferroni'
  datatag <- 'bc'
  bc <- bc %>%
    dplyr::filter(logFC>=0.5) 
  obs_genes_symb <- bc$symbol
  pathway_fn <- '/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/metastasis_material/materials/bulkRNAseq/pathway_set/h.all.v7.0.symbols.gmt'
  pw_df <- get_gprofiler_pathways_obsgenes(obs_genes_symb, save_dir, datatag, 
                                           custom_id="gp__WJMi_fsQ3_9pc", pathway_fn=pathway_fn, 
                                           save_data=T, correction_func='gSCS')
  pw_df
  
  signif_genes_mitotic <- unlist(strsplit(pw_df$signif_genes[2], ","))
  sum(signif_genes_mitotic %in% cnv$gene_symbol)
  signif_genes_mitotic[signif_genes_mitotic %in% cnv$gene_symbol]
  signif_genes_myogenesis <- unlist(strsplit(pw_df$signif_genes[1], ","))
  sum(signif_genes_myogenesis %in% cnv$gene_symbol)
  signif_genes_myogenesis[signif_genes_myogenesis %in% cnv$gene_symbol]
  
  
  ref_set <- fgsea::gmtPathways(pathway_fn)
  ref_set$HALLMARK_MITOTIC_SPINDLE
  names(ref_set)
  for(r in names(ref_set)){
    sz <- length(intersect(obs_genes, ref_set[[r]]))
    if(sz>2){
      print(r)
      print(sz)
    }
  }
  ref_dir <- '/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/metastasis_material/materials/bulkRNAseq/pathway_set/'
  cosmic <- data.table::fread(paste0(ref_dir,'cosmic.csv.gz'))
  corefitness <- data.table::fread(paste0(ref_dir,'CoreFitness_corrected.csv.gz'))
  dim(cosmic)
  dim(corefitness)
  head(cosmic)
  intersect(obs_genes,corefitness$gene_symb)
}

## For each pathway, get all genes, and create an unweighted connection between every two of them
## Check cis genes, get how many trans genes in the list. And calculate the percentage of cis-effect
##
##
##
##
##
##
##
get_connected_genes <- function(){
  pathway_fn <- '/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/metastasis_material/materials/bulkRNAseq/pathway_set/h.all.v7.0.symbols.gmt'
  ref_set <- fgsea::gmtPathways(pathway_fn)
  names(ref_set)
  total_df <- tibble::tibble()
  for(r in names(ref_set)){
    # print(r)
    df <- as.data.frame(t(combn(unlist(ref_set[[r]]),2)))
    df$ref_set <- as.character(r)
    total_df <- dplyr::bind_rows(total_df, df)
  }
  dim(total_df)
  dim(df)
  head(total_df)
  input_dir <- "/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/metastasis_material/materials/bulkRNAseq/SA919_DE_analysis_DESeq2_Hoa_09April2024/"
  save_dir <- paste0(input_dir, 'test_09Nov2024/')
  
  cnv <- data.table::fread(paste0(save_dir, 'cnv_variance_CB.csv.gz'))
  dim(cnv)
  bb <- data.table::fread(paste0(input_dir, 'SA919_Bmet_Bpri_DESeq2/DE_signif_genes.csv.gz'))
  bc <- data.table::fread(paste0(input_dir, 'SA919_CMetMainMix_Bpri_DESeq2/DE_signif_genes.csv.gz'))
  
  cis_genes <- intersect(bc$symbol, cnv$gene_symbol)
  length(cis_genes)
  trans_genes <- bc$symbol[!bc$symbol %in% cis_genes]
  total_df1 <- total_df %>%
    dplyr::filter(V1 %in% cis_genes & V2 %in% trans_genes)
  total_df2 <- total_df %>%
    dplyr::filter(V1 %in% trans_genes & V2 %in% cis_genes)
  # combn(c('a','b','c'),2)
  dim(total_df1)
  dim(total_df2)
  obs_genes_all1 <- c(total_df1$V1, total_df1$V2) 
  obs_genes_all2 <- c(total_df2$V1, total_df2$V2)
  obs_genes <- c(obs_genes_all1, obs_genes_all2)
  obs_genes <- unique(obs_genes)
  length(obs_genes)
  sum(obs_genes %in% cis_genes)
  obs_genes[obs_genes %in% cis_genes] #"PTGER4" "PFKFB3" "HMGCS1" "TRIO" "CEP72"  "IL2RA"  "ANKH"   
  sum(obs_genes %in% trans_genes)
  dim(bc)
  426/length(trans_genes) * 100 #4.6% trans genes affected by cis genes
  datatag <- 'SA919_cis_test'
  pw_df <- get_gprofiler_pathways_obsgenes_(obs_genes, save_dir, datatag, 
                                           custom_id=NULL, pathway_fn=pathway_fn, 
                                           save_data=F, correction_func='gSCS')
  dim(pw_df)
  head(pw_df)
  
  mitotic <- 'ERBB3,SPEG,MYH7,ABLIM1,MYH9,PVALB,TNNC2,MYOM1,CKM,BAG1,MYH1,MYH3,TNNC1,TNNT2,SORBS3,SPDEF,TNNT3,TNNI2,ITGA7,CHRNA1,KCNH1,PSEN2,TNNI1,NAV2,TCAP,SMTN,MAPK12,MYL6B,SVIL,SPTAN1,ACTN3,MYH4'
  
  mitotic_genes <- unlist(strsplit(mitotic, ','))
  sum(obs_genes %in% mitotic_genes)
  sum(cis_genes %in% mitotic_genes)
  cis_genes[cis_genes %in% mitotic_genes]
  myo <- 'VCL,TRIO,SMC1A,CLASP1,KIF3C,CTTN,KIF4A,SUN2,MYH9,NIN,KIF3B,CEP72,LMNB1,SPTBN1,ALMS1,ARHGEF2,PREX1,KLC1,MYH10,ESPL1,FLNB,CSNK1D,KIF2C,INCENP,FARP1,TIAM1,WASF2,PCNT,ARF6,ATG4B,CKAP5,BCR,FLNA,SPTAN1,PRC1'
  myo_genes <- unlist(strsplit(myo, ','))
  sum(obs_genes %in% myo_genes)
  length(myo_genes)
}
# 

# to get 1st degree and 2nd degree in trans, you can use the attached influence graph.  
# If a a gene is DE (FDR < 0.01) trans gene  
# if it is directly connected with a DE cis gene, label it in trans 1st deg with FDR < 0.01.  
# else if it is directly connected with a DE 1st deg trans gene, label it trans 2nd deg with FDR < 0.01.  
# else label it in trans other  
# And so on.

## to do: check cis genes: are in sanger genes set, or cosmic gene set?
get_connected_genes_network <- function(){
  ## First, for clone B, C. 
  input_dir <- "/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/metastasis_material/materials/bulkRNAseq/SA919_DE_analysis_DESeq2_Hoa_09April2024/"
  save_dir <- paste0(input_dir, 'test_09Nov2024/')
  
  # cnv <- data.table::fread(paste0(save_dir, 'cnv_variance_CB.csv.gz'))
  # dim(cnv)
  bb <- data.table::fread(paste0(input_dir, 'SA919_Bmet_Bpri_DESeq2/DE_signif_genes.csv.gz'))
  bc <- data.table::fread(paste0(input_dir, 'SA919_CMetMainMix_Bpri_DESeq2/DE_signif_genes.csv.gz'))
  # bc <- data.table::fread(paste0(input_dir, 'SA919_CMetMain_Bpri_DESeq2/DE_signif_genes.csv.gz'))
  dim(bb)
  dim(bc)
  # dir.create(save_dir)
  
  dif_genes <- bc$symbol[!bc$symbol %in% bb$symbol]
  length(dif_genes)
  sum(dif_genes %in% cna$gene_symbol) # 73 cis genes and 66 cis genes are only in C vs B, 7 genes are cis genes in CB, and trans genes in BB
  
  ## what is the 66 cis genes? 
  
  obs_genes_symb <- intersect(bb$symbol, bc$symbol)
  # obs_genes_symb <- intersect(bb$g_ID, bc$g_ID)
  length(obs_genes_symb) #447 genes in common between bb and bc, metastasis genes
  
  cna <- data.table::fread(paste0(save_dir, 'cnv_variance_CB.csv.gz'))
  dim(cna)
  head(cna)
  sum(obs_genes_symb %in% cna$gene_symbol)
  sum(bc$symbol %in% cna$gene_symbol) # 73 cis genes
  # cna <- cna %>%
  #   dplyr::filter(chr %in% c(5, 7, 10))
  dim(cna)
  
 
  var_genes <- cna$gene_symbol
  length(var_genes)
  bc <- bc %>%
    # dplyr::filter(!Gene.name %in% obs_genes_symb) %>%
    dplyr::rename(gene_symbol=symbol)
  
  cis_genes <- intersect(bc$gene_symbol, cna$gene_symbol)
  length(cis_genes) # 73 cis genes
  
  trans_genes <- bc$gene_symbol[!bc$gene_symbol %in% cis_genes]
  length(trans_genes)
  
  input_dir <- "/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/metastasis_material/materials/bulkRNAseq/"
  inf_graph <- read.delim(paste0(input_dir, 'old_analyses/SA919_DE_analysis/influence_graph.txt'), header = TRUE, sep = "\t", dec = ".")
  dim(inf_graph)
  head(inf_graph)
  # to get 1st degree and 2nd degree in trans, you can use the attached influence graph.  
  # If a a gene is DE trans gene  
  # if it is directly connected with a DE cis gene, label it in trans 1st deg.  
  
  inf_graph <- inf_graph %>%
    dplyr::mutate(desc=paste0(gene_a, '_',gene_b,'_',weight))
  
  direct_df <- inf_graph %>%
    dplyr::filter(gene_a %in% trans_genes & gene_b %in% cis_genes)
  dim(direct_df)
  
  direct_df1 <- direct_df %>%
    dplyr::mutate(connect_gene=paste0(gene_a, '_',gene_b)) %>%
    dplyr::group_by(connect_gene, gene_a, gene_b) %>%
    dplyr::summarise(weight=max(weight)) %>%
    dplyr::mutate(desc=paste0(gene_a, '_',gene_b,'_',weight))
  direct_df <- direct_df %>%
    dplyr::filter(desc %in% direct_df1$desc)
  dim(direct_df)
  head(direct_df)
  
  direct_df2 <- direct_df %>%
    dplyr::select(gene_b, desc) %>%
    dplyr::mutate(desc=paste0(desc, ': direct_1st_layer')) %>%
    dplyr::rename(C_vs_B_cis_gene=gene_b, connect_to_cis_gene=desc)
  head(direct_df2)
  dim(direct_df2)
  # View(direct_df2)
  remain_genes <- trans_genes[!trans_genes %in% direct_df$gene_a]
  length(remain_genes)
  head(direct_df)
  length(unique(direct_df$gene_b))
  indirect_df <- inf_graph %>%
    dplyr::filter(gene_a %in% trans_genes & !desc %in% direct_df$desc) %>%
    dplyr::select(desc, gene_b) %>%
    dplyr::rename(desc1=desc, gene_1st=gene_b)
  
  indirect_df2 <- inf_graph %>%
    dplyr::filter(gene_a %in% remain_genes 
                  & gene_b %in% unique(direct_df$gene_a))
  dim(indirect_df2)
  length(unique(indirect_df2$gene_a))
  length(unique(direct_df$gene_a))
  colnames(direct_df)
  colnames(indirect_df2)
  
  # dplyr::full_join(indirect_df, by=c('gene_a'='gene_1st'), relationship = "many-to-many") %>%
  # indirect_df3 <- indirect_df2 %>%
  #   dplyr::group_by(gene_b) %>%
  #   dplyr::summarise(desc2=max(weight)) %>%
  #   dplyr::mutate(desc2=paste0(gene_b, '_',desc2))
  # indirect_df2 <- indirect_df2 %>%
  #   dplyr::mutate(desc2=paste0(gene_b, '_',weight)) %>%
  #   dplyr::filter(desc2 %in% indirect_df3$desc2)
  
  dim(indirect_df2)
  length(unique(indirect_df2$C_vs_B_cis_gene))
  head(indirect_df2)
  indirect_df2 <- indirect_df2 %>%
    dplyr::select(gene_b, desc) %>%
    dplyr::mutate(desc=paste0(desc, ': indirect_2nd_layer')) %>%
    dplyr::rename(C_vs_B_DE_gene=gene_b, connect_to_cis_gene=desc)
  
  indirect_df2 <- indirect_df2[!duplicated(indirect_df2$C_vs_B_DE_gene),]
  direct_df2 <- direct_df2[!duplicated(direct_df2$C_vs_B_DE_gene),]
  total_df <- dplyr::bind_rows(direct_df2, indirect_df2)
  dim(total_df)
  head(total_df)
  dim(direct_df2)
  dim(indirect_df2)
  data.table::fwrite(total_df, paste0(save_dir, 'cloneC_vs_cloneB_DE_genes_cis_connected.csv'))
  
  input_dir <- "/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/metastasis_material/materials/bulkRNAseq/SA919_DE_analysis_DESeq2_Hoa_09April2024/"
  save_dir <- paste0(input_dir, 'test_09Nov2024/')
  
  pw_df <- data.table::fread(paste0(save_dir, 'pathways_bc.csv.gz'))
  dim(pw_df)
  colnames(pw_df)
  pw_df <- as.data.frame(pw_df)
  pw_df$reference_set <- ifelse(pw_df$reference_set=="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_EMT",pw_df$reference_set)
  rownames(pw_df) <- pw_df$reference_set
  pattern_use <- '^hallmark_'
  # geneSets <- list()
  ref_genes <- tibble::tibble()
  for(p in rownames(pw_df)){
    
    gene_ls <- strsplit(as.character(pw_df[p,'signif_genes']),',')
    gene_ls <- unlist(gene_ls, use.names=FALSE)
    p <- gsub(pattern_use,'',tolower(p))
    # p <- gsub(pattern_use,'',toupper(p))
    tmp <- tibble::tibble(gene=gene_ls, significant_pathway=p)
    ref_genes <- dplyr::bind_rows(ref_genes, tmp)
    # geneSets[[p]] <- gene_ls
    # print(length(gene_ls))
  }
  dim(ref_genes)
  head(ref_genes)
  ref_genes <- ref_genes %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(significant_pathway=paste(significant_pathway, collapse = '_'))
  dim(ref_genes)
  
  ref_genes_connected <- ref_genes %>%
    dplyr::filter(gene %in% inf_graph$gene_a)
  ref_genes_connected <- ref_genes %>%
    dplyr::filter(gene %in% direct_df$gene_a)
  dim(ref_genes_connected) 
  cis_genes
  ref_genes$gene
  head(ref_genes)
  total_df1 <- total_df
  head(total_df1)
  dim(total_df1)
  total_df1 <- total_df1 %>%
    left_join(ref_genes, by=c('C_vs_B_DE_gene'='gene'))
  View(total_df1)
  sum(!is.na(total_df1$significant_pathway))
  data.table::fwrite(total_df1, paste0(save_dir, 'cloneC_vs_cloneB_DE_genes_cis_connected.csv'))
}
