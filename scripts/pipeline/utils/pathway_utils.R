suppressPackageStartupMessages({
  require(RColorBrewer)
  require(fgsea)
  require(ggplot2)
  require(data.table)
  require(dplyr)
  require(ggrepel)
  options(dplyr.summarise.inform = FALSE)
  options(tidyverse.quiet = TRUE)
  # require(DOSE)
})
library(extrafont)
font_import(prompt=F, paths ='/usr/share/fonts/truetype/myfonts/') # import Helvetica font
fonts()

# A vector of genes as input
get_gprofiler_pathways_obsgenes <- function(obs_genes_symb, save_dir, datatag, 
                                   custom_id=NULL, pathway_fn=NULL, save_data=F){
  library(gprofiler2)
  if(is.null(pathway_fn)){
    pathway_fn = '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/pathway_set/h.all.v7.0.symbols.gmt'  
  }
  
  ref_set <- fgsea::gmtPathways(pathway_fn)
  # for(s in names(ref_set)){ ## just quick check the size of each set
  #   print(s)
  #   print(length(ref_set[[s]]))
  # }
  if(is.null(custom_id)){
    custom_id <- gprofiler2::upload_GMT_file(pathway_fn)  
  }
  stat <- NULL
  ## correction_method: one of 'fdr', 'gSCS', 'bonferroni'
  gostres <- gprofiler2::gost(list(obs_genes_symb), organism = custom_id, correction_method='gSCS')
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
      added_time <- gsub(':','',format(Sys.time(), "%Y%b%d_%X"))
      data.table::fwrite(stat, paste0(save_dir, 'pathways_',added_time,'.csv.gz'))  
    }
  }  
  return(stat)
  
}  
get_gprofiler_pathways <- function(genes_df, save_dir, datatag, 
                                   custom_id=NULL, pathway_fn=NULL, save_data=F){
  library(gprofiler2)
  if(is.null(pathway_fn)){
    pathway_fn = '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/pathway_set/c2.cp.kegg.v7.1.symbols.gmt'  
  }
  
  ref_set <- fgsea::gmtPathways(pathway_fn)
  if(is.null(custom_id)){
    custom_id <- gprofiler2::upload_GMT_file(pathway_fn)  
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
get_custom_pathway_results <- function(deg_df,                      # named vector of statistical significance 
                                       desc='', base_name = '',  
                                       pathway_name=c('custom_pathways','hallmark'), #'cosmic' or 'cisplatin_resistance', or 'metastasis'
                                       groups_use=c("UT","UU"),    # vector of 2 elements, 2 group name used for DE analysis
                                       save_dir = "/home/htran/",
                                       reference_genes_set=NULL,
                                       n_top=20,
                                       ref_dif=NULL){
  # library(fgsea)
  save_dir_log <- save_dir
  save_dir <- paste0(save_dir, pathway_name, '/')
  if(!dir.exists(save_dir)){
    dir.create(save_dir)
  }
  print(save_dir)
  # save_fn_txt = paste0(base_name,"_",desc,".txt")
  # print(save_fn_txt)
  # Can use deframe: first col as name, second col as value, deg_stat <- deframe(markers_ls_output)
  # colFC <- colnames(deg_df)[grepl('FC',colnames(deg_df))]
  # if(!is.null(colFC) & colFC!='logFC'){
  #   colnames(deg_df)[which(colnames(deg_df)==colFC)] <- 'logFC'
  # }else{
  #   stop('Do not exist logFC column, double check input data!!!')
  # }
  # gs <- colnames(deg_df)[grepl('symb',colnames(deg_df))]
  # if(!is.null(gs) & gs!='gene_symb'){
  #   colnames(deg_df)[which(colnames(deg_df)==gs)] <- 'gene_symb'
  # }else{
  #   stop('Do not exist gene symbol column, double check input data!!!')
  # }
  
  # Deal with overlapping and NA gene symbols
  deg_df <- deg_df %>%
    dplyr::select(gene_symbol, logFC) %>%
    na.omit() %>%
    distinct() %>%
    group_by(gene_symbol) %>%
    summarize(logFC=mean(logFC))
  
  deg_stat <- deg_df$logFC
  names(deg_stat) <- deg_df$gene_symbol   # TO DO: check duplicated genes, remove duplicated columns cause gene symbol are not unique in some case
  # head(deg_stat)
  
  if(is.null(ref_dif)){
    ref_dif <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/'  
  }
  
  # ref_genes_ls <- ref_genes_ls[[pathway_name]]
  # pathways: List of gene sets to check, ex in data(examplePathways)
  # stats: Named vector of gene-level stats. Names should be the same as in ’pathways’
  # nperm: Number of permutations to do. Minimial possible nominal p-value is about 1/nperm
  print(pathway_name)
  if(pathway_name=='custom_pathways'){
    ref_genes_ls <- load_reference_genes_set(ref_dif)
    if(is.null(reference_genes_set)){ # select all cancer reference sets
      ref_genes_ls <- ref_genes_ls[!names(ref_genes_ls) %in% c('METASTASIS')] #,'BS_ESSENTIAL_CANCER'
      # ref_genes_ls <- list(ref_genes_ls[['COSMIC']])
      # names(ref_genes_ls) <- 'COSMIC'
    }else{
      ref_genes_ls <- list(ref_genes_ls[['COSMIC']])
      names(ref_genes_ls) <- 'COSMIC'  # more than 500 reference genes
    }
            # maxsz <- length(ref_genes_ls1)
    # gsea_out1 <- fgsea(pathways=ref_genes_ls1, stats=deg_stat, nperm=10000, maxSize = maxsz)  
    # print(dim(gsea_out1))
    # gsea_out <- fgsea(pathways=ref_genes_ls, stats=deg_stat, nperm=10000) 
    gsea_out <- fgsea(pathways=ref_genes_ls, stats=deg_stat) 
    print(dim(gsea_out))
  } else if(pathway_name=='hallmark'){
    gmt_fn <- paste0(ref_dif,'pathway_set/h.all.v7.0.symbols.gmt')
    ref_set <- fgsea::gmtPathways(gmt_fn)
    # class(Hs.H) # a list 
    # class(Hs.H$HALLMARK_HYPOXIA) # a vector of characters
    # gsea_out <- fgsea(pathways=Hs.H, stats=deg_stat, nperm=10000)
    gsea_out <- fgsea(pathways=ref_set, stats=deg_stat)  #, scoreType = "pos"  
    # if(grepl('_IU', base_name)){
    #   gsea_out <- fgsea(pathways=ref_set, stats=deg_stat, scoreType = "pos")  #in case all upregulated genes
    # }else{
    #   gsea_out <- fgsea(pathways=ref_set, stats=deg_stat)  #, scoreType = "pos"  
    # }
    
  }  else if(pathway_name=='kegg'){  # TO DO
    gmt_fn <- paste0(ref_dif,'pathway_set/c2.cp.kegg.v7.1.symbols.gmt')
    ref_set <- fgsea::gmtPathways(gmt_fn)
    # class(Hs.H) # a list 
    # class(Hs.H$HALLMARK_HYPOXIA) # a vector of characters
    # gsea_out <- fgsea(pathways=Hs.H, stats=deg_stat, nperm=10000)
    gsea_out <- fgsea(pathways=ref_set, stats=deg_stat) 
    # if(grepl('_IU', base_name)){
    #   gsea_out <- fgsea(pathways=Hs.H, stats=deg_stat, scoreType = "pos")  #in case all upregulated genes
    # }else{
    #   gsea_out <- fgsea(pathways=Hs.H, stats=deg_stat)  #, scoreType = "pos"  
    # }
    # gsea_out <- gsea_out %>%
    #   filter(pval<0.05)
    # gsea_out$pathway
  } else{
    stop('Check input parameters')
  }
  
  plottitle <- paste0(base_name,"_",paste(groups_use, collapse='_'))
  # plottitle <- 'positive_increase_late_time_points_ref_sets'
    # gsea_out$reference_genes_set <- reference_genes_set  
  gsea_out$datatag <- base_name
  gsea_out$desc <- desc
  gsea_out$signf_genes <- ''
  for(i in rep(1:length(gsea_out$pathway), 1)){
    signf_genes <- unlist(gsea_out$leadingEdge[[i]])
    gsea_out$nb_signf_genes[i] <- length(unlist(signf_genes))
    gsea_out$signf_genes[i] <- paste(signf_genes, collapse=',')  
  }
  gsea_out$leadingEdge <- NULL
  
  
  
  # data.table::fwrite(gsea_out, file=paste0(save_dir, reference_genes_set, plottitle,'.csv'))
  
  # topPathways5Up <- gsea_out[(ES > 0) & (padj < 0.05)][head(order(padj), n=n_top), pathway]
  # topPathways5Down <- gsea_out[ES < 0 & (padj < 0.05)][head(order(padj), n=n_top), pathway]
  
  # using p value threshold
  topPathways5Up <- gsea_out[(ES > 0) & (pval < 0.05)][head(order(padj), n=n_top), pathway]
  topPathways5Down <- gsea_out[ES < 0 & (pval < 0.05)][head(order(padj), n=n_top), pathway]
  topPathways <- c(topPathways5Up, rev(topPathways5Down))
 
  # save_fn = paste0("summary_pathway_cls_",base_name,"_",paste(groups_use, collapse='_'),".rds")
  
  if(length(topPathways5Up)>0){
    topPathwaysUp_df <- gsea_out[match(topPathways5Up, pathway)]
    pathway_barplot(topPathwaysUp_df, plottitle, xtitle="padj", save_dir, typePth='upreg', x="Count", color='pval', showCategory=n_top)
    # fwrite(gsea_hallmark[match(topPathways5Up, pathway)], file=paste0(save_dir, "top_up_",save_fn_txt), sep="\t", sep2=c("", " ", ""))
    # top_pathway_ext <- topPathwaysUp_df[topPathwaysUp_df$padj < 0.05,]
    # pw_df_up <- get_hm_data(top_pathway_ext, markers_ls_output, save_dir, base_name, 'up_pathway')
  }
  
  # cat(paste0("\n Number of down pathways: ",length(topPathways5Down)), file = paste0(save_dir_log,"de_analysis_log.txt"), append = TRUE)
  
  if(length(topPathways5Down)>0){
    topPathwaysDown_df <- gsea_out[match(topPathways5Down, pathway)]
    pathway_barplot(topPathwaysDown_df, plottitle, xtitle="padj", save_dir, typePth='downreg', x="Count", color='pval', showCategory=n_top)
    # top_pathway_ext <- topPathwaysDown_df[topPathwaysDown_df$padj < 0.05,]
    # print(dim(top_pathway_ext))
    # pw_df_down <- get_hm_data(top_pathway_ext, markers_ls_output, save_dir, base_name, 'down_pathway')
    
  }  
  if(length(topPathways)>0){
    # save_fn_txt = paste0(base_name,"_",desc,".txt")
    # data.table::fwrite(gsea_out[match(topPathways, pathway)], file=paste0(save_dir, "pathways_",save_fn_txt), sep="\t", sep2=c("", " ", ""))
    gsea_out <- as.data.frame(gsea_out)
    # gsea_out <- gsea_out %>%
    #   filter(padj<0.05)
    
    gsea_out <- gsea_out %>%
      dplyr::filter(pval<0.05)
    # rownames(gsea_out) <- gsea_out$pathway
    # gsea_out$nb_signf_genes <- 0 
    # for(pw in gsea_out$pathway){
    #   genes <- unique(unlist(strsplit(unlist(gsea_out[pw,'leadingEdge']),',')))
    #   gsea_out[pw,'nb_signf_genes'] <- length(genes) 
    # }
    # print(gsea_out$nb_signf_genes)
    data.table::fwrite(gsea_out,paste0(save_dir, "signf_pathways_",base_name,"_",desc,".csv"))
    return(gsea_out) #[padj < 0.05]
  } else{
    print("Do not exist any pathway!!!")
    return(NULL)
  }
  
}


get_pathway_results <- function(markers_ls_output,                      # named vector of statistical significance 
                                base_name = '',    
                                gmt_fn="biodatabase/h.all.v7.0.symbols.gmt",
                                pathway_name='hallmark',
                                groups_use=c("UT","UU"),    # vector of 2 elements, 2 group name used for DE analysis
                                save_dir = "/home/htran/",
                                n_top=20){
  
  
  # library(fgsea)
  save_dir_log <- save_dir
  save_dir <- paste0(save_dir,pathway_name,'/')
  if(!dir.exists(save_dir)){
    dir.create(save_dir)
  }
  print(save_dir)
  # Can use deframe: first col as name, second col as value, deg_stat <- deframe(markers_ls_output)
  deg_stat <- markers_ls_output$logFC
  names(deg_stat) <- markers_ls_output$gene_symb   # TO DO: check duplicated genes, remove duplicated columns cause gene symbol are not unique in some case
  Hs.H <- gmtPathways(gmt_fn) 
  
  # Hs.c5 <- gmtPathways(paste0(gmt_dir, gmt_ls[2])) 
  # Hs.c6 <- gmtPathways(paste0(gmt_dir, gmt_ls[3])) 
  
  # pathways: List of gene sets to check, ex in data(examplePathways)
  # stats: Named vector of gene-level stats. Names should be the same as in ’pathways’
  # nperm: Number of permutations to do. Minimial possible nominal p-value is about 1/nperm
  gsea_hallmark <- fgsea(pathways=Hs.H, stats=deg_stat, nperm=10000)
  # gsea_go <- fgsea(pathways=Hs.c5, stats=deg_stat, nperm=10000)
  # gsea_oncogenic <- fgsea(pathways=Hs.c6, stats=deg_stat, nperm=10000)
  
  plottitle <- paste0(base_name,"_",groups_use[1],"_vs_",groups_use[2])
  
  
  topPathways5Up <- gsea_hallmark[(ES > 0) & (padj < 0.05)][head(order(padj), n=n_top), pathway]
  topPathways5Down <- gsea_hallmark[ES < 0 & (padj < 0.05)][head(order(padj), n=n_top), pathway]
  topPathways <- c(topPathways5Up, rev(topPathways5Down))
  save_fn = paste0("summary_pathway_cls_",base_name,"_",groups_use[1],"_",groups_use[2],".rds")
  save_fn_txt = paste0("pathway_cls_",base_name,"_",groups_use[1],"_",groups_use[2],".txt")
  
  cat(paste0("\n Number of up pathways: ",length(topPathways5Up)), file = paste0(save_dir_log,"de_analysis_log.txt"), append = TRUE)
  
  if(length(topPathways5Up)>0){
    
    # png(paste0(save_dir,"up_pathway_cls_",base_name,"_",groups_use[1],"_",groups_use[2],".png"), height = 2*500, width=2*900,res = 2*72)
    #   plotGseaTable(Hs.H[topPathways5Up], deg_stat, gsea_hallmark, gseaParam = 0.5)
    # dev.off()
    
    topPathwaysUp_df <- gsea_hallmark[match(topPathways5Up, pathway)]
    pathway_barplot(topPathwaysUp_df, plottitle, xtitle="padj", save_dir, typePth='upreg', x="Count", color='pval', showCategory=n_top)
    # fwrite(gsea_hallmark[match(topPathways5Up, pathway)], file=paste0(save_dir, "top_up_",save_fn_txt), sep="\t", sep2=c("", " ", ""))
    
    top_pathway_ext <- topPathwaysUp_df[topPathwaysUp_df$padj < 0.05,]
    # pw_df_up <- get_hm_data(top_pathway_ext, markers_ls_output, save_dir, base_name, 'up_pathway')
    
  }
  
  
  cat(paste0("\n Number of down pathways: ",length(topPathways5Down)), file = paste0(save_dir_log,"de_analysis_log.txt"), append = TRUE)
  
  if(length(topPathways5Down)>0){
    # print("Contain down-pathway")
    # print(length(topPathways5Down))
    # png(paste0(save_dir,"down_pathway_cls_",base_name,"_",groups_use[1],"_",groups_use[2],".png"), height = 2*500, width=2*900,res = 2*72)
    #   plotGseaTable(Hs.H[topPathways5Down], deg_stat, gsea_hallmark, gseaParam = 0.5)
    # dev.off()
    topPathwaysDown_df <- gsea_hallmark[match(topPathways5Down, pathway)]
    pathway_barplot(topPathwaysDown_df, plottitle, xtitle="padj", save_dir, typePth='downreg', x="Count", color='pval', showCategory=n_top)
    # fwrite(gsea_hallmark[match(topPathways5Down, pathway)], file=paste0(save_dir, "top_down_",save_fn_txt), sep="\t", sep2=c("", " ", ""))  
    
    top_pathway_ext <- topPathwaysDown_df[topPathwaysDown_df$padj < 0.05,]
    print(dim(top_pathway_ext))
    # pw_df_down <- get_hm_data(top_pathway_ext, markers_ls_output, save_dir, base_name, 'down_pathway')
    
  }  
  if(length(topPathways)>0){
    # png(paste0(save_dir,"summary_pathway_cls_",base_name,"_",groups_use[1],"_",groups_use[2],".png"), height = 2*800, width=2*900,res = 2*72)
    #   plotGseaTable(Hs.H[topPathways], deg_stat, gsea_hallmark, gseaParam = 0.5)
    # dev.off()
    
    # collapsedPathways <- collapsePathways(gsea_hallmark[order(pval)][padj < 1],   #padj < 0.01
    #                                       Hs.H, deg_stat)
    # mainPathways <- gsea_hallmark[pathway %in% collapsedPathways$mainPathways][
    #   order(-NES), pathway]
    # if(length(mainPathways)>0){
    #   png(paste0(save_dir,"summary_collapse_pathway_cls_",base_name,"_",groups_use[1],"_",groups_use[2],".png"), height = 2*800, width=2*900,res = 2*72)
    #   plotGseaTable(Hs.H[mainPathways], deg_stat, gsea_hallmark, 
    #                 gseaParam = 0.5)
    #   dev.off()  
    #   fgseaResMain <- gsea_hallmark[match(mainPathways, pathway)]
    #   # fgseaResMain[, leadingEdge := lapply(leadingEdge, mapIds, x=org.Hs.eg.db, keytype="ENTREZID", column="SYMBOL")]
    #   fwrite(fgseaResMain, file=paste0(save_dir, "main_",save_fn_txt), sep="\t", sep2=c("", " ", ""))
    # }
    
    topPathways_df <- gsea_hallmark[match(topPathways, pathway)]
    
    
    # pathway_barplot(topPathways_df, plottitle, xtitle="padj", save_dir, typePth='total', x="Count", color='pval', showCategory=n_top)
    
    pathway_ls <- list(
      stat = deg_stat,
      topPathwaysUp = topPathways5Up,
      topPathwaysDown = topPathways5Down,
      gsea_hallmark = gsea_hallmark,
      base_name = base_name, 
      groups = groups_use
    )
    
    # library(data.table)
    fwrite(gsea_hallmark, file=paste0(save_dir, save_fn_txt), sep="\t", sep2=c("", " ", ""))
    # Readable version
    # library(org.Hs.eg.db)
    fwrite(gsea_hallmark[match(topPathways, pathway)], file=paste0(save_dir, "summary_topPathways_",save_fn_txt), sep="\t", sep2=c("", " ", ""))
    
    saveRDS(pathway_ls, file = paste0(save_dir, save_fn))
    # print(paste0("Save pathway ",save_fn," in the directory ",save_dir))
    return(pathway_ls)
    
  } else{
    print("Do not exist any pathway!!!")
    return(FALSE)
  }
  
}

# deg_stat <- markers_ls_output$avg_log2FC
# names(deg_stat) <- markers_ls_output$gene_symb
load_pathway <- function(deg_df, save_dir, base_name='SA', 
                         groups_use='', gmt_dir=NULL, gmt_ls=NULL, pathway_names=NULL){
  
  stat_col <- colnames(deg_df)[grepl('FC',colnames(deg_df))]
  if(!is.null(stat_col) & stat_col!='logFC'){
    colnames(deg_df)[which(colnames(deg_df)==stat_col)] <- 'logFC'
  }
  sub_dir <- paste0(base_name,"_",paste(groups_use, collapse = '_'))
  save_dir_pw <- paste0(save_dir,sub_dir,"/")
  if (!file.exists(save_dir_pw)){
    dir.create(save_dir_pw)
  }
  
  print("Pathway analysis....")
  if(is.null(gmt_dir)){
    gmt_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/"  
  }
  if(is.null(gmt_ls)){
    gmt_ls <- c("h.all.v7.0.symbols.gmt") #,"c2.cp.kegg.v7.1.symbols.gmt","GO_c5.all.v7.1.symbols.gmt"
    # pathway_names <- c("hallmark","kegg") #,"go"
    pathway_names <- c("hallmark")
  }  
  
  # gmtfile <- paste0(gmt_dir, gmt_ls[1])
  print("Get pathway")
  # Deal with overlapping and NA gene symbols, avg_log2FC
  deg_df <- deg_df %>%
    dplyr::select(gene_symb, logFC) %>%
    na.omit() %>%
    distinct() %>%
    group_by(gene_symb) %>%
    summarize(logFC=mean(logFC))
  
  
  for(i in rep(1:length(gmt_ls),1)){
    # gmt_fn <- paste0(gmt_dir, gmt_ls[i])
    # print(gmt_ls[i])
    print(paste0("\n Pathway analysis, and using gene sets: ",
                 gmt_ls[i]," pathway is: ",pathway_names[i]))
    # print(pathway_names[i])
    pathway_ls <- get_custom_pathway_results(deg_df, 
                                      base_name, paste0(gmt_dir, gmt_ls[i]), 
                                      pathway_names[i],
                                      groups_use,
                                      save_dir_pw, 30)
    pathway_ls <- get_pathway_results(deg_df, 
                                      base_name, paste0(gmt_dir, gmt_ls[i]), 
                                      pathway_names[i],
                                      groups_use,
                                      save_dir_pw, 30)
  }  
}


pathway_barplot <- function(df, plottitle="Pathway", xtitle="pval", save_dir="",
                            typePth='upreg', x="Count", color='padj', showCategory=10) {
  # library(enrichplot)
  # library(DOSE)
  ## source code based on barplot function in enrichplot package 
  ## See more at https://bioconductor.org/packages/release/bioc/html/enrichplot.html
  if(nrow(df)>showCategory){
    df <- df[1:showCategory,]  
  }
  
  
  colorBy <- match.arg(color, c("pval", "p.adjust", "qvalue","padj"))
  if (x == "geneRatio" || x == "GeneRatio") {
    x <- "GeneRatio"
  } else if (x == "count" || x == "Count") {
    x <- "Count"
  } else{
    x <- "Count"
  }
  
  # change x to log fold change here 
  if(typePth=="upreg"){
    lowC <- "#006400"
    highC <- "#90EE90"
    # plottitle <- paste0(plottitle,"_up-regulated")
    plottitle <- 'Up-regulated pathways'
  } else if(typePth=="downreg"){
    lowC <- "#4B0082"
    highC <- "#EE82EE"
    # plottitle <- paste0(plottitle,"_down-regulated")
    plottitle <- 'Down-regulated pathways'
  } else{
    lowC <- "red"
    highC <- "blue"
    plottitle <- paste0(plottitle,"_total")
  }
  # df <- fortify(object, showCategory=showCategory, by=x, ...)
  if(!is.data.frame(df)){
    df <- as.data.frame(df)
  }
  # remove pathway_ str at the beginning of each description, change the text from capital to lower character
  # pathway_ls <- df$pathway
  # pathway_ls2 <- c()
  # for(p in pathway_ls){
  #   p <- str_sub(p, 10, str_length(p))
  #   pathway_ls2 <- c(pathway_ls2, p)
  # }
  # df$pathways <- pathway_ls2
  
  
  if(colorBy %in% colnames(df)) {
    p <- ggplot(df, aes(x = reorder(pathway,-padj), y = padj, fill = pval)) +
      DOSE::theme_dose(12) +
      # scale_fill_continuous(type = "gradient", name = color, guide=guide_colorbar(reverse=TRUE))
      scale_fill_continuous(low=lowC, high=highC, name = color, guide=guide_colorbar(reverse=TRUE)) 
  } else {
    p <- ggplot(df, aes_string(x = "pathway", y = xtitle)) + theme_dose(12) 
  }
  p <- p + geom_bar(stat = "identity", width=0.4) + coord_flip() 
  # ggtitle(title) + xlab(NULL) + ylab(xtitle)
  p <- p + labs(x="", y=xtitle, title=plottitle)
  p <- p + theme(plot.title = element_text(color="black", size=11, hjust = 0.5),
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                 axis.text.y = element_text(color="black", size=9),
                 axis.text.x = element_text(color="black", size=7, angle=90),
                 axis.title = element_text(color="black", size=15))
  
  plottitle <- gsub(' ','',plottitle)
  png(paste0(save_dir,"barplot_",plottitle,".png"), height = 2*25*nrow(df)+100, width=2*600,res = 2*72)
  print(p)
  dev.off()
  
  return(p)
}

load_reference_genes_set <- function(ref_dif){
  pathway_name=c('cosmic','cisplatin_resistance','broadsanger','core_fitness','metastasis') #'cosmic' or 'cisplatin_resistance', or 'metastasis'
  # ref_dif <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/'
  if(is.null(ref_dif)){
    ref_dif <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/'  
  }
  ref_fn <- paste0(ref_dif, 'reference_cancer_genes_set.rds')
  if(file.exists(ref_fn)){
    ref_genes_ls <- readRDS(ref_fn)
    # print(names(ref_genes_ls))
    return(ref_genes_ls)
  }else{
    ref_genes_ls <- list()
    if('cosmic' %in% pathway_name){
      ref_genes <- read.table(paste0(ref_dif, 'oncogene_cosmic.txt'),sep='\t',header = T, check.names = F, stringsAsFactors = F)
      print(length(ref_genes$Gene_Symbol))  
      ref_genes_ls[['COSMIC']] <- toupper(ref_genes$Gene_Symbol)
    }
    if('cisplatin_resistance' %in% pathway_name){
      ref_genes <- read.table(paste0(ref_dif, 'cisplatin_resistance_genes.txt'),sep='\t',header = T, check.names = F, stringsAsFactors = F)
      print(length(ref_genes$gene_symbol))  
      ref_genes_ls[['CISPLATIN_RESISTANCE']] <- toupper(ref_genes$gene_symbol)
    }
    if('metastasis' %in% pathway_name){
      ref_genes <- data.table::fread(paste0(ref_dif, 'metastasis_genes/Metastasis_ref_genes_HK.csv'), header=T) %>% as.data.frame()
      print(length(ref_genes$gene_symbol))  
      
      ref_genes_ls[['METASTASIS']] <- toupper(ref_genes$gene_symbol)
    }
    # "/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/"
    if('broadsanger' %in% pathway_name){
      ref_genes <- data.table::fread(paste0(ref_dif, 'integrated_broad_sanger_essential_983_genes.csv')) %>% as.data.frame()
      print(length(ref_genes$gene_symbol))  
      ref_genes_ls[['BS_ESSENTIAL_CANCER']] <- toupper(ref_genes$gene_symbol)
    }
    if('core_fitness' %in% pathway_name){
      ref_genes <- data.table::fread(paste0(ref_dif, 'Behan_CFgenes.csv')) %>% as.data.frame()
      colnames(ref_genes)[which(colnames(ref_genes) == "ADAM PanCancer Core-Fitness genes")] <- "gene_symbol"
      print(length(ref_genes$gene_symbol))  
      ref_genes_ls[['CORE_FITNESS']] <- toupper(ref_genes$gene_symbol)
    }
    print(names(ref_genes_ls))
    saveRDS(ref_genes_ls, ref_fn)
    return(ref_genes_ls)
  }
  
}




get_confident_interval_genes <- function(obs_genes, ref_genes, datatag='', used_gene_set=' '){
  print(used_gene_set)
  ref_genes <- unique(ref_genes) # genes symbols from reference sets that need to test significance level
  print(length(ref_genes))
  # CI_out <- get_bootstrap_stat(obs_genes, ref_genes, genome_genes=NULL, nsamples=10000) #nsamples=1000
  # CI_out[['reference_set']] <- used_gene_set
  # CI_out[['datatag']] <- datatag
  ref_dif <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/'
  genome_genes_df <- read.csv(paste0(ref_dif, 'Symbol_ensembl.csv'), check.names = F, stringsAsFactors = F)  
  # dim(genome_genes_df)
  genome_genes <- unique(genome_genes_df$Symbol) # entire genes set
  nb_rand <- 100
  # pval_ls <- c()
  # for(i in seq_len(nb_rand)){
  #   pval <- get_bootstrap_stat(obs_genes, ref_genes, genome_genes, nsamples=1000)
  #   pval_ls <- c(pval_ls, pval)
  # }
  pval_ls <- parallel::mclapply(seq_len(nb_rand), function(x) {
    pval <- get_bootstrap_stat(obs_genes, ref_genes, genome_genes, nsamples=10000)
    return(pval)
  }, mc.cores = 4)
  
  padj <- round(p.adjust(pval_ls, "BH"), 3)
  padj <- padj[1]
  CI_out <- data.frame(reference_set=used_gene_set, datatag=datatag, padj=padj)
  return(CI_out)
}



# round(p.adjust(p, "BH"), 3)
get_bootstrap_stat <- function(our_obs, ref_genes, genome_genes=NULL, nsamples=1000){
  # set.seed(42)
  if(is.null(genome_genes)){
    ref_dif <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/'
    genome_genes_df <- read.csv(paste0(ref_dif, 'Symbol_ensembl.csv'), check.names = F, stringsAsFactors = F)  
    # dim(genome_genes_df)
    genome_genes <- unique(genome_genes_df$Symbol) # entire genes set
  }
  
  # length(genome_genes)
  resamples <- lapply(1:nsamples, function(i) sample(genome_genes, size=length(our_obs), replace = T))
  count_occurance <- function(s){
    return(sum(s %in% ref_genes))
  }
  our <- sum(our_obs %in% ref_genes)
  # print('Ours: ')
  # print(our)
  occurs <- sapply(resamples, count_occurance)
  names(occurs) <- paste0('R',seq(1:length(occurs)))
  # print('Randomized: ')
  # print(summary(occurs))
  occurs['our'] <- our
  # length(occurs)
  r <- rank(occurs) #,ties.method = "max"
  pval <- (nsamples + 1 - r['our'])/nsamples
  # return(list(CI=as.numeric(r['our']),pval=as.numeric(pval)))
  return(as.numeric(pval))
}

get_bootstrap_stat_each_series <- function(deg_df, save_dir, datatag=''){
  if(!dir.exists(save_dir)){ 
    dir.create(save_dir)
  }
  deg_df <- deg_df[!duplicated(deg_df$gene_symbol),]
  dim(deg_df)
  meta_type <- as.data.frame(table(deg_df$gene_type, deg_df$gene_type_module)) %>%
      dplyr::filter(Freq>0) %>%
      dplyr::select(-Freq) %>%
      dplyr::rename(gene_type=Var1, gene_type_module=Var2)
  head(meta_type)
  ref_dif <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/'
  gmt_fn <- paste0(ref_dif,'pathway_set/h.all.v7.0.symbols.gmt')
  ref_set <- fgsea::gmtPathways(gmt_fn)
  
  genes_set <- names(ref_set)
  total_stat <- tibble::tibble()
  for(gt in unique(deg_df$gene_type)){
    stat <- tibble::tibble()
    obs_genes <- deg_df %>% 
      dplyr::filter(gene_type==gt) %>% 
      dplyr::pull(gene_symbol)
    for(used_gene_set in genes_set){
      ref_genes <- ref_set[[used_gene_set]]
      # obs_genes <- deg_df$gene_symbol
      intersect_genes <- intersect(obs_genes, ref_genes)
      thrs <- length(obs_genes)/10
      if(thrs<10){
        thrs <- 10
      }
      if(length(intersect_genes)>=thrs){ #
        ci_df <- get_confident_interval_genes(obs_genes, ref_genes, datatag, used_gene_set)
        ci_df$gene_type <- gt
        ci_df$signif_genes <- paste(intersect_genes, collapse=',')
        ci_df$nb_signf_genes <- length(intersect_genes)
        stat <- dplyr::bind_rows(stat, ci_df)  
      }
      
    }
    total_stat <- dplyr::bind_rows(total_stat, stat)  
  }
  total_stat <- total_stat %>% inner_join(meta_type, by='gene_type')
  data.table::fwrite(total_stat, paste0(save_dir, datatag,'_genes_modules_pathways.csv'))
  return(total_stat)
}

# https://github.com/ctlab/fgsea/issues/87
get_pathway_trajectory_genes <- function(deg_df, save_dir,
                                       desc='', base_name = '',  
                                       pathway_name=c('hallmark','kegg','go'),
                                       reference_genes_set=NULL,
                                       ref_dif=NULL){
  # library(fgsea)
  set.seed(42)
  # save_dir_pw <- paste0(save_dir, pathway_name, '/')
  if(!dir.exists(save_dir)){
    dir.create(save_dir)
  }
  print(save_dir)
  
  # Deal with overlapping and NA gene symbols
  # stat_col <- colnames(deg_df)[grepl(pattern = 'Stat',colnames(deg_df))]
  stat_col <- colnames(deg_df)[grepl(pattern = 'FC',colnames(deg_df))]
  # stat_col <- 'pvalue'
  print(stat_col)
  deg_df <- deg_df %>%
    dplyr::rename(stat=!!sym(stat_col))
  deg_df <- deg_df %>%
    dplyr::select(gene_symbol, stat) %>%
    na.omit() %>%
    distinct() %>%
    group_by(gene_symbol) %>%
    summarize(stat=mean(stat))
  
  deg_stat <- deg_df$stat
  names(deg_stat) <- deg_df$gene_symbol   # TO DO: check duplicated genes, remove duplicated columns cause gene symbol are not unique in some case
  # head(deg_stat)
  print(length(deg_stat))
  
  if(is.null(ref_dif)){
    ref_dif <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/'  
  }
  
  # ref_genes_ls <- ref_genes_ls[[pathway_name]]
  # pathways: List of gene sets to check, ex in data(examplePathways)
  # stats: Named vector of gene-level stats. Names should be the same as in ’pathways’
  # nperm: Number of permutations to do. Minimial possible nominal p-value is about 1/nperm
  print(pathway_name)
  if(pathway_name=='hallmark'){
    gmt_fn <- paste0(ref_dif,'pathway_set/h.all.v7.0.symbols.gmt')
    ref_set <- fgsea::gmtPathways(gmt_fn)
    # gsea_out <- fgsea(pathways=Hs.H, stats=deg_stat, nperm=10000)
    gsea_out <- fgsea(pathways=ref_set, stats=deg_stat, eps=0.0)  # , scoreType = "pos"
  } else if(pathway_name=='kegg'){
    gmt_fn <- paste0(ref_dif,'pathway_set/c2.cp.kegg.v7.1.symbols.gmt')
    ref_set <- fgsea::gmtPathways(gmt_fn)
    gsea_out <- fgsea(pathways=ref_set, stats=deg_stat, eps=0.0) #, scoreType = "pos"
  } else if(pathway_name=='go'){
    gmt_fn <- paste0(ref_dif,'pathway_set/GO_c5.all.v7.1.symbols.gmt')
    ref_set <- fgsea::gmtPathways(gmt_fn)
    gsea_out <- fgsea(pathways=ref_set, stats=deg_stat, eps=0.0) 
  } else if(pathway_name=='reactome'){
    # BiocManager::install('reactome.db')
    stop('Need to convert ens gene id to entrez gene id in this case')
    pathways <- reactomePathways(names(deg_stat)) # TO DO: convert to entrez id names in this case
    fgseaRes <- fgsea(pathways, deg_stat, maxSize=500)
    print(sum(fgseaRes$padj<0.05))
  }else{
    stop('Check input parameters')
  }
  # plottitle <- 'positive_increase_late_time_points_ref_sets'
  # gsea_out$reference_genes_set <- reference_genes_set  
  gsea_out$datatag <- base_name
  gsea_out$desc <- desc
  gsea_out$signf_genes <- ''
  for(i in rep(1:length(gsea_out$pathway), 1)){
    signf_genes <- unlist(gsea_out$leadingEdge[[i]])
    gsea_out$nb_signf_genes[i] <- length(unlist(signf_genes))
    gsea_out$signf_genes[i] <- paste(signf_genes, collapse=',')  
  }
  gsea_out$leadingEdge <- NULL
  
  
  
  # data.table::fwrite(gsea_out, file=paste0(save_dir, reference_genes_set, plottitle,'.csv'))
  
  # topPathways5Up <- gsea_out[(ES > 0) & (padj < 0.05)][head(order(padj), n=n_top), pathway]
  # topPathways5Down <- gsea_out[ES < 0 & (padj < 0.05)][head(order(padj), n=n_top), pathway]
  
  # using p value threshold
  n_top=30
  topPathwaysUp <- gsea_out[(ES > 0) & (pval < 0.05)][head(order(padj), n=n_top), pathway]
  topPathwaysDown <- gsea_out[ES < 0 & (pval < 0.05)][head(order(padj), n=n_top), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  
  # save_fn = paste0("summary_pathway_cls_",base_name,"_",paste(groups_use, collapse='_'),".rds")
  
  
  if(length(topPathways)>0){
    gsea_out <- as.data.frame(gsea_out)
    gsea_out <- gsea_out %>%
      dplyr::filter(padj<0.05)
    # dim(gsea_out)
    data.table::fwrite(gsea_out,paste0(save_dir, pathway_name, "_signf_pathways_",base_name,"_",desc,".csv"))
    print(gsea_out$pathway)
    return(gsea_out) #[padj < 0.05]
  } else{
    print("Do not exist any pathway!!!")
    return(NULL)
  }
  
}

# nb genes in col, genes modules in rows, pathway as annotated labels
viz_pathways_stat <- function(pathway_stat, save_dir){
  my_font <- "Helvetica"
  pathway_stat$datatag <- factor(pathway_stat$datatag, levels=unique(pathway_stat$datatag))
  
  p <- ggplot(pathway_stat, aes(x=pathway, y=datatag, color=-log10(padj))) +
    geom_point(aes(size=log2(nb_signf_genes))) +#
    viridis::scale_color_viridis(discrete=F, alpha=0.8)  + 
    # facet_grid(PDX ~ type, scales="free_y", space='free',drop=T) + # PDX ~ ref_gene, . ~ PDX
    theme_bw() + 
    theme(strip.text = element_text(size=9, color="black", family=my_font),
          strip.background = element_blank(),
          legend.position = "bottom",
          legend.box="vertical", 
          legend.margin=margin(),
          legend.text = element_text(size=8, hjust = 0.5, family=my_font),
          # panel.grid.major.x = element_blank(),
          # panel.grid.minor.x = element_blank(),
          # axis.ticks.y = element_blank(),
          axis.text.y = element_text(size=10, hjust = 0.5, family=my_font),
          # axis.text.y = element_text(size=9, hjust = 0.5, family=my_font, angle = 90),
          # text = element_text(size = 7, hjust = 0.5, family=my_font),
          axis.text.x = element_text(size=10, hjust = 0.5, family=my_font, angle = 90),  #, angle = 90
          # axis.text.y = element_text(size=9, hjust = 0.5, family=my_font),
          plot.title = element_text(size=13, hjust=0.5, family=my_font),# face="bold"
          axis.title = element_blank()) 
  p <- p + labs(title='Significant Pathways')
  p <- p + guides(size=guide_legend(title="log2(#genes)"))#color=guide_legend(title="P-Adj"), 
  # p
  # save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
  png(paste0(save_dir,"trajectory_pathways.png"), height = 2*500, width=2*200, res = 2*72)
  print(p)
  dev.off()
  return(p)
}  
viz_chromatin_barplot <- function(df, datatag, save_dir){
  my_font <- "Helvetica"
  df$analysis <- 'Chromatin Status'
  cols_use <- c('red','green','blue')
  names(cols_use) <- c("Active","Bivalent","Repressed")
  df$gene_type_module <- df$Module
  p <- ggplot(df, aes(x=gene_type_module, y=1, fill=chromatin_status)) + #, color=-log10(padj)
    geom_bar(stat="identity", alpha=0.6) +#4.5
    # viridis::scale_color_viridis(discrete=F, alpha=0.8)  + 
    scale_fill_manual(values = cols_use) + 
    facet_grid(~ patient,
               # strip.position = "bottom",
               switch = "x", space='free',scales = "free_x"
               # nrow = 1
    ) +
    theme_bw()+
    theme(strip.text.x = element_blank(), #element_text(size=11, color="black", family=my_font)
          strip.background = element_blank(),
          legend.position = "bottom",
          legend.box="vertical", 
          legend.margin=margin(),
          legend.text = element_text(size=9, hjust = 0.5, family=my_font),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_text(size=11, hjust = 0.5, family=my_font, color="black"),
          axis.text.y = element_blank(),
          panel.border = element_blank(),
          # axis.text.y = element_text(size=9, hjust = 0.5, family=my_font, angle = 90),
          # text = element_text(size = 7, hjust = 0.5, family=my_font),
          axis.text.x = element_text(size=12, hjust = 0.5, family=my_font, color="black"),  #, angle = 90
          # axis.text.y = element_text(size=9, hjust = 0.5, family=my_font),
          plot.title = element_text(size=13, hjust=0.5, family=my_font))# face="bold") 
  p <- p + labs(title=NULL, y=NULL, x=NULL) #, x='Chromatin status of gene modules'
  p <- p + guides(fill=guide_legend(title="Chromatin St", nrow = 3, override.aes = list(size=0.06)))
  plg <- cowplot::get_legend(p)
  p <- p + theme(legend.position = "none")
  # p <- p + geom_text_repel(data = pathway_stat, aes(label = pathway), size = 3.7)
  # p <- p + annotate("text", x = pathway_stat$datatag, y = pathway_stat$nb_signf_genes+2, 
  #                   label = pathway_stat$pathway, size=3.5)
  # png(paste0(save_dir, datatag, "_trajectory_pathways.png"), height = 2*500, width=2*400, res = 2*72)
  # print(p)
  # dev.off()
  
  # saveRDS(p, paste0(save_dir, datatag, "_chromatin_plt.rds"))
  # saveRDS(plg, paste0(save_dir, datatag, "_chromatin_plt_lg.rds"))
  return(list(p=p, plg=plg))
  
}

viz_pathways_combined_dotplot <- function(pathway_stat, datatag, save_dir){
  
  patient_names <- c('Pt4','Pt5','Pt6')
  names(patient_names) <- c('SA609','SA535','SA1035')
  pathway_stat$patient <- patient_names[pathway_stat$datatag]
  # pathway_stat$pathway <- ifelse(pathway_stat$pathway=='epithelial_mesenchymal_transition','epithelial_mesenchymal_trans',pathway_stat$pathway)
  # pathway_stat$pathway <- gsub('response','res',pathway_stat$pathway)
  my_font <- "Helvetica"
  p <- ggplot(pathway_stat, aes(x=gene_type_module, y=pathway)) + #, color=-log10(padj)
    geom_point(aes(size=log2(nb_signf_genes)+0.2), color='red', alpha=0.7) +#4.5
    # viridis::scale_color_viridis(discrete=F, alpha=0.8)  + 
    facet_grid(~ patient,
               # strip.position = "bottom", #switch = "x", 
               space='free',scales = "free_x", drop = F 
               # nrow = 1
    ) +
    theme_bw()+
    theme(strip.text = element_text(size=11, color="black", family=my_font, face='bold'),
          strip.background = element_blank(),
          legend.position = "bottom",
          legend.box="vertical", 
          legend.margin=margin(),
          legend.text = element_text(size=8, hjust = 0.5, family=my_font),
          # panel.grid.major.x = element_blank(),
          # panel.grid.minor.x = element_blank(),
          # axis.ticks.y = element_blank(),
          axis.title.y = element_text(size=10, hjust = 0.5, family=my_font),
          axis.text.y = element_text(size=11, hjust = 0.5, family=my_font, color="black"),
          # axis.text.y = element_text(size=9, hjust = 0.5, family=my_font, angle = 90),
          # text = element_text(size = 7, hjust = 0.5, family=my_font),
          # axis.text.x = element_text(size=12, hjust = 0.5, family=my_font, color="black"),  #, angle = 90
          axis.text.x = element_blank(),
          # axis.text.y = element_text(size=9, hjust = 0.5, family=my_font),
          plot.title = element_text(size=13, hjust=0.5, family=my_font))# face="bold") 
  p <- p + labs(title=NULL, y=NULL, x=NULL) #x='Gene modules'
  p <- p + guides(size=guide_legend(title="log2\n(#hallmark genes)", nrow = 3))
  plg <- cowplot::get_legend(p)
  p <- p + theme(legend.position = "none")
  
  # p <- p + geom_text_repel(data = pathway_stat, aes(label = pathway), size = 3.7)
  # p <- p + annotate("text", x = pathway_stat$datatag, y = pathway_stat$nb_signf_genes+2, 
  #                   label = pathway_stat$pathway, size=3.5)
  # png(paste0(save_dir, datatag, "_trajectory_pathways.png"), height = 2*500, width=2*400, res = 2*72)
  # print(p)
  # dev.off()
  # saveRDS(p, paste0(save_dir, datatag, "_trajectory_pathways.rds"))
  return(list(p=p, plg=plg))
  
}

# pathway_stat <- pathway_stat1
# nb genes in col, genes modules in rows, pathway as annotated labels
viz_pathways_dotplot <- function(pathway_stat, datatag, save_dir){
  
  pathway_stat$pathway <- ifelse(pathway_stat$pathway=='epithelial_mesenchymal_transition','EMT',pathway_stat$pathway)
  pathway_stat$pathway <- gsub('response','res',pathway_stat$pathway)
  my_font <- "Helvetica"
  p <- ggplot(pathway_stat, aes(x=datatag, y=pathway)) + #, color=-log10(padj)
    geom_point(aes(size=log2(nb_signf_genes)+0.2), color='red', alpha=0.7) +#4.5
    # viridis::scale_color_viridis(discrete=F, alpha=0.8)  + 
    theme_bw()+
    theme(strip.text = element_text(size=9, color="black", family=my_font),
          strip.background = element_blank(),
          legend.position = "bottom",
          legend.box="vertical", 
          legend.margin=margin(),
          legend.text = element_text(size=8, hjust = 0.5, family=my_font),
          # panel.grid.major.x = element_blank(),
          # panel.grid.minor.x = element_blank(),
          # axis.ticks.y = element_blank(),
          axis.title.y = element_text(size=10, hjust = 0.5, family=my_font),
          axis.text.y = element_text(size=11, hjust = 0.5, family=my_font, color="black"),
          # axis.text.y = element_text(size=9, hjust = 0.5, family=my_font, angle = 90),
          # text = element_text(size = 7, hjust = 0.5, family=my_font),
          axis.text.x = element_text(size=12, hjust = 0.5, family=my_font, angle = 90, color="black"),  #, angle = 90
          # axis.text.y = element_text(size=9, hjust = 0.5, family=my_font),
          plot.title = element_text(size=13, hjust=0.5, family=my_font))# face="bold") 
  p <- p + labs(title=NULL, y=NULL, x='Gene modules') 
  p <- p + guides(size=guide_legend(title="log2 (#hallmark genes)"))
  # p  
  
  # p <- p + geom_text_repel(data = pathway_stat, aes(label = pathway), size = 3.7)
  # p <- p + annotate("text", x = pathway_stat$datatag, y = pathway_stat$nb_signf_genes+2, 
  #                   label = pathway_stat$pathway, size=3.5)
  png(paste0(save_dir, datatag, "_trajectory_pathways.png"), height = 2*500, width=2*400, res = 2*72)
  print(p)
  dev.off()
  saveRDS(p, paste0(save_dir, datatag, "_trajectory_pathways.rds"))
  return(p)
  
}


# nb genes in col, genes modules in rows, pathway as annotated labels
viz_pathways <- function(pathway_stat, datatag, save_dir){
  require(ggrepel)
  my_font <- "Helvetica"
  p <- ggplot(pathway_stat, aes(x=datatag, y=nb_signf_genes)) + #, color=-log10(padj)
    geom_point(aes(size=log2(nb_signf_genes)+0.3), color='red', alpha=0.7) +#4.5
    # viridis::scale_color_viridis(discrete=F, alpha=0.8)  + 
    theme_bw()+
    theme(strip.text = element_text(size=9, color="black", family=my_font),
          strip.background = element_blank(),
          legend.position = "none",
          legend.box="vertical", 
          legend.margin=margin(),
          legend.text = element_text(size=8, hjust = 0.5, family=my_font),
          # panel.grid.major.x = element_blank(),
          # panel.grid.minor.x = element_blank(),
          # axis.ticks.y = element_blank(),
          axis.title.y = element_text(size=10, hjust = 0.5, family=my_font),
          axis.text.y = element_text(size=9, hjust = 0.5, family=my_font),
          # axis.text.y = element_text(size=9, hjust = 0.5, family=my_font, angle = 90),
          # text = element_text(size = 7, hjust = 0.5, family=my_font),
          axis.text.x = element_text(size=11, hjust = 0.5, family=my_font, angle = 90),  #, angle = 90
          # axis.text.y = element_text(size=9, hjust = 0.5, family=my_font),
          plot.title = element_text(size=13, hjust=0.5, family=my_font))# face="bold") 
  p <- p + labs(title=NULL, x=NULL, y='# hallmark genes') + 
    ylim(min(pathway_stat$nb_signf_genes)-5, max(pathway_stat$nb_signf_genes)+5)
  
  p <- p + geom_text_repel(data = pathway_stat, aes(label = pathway), size = 3.7)
  # p <- p + annotate("text", x = pathway_stat$datatag, y = pathway_stat$nb_signf_genes+2, 
  #                   label = pathway_stat$pathway, size=3.5)
  png(paste0(save_dir, datatag, "_trajectory_pathways.png"), height = 2*500, width=2*400, res = 2*72)
  print(p)
  dev.off()
  saveRDS(p, paste0(save_dir, datatag, "_trajectory_pathways.rds"))
  return(p)
  
}

# genes_set is a list of 2, 3 vectors
# SNP_pop_1=paste(rep("SNP_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
# SNP_pop_2=paste(rep("SNP_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
# SNP_pop_3=paste(rep("SNP_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
# ex: genes_set <- list(SNP_pop_1 , SNP_pop_2 , SNP_pop_3)
# ex: genes_set_names <- c("SNP pop 1" , "SNP pop 2 " , "SNP pop 3")
# ex: plottitle <- "Four sets Venn Diagram \ngenerated by "
viz_vennDiagram_ggplot <- function(genes_set, save_dir, datatag, 
                                   plottitle='', genes_set_names=NULL,
                                   ht=500, wd=600){
  # names(genes_set) <- genes_set_names
  library("ggVennDiagram")
  set.seed(42)
  if(is.null(genes_set_names)){
    genes_set_names <- names(genes_set)  
  }
  
  if(class(genes_set)!='list'){
    stop('Check input data first')
  }
  p <- ggVennDiagram(genes_set, label_alpha = 0,
                     category.names = genes_set_names,
                     edge_size = 0.1,
  ) +
    ggplot2::scale_fill_distiller(palette = "RdBu") + 
    scale_x_continuous(expand = expansion(mult = .2))
  # ggplot2::scale_fill_manual(values = cols_use)
  
  p <- p + labs(title = plottitle)
  png(paste0(save_dir,datatag,'_',plottitle,"_venn_diagramm_ggplot.png"), 
      height = 2*ht, width=2*wd, res = 2*72)
  print(p)
  dev.off()
  return(p)
}

viz_pathways_barplot <- function(df, xplt='pathway', yplt='nb_signf_genes'){
  df$pathway <- ifelse(df$pathway=='epithelial_mesenchymal_transition','EMT',df$pathway)
  df$pathway <- gsub('response','res',df$pathway)
  
  colnames(df)
  my_font <- "Helvetica"
  p <- ggplot(df) + 
    geom_bar(stat="identity", aes_string(x=xplt, y=yplt),width=0.4, color='grey', fill='grey') + 
    coord_flip() + 
    facet_grid(rows = vars(gene_type_module), scales = "free", space = "free") + 
    theme_bw() + 
    theme(#strip.text.y = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      strip.background = element_rect(color='white', fill='white'),
      axis.text.x = element_text(color="black", size=10, family=my_font),
      axis.text.y = element_text(color="black", size=10, family=my_font)) + 
    labs(x=NULL, y='#signf genes', title = NULL)
  # p
  return(p)
}

# It's not a parameter from the original GSEA paper. It was developed by request for one-tailed tests, when you are interested in either only positive enrichment ("pos") or negateive ("neg"). You can check out discussion here: #27
# 
# There are two particular use cases:
# 
# When your stats vector is positive (not signed, as usual), for example, the absolute value of logFC, and you're interested in deregulated pathways that has more positive logFC compared to a random pathway
# When your stats vector is not balanced, for example, when you do differeitntial expression in scRNA-seq between a single cluster and all other clusters. Usually, there are some specific genes (with high positive metric), but no negative specific genes (so there are no highly negative values, only somewhat negative). Similarly, you're looking for a pathways "positively" specific to the cluster, and there are no "negatively" specific pathways.
# Some people would fine it useful, but we didn't do any analysis to show what are properties of this non-standard score types.