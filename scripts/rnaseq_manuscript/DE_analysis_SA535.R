
suppressPackageStartupMessages({
  require(RColorBrewer)
  require(glue)
  require(ggrepel)
  require(tidyr)
  require(fgsea)
  # require(scales)
  require(ComplexHeatmap)
  # require(DOSE)
  # require(enrichplot)
  # require(ggraph)
  require(ggplot2)
  # require(igraph)
  require(tidyverse)
  require(data.table)
  require(Seurat)
  require(dplyr)
})

## DE analysis
main(){
  test_use='wilcox'
  pAdjustThrs=0.05
  minLogFC=0.5
  nbtopup = 30
  nbtopdown = 30
  save_data=T
  
  
  base_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/'
  output_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/normalized_introns/'
  datatag <- 'SA535' 
  sce_tmp <- readRDS(paste0(base_dir,'filtered_introns/SCRNA10X_SA_CHIP0077_004.rds'))
  genes_map <- as.data.frame(rowData(sce_tmp))
  dim(genes_map)
  colnames(genes_map) <- c('ens_gene_id','gene_symb')
  rownames(genes_map) <- genes_map$ens_gene_id
  
  save_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/DE_analysis_introns/'
  dir.create(save_dir)
  input_dir=""
  # base_name <- datatag
  # feature_use="Grouping"
  # 
  srt <- readRDS(paste0(output_dir,'SA535_scTransform_srt.rds'))
  meta_data <- srt@meta.data
  rownames(meta_data)[1]
  meta_data$cell_id <- paste0(meta_data$library_id,'_',meta_data$Barcode)
  umap_df <- data.table::fread(paste0(output_dir, datatag, "_scTransform_umap_filtered_cells_clones.csv.gz"))
  dim(umap_df)  
  dim(srt)
  colnames(umap_df)
  summary(as.factor(umap_df$clone_id))
  umap_df <- umap_df %>%
    dplyr::mutate(clone_id=
                    case_when(
                      clone_id=='None' ~ 'C',
                      TRUE ~ clone_id
                    ))
  umap_df <- umap_df %>%
    dplyr::select(Site_origin,Grouping, pdxid, clone_id, mouse_id, cell_id)
  meta_data <- meta_data %>%
    dplyr::left_join(umap_df, by='cell_id')
  dim(meta_data)
  rownames(meta_data) <- meta_data$cell_id
  
  srt <- AddMetaData(object = srt, metadata = as.data.frame(meta_data)) #, col.name =colnames(meta_data)
  sum(colnames(srt)==srt$cell_id)
  
  # meta_data <- srt@meta.data
  
  
  
  groups_use_ls <- list(T1=c('H','G'))
  # site_ls <- list(G=c("Primary"),H=c("Metastasis"))
  site_ls <- c("Metastasis","Primary")
  
  groups_use_ls <- list(T1=c('H','H'))
  # site_ls <- list(G=c("Primary"),H=c("Metastasis"))
  site_ls <- c("Metastasis","Primary")
  
  groups_use_ls <- list(T1=c('M','M'))
  # site_ls <- list(M=c("Primary"),M=c("Metastasis"))
  site_ls <- c("Metastasis","Primary")
  
  groups_use_ls <- list(T1=c('J','M'))
  # site_ls <- list(M=c("Primary"),M=c("Metastasis"))
  site_ls <- c("Metastasis","Primary")
  
  groups_use_ls <- list(T1=c('P','M'))
  # site_ls <- list(M=c("Primary"),M=c("Metastasis"))
  site_ls <- c("Metastasis","Primary")
  
  groups_use_ls <- list(T1=c('O','M'))
  # site_ls <- list(M=c("Primary"),M=c("Metastasis"))
  site_ls <- c("Metastasis","Primary")
  # srt$clone_site <- paste0(srt$clone_id,'_',srt$Grouping)
  # Idents(object = srt) <- "clone_site"
  # unique(srt$clone_site)
  # meta_data$
  groups_use <- groups_use_ls$T1
  for(groups_use in groups_use_ls){
    
    # cells_use_g1 <- meta_data %>%
    #   dplyr::filter(clone_id==groups_use[1] & Grouping==site_ls[[groups_use[1]]]) %>%
    #   dplyr::pull(cell_id)
    # 
    # cells_use_g2 <- meta_data %>%
    #   dplyr::filter(clone_id==groups_use[2] & Grouping==site_ls[[groups_use[2]]]) %>%
    #   dplyr::pull(cell_id)
    
    cells_use_g1 <- meta_data %>%
      dplyr::filter(clone_id==groups_use[1] & Grouping==site_ls[1]) %>%
      dplyr::pull(cell_id)
    
    cells_use_g2 <- meta_data %>%
      dplyr::filter(clone_id==groups_use[2] & Grouping==site_ls[2]) %>%
      dplyr::pull(cell_id)
    print(length(cells_use_g1))
    print(length(cells_use_g2))
    if(length(cells_use_g1) < 40 || length(cells_use_g2) < 40){
      print(groups_use)
      print("There are no cells or small nb cells 
            only which satisfy the input condition ")
      print(paste0("\n Observed clones: ",groups_use[1],' vs ',groups_use[2], 
                   "  nb cells in ",groups_use[1], ":",length(cells_use_g1),
                   "  nb cells in ",groups_use[2], ":",length(cells_use_g2)))
      
    } else{
      # base_name <- paste0(datatag,'_', groups_use[1],'_',site_ls[[groups_use[1]]],'_', 
      #                     groups_use[2],'_',site_ls[[groups_use[2]]])
      
      base_name <- paste0(datatag,'_', groups_use[1],'_',site_ls[1],'_', 
                          groups_use[2],'_',site_ls[2])
      groups_use_desc <- c(paste0(groups_use[1],'_',site_ls[1]),
                           paste0(groups_use[2],'_',site_ls[2]))
      pathway_ls <- calculate_DE_analysis_v2(base_name=base_name, meta_data=meta_data, 
                                             cells_use_g1, cells_use_g2,
                                             feature_use="clone_site",
                                             srt, genes_map, groups_use = groups_use_desc,
                                             test_use, save_dir, input_dir,
                                             pAdjustThrs=pAdjustThrs, minLogFC=minLogFC,
                                             nbtopup = 30, nbtopdown = 30, 
                                             save_data=T, viz=F)
      
      
    }  
    
    
  }
  
  df1 <- data.table::fread(paste0(save_dir,'SA535_H_Metastasis_G_Primary/','de_significant_genes.csv.gz'))
  df2 <- data.table::fread(paste0(save_dir,'SA535_H_Metastasis_H_Primary/','de_significant_genes.csv.gz'))
  
  df1 <- data.table::fread(paste0(save_dir,'SA535_M_Metastasis_M_Primary/','de_significant_genes.csv.gz'))
  df2 <- data.table::fread(paste0(save_dir,'SA535_O_Metastasis_M_Primary/','de_significant_genes.csv.gz'))
  
  dim(df1)
  dim(df2)  
  sum(df1$avg_log2FC>0)
  sum(df2$avg_log2FC>0)
}
calculate_DE_analysis_v2 <- function(base_name, meta_data, cells_use_g1, cells_use_g2, 
                                     feature_use="Grouping",
                                     srt, genes_map, groups_use, test_use='wilcox', 
                                     save_dir="", input_dir="",
                                     pAdjustThrs=0.05, minLogFC=0.25,
                                     nbtopup = 30, nbtopdown = 30, save_data=T, viz=F){
  
  # print(head(meta_data))
  print(dim(meta_data))
  # sub_dir <- paste0(base_name,"_",groups_use[1],"_",groups_use[2])
  sub_dir <- base_name
  save_dir_pw <- paste0(save_dir,sub_dir,"/")
  if (!file.exists(save_dir_pw)){
    dir.create(save_dir_pw)
  }
  cells_use <- c(cells_use_g1, cells_use_g2)
  # cells_use <- rownames(meta_data)
  cat("DE analysis ", file = paste0(save_dir_pw,"de_analysis_log.txt"))
  print(length(cells_use))
  cat(paste0("\n Observed clones: ",groups_use[1],' vs ',groups_use[2]," Nb obs cells: ",length(cells_use), 
             "  nb cells in ",groups_use[1], ":",length(cells_use_g1),
             "  nb cells in ",groups_use[2], ":",length(cells_use_g2)), file = paste0(save_dir_pw,"de_analysis_log.txt"), append = TRUE)
  
  
  # length(cells_use)
  srt2 <- srt[,as.character(cells_use)]  # subset cells
  dim(srt2)
  Idents(object = srt2) <- feature_use
  # unique(srt2$library_id)
  # print(Idents(object = srt2))
  # rownames(genes_map) <- genes_map$gene_ens
  # Differential Expression Analysis 
  
  markers_ls <- Seurat::FindMarkers(srt2, slot="data",  # using normalized data   
                                    ident.1=groups_use[1], ident.2=groups_use[2],    #Vector of cell names belonging to group 1, group 2
                                    logfc.threshold=minLogFC, test.use=test_use)
  
  print(paste0("Total DE marker genes: ",nrow(markers_ls)))
  cat(paste0("\n Total DE marker genes: ",nrow(markers_ls)), file = paste0(save_dir_pw,"de_analysis_log.txt"), append = TRUE)
  markers_ls$gene_symb <- as.character(genes_map[rownames(markers_ls),"gene_symb"])
  markers_ls$GENEID <- rownames(markers_ls)
  # View(head(markers_ls))
  # markers_ls$avg_log2FC
  # Get top marker genes
  markers_ls_output <- markers_ls[markers_ls$p_val_adj<pAdjustThrs & abs(markers_ls$avg_log2FC)>minLogFC,]
  markers_ls_output <- markers_ls_output[order(markers_ls_output$avg_log2FC,decreasing = T),]  
  print(paste0("Nb selected marker genes: ",nrow(markers_ls_output)))
  cat(paste0("\n Nb selected marker genes: ",nrow(markers_ls_output)), file = paste0(save_dir_pw,"de_analysis_log.txt"), append = TRUE)
  # Get top up-regulated genes
  markers_ls_upreg <- markers_ls_output[markers_ls_output$avg_log2FC>minLogFC,]
  markers_ls_upreg <- markers_ls_upreg[order(markers_ls_upreg$avg_log2FC,decreasing = T),] 
  
  
  # dim(markers_ls_upreg)
  if(nrow(markers_ls_upreg) < nbtopup){
    nbtopup <- nrow(markers_ls_upreg)
  }
  markers_ls_upreg_top <- markers_ls_upreg[1:nbtopup,]
  
  # Get top down-regulated genes
  markers_ls_downreg <- markers_ls_output[markers_ls_output$avg_log2FC<(-minLogFC),]
  markers_ls_downreg <- markers_ls_downreg[order(markers_ls_downreg$avg_log2FC,decreasing = F),] 
  # dim(markers_ls_downreg)
  if(nrow(markers_ls_downreg) < nbtopdown){
    nbtopdown <- nrow(markers_ls_downreg)
  }
  markers_ls_downreg_top <- markers_ls_downreg[1:nbtopdown,]
  
  markers_ls_upreg_top$genes_deg <- 'up-regulated'
  markers_ls_downreg_top$genes_deg <- 'down-regulated'
  genes_deg_df <- rbind(markers_ls_upreg_top, markers_ls_downreg_top)
  
  if(save_data){
    # write.table(markers_ls_output,paste0(save_dir_pw,"de_significant_genes.txt"),sep="\t",col.names=NA,quote=F)
    data.table::fwrite(markers_ls_output,paste0(save_dir_pw,"de_significant_genes.csv.gz"))
    data.table::fwrite(markers_ls,paste0(save_dir_pw,"total_markers.csv.gz"))
    # write.table(markers_ls,paste0(save_dir_pw,"total_markers.txt"),sep="\t",col.names=NA,quote=F)
    # ,paste0(save_dir_pw,"de_significant_genes.txt")convert_de_csv(sub_dir, save_dir_pw)
    
    # write.table(genes_deg_df, paste0(save_dir_pw,"de_topup_",nbtopup,"_topdown_",nbtopdown,".txt"),sep="\t",col.names=NA,quote=F)
  }
  
  print("Plot DE genes")
  # plttitle <- paste0(base_name,"; ",groups_use[1]," vs. ",groups_use[2])
  plttitle <- base_name
  topGenes <- genes_deg_df$gene_symb
  # markers_ls_tmp <- markers_ls
  # rownames(markers_ls_tmp) <- markers_ls_tmp$gene_symb
  # plot_DE_genes(markers_ls_tmp, topGenes, nrow(markers_ls_output), 
  #               FCcutoff=minLogFC, pCutoff=pAdjustThrs, plttitle, save_dir_pw)
  # 
  # summary_genes_ls <- list(markers_ls=markers_ls,
  # markers_ls_output=markers_ls_output,
  # top_genes=genes_deg_df,
  # markers_ls_upreg=markers_ls_upreg,
  # markers_ls_downreg=markers_ls_downreg)
  # if(save_data){
  #   saveRDS(summary_genes_ls,file = paste0(save_dir_pw,'summary_genes_ls.rds'))
  # }
  
  p <- plot_DE_genes_ggplot(markers_ls, topGenes, capstr="", 
                            FDRcutoff=0.01, logFCcutoff=0.25, pValuecutoff=0.05,
                            plttitle, save_dir_pw, legendVisible=F,
                            iscaption=TRUE, legend_verbose='none', save_plot=TRUE, xl=NULL)
  # p  
  # Plot heatmap
  # norm_df <- GetAssayData(object = srt2, slot = "data")
  # ext <- rownames(genes_deg_df)[rownames(genes_deg_df) %in% rownames(norm_df)]
  # norm_df <- norm_df[ext,] # get correct order
  # rownames(norm_df) <- genes_deg_df[rownames(norm_df),'gene_symb']
  # clusters <- srt2@meta.data[colnames(norm_df), feature_use]
  # genes_cls <- genes_deg_df[ext,]$genes_deg
  # print("Plot heatmap")
  # p <- Heatmap(as.matrix(norm_df),show_column_names=FALSE, 
  #              cluster_rows=F,cluster_columns=F,
  #              column_split=clusters, row_split = genes_cls) 
  # plttitle_hm <- paste0(base_name,"_",groups_use[1],"_vs_",groups_use[2])
  # png(paste0(save_dir_pw,"heatmap_",plttitle_hm,".png"), height = 2*500, width=2*700,res = 2*72)
  # print(p)
  # dev.off()
  
  # Plot the top 20 reference genes cisplatin
  # gmt_dir <- paste0(input_dir,"biodatabase/")
  # cis_genes <- read.csv(paste0(gmt_dir,'cisplatin_20_genes.csv'), header = T, sep=',', check.names=F, stringsAsFactors=F)
  
  # markers_ls_output <- read.table(paste0(save_dir,'SA1035_E_D/de_significant_genes.txt'), header=T,sep='\t',row.names = 1)
  # head(markers_ls_output)
  
  
  # if(length(unique(markers_ls_output$gene_symb))!=nrow(markers_ls_output)){
  #   print("Duplicate gene symbol, remove duplicated genes from the mtx")
  #   deg_mtx <- markers_ls_output[!duplicated(markers_ls_output$gene_symb),]
  # } else{
  #   deg_mtx <- markers_ls_output
  # }
  # c_use <- markers_ls$gene_symb %in% cis_genes$gene
  # if(sum(c_use)>0){
  #   cat(paste0("\n Nb cells exist in top20 ref genes: ",sum(c_use)), file = paste0(save_dir_pw,"de_analysis_log.txt"), append = TRUE)
  #   norm_df <- GetAssayData(object = srt2, slot = "data")
  #   # dim(norm_df)
  #   norm_df <- norm_df[rownames(markers_ls),]
  #   norm_df <- norm_df[c_use,]
  #   rownames(norm_df) <- markers_ls[rownames(norm_df),'gene_symb']
  # 
  #   # norm_df <- norm_df[rownames(norm_df) %in% cis_genes$gene,]
  #   print(dim(norm_df))
  #   cell_clusters <- srt2@meta.data[colnames(norm_df), feature_use]
  # # # genes_cls <- rownames(norm_df) %in% cis_genes$gene
  # # # names(genes_cls) <- rownames(norm_df)
  #   # markers_ls$is_included_top20_ref <- markers_ls_output$gene_symb %in% cis_genes$gene
  #   
  #   # print("Plot heatmap")
  #   # p <- Heatmap(as.matrix(norm_df),show_column_names=F, show_row_names =T,
  #   #              cluster_rows=F,cluster_columns=F, column_split = cell_clusters)  #,row_split = genes_cls
  #   # ht <- nrow(norm_df)
  #   # png(paste0(save_dir_pw,"hm_ref_cisplatin_genes.png"), height = 2*40*ht+40, width=2*700,res = 2*72)
  #   # print(p)
  #   # dev.off()
  # } else{
  #   print("No gene match the reference 20 gene list")
  # }
  
  
  
  
  # Pathway analysis
  print("Pathway analysis....")
  # gmt_dir <- paste0(input_dir,"biodatabase/")
  gmt_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/pathway_set/'
  deg_stat <- markers_ls_output$avg_log2FC
  names(deg_stat) <- markers_ls_output$gene_symb
  # # print(head(deg_stat))
  
  # comment back
  method_use <- paste0("seurat_DEG_",test_use)
  gmt_ls <- c("h.all.v7.0.symbols.gmt") #,"c2.cp.kegg.v7.1.symbols.gmt","GO_c5.all.v7.1.symbols.gmt"
  # pathway_names <- c("hallmark","kegg") #,"go"
  pathway_names <- c("hallmark")
  
  gmtfile <- paste0(gmt_dir, gmt_ls[1])
  
  # de_genes <- markers_ls_upreg
  # if(viz){
  #   viz_pathway(gmt_dir, markers_ls_upreg, save_dir_pw, base_name, tag='up_pathway')
  #   viz_pathway(gmt_dir, markers_ls_downreg, save_dir_pw, base_name, tag='down_pathway')
  # }
  
  print("Get pathway")
  # if(!viz){
  # Deal with overlapping and NA gene symbols
  markers_ls <- markers_ls %>%
    dplyr::select(gene_symb, avg_log2FC) %>%
    na.omit() %>%
    distinct() %>%
    group_by(gene_symb) %>%
    summarize(avg_log2FC=mean(avg_log2FC))
  
  for(i in rep(1:length(gmt_ls),1)){
    # gmt_fn <- paste0(gmt_dir, gmt_ls[i])
    print(gmt_ls[i])
    cat(paste0("\n Pathway analysis, and using gene sets: ",gmt_ls[i]," pathway is: ",pathway_names[i]), file = paste0(save_dir_pw,"de_analysis_log.txt"), append = TRUE)
    print(pathway_names[i])
    markers_ls$avg_logFC <- markers_ls$avg_log2FC
    pathway_ls <- get_pathway_results(markers_ls, method_use,
                                      base_name, paste0(gmt_dir, gmt_ls[i]), 
                                      pathway_names[i],
                                      groups_use,
                                      save_dir_pw, 30)
  }
  
  # # return(pathway_ls)
  # }
  # print("Generate genes network")
  # generate_genes_network(save_dir_pw,
  #                        datatag=sub_dir,
  #                        de_desc=sub_dir,
  #                        genes_set='hallmark',
  #                        pattern_use='^hallmark_'
  # )
}
