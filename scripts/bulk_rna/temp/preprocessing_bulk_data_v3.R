suppressPackageStartupMessages({
  library(tidyverse)
  library(annotables)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  # library(scales)
  # library(ggrepel)
  library(stringr)
  # library(scran)
  library(SingleCellExperiment)
  # library(tximport)
  
})


script_dir <- '/home/htran/Projects/hakwoo_project/metastasis_material/scripts/bulk_rna/'
source(paste0(script_dir, 'bulk_utils.R'))

normalizing_cloneA <- function(){
  input_dir <- '/home/htran/storage/rnaseq_datasets/bulk_metastasis/mixing_exp_SA919_bulkRNAseq/'  
  datatag <- 'SA919_cloneA'
  save_dir <- '/home/htran/storage/datasets/metastasis_results/bulk_SA919/mixing_SA919/'
  bulk_sids_A_main <- c('SA919X4XB40503','SA919X4XB09563') #main_A_met, main_A_primary
  bulk_sids_A_mixingM42 <- c('AT24180')
  
  df_counts_A_main <- load_raw_counts_kallisto(input_dir, datatag, save_dir, sample_ids=bulk_sids_A_main)
  
  df_counts_A_mixingM42 <- load_raw_counts_kallisto(input_dir, datatag, save_dir, sample_ids=bulk_sids_A_mixingM42)
  class(df_counts_A_mixingM42)
  df_counts_A_mixingM42 <- df_counts_A_mixingM42 %>%
    dplyr::filter(AT24180>0)
  
  df_counts_A_mixingM42 <- df_counts_A_mixingM42 %>%
    tidyr::pivot_longer(!ens_gene_id, names_to = 'sample_id', values_to = 'raw_counts')
  df_counts_A_main <- df_counts_A_main %>%
    tidyr::pivot_longer(!ens_gene_id, names_to = 'sample_id', values_to = 'raw_counts')
  head(df_counts_A_mixingM42)
  head(df_counts_A_main)
  
  df <- dplyr::bind_rows(df_counts_A_main, df_counts_A_mixingM42)
  dim(df)
  df <- df %>%
    tidyr::pivot_wider(values_from = 'raw_counts', names_from = 'sample_id')
  head(df)
  df1 <- df %>%
    na.omit()
  dim(df1)
  
  colnames(df_counts) #bulk_sids_A_42
  ## Loading counts - output of above function, and normalize data
  datatag <- 'SA919_cloneA'
  df_counts_fn <- paste0(save_dir, datatag,'_total_raw_counts.csv.gz')
  data.table::fwrite(df1, df_counts_fn)
  df_normalized <- normalize_by_size_factor(df_counts_fn, datatag, save_dir)
  
  ## To Do
  ## Get DE genes for 2 comparisons, and check the values here
  c1 <- '42_P_A_vs_origP_A_DESeq_All.csv'
  c1 <- '42_P_A_vs_origM_A_DESeq_All.csv'
  tmp_de <- data.table::fread(paste0(de_results_dir, c1))
  
  ## Only up-regulated genes
  tmp_de <- tmp_de %>%
    dplyr::select(all_of(c('g_ID','Gene.name','log2FoldChange'))) %>%
    dplyr::filter(log2FoldChange>=0.5) %>%
    dplyr::mutate(log2FoldChange=round(as.numeric(log2FoldChange), 2),
                  comp_DE=c1)
  # tmp_de$g_ID[1:3]
  dim(tmp_de)
  df_normalized_c1 <- df_normalized %>%
    dplyr::select(SA919X4XB09563, AT24180, ens_gene_id) %>%
    dplyr::filter(ens_gene_id %in% tmp_de$g_ID)
  # TP: 17/23, FP: 6/23 genes
  
  df_normalized_c1 <- df_normalized %>%
    dplyr::select(SA919X4XB40503, AT24180, ens_gene_id) %>%
    dplyr::filter(ens_gene_id %in% tmp_de$g_ID)
  
  dim(df_normalized_c1)
  # View(df_normalized_c1)
  sum(df_normalized_c1$AT24180>df_normalized_c1$SA919X4XB40503)
  # TP: 129/222, FP: 93/222 genes
  
  cor.test(df_normalized$AT24180, df_normalized$SA919X4XB09563, method="pearson")
  cor.test(df_normalized$AT24180, df_normalized$SA919X4XB40503, method="pearson")
  
  ## To Do: redo analysis here using DESeq2, 
  ## Is this possible to use only 2 samples for comparisons? 
}
main <- function(){
  input_dir <- '/home/htran/storage/rnaseq_datasets/bulk_metastasis/mixing_exp_SA919_bulkRNAseq/'  
  datatag <- 'SA919_mixing'
  save_dir <- '/home/htran/storage/datasets/metastasis_results/bulk_SA919/mixing_SA919/'
  df_counts <- load_raw_counts_kallisto(input_dir, datatag, save_dir, sample_ids=NULL)
  
  ## Loading counts - output of above function, and normalize data
  df_counts_fn <- paste0(save_dir, datatag,'_total_raw_counts.csv.gz')
  df_normalized <- normalize_by_size_factor(df_counts_fn, datatag, save_dir)
  
  ## Meta data
  save_dir <- '/home/htran/storage/datasets/metastasis_results/bulk_SA919/mixing_SA919/'
  meta_samples_fn <- paste0(save_dir,'Samples-Metastasis_Hakwoo_bulkRNA_mixing_exp.csv')
  get_meta_samples_SA919_mixing_exp(meta_samples_fn, datatag, save_dir)
  
  ## Normalize data for GS
  # Noted: in total of 20738 genes, there are 1252 genes with NA values in mixing exp samples
  # but not NA in main exp. The reason maybe due to reference mapping from transcript ids
  # into gene id. I will take a look at this later after tmr meeting but you can use this file for a draft calculation.
  datatag <- 'Mixed_Ex3_GS'
  df_counts_fn <- paste0(save_dir, datatag, '_raw_counts.csv.gz')
  df_normalized <- normalize_by_size_factor(df_counts_fn, datatag, save_dir)
  
}

# main()


process_metadata <- function(){
  save_dir <- '/home/htran/storage/datasets/metastasis_results/bulk_SA919/mixing_SA919/'
  datatag <- 'Mixed_Ex3_GS'
  
  meta_samples_fn <- paste0(save_dir,'Samples-Metastasis_Hakwoo_bulkRNA_mixing_exp.csv')
  metasamples <- data.table::fread(meta_samples_fn)
  
  # head(metasamples)
  colnames(metasamples) <- gsub(' ','_', colnames(metasamples))
  metasamples$Bulk_RNA_ATID
  # View(metasamples)
  metasamples <- metasamples %>%
    dplyr::filter(Bulk_RNA_ATID!='')
  dim(metasamples) ## 13 samples
  # View(metasamples)
  metasamples$library_id <- sapply(strsplit(metasamples$`Library_ID(s)`,','), function(x){
    return(as.character(x[1]))
  })
  # sum(metasamples$library_id %in% mixing_meta_df$library_id)
  # metasamples$`Library_ID(s)`
  mixing_meta_df <- data.table::fread('/home/htran/Projects/hakwoo_project/metastasis_material/materials/dlp_trees/SA919_mixing_experiment/Metastasis_Hakwoo_mixing_exp_SA919_results.csv')
  dim(mixing_meta_df)
  # mixing_ids_df <- data.table::fread('/home/htran/Projects/hakwoo_project/metastasis_material/materials/dlp_trees/SA919_mixing_experiment/Metastasis_Hakwoo_mixing_exp_SA919_libraryIds.csv')
  # # View(mixing_meta_df)
  # colnames(mixing_meta_df)
  mixing_meta_df <- mixing_meta_df %>%
    dplyr::select(library_id, main_clone, mainsite)
  metasamples <- metasamples %>%
    dplyr::left_join(mixing_meta_df, by='library_id')
  dim(metasamples)
  metasamples <- as.data.frame(metasamples)
  rownames(metasamples) <- metasamples$Bulk_RNA_ATID
  
  # View(metasamples)
  norm_df <- normalized_df %>%
    tibble::column_to_rownames('ens_gene_id')
  cols_use <- colnames(norm_df)[colnames(norm_df) %in% metasamples$Bulk_RNA_ATID]
  norm_df <- norm_df[ ,cols_use]
  dim(norm_df)
  metasamples$mouse_id <- paste0('M',stringr::str_sub(metasamples$Sample_Name, 4,5))
  colnames(metasamples)
  
  metasamples <- metasamples %>% 
    dplyr::rename(pdxid=SA_ID, origin=Anatomical_Site, 
                  sample_id=AT_ID, bulk_sid=Bulk_RNA_ATID) %>%
    dplyr::mutate(experiment='mixing_exp')
  
  main_meta_df <- data.table::fread('/home/htran/Projects/hakwoo_project/metastasis_material/materials/bulkRNAseq/SA919/library_groupings_bulk_SA919_cloneIds.csv')
  dim(main_meta_df)
  colnames(main_meta_df)
  main_meta_df <- main_meta_df %>% 
    dplyr::rename(main_clone=clone_id) %>% 
    dplyr::select(-nb_cells)  %>%
    dplyr::mutate(experiment='main_exp')
  
  total_meta <- dplyr::bind_rows(main_meta_df, metasamples)
  dim(total_meta)
  # View(total_meta)
  
  data.table::fwrite(total_meta, paste0(save_dir, 'metadata_Hakwoo_bulkRNA_mixing_main_exp.csv'))
  total_meta <- data.table::fread(paste0(save_dir, 'metadata_Hakwoo_bulkRNA_mixing_main_exp.csv'))
  
  total_meta$transplanted_mouse_id <- ifelse(total_meta$experiment=='main_exp','',
                                             stringr::str_sub(total_meta$Sample_Name,4,5))
  data.table::fwrite(total_meta, paste0(save_dir, 'metadata_Hakwoo_bulkRNA_mixing_main_exp.csv'))
  script_dir <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/bulkRNAseq/'
  data.table::fwrite(total_meta, paste0(script_dir, 'metadata_Hakwoo_bulkRNA_mixing_main_exp.csv'))
  
}

process_metadata_comparisons <- function(){
  save_dir <- '/home/htran/storage/datasets/metastasis_results/bulk_SA919/mixing_SA919/'
  script_dir <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/bulkRNAseq/'
  datatag <- 'Mixed_Ex3_GS'
  total_meta <- data.table::fread(paste0(save_dir, 'metadata_Hakwoo_bulkRNA_mixing_main_exp.csv'))
  dim(total_meta)
  total_meta$comp_DE <- paste0(total_meta$transplanted_mouse_id,'_',
                              stringr::str_sub(total_meta$mainsite, 1,1),'_',
                              total_meta$main_clone,'_vs_','orig')
  colnames(total_meta)
  t <- total_meta %>%
    dplyr::filter(experiment=='main_exp' & main_clone=='C')

  
  comp_C <- total_meta %>%
    dplyr::filter(experiment=='mixing_exp' & main_clone=='C')
  comp_C <- comp_C %>%
    dplyr::mutate(versus_group2='SA919X7XB05604', group2_origin='SupraSpinal', 
                  group2_main_clone='C', comp_DE=paste0(comp_DE,'M','_C','_DESeq_All.csv'))
  # comp_C$comp_DE
  dim(comp_C)  
  ## Just checking about naming DE results files 
  de_results_dir <- paste0(script_dir, 'Bulk_RNA_Seq_GS/')
  # for(c in comp_C$comp_DE){
  #   if(file.exists(paste0(de_results_dir, c))){
  #     print('ok')
  #   }else{
  #     print(c)
  #   }
  # }
  comp_A <- total_meta %>%
    dplyr::filter(experiment=='mixing_exp' & main_clone=='A')
  # t <- total_meta %>%
  #   dplyr::filter(experiment=='main_exp' & main_clone=='A')
  # View(t)
  comp_A1 <- comp_A %>%
    dplyr::mutate(versus_group2='SA919X4XB09563', group2_origin='Primary', 
                  group2_main_clone='A', comp_DE=paste0(comp_DE,'P','_A','_DESeq_All.csv'))
  
  comp_A2 <- comp_A %>%
    dplyr::mutate(versus_group2='SA919X4XB40503', group2_origin='Axillary', 
                  group2_main_clone='A', comp_DE=paste0(comp_DE,'M','_A','_DESeq_All.csv'))
  comp_A <- dplyr::bind_rows(comp_A1, comp_A2)
  dim(comp_A1)
  dim(comp_A2)
  dim(comp_A)
  # for(c in comp_A$comp_DE){
  #   if(file.exists(paste0(de_results_dir, c))){
  #     print('ok')
  #     # print(c)
  #   }else{
  #     print('!!')
  #     print(c)
  #   }
  # }
  comp_all <- dplyr::bind_rows(comp_C, comp_A1, comp_A2)
  data.table::fwrite(comp_all, paste0(save_dir, 'list_DE_comparisons.csv'))
  script_dir <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/bulkRNAseq/'
  data.table::fwrite(comp_all, paste0(script_dir, 'list_DE_comparisons.csv'))
  
}



load_DE_genes_v2 <- function(){
  script_dir <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/bulkRNAseq/'
  de_results_dir <- paste0(script_dir, 'Bulk_RNA_Seq_GS/')
  comp_all <- data.table::fread(paste0(script_dir, 'list_DE_comparisons.csv'))
  comp_C <- comp_all %>%
    dplyr::filter(group2_main_clone=='C' & transplanted_mouse_id=='10')
  dim(comp_C)
  
  # View(comp_C)
  comp_used <- unique(comp_C$comp_DE)
  # comp_used <- comp_used[!comp_used %in% c("45_M_C_vs_origM_C_DESeq_All.csv", "42_M_CC1_vs_origM_C_DESeq_All.csv")]
  comp_used
  de_all <- tibble::tibble()
  nb_DE_genes_max <- 0
  for(c in comp_used){
    if(file.exists(paste0(de_results_dir, c))){
      print('Read file')
      print(c)
      tmp_de <- data.table::fread(paste0(de_results_dir, c))
      # print(colnames(tmp_de))
      # sum(as.numeric(tmp_de$log2FoldChange))
      
      tmp_de <- tmp_de %>%
        dplyr::select(all_of(c('g_ID','Gene.name','log2FoldChange'))) %>%
        dplyr::mutate(log2FoldChange=round(as.numeric(log2FoldChange), 2),
                      comp_DE=c) %>%
        dplyr::filter(log2FoldChange>=0.5)
        # dplyr::filter(abs(log2FoldChange)>=1)%>%
        # dplyr::mutate(
        #   logFC_direction=case_when(
        #     log2FoldChange>0 ~ 1,
        #     log2FoldChange<0 ~ -1,
        #     TRUE ~ 0
        #   )
        # )
      if(dim(tmp_de)[1]>nb_DE_genes_max){
        nb_DE_genes_max <- dim(tmp_de)[1]
      }
      print(dim(tmp_de))
      # head(tmp_de)
      de_all <- dplyr::bind_rows(de_all, tmp_de)
      
    }else{
      print('Issue with read file, double check!!!')
      print(c)
    }
  }
  
  ## Get statistical values
  stat <- tibble::tibble()
  for(c in comp_used){
    tmp_de <- de_all %>%
      dplyr::filter(comp_DE==c)
    t <- tibble::tibble(pct_regulated_genes=round(dim(tmp_de)[1]/nb_DE_genes_max, 2),
                        comp_DE=c,
                        avg_logFC=round(mean(abs(tmp_de$log2FoldChange)), 2))
    stat <- dplyr::bind_rows(stat, t)
  }  
  stat
  dim(comp_C)
  stat <- stat %>% 
    left_join(comp_C, by='comp_DE')
  dim(stat)
  data.table::fwrite(stat, paste0(save_dir, 'SA919_mixing_DE_analysis_stat_cloneC.csv.gz'))
  
  
  ## For clone A
  comp_A <- comp_all %>%
    dplyr::filter(group2_main_clone=='A' & transplanted_mouse_id=='45')
  comp_A <- comp_all %>%
    dplyr::filter(group2_main_clone=='A' & transplanted_mouse_id=='33')
  comp_A <- comp_all %>%
    dplyr::filter(group2_main_clone=='A' & transplanted_mouse_id=='34' & origin=='Primary')
  comp_A <- comp_all %>%
    dplyr::filter(group2_main_clone=='A' & transplanted_mouse_id=='34' & origin!='Primary')
  dim(comp_A)
  View(comp_A)
  comp_used <- unique(comp_A$comp_DE)
  # comp_used <- comp_used[!comp_used %in% c("45_M_C_vs_origM_C_DESeq_All.csv", "42_M_CC1_vs_origM_C_DESeq_All.csv")]
  comp_used
  
  de_all <- tibble::tibble()
  # nb_DE_genes <- 0
  for(c in comp_used){
    if(file.exists(paste0(de_results_dir, c))){
      print('Read file')
      print(c)
      tmp_de <- data.table::fread(paste0(de_results_dir, c))
      # print(colnames(tmp_de))
      # sum(as.numeric(tmp_de$log2FoldChange))
      
      tmp_de <- tmp_de %>%
        dplyr::select(all_of(c('g_ID','Gene.name','log2FoldChange'))) %>%
        dplyr::mutate(log2FoldChange=round(as.numeric(log2FoldChange), 2),
                      comp_DE=c) %>%
        dplyr::filter(abs(log2FoldChange)>=0.5)
        # dplyr::filter(abs(log2FoldChange)>=1)%>%
        # dplyr::mutate(
        #   logFC_direction=case_when(
        #     log2FoldChange>0 ~ 1,
        #     log2FoldChange<0 ~ -1,
        #     TRUE ~ 0
        #   )
        # )
      # nb_DE_genes <- nb_DE_genes + dim(tmp_de)[1]
      # if(dim(tmp_de)[1]>nb_DE_genes_max){
      #   nb_DE_genes_max <- dim(tmp_de)[1]
      # }
      print(dim(tmp_de))
      # head(tmp_de)
      de_all <- dplyr::bind_rows(de_all, tmp_de)
      
    }else{
      print('Issue with read file, double check!!!')
      print(c)
    }
  }
  
  nb_DE_genes <- dim(de_all)[1]
  nb_DE_genes
  ## Get statistical values
  stat <- tibble::tibble()
  for(c in comp_used){
    print(c)
    tmp_de <- de_all %>%
      dplyr::filter(comp_DE==c)
    t <- tibble::tibble(pct_regulated_genes=round(dim(tmp_de)[1]/nb_DE_genes,2),
                        comp_DE=c,
                        avg_logFC=round(mean(abs(tmp_de$log2FoldChange)),2))
    stat <- dplyr::bind_rows(stat, t)
  }  
  stat
  # dim(comp_A)
  stat <- stat %>% 
    left_join(comp_A, by='comp_DE')
  dim(stat)
  data.table::fwrite(stat, paste0(save_dir, 'SA919_mixing_DE_analysis_stat_cloneA.csv.gz'))
  
  
  genes_used <- unique(de_all$g_ID)
  length(genes_used)  
  unique(de_all$comp_DE)
  de_all1 <- de_all %>%
    dplyr::select(all_of(c('g_ID','logFC_direction','comp_DE'))) %>%
    dplyr::mutate(comp_DE=gsub('_DESeq_All.csv','',comp_DE))%>%
    dplyr::group_by(g_ID, comp_DE) %>%
    dplyr::summarise(logFC_direction=mean(logFC_direction))%>%
    # dplyr::rename(gene_symbol=`Gene.name`)%>%
    tidyr::pivot_wider(names_from = 'comp_DE', values_from = 'logFC_direction', values_fill = 0) %>%
    tibble::column_to_rownames('g_ID')
  
  dim(de_all1)    
  head(de_all1)
  # de_all1 <- de_all1[,1:4]
  p <- ComplexHeatmap::Heatmap(as.matrix(de_all1), na_col = "white",
                               # col = col_fun,
                               show_column_names=T,
                               show_row_names = F,
                               # cluster_rows=T,
                               # cluster_columns=T,
                               # clustering_distance_columns = "pearson",
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
                               # top_annotation=top_anno,
                               # left_annotation = left_anno,
                               # cell_fun = cell_func,
                               row_dend_reorder=T
  )
  # p
  
}

get_DE_genes_cloneA <- function(){
  script_dir <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/bulkRNAseq/'
  de_results_dir <- paste0(script_dir, 'Bulk_RNA_Seq_GS/')
  comp_all <- data.table::fread(paste0(script_dir, 'list_DE_comparisons.csv'))
  ## For clone A
  comp_A <- comp_all %>%
    dplyr::filter(group2_main_clone=='A' & transplanted_mouse_id=='42')
  # comp_A <- comp_all %>%
  #   dplyr::filter(group2_main_clone=='A' & transplanted_mouse_id=='45')
  
  comp_used <- unique(comp_A$comp_DE)
  # comp_used <- comp_used[!comp_used %in% c("45_M_C_vs_origM_C_DESeq_All.csv", "42_M_CC1_vs_origM_C_DESeq_All.csv")]
  comp_used
  # c <- comp_used[1]
  # tmp_de <- data.table::fread(paste0(de_results_dir, c))
  
  de_all <- tibble::tibble()
  # nb_DE_genes <- 0
  for(c in comp_used){
    if(file.exists(paste0(de_results_dir, c))){
      print('Read file')
      print(c)
      tmp_de <- data.table::fread(paste0(de_results_dir, c))
      # print(colnames(tmp_de))
      # sum(as.numeric(tmp_de$log2FoldChange))
      
      # tmp_de <- tmp_de %>%
      #   dplyr::select(all_of(c('g_ID','Gene.name','log2FoldChange'))) %>%
      #   dplyr::mutate(log2FoldChange=round(as.numeric(log2FoldChange), 2),
      #                 comp_DE=c) %>%
      #   dplyr::filter(abs(log2FoldChange)>=0.5)
      tmp_de <- tmp_de %>%
        dplyr::select(all_of(c('g_ID','Gene.name','log2FoldChange'))) %>%
        dplyr::filter(log2FoldChange>=1) %>%
        dplyr::mutate(log2FoldChange=round(as.numeric(log2FoldChange), 2),
                      comp_DE=c)
      # dplyr::filter(abs(log2FoldChange)>=1)%>%
      # dplyr::mutate(
      #   logFC_direction=case_when(
      #     log2FoldChange>0 ~ 1,
      #     log2FoldChange<0 ~ -1,
      #     TRUE ~ 0
      #   )
      # )
      # nb_DE_genes <- nb_DE_genes + dim(tmp_de)[1]
      # if(dim(tmp_de)[1]>nb_DE_genes_max){
      #   nb_DE_genes_max <- dim(tmp_de)[1]
      # }
      print(dim(tmp_de))
      # head(tmp_de)
      de_all <- dplyr::bind_rows(de_all, tmp_de)
      
    }else{
      print('Issue with read file, double check!!!')
      print(c)
    }
  }
  de_all$Gene.name[1]
  de_all <- de_all %>%
    dplyr::filter(!grepl('MT-', Gene.name))
  nb_DE_genes <- dim(de_all)[1]
  nb_DE_genes
  genes_used <- unique(de_all$g_ID)
  genes_used[1:5]
  length(genes_used)
}

get_hierachical_clustering_cloneA <- function(){
  save_dir <- '/home/htran/storage/datasets/metastasis_results/bulk_SA919/mixing_SA919/'
  datatag <- 'Mixed_Ex3_GS'
  normalized_df <- data.table::fread(paste0(save_dir, datatag, '_sizefactor_normalized.csv.gz'))
  dim(normalized_df)  
  colnames(normalized_df)
  # t <- normalized_df
  # t$ens_gene_id <- NULL
  # colSums(t)
  ## See process_metadata() above
  total_meta <- data.table::fread(paste0(save_dir, 'metadata_Hakwoo_bulkRNA_mixing_main_exp.csv'))
  # total_meta <- total_meta %>%
  #   dplyr::filter(main_clone=='C')
  # View(total_meta)
  
  
  fpkm_df <- get_normalized_data(normalized_df)
  fpkm_df <- data.table::fread(paste0(save_dir, 'SA919_mixing_normalized_fpkm.csv.gz'))
  dim(fpkm_df)
  colnames(fpkm_df)
  normalized_df <- fpkm_df
  
  main_A_primary <- total_meta %>%
    dplyr::filter(main_clone=='A' & experiment=='main_exp' & origin=='Primary' & pdxid=='X08472164') %>%
    dplyr::pull(bulk_sid)
  main_A_met <- total_meta %>%
    dplyr::filter(main_clone=='A' & experiment=='main_exp' & origin=='Axillary' & pdxid=='X08472164') %>%
    dplyr::pull(bulk_sid)
  
  
  script_dir <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/bulkRNAseq/'
  de_results_dir <- paste0(script_dir, 'Bulk_RNA_Seq_GS/')
  comp_all <- data.table::fread(paste0(script_dir, 'list_DE_comparisons.csv'))
  comp_A_42 <- comp_all %>%
    dplyr::filter(group2_main_clone=='A' & mainsite=='Primary' & transplanted_mouse_id =='42')
  # comp_A_45 <- comp_all %>%
  #   dplyr::filter(group2_main_clone=='A' & mainsite=='Primary' & transplanted_mouse_id =='45')
  # dim(comp_A)
  # colnames(comp_A)
  # View(comp_A)
  
  bulk_sids_A_42 <- unique(c(main_A_met, main_A_primary, comp_A_42$bulk_sid))
  # bulk_sids_A_45 <- unique(c(main_A_met, main_A_primary, comp_A_45$bulk_sid))
  
  # norm_A <- normalized_df %>%
  #   dplyr::select(all_of(c(bulk_sids_A,'ens_gene_id'))) %>%
  #   dplyr::filter(ens_gene_id %in% genes_used) %>%
  #   tibble::column_to_rownames('ens_gene_id') %>%
  #   as.data.frame()
  # norm_A <- norm_A[genes_used, ]
  # 
  # fpkm_A <- fpkm_df %>%
  #   dplyr::select(all_of(c(bulk_sids_A,'ens_gene_id'))) %>%
  #   dplyr::filter(ens_gene_id %in% genes_used) %>%
  #   tibble::column_to_rownames('ens_gene_id')  %>%
  #   as.data.frame()
  # fpkm_A <- fpkm_A[genes_used, ]
  
  norm_A <- normalized_df %>%
      dplyr::select(all_of(c(bulk_sids_A_42,'ens_gene_id')))
  # norm_A <- normalized_df %>%
  #   dplyr::select(all_of(c(bulk_sids_A_45,'ens_gene_id')))
  dim(norm_A)  
  colnames(norm_A)
  total_meta$bulk_sid
  metasamples <- total_meta
  # View(metasamples)
  colSums(norm_A[,1:3])
  
  norm_df <- norm_A
  norm_df <- as.data.frame(norm_df)
  rownames(norm_df) <- norm_df$ens_gene_id
  norm_df$ens_gene_id <- NULL
  rownames(norm_df)[1]
  # genes_used <- rownames(norm_df)[rownames(norm_df) %in% genes_used]
  genes_used2 <- intersect(rownames(norm_df), genes_used)
  length(genes_used2)
  norm_df <- norm_df[genes_used2,]
  dim(norm_df)
  
  # For only DE genes
  # For all genes --> similar results. 
  ## To Do: double check again DE analysis, using DESeq2
  cor.test(norm_df$AT24180, norm_df$SA919X4XB09563,method="pearson") #pri M42 and pri main A: 0.83
  cor.test(norm_df$AT24180, norm_df$SA919X4XB40503,method="pearson") #pri M42 and met main A: 0.93
  cor.test(norm_df$SA919X4XB09563, norm_df$SA919X4XB40503,method="pearson") #main exp: pri A and met A: 0.96
  
  
  metasamples <- metasamples %>%
    as.data.frame() %>%
    dplyr::filter(bulk_sid %in% colnames(norm_df))
  rownames(metasamples) <- metasamples$bulk_sid
  # metasamples$transplanted_mouse_id <- ifelse(metasamples$experiment=='main_exp',1,metasamples$transplanted_mouse_id)
  unique(metasamples$origin)
  predefined_cols <- c('#EE220C','#004D80')
  names(predefined_cols) <- c('Metastasis','Primary')
  # colors_use <- predefined_cols[metasamples[colnames(norm_df),'mainsite']]
  # cols_use=c("#66A61E", "#66A61E", "#666666") # green and grey
  # names(colors_use) <- colnames(norm_df)
  # unique(metasamples[colnames(norm_df),'main_clone'])
  clones_color <- c('A'='#66C2A5','B'='#FC8D62','C'='#8DA0CB')
  # colors_clone_use <- clones_color[metasamples[colnames(norm_df),'main_clone']]
  # cols_use=c("#66A61E", "#66A61E", "#666666") # green and grey
  # names(colors_clone_use) <- colnames(norm_df)
  mouse_cols <- brewer.pal(8, "Accent")[1:length(unique(metasamples$transplanted_mouse_id))]
  names(mouse_cols) <- unique(metasamples$transplanted_mouse_id)
  
  origin_cols <- brewer.pal(8, "Dark2")[1:length(unique(metasamples$origin))]
  names(origin_cols) <- unique(metasamples$origin)
  
  
  top_anno = ComplexHeatmap::HeatmapAnnotation(MainSite = factor(metasamples[colnames(norm_df),'mainsite']),
                                               MainClone=factor(metasamples[colnames(norm_df),'main_clone']),
                                               # MouseId=factor(metasamples[colnames(norm_df),'transplanted_mouse_id']),
                                               # Origin=factor(metasamples[colnames(norm_df),'origin']),
                                               col = list(MainSite=predefined_cols, 
                                                          MainClone=clones_color#,
                                                          # MouseId=mouse_cols#,
                                                          # Origin=origin_cols
                                               )
  ) #, MainClone=colors_use
  # library(ComplexHeatmap)
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
  p
  
  # colSums(norm_df)
  png(paste0(save_dir,datatag,"_hierarchial_clusters.png"), height = 2*800, width=2*850, res = 2*72)
  print(p)
  dev.off()  
  
}
load_DE_genes <- function(){
  script_dir <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/bulkRNAseq/'
  de_results_dir <- paste0(script_dir, 'Bulk_RNA_Seq_GS/')
  comp_all <- data.table::fread(paste0(script_dir, 'list_DE_comparisons.csv'))
  comp_C <- comp_all %>%
    dplyr::filter(group2_main_clone=='C')
  dim(comp_C)
  unique(comp_C$comp_DE)
  # View(comp_C)
  comp_used <- unique(comp_C$comp_DE)
  # comp_used <- comp_used[!comp_used %in% c("45_M_C_vs_origM_C_DESeq_All.csv", "42_M_CC1_vs_origM_C_DESeq_All.csv")]
  
  de_all <- tibble::tibble()
  nb_DE_genes_max <- 0
  for(c in comp_used){
    if(file.exists(paste0(de_results_dir, c))){
      print('Read file')
      print(c)
      tmp_de <- data.table::fread(paste0(de_results_dir, c))
      # print(colnames(tmp_de))
      # sum(as.numeric(tmp_de$log2FoldChange))
      
      tmp_de <- tmp_de %>%
        dplyr::select(all_of(c('g_ID','Gene.name','log2FoldChange'))) %>%
        dplyr::mutate(log2FoldChange=round(as.numeric(log2FoldChange), 2),
                      comp_DE=c) %>%
        dplyr::filter(abs(log2FoldChange)>=1)%>%
        dplyr::mutate(
          logFC_direction=case_when(
            log2FoldChange>0 ~ 1,
            log2FoldChange<0 ~ -1,
            TRUE ~ 0
          )
        )
      if(dim(tmp_de)[1]>nb_DE_genes_max){
        nb_DE_genes_max <- dim(tmp_de)[1]
      }
      print(dim(tmp_de))
      # head(tmp_de)
      de_all <- dplyr::bind_rows(de_all, tmp_de)
      
    }else{
      print('Issue with read file, double check!!!')
      print(c)
    }
  }
  dim(de_all)
  nb_DE_genes_max
  
  ## Get statistical values
  stat <- tibble::tibble()
  for(c in comp_used){
    tmp_de <- de_all %>%
      dplyr::filter(comp_DE==c)
    t <- tibble::tibble(pct_regulated_genes=round(100*dim(tmp_de)[1]/nb_DE_genes_max, 2),
                        comp_DE=c,
                        avg_logFC=round(mean(abs(tmp_de$log2FoldChange)), 2))
    stat <- dplyr::bind_rows(stat, t)
  }  
  # stat
  dim(comp_C)
  stat <- stat %>% 
    left_join(comp_C, by='comp_DE')
  dim(stat)
  data.table::fwrite(stat, paste0(save_dir, 'SA919_mixing_DE_analysis_stat_cloneC.csv.gz'))
  
  
  ## For clone A
  comp_A <- comp_all %>%
    dplyr::filter(group2_main_clone=='A')
  dim(comp_A)
  
  comp_used <- unique(comp_A$comp_DE)
  # comp_used <- comp_used[!comp_used %in% c("45_M_C_vs_origM_C_DESeq_All.csv", "42_M_CC1_vs_origM_C_DESeq_All.csv")]
  
  
  de_all <- tibble::tibble()
  nb_DE_genes_max <- 0
  for(c in comp_used){
    if(file.exists(paste0(de_results_dir, c))){
      print('Read file')
      print(c)
      tmp_de <- data.table::fread(paste0(de_results_dir, c))
      # print(colnames(tmp_de))
      # sum(as.numeric(tmp_de$log2FoldChange))
      
      tmp_de <- tmp_de %>%
        dplyr::select(all_of(c('g_ID','Gene.name','log2FoldChange'))) %>%
        dplyr::mutate(log2FoldChange=round(as.numeric(log2FoldChange), 2),
                      comp_DE=c) %>%
        dplyr::filter(abs(log2FoldChange)>=1)%>%
        dplyr::mutate(
          logFC_direction=case_when(
            log2FoldChange>0 ~ 1,
            log2FoldChange<0 ~ -1,
            TRUE ~ 0
          )
        )
      if(dim(tmp_de)[1]>nb_DE_genes_max){
        nb_DE_genes_max <- dim(tmp_de)[1]
      }
      print(dim(tmp_de))
      # head(tmp_de)
      de_all <- dplyr::bind_rows(de_all, tmp_de)
      
    }else{
      print('Issue with read file, double check!!!')
      print(c)
    }
  }
  dim(de_all)
  nb_DE_genes_max
  
  ## Get statistical values
  stat <- tibble::tibble()
  for(c in comp_used){
    print(c)
    tmp_de <- de_all %>%
      dplyr::filter(comp_DE==c)
    t <- tibble::tibble(pct_regulated_genes=round(100*dim(tmp_de)[1]/nb_DE_genes_max,2),
                        comp_DE=c,
                        avg_logFC=round(mean(abs(tmp_de$log2FoldChange)),2))
    stat <- dplyr::bind_rows(stat, t)
  }  
  # stat
  dim(comp_A)
  stat <- stat %>% 
    left_join(comp_A, by='comp_DE')
  dim(stat)
  data.table::fwrite(stat, paste0(save_dir, 'SA919_mixing_DE_analysis_stat_cloneA.csv.gz'))
  
  
  genes_used <- unique(de_all$g_ID)
  length(genes_used)  
  unique(de_all$comp_DE)
  de_all1 <- de_all %>%
    dplyr::select(all_of(c('g_ID','logFC_direction','comp_DE'))) %>%
    dplyr::mutate(comp_DE=gsub('_DESeq_All.csv','',comp_DE))%>%
    dplyr::group_by(g_ID, comp_DE) %>%
    dplyr::summarise(logFC_direction=mean(logFC_direction))%>%
    # dplyr::rename(gene_symbol=`Gene.name`)%>%
    tidyr::pivot_wider(names_from = 'comp_DE', values_from = 'logFC_direction', values_fill = 0) %>%
    tibble::column_to_rownames('g_ID')
    
  dim(de_all1)    
  head(de_all1)
  # de_all1 <- de_all1[,1:4]
  p <- ComplexHeatmap::Heatmap(as.matrix(de_all1), na_col = "white",
                               # col = col_fun,
                               show_column_names=T,
                               show_row_names = F,
                               # cluster_rows=T,
                               # cluster_columns=T,
                               # clustering_distance_columns = "pearson",
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
                               # top_annotation=top_anno,
                               # left_annotation = left_anno,
                               # cell_fun = cell_func,
                               row_dend_reorder=T
  )
  p
  
}

get_normalized_data <- function(norm_df){
  ## Get gene length and normalize data
  save_dir <- '/home/htran/storage/datasets/metastasis_results/bulk_SA919/mixing_SA919/'
  genes_length <- data.table::fread(paste0(save_dir, 'SA919_mixing_total_raw_length.csv.gz'))
  dim(genes_length)
  head(genes_length)
  colnames(genes_length) <- c('gene_length','ens_gene_id_with_version')
  genes_length$ens_gene_id <- get_geneId_without_version(genes_length$ens_gene_id_with_version)
  # sum(genes_length$ens_gene_id %in% rownames(norm_df))
  genes_length <- genes_length %>%
    dplyr::filter(ens_gene_id %in% norm_df$ens_gene_id) %>%
    dplyr::select(-ens_gene_id_with_version)
    # tibble::column_to_rownames('ens_gene_id')
  dim(genes_length)
  dim(norm_df)
  norm_df <- norm_df %>%
    dplyr::left_join(genes_length, by='ens_gene_id') %>%
    dplyr::filter(!is.na(gene_length))
  colnames(norm_df)
  norm_df <- as.data.frame(norm_df)
  fpkm_df <- tibble::tibble()
  used_sids <- colnames(norm_df)[!colnames(norm_df) %in% c('ens_gene_id','gene_length')]
  for(s in used_sids){
    fpkm_val <- get_normalized_FPKM(norm_df[,s], norm_df$gene_length)
    print(length(fpkm_val))
    tmp <- tibble::tibble(bulk_sid=s, fpkm=fpkm_val, ens_gene_id=norm_df$ens_gene_id)
    fpkm_df <- dplyr::bind_rows(fpkm_df, tmp)
  }
  # head(fpkm_df)
  fpkm_df <- fpkm_df %>%
    tidyr::pivot_wider(names_from = 'bulk_sid', values_from = 'fpkm')
  # head(fpkm_df)
  data.table::fwrite(fpkm_df, paste0(save_dir, 'SA919_mixing_normalized_fpkm.csv.gz'))
  # length(norm_df$SA919X4XB40503)
  # length(norm_df$gene_length)
  # norm_df$tpm_SA919X4XB40503 <- get_normalized_TPM(norm_df$SA919X4XB40503, norm_df$gene_length)
  # norm_df$tpm_SA919X4XB09563 <- get_normalized_TPM(norm_df$SA919X4XB09563, norm_df$gene_length)
  # norm_df$tpm_AT24180 <- get_normalized_TPM(norm_df$AT24180, norm_df$gene_length)
  # norm_df$fpkm_AT24180 <- get_normalized_FPKM(norm_df$AT24180, norm_df$gene_length)
  # norm_df$fpkm_SA919X4XB40503 <- get_normalized_FPKM(norm_df$SA919X4XB40503, norm_df$gene_length)
  # norm_df$fpkm_SA919X4XB09563 <- get_normalized_FPKM(norm_df$SA919X4XB09563, norm_df$gene_length)
}
## To Do:   
## Get normalized data, convert it to cpm/tpm/fpkm format if needed
## Selecting only clone C
## Do clustering using ComplexHeatmap with Pearson correlation as distance
## Visualize output
get_clustering <- function(){
  
  
  save_dir <- '/home/htran/storage/datasets/metastasis_results/bulk_SA919/mixing_SA919/'
  datatag <- 'Mixed_Ex3_GS'
  normalized_df <- data.table::fread(paste0(save_dir, datatag, '_sizefactor_normalized.csv.gz'))
  dim(normalized_df)  
  colnames(normalized_df)
  
    
  ## See process_metadata() above
  total_meta <- data.table::fread(paste0(save_dir, 'metadata_Hakwoo_bulkRNA_mixing_main_exp.csv'))
  # total_meta <- total_meta %>%
  #   dplyr::filter(main_clone=='C')
  # View(total_meta)
  
  
  fpkm_df <- get_normalized_data(normalized_df)
  # fpkm_df <- data.table::fread(paste0(save_dir, 'SA919_mixing_normalized_fpkm.csv.gz'))
  
  
  main_C <- total_meta %>%
    dplyr::filter(main_clone=='C' & experiment=='main_exp') %>%
    dplyr::pull(bulk_sid)
  main_A_primary <- total_meta %>%
    dplyr::filter(main_clone=='A' & experiment=='main_exp' & origin=='Primary' & pdxid=='X08472164') %>%
    dplyr::pull(bulk_sid)
  main_A_met <- total_meta %>%
    dplyr::filter(main_clone=='A' & experiment=='main_exp' & origin=='Axillary' & pdxid=='X08472164') %>%
    dplyr::pull(bulk_sid)
  
    
  script_dir <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/bulkRNAseq/'
  de_results_dir <- paste0(script_dir, 'Bulk_RNA_Seq_GS/')
  comp_all <- data.table::fread(paste0(script_dir, 'list_DE_comparisons.csv'))
  comp_C <- comp_all %>%
    dplyr::filter(group2_main_clone=='C')
  comp_A <- comp_all %>%
    dplyr::filter(group2_main_clone=='A')
  dim(comp_C)
  dim(comp_A)
  unique(comp_C$comp_DE)
  
  
  
  
  
  genes_used <- intersect(fpkm_df$ens_gene_id, normalized_df$ens_gene_id)
  
  norm_C <- normalized_df %>%
    dplyr::select(all_of(c(bulk_sids_C,'ens_gene_id'))) %>%
    dplyr::filter(ens_gene_id %in% genes_used) %>%
    tibble::column_to_rownames('ens_gene_id') %>%
    as.data.frame()
  norm_C <- norm_C[genes_used, ]
  bulk_sids_C <- c(main_C, comp_C$bulk_sid)
  fpkm_C <- fpkm_df %>%
    dplyr::select(all_of(c(bulk_sids_C,'ens_gene_id'))) %>%
    dplyr::filter(ens_gene_id %in% genes_used) %>%
    tibble::column_to_rownames('ens_gene_id')  %>%
    as.data.frame()
  fpkm_C <- fpkm_C[genes_used, ]
  dim(norm_C)
  dim(fpkm_C)
  colnames(norm_C)
  # normalized_df$ens_gene_id
  # backup_norm <- norm_df
  
  ## Get correlation
  cor_df <- tibble::tibble()
  for(sc in comp_C$bulk_sid){
    c <- cor.test(norm_C[,main_C], norm_C[,sc],method="pearson")
    cp <- cor.test(fpkm_C[,main_C], fpkm_C[,sc],method="pearson")
    ctmp <- tibble::tibble(cor_norm_counts=round(c$estimate,2),
                           cor_norm_counts_pval=c$p.value,
                           cor_norm_fpkm=round(cp$estimate,2),
                           cor_norm_fpkm_pval=cp$p.value,
                           bulk_sid=sc, versus_group2=main_C)
    cor_df <- dplyr::bind_rows(cor_df, ctmp)
  }
  dim(cor_df)
  
  cor_df <- cor_df %>%
    dplyr::left_join(comp_C, by=c('bulk_sid','versus_group2'))
  data.table::fwrite(cor_df, paste0(save_dir, 'SA919_mixing_correlation_DE_analysis_cloneC.csv.gz'))
  
  
  norm_A <- normalized_df %>%
    dplyr::select(all_of(c(bulk_sids_A,'ens_gene_id'))) %>%
    dplyr::filter(ens_gene_id %in% genes_used) %>%
    tibble::column_to_rownames('ens_gene_id') %>%
    as.data.frame()
  norm_A <- norm_A[genes_used, ]
  bulk_sids_A <- c(main_A_met, main_A_primary, comp_A$bulk_sid)
  fpkm_A <- fpkm_df %>%
    dplyr::select(all_of(c(bulk_sids_A,'ens_gene_id'))) %>%
    dplyr::filter(ens_gene_id %in% genes_used) %>%
    tibble::column_to_rownames('ens_gene_id')  %>%
    as.data.frame()
  fpkm_A <- fpkm_A[genes_used, ]
  # normalized_df$ens_gene_id
  # backup_norm <- norm_df
  
  ## Get correlation
  cor_df <- tibble::tibble()
  for(sa in comp_A$bulk_sid){
    cpri <- cor.test(norm_A[,main_A_primary], norm_A[,sa],method="pearson")
    cmet <- cor.test(norm_A[,main_A_met], norm_A[,sa],method="pearson")
    cp_pri <- cor.test(fpkm_A[,main_A_primary], fpkm_A[,sa],method="pearson")
    cp_met <- cor.test(fpkm_A[,main_A_met], fpkm_A[,sa],method="pearson")
    ctmp_pri <- tibble::tibble(cor_norm_counts=round(cpri$estimate,2),
                           cor_norm_counts_pval=cpri$p.value,
                           cor_norm_fpkm=round(cp_pri$estimate,2),
                           cor_norm_fpkm_pval=cp_pri$p.value,
                           bulk_sid=sa, versus_group2=main_A_primary)
    ctmp_met <- tibble::tibble(cor_norm_counts=round(cmet$estimate,2),
                           cor_norm_counts_pval=cmet$p.value,
                           cor_norm_fpkm=round(cp_met$estimate,2),
                           cor_norm_fpkm_pval=cp_met$p.value,
                           bulk_sid=sa, versus_group2=main_A_met)
    cor_df <- dplyr::bind_rows(cor_df, ctmp_pri)
    cor_df <- dplyr::bind_rows(cor_df, ctmp_met)
  }
  dim(cor_df)
  dim(comp_A)
  cor_df <- cor_df %>%
    dplyr::left_join(comp_A, by=c('bulk_sid','versus_group2'))
  data.table::fwrite(cor_df, paste0(save_dir, 'SA919_mixing_correlation_DE_analysis_cloneA.csv.gz'))
  # View(cor_df)
  # colnames(norm_df)
  ## Just testing
  # main_A_met <- 'SA919X4XB40503'
  # main_A_pri <- 'SA919X4XB09563'
  # mixing_A_pri <- 'AT24180'
  # obs_cols <- c(main_A_met, main_A_pri, mixing_A_pri)
  # norm_df <- normalized_df %>%
  #   tibble::column_to_rownames('ens_gene_id') %>%
  #   dplyr::select(all_of(obs_cols))
  
  # norm_df <- norm_df[,obs_cols]
  # dim(norm_df)
  # ## No potential to use here
  # cor.test(norm_df$AT24180,norm_df$SA919X4XB40503,method="spearman")
  # cor.test(norm_df$AT24180,norm_df$SA919X4XB09563,method="spearman")
  # rownames(norm_df)[1:5]
  
  
  ## Note: here we used all genes, TODO: selecting only DE genes, and redo this calculation
  # cor.test(norm_df$tpm_AT24180,norm_df$tpm_SA919X4XB40503,method="pearson")
  # cor.test(norm_df$tpm_AT24180,norm_df$tpm_SA919X4XB09563,method="pearson")
  # 
  # cor.test(norm_df$fpkm_AT24180,norm_df$fpkm_SA919X4XB40503,method="pearson")
  # cor.test(norm_df$fpkm_AT24180,norm_df$fpkm_SA919X4XB09563,method="pearson")
  # 
  # cor.test(norm_df$AT24180,norm_df$SA919X4XB40503,method="pearson")
  # cor.test(norm_df$AT24180,norm_df$SA919X4XB09563,method="pearson")
  
  
  
  
  
  # View(metasamples)
  total_meta$bulk_sid
  metasamples <- total_meta
  metasamples <- metasamples %>%
    as.data.frame() %>%
    dplyr::filter(bulk_sid %in% colnames(norm_df))
  rownames(metasamples) <- metasamples$bulk_sid
  # metasamples$transplanted_mouse_id <- ifelse(metasamples$experiment=='main_exp',1,metasamples$transplanted_mouse_id)
  unique(metasamples$origin)
  predefined_cols <- c('#EE220C','#004D80')
  names(predefined_cols) <- c('Metastasis','Primary')
  # colors_use <- predefined_cols[metasamples[colnames(norm_df),'mainsite']]
  # cols_use=c("#66A61E", "#66A61E", "#666666") # green and grey
  # names(colors_use) <- colnames(norm_df)
  # unique(metasamples[colnames(norm_df),'main_clone'])
  clones_color <- c('A'='#66C2A5','B'='#FC8D62','C'='#8DA0CB')
  # colors_clone_use <- clones_color[metasamples[colnames(norm_df),'main_clone']]
  # cols_use=c("#66A61E", "#66A61E", "#666666") # green and grey
  # names(colors_clone_use) <- colnames(norm_df)
  mouse_cols <- brewer.pal(8, "Accent")[1:length(unique(metasamples$transplanted_mouse_id))]
  names(mouse_cols) <- unique(metasamples$transplanted_mouse_id)
  
  origin_cols <- brewer.pal(8, "Dark2")[1:length(unique(metasamples$origin))]
  names(origin_cols) <- unique(metasamples$origin)
  
  
  top_anno = ComplexHeatmap::HeatmapAnnotation(MainSite = factor(metasamples[colnames(norm_df),'mainsite']),
                                               MainClone=factor(metasamples[colnames(norm_df),'main_clone']),
                                               # MouseId=factor(metasamples[colnames(norm_df),'transplanted_mouse_id']),
                                               # Origin=factor(metasamples[colnames(norm_df),'origin']),
                                               col = list(MainSite=predefined_cols, 
                                                          MainClone=clones_color#,
                                                          # MouseId=mouse_cols#,
                                                          # Origin=origin_cols
                                                          )
                                               ) #, MainClone=colors_use
  # library(ComplexHeatmap)
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
  p
  
  
  png(paste0(save_dir,datatag,"_hierarchial_clusters.png"), height = 2*800, width=2*850, res = 2*72)
  print(p)
  dev.off()  
  
}  