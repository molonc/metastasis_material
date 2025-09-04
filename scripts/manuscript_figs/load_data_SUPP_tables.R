
suppressPackageStartupMessages({
  library("dplyr")
  library("ggplot2")
})


# Get supp tables for manuscript submission

# supp table 2: ihc tnbc markers
input_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/IHC_tumour_markers/'
df1 <- data.table::fread(paste0(input_dir, 'IHC_TMA_scores_metastasis_project/IHC_scoring_19-012_Hakwoo_PDX_BrCa.csv'))
df2 <- data.table::fread(paste0(input_dir, 'IHC_TMA_scores_metastasis_project/IHC_scoring_19-006_HakWoo_PDX_BrCa__MASTERFILE2.csv'))
df1$Captured_ID <- 'IHC_scoring_19-012'
df2$Captured_ID <- 'IHC_scoring_19-006'


dim(df1)
dim(df2)
colnames(df1)
colnames(df2)
View(df2[1:3,1:6])
View(df1[1:3,1:6])

colnames(df1) <- gsub(' ', '_', colnames(df1))
df1$Core_ID <- as.character(df1$Core_ID)
df1 <- as.data.frame(df1)
for(cn in colnames(df1)){
  df1[,cn] <- as.character(df1[,cn])
}
  
  
colnames(df2) <- gsub(' ', '_', colnames(df2))
df2$V1 <- NULL
df2$note <- df2$others
df2$others <- NULL
df2$ER_intensity <- as.character(df2$ER_intensity)
df2 <- as.data.frame(df2)
for(cn in colnames(df2)){
  df2[,cn] <- as.character(df2[,cn])
}


sum(colnames(df2) %in% colnames(df1))
colnames(df2)[!colnames(df2) %in% colnames(df1)]
# Capture
# Block_ID = Accession #
# bio_features <- c('Captured_ID','Core_ID','Block_ID','PDX_ID','SA_ID')
typeof(df1$Core_ID)
typeof(df2$Core_ID)
df <- dplyr::bind_rows(df1, df2)
dim(df)
View(head(df))

df <- df %>%
  dplyr::select(Captured_ID, everything())
df$Captured_ID <- gsub('IHC_scoring_', '', df$Captured_ID)

df <- df %>%
  dplyr::rename(FileName=Captured_ID)

data.table::fwrite(df, paste0(input_dir, 'SUPP_Table2_IHC_scoring_19-012_19-006_Hakwoo_PDX_BrCa.csv'))





get_expression_data_table7 <- function(obs_clones, subtag){
  # data_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/preprocessed_09April2024/SA919_DE_analysis_DESeq2_Hoa_09April2024/SA919_CMetMainMix_Bpri_DESeq2/'
  # de_df <- data.table::fread(paste0(data_dir, 'DE_signif_genes.csv.gz'))  
  # dim(de_df)
  sids <- paste0('clone', obs_clones)
  data_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/results_bulkRNAseq/SA919_full/'
  save_fig_dir <- paste0(data_dir, subtag, '/')
  exp_df <- data.table::fread(paste0(save_fig_dir, subtag, '_DE_genes.csv.gz'))
  exp_df$DE_comparison <- subtag
  print(subtag)
  print(dim(exp_df))
  
  if(!'symbol' %in% colnames(exp_df)){
    ref <- annotables::grch38 %>%
      dplyr::select(ensgene, symbol, chr) %>%
      dplyr::rename(ensembl_gene_id=ensgene, gene_symbol=symbol)
    ref <- ref[!duplicated(ref$ensembl_gene_id),]
    exp_df <- exp_df %>%
      dplyr::left_join(ref, by='ensembl_gene_id')
  }
  print(colnames(exp_df))
  return(exp_df)
}

load_supp_table7 <- function(){
  input_dir <- "/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/preprocessed_09April2024/"
  save_dir <- paste0(input_dir, 'drivernet_demo/')
  # 
  # obs_clones <- c('A','B')
  
  de_comps <- list(c('A','B'), c('A','C'))
  total_de <- tibble::tibble()
  for(obs_clones in de_comps){
    # print(obs_clones)
    subtag <- paste0(obs_clones[2],'met_',obs_clones[1],'pri')
    print(subtag)
    de_df <- get_expression_data_table7(obs_clones, subtag)
    total_de <- dplyr::bind_rows(total_de, de_df)
  }
  
  
  obs_clones <- c('B','C')
  subtag <- paste0(obs_clones[2],'met_',obs_clones[1],'pri')
  print(subtag)
  data_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/results_bulkRNAseq/SA919_full/'
  save_fig_dir <- paste0(data_dir, subtag, '/')
  exp_df <- data.table::fread(paste0(save_fig_dir, subtag, '_DE_genes.csv.gz'))
  exp_df$DE_comparison <- subtag
  print(subtag)
  print(dim(exp_df))
  print(colnames(exp_df))
  sum(exp_df$log2FoldChange==exp_df$logFC)
  exp_df$logFC <- NULL
  exp_df <- exp_df %>%
    dplyr::rename(ensembl_gene_id=ens_gene_id, 
                  gene_symbol = symbol)
  
  colnames(exp_df)[!colnames(exp_df) %in% colnames(total_de)]
  total_de <- dplyr::bind_rows(total_de, exp_df)
  dim(total_de)
  View(head(total_de))
  colnames(total_de)
  total_de <- total_de %>%
    dplyr::mutate(baseMean=round(baseMean, 2),
                  log2FoldChange=round(log2FoldChange, 2),
                  lfcSE=round(lfcSE, 2),
                  stat=round(stat, 2), 
                  pvalue=round(pvalue, 5),
                  padj=round(padj, 5))
  
  save_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/supp_tables/'
  data.table::fwrite(total_de, paste0(save_dir, 'SuppTable7_DE_analysis_bulkRNAseq_part1.csv'))
  
  
  input_dir <- "/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/results_bulkRNAseq/SA919_full/"
  de_comps <- list(c('A','A'), c('B','B'))
  total_de <- tibble::tibble()
  for(obs_clones in de_comps){
    # print(obs_clones)
    subtag <- paste0(obs_clones[2],'met_',obs_clones[1],'pri')
    print(subtag)
    de_df <- get_expression_data_table7(obs_clones, subtag)
    total_de <- dplyr::bind_rows(total_de, de_df)
  }
  total_de <- total_de %>%
    dplyr::rename(gene_symbol = symbol) %>%
    dplyr::mutate(baseMean=round(baseMean, 2),
                  log2FoldChange=round(log2FoldChange, 2),
                  lfcSE=round(lfcSE, 2),
                  stat=round(stat, 2), 
                  pvalue=round(pvalue, 5),
                  padj=round(padj, 5))
  colnames(total_de)
  save_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/supp_tables/'
  data.table::fwrite(total_de, paste0(save_dir, 'SuppTable7_DE_analysis_bulkRNAseq_part2.csv'))
  
  df1 <- data.table::fread(paste0(save_dir, 'SuppTable7_DE_analysis_bulkRNAseq_part1.csv'))
  df2 <- data.table::fread(paste0(save_dir, 'SuppTable7_DE_analysis_bulkRNAseq_part2.csv'))
  dim(df1)
  dim(df2)
  
  total_df <- dplyr::bind_rows(df1, df2)
  data.table::fwrite(total_df, paste0(save_dir, 'SuppTable7_DE_analysis_bulkRNAseq_3Sept2025.csv'))
  dim(total_df)
  
  
  
}



load_supp_table8 <- function(){
  input_dir <- "/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/preprocessed_09April2024/"
  save_dir <- paste0(input_dir, 'drivernet_demo/drivernet_output_metastasis_proj/')
  # obs_clones <- c('B','C')
  # obs_clones <- c('A','B')
  obs_clones <- c('A','C')
  subtag <- paste0(obs_clones[2],'met_',obs_clones[1],'pri')
  ac_df <- data.table::fread(paste0(save_dir, subtag,'_significant_genes_DriverNet.csv'))
  dim(ac_df)
  ac_df$DE_comparison <- subtag
  typeof(ac_df$connected.drivers)
  ac_df <- ac_df %>%
    dplyr::mutate(connected.drivers=as.character(connected.drivers))
  obs_clones <- c('A','B')
  subtag <- paste0(obs_clones[2],'met_',obs_clones[1],'pri')
  subtag
  ab_df <- data.table::fread(paste0(save_dir, subtag,'_significant_genes_DriverNet.csv'))
  dim(ab_df)
  colnames(ab_df)
  ab_df$DE_comparison <- subtag
  typeof(ab_df$connected.drivers)
  
  
  abt <- ab_df %>%
    dplyr::filter(`p-value`<0.05)
  
  obs_clones <- c('B','C')
  subtag <- paste0(obs_clones[2],'met_',obs_clones[1],'pri')
  bc_df <- data.table::fread(paste0(save_dir, subtag,'_significant_genes_DriverNet.csv'))
  dim(bc_df)
  bc_df$DE_comparison <- subtag
  typeof(bc_df$connected.drivers)
  bc_df <- bc_df %>%
    dplyr::mutate(connected.drivers=as.character(connected.drivers))
  
  total_df <- dplyr::bind_rows(ab_df, ac_df)
  total_df <- dplyr::bind_rows(total_df, bc_df)
  dim(total_df)  
  colnames(total_df)
  View(head(total_df))
  total_df <- total_df %>%
    dplyr::rename(median_copy_number_cloneA=A, 
                  median_copy_number_cloneB=B,
                  median_copy_number_cloneC=C) %>%
    dplyr::select(DE_comparison, gene_symbol, gene_type, everything())
  
  save_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/supp_tables/'
  data.table::fwrite(total_df, paste0(save_dir, 'SuppTable8_DE_analysis_DriverNet_output.csv'))
  
  
  
  
}  