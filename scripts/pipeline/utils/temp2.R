library(dplyr)
library(tidyverse)
pt_use <- 'hallmark/'
pattern_use <- 'HALLMARK_'

# pt_use <- 'kegg/'
# pattern_use <- 'KEGG_'

# obs_pw <- 'top_down_pathway_cls_'
obs_pw <- 'top_up_pathway_cls_'

input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_total/deg_analysis/'
tag <- 'SA535_GF_D'
up_SA535_GF_D <- read.table(paste0(input_dir,'SA535_GF_drugholiday_D/',pt_use, obs_pw, tag,'.txt'),
                            sep = '\t', header=T, check.names=F)
up_SA535_GF_D <- up_SA535_GF_D %>% dplyr::select(pathway, padj)
colnames(up_SA535_GF_D)[which(colnames(up_SA535_GF_D)=='padj')] <- tag

tag <- 'SA535_E_D'
up_SA535_E_D <- read.table(paste0(input_dir,'SA535_E_drugholiday_D/',pt_use,obs_pw,tag,'.txt'),
                           sep = '\t', header=T, check.names=F)
up_SA535_E_D <- up_SA535_E_D %>% dplyr::select(pathway, padj)
colnames(up_SA535_E_D)[which(colnames(up_SA535_E_D)=='padj')] <- tag

tag <- 'SA535_E_C' #no kegg down pathway
up_SA535_E_C <- read.table(paste0(input_dir,'SA535_E_drugholiday_C/',pt_use,obs_pw,tag,'.txt'),
                           sep = '\t', header=T, check.names=F)
up_SA535_E_C <- up_SA535_E_C %>% dplyr::select(pathway, padj)
colnames(up_SA535_E_C)[which(colnames(up_SA535_E_C)=='padj')] <- tag


#SA1035 pathway 
input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_rna/deg_analysis/'
tag <- 'SA1035_E_D' # no drug holiday down
up_SA1035_E_D <- read.table(paste0(input_dir,'SA1035_E_drugholiday_D/',pt_use,obs_pw,tag,'.txt'),
                            sep = '\t', header=T, check.names=F)
up_SA1035_E_D <- up_SA1035_E_D %>% dplyr::select(pathway, padj)
colnames(up_SA1035_E_D)[which(colnames(up_SA1035_E_D)=='padj')] <- tag


tag <- 'SA1035_B_D'  # no down-pathway
up_SA1035_B_D <- read.table(paste0(input_dir,'SA1035_B_drugholiday_D/',pt_use,obs_pw,tag,'.txt'),
                            sep = '\t', header=T, check.names=F)
up_SA1035_B_D <- up_SA1035_B_D %>% dplyr::select(pathway, padj)
colnames(up_SA1035_B_D)[which(colnames(up_SA1035_B_D)=='padj')] <- tag

tag <- 'SA1035_E_A'
up_SA1035_E_A <- read.table(paste0(input_dir,'SA1035_E_drugholiday_A/',pt_use,obs_pw,tag,'.txt'),
                            sep = '\t', header=T, check.names=F)
up_SA1035_E_A <- up_SA1035_E_A %>% dplyr::select(pathway, padj)
colnames(up_SA1035_E_A)[which(colnames(up_SA1035_E_A)=='padj')] <- tag

tag <- 'SA1035_B_A' 
up_SA1035_B_A <- read.table(paste0(input_dir,'SA1035_B_drugholiday_A/',pt_use,obs_pw,tag,'.txt'),
                            sep = '\t', header=T, check.names=F)
up_SA1035_B_A <- up_SA1035_B_A %>% dplyr::select(pathway, padj)
colnames(up_SA1035_B_A)[which(colnames(up_SA1035_B_A)=='padj')] <- tag


#SA609
input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA609_rna/deg_analysis/'
tag <- 'SA609_R_UTTTT_E'
tag1 <- 'SA609_R_E'
upSA609_R_UTTTT_E <- read.table(paste0(input_dir,'SA609_R_UTTTT_E/',pt_use,obs_pw,tag1,'.txt'),
                                sep = '\t', header=T, check.names=F)
upSA609_R_UTTTT_E <- upSA609_R_UTTTT_E %>% dplyr::select(pathway, padj)
colnames(upSA609_R_UTTTT_E)[which(colnames(upSA609_R_UTTTT_E)=='padj')] <- tag


tag <- 'SA609_R_UTTT_E'
upSA609_R_UTTT_E <- read.table(paste0(input_dir,'SA609_R_UTTT_E/',pt_use,obs_pw,tag1,'.txt'),
                               sep = '\t', header=T, check.names=F)
upSA609_R_UTTT_E <- upSA609_R_UTTT_E %>% dplyr::select(pathway, padj)
colnames(upSA609_R_UTTT_E)[which(colnames(upSA609_R_UTTT_E)=='padj')] <- tag

tag <- 'SA609_R_UTT_E'
upSA609_R_UTT_E <- read.table(paste0(input_dir,'SA609_R_UTT_E/',pt_use,obs_pw,tag1,'.txt'),
                              sep = '\t', header=T, check.names=F)
upSA609_R_UTT_E <- upSA609_R_UTT_E %>% dplyr::select(pathway, padj)
colnames(upSA609_R_UTT_E)[which(colnames(upSA609_R_UTT_E)=='padj')] <- tag





df_ls <- list(up_SA535_GF_D,up_SA535_E_D,up_SA535_E_C,
              up_SA1035_B_D,up_SA1035_E_A,up_SA1035_B_A)

pathway_ls <- list()
c <- 0
for(df in df_ls){
  c <- c+1
  pathway_ls[[c]] <- df$pathway
}
length(pathway_ls)
pathway_ls <- unique(unlist(pathway_ls))
pathway_stat <- data.frame(pathway=pathway_ls, row.names = pathway_ls)

# Up pathway
# clusters <- c(rep('SA535',3),rep('SA1035',4),rep('SA609',3))
# for(df in df_ls){
#   pathway_stat <- pathway_stat %>% left_join(df, by = 'pathway')
# }

clusters <- c(rep('SA535',3), rep('SA1035',3))
for(df in df_ls){
  pathway_stat <- pathway_stat %>% left_join(df, by = 'pathway')
}
# pathway_stat$pathway <- rownames(pathway_stat)
save_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_total/deg_analysis/summary/'
if(!file.exists(save_dir)){
  dir.create(save_dir)
}
pathway_stat$pathway <- gsub(pattern_use,'',pathway_stat$pathway)
write.csv(pathway_stat,file = paste0(save_dir,pattern_use,'_summary_SA535_SA1035_SA609_drugholiday_down_pathways.csv'), row.names=F, quote=F)
dim(pathway_stat)

pathway_stat <- column_to_rownames(pathway_stat, 'pathway')


hm <- ComplexHeatmap::Heatmap(as.matrix(pathway_stat), na_col = "white",
                              show_column_names=T, 
                              show_row_names = T,
                              cluster_rows=F,cluster_columns=F,
                              column_split = clusters,
                              name = 'P_adj', 
                              column_names_gp = grid::gpar(fontsize = 10), 
                              row_names_gp = grid::gpar(fontsize = 10)) 


png(paste0(save_dir,pattern_use,'_summary_SA535_SA1035_SA609_drugholiday_down_pathways.png'), height = 2*30*dim(pathway_stat)[1]+100, width=2*1100,res = 2*72)
print(hm)
dev.off()

