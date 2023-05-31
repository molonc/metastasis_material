

results_10x_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_rna/deg_analysis/'


SA1035_B_D <- read.table(paste0(results_10x_dir,sub_fn[1],'/top_up_pathway_cls_',sub_fn[1],'.txt'),
                         sep = '\t', header=T, check.names=F)
SA1035_E_D <- read.table(paste0(results_10x_dir,sub_fn[2],'/top_up_pathway_cls_',sub_fn[2],'.txt'),
                         sep = '\t', header=T, check.names=F)
View(head(SA1035_B_D))


results_10x_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_total/deg_analysis/'
sub_fn <- c('SA535_GF_D')


results_10x_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA609_rna/deg_analysis/'



# markers_ls_output <- read.table(paste0(save_dir,'SA1035_E_D/de_significant_genes.txt'), header=T,sep='\t',row.names = 1)
save_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_total/deg_analysis/SA535_GF_D/'

up_v2 <- read.table(paste0(save_dir,'top_up_pathway_cls_',sub_fn[1],'.txt'),
                         sep = '\t', header=T, check.names=F)
down_v2 <- read.table(paste0(save_dir,'top_down_pathway_cls_',sub_fn[1],'.txt'),
                    sep = '\t', header=T, check.names=F)
summary(up_v2$padj)  # [0.002, 0.16]
summary(up_v2$pval)  # [0.0001, 0.07]
up_v2 <- up_v2[up_v2$padj<0.05,] # 12 significant pathways
dim(up_v2) 

up_v2 <- up_v2[up_v2$pval<0.05,]
dim(up_v2)  # 16 significant pathways


down_v2 <- down_v2[down_v2$padj<0.05,] # 1 significant pathways
dim(down_v2)

down_v2 <- down_v2[down_v2$pval<0.05,] # 2 significant pathways
dim(down_v2)


# 10  8
# 1 8

save_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_total/deg_analysis/SA535_GF_D_v1/'
up_v1 <- read.table(paste0(save_dir,'top_up_pathway_cls_',sub_fn[1],'.txt'),
                    sep = '\t', header=T, check.names=F)
down_v1 <- read.table(paste0(save_dir,'top_down_pathway_cls_',sub_fn[1],'.txt'),
                      sep = '\t', header=T, check.names=F)

up_v1 <- up_v1[up_v1$padj<0.05,] # 10 significant pathways
dim(up_v1) 
summary(up_v1$padj)  # [0.002, 0.294]
summary(up_v1$pval)  # [0.0001, 0.141]

up_v1 <- up_v1[up_v1$pval<0.05,] # 11 significant pathways
dim(up_v1)  # 11 significant pathways


down_v1 <- down_v1[down_v1$padj<0.05,] # 1 significant pathways
dim(down_v1)

down_v1 <- down_v1[down_v1$pval<0.05,] # 1 significant pathways
dim(down_v1)


markers_ls <- read.table(paste0(save_dir,'total_markers.txt'), header=T,sep='\t',row.names = 1)
gmt_dir <- paste0(input_dir,"/biodatabase/")



#SA535 pathway
library(dplyr)
library(tidyverse)
pt_use <- 'hallmark/'
pattern_use <- 'HALLMARK_'

pt_use <- 'kegg/'
pattern_use <- 'KEGG_'

# obs_pw <- 'top_down_pathway_cls_'
obs_pw <- 'top_up_pathway_cls_'

input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_total/deg_analysis/'
tag <- 'SA535_GF_D'
up_SA535_GF_D <- read.table(paste0(input_dir,'SA535_GF_treated_D/',pt_use, obs_pw, tag,'.txt'),
                    sep = '\t', header=T, check.names=F)
up_SA535_GF_D <- up_SA535_GF_D %>% dplyr::select(pathway, padj)
colnames(up_SA535_GF_D)[which(colnames(up_SA535_GF_D)=='padj')] <- tag

tag <- 'SA535_E_D'
up_SA535_E_D <- read.table(paste0(input_dir,'SA535_E_treated_D/',pt_use,obs_pw,tag,'.txt'),
                            sep = '\t', header=T, check.names=F)
up_SA535_E_D <- up_SA535_E_D %>% dplyr::select(pathway, padj)
colnames(up_SA535_E_D)[which(colnames(up_SA535_E_D)=='padj')] <- tag

tag <- 'SA535_E_C' #no kegg down pathway
up_SA535_E_C <- read.table(paste0(input_dir,'SA535_E_treated_C/',pt_use,obs_pw,tag,'.txt'),
                           sep = '\t', header=T, check.names=F)
up_SA535_E_C <- up_SA535_E_C %>% dplyr::select(pathway, padj)
colnames(up_SA535_E_C)[which(colnames(up_SA535_E_C)=='padj')] <- tag


#SA1035 pathway 
input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_rna/deg_analysis/'
tag <- 'SA1035_E_D'
up_SA1035_E_D <- read.table(paste0(input_dir,'SA1035_E_treated_D/',pt_use,obs_pw,tag,'.txt'),
                            sep = '\t', header=T, check.names=F)
up_SA1035_E_D <- up_SA1035_E_D %>% dplyr::select(pathway, padj)
colnames(up_SA1035_E_D)[which(colnames(up_SA1035_E_D)=='padj')] <- tag


tag <- 'SA1035_B_D'  # remove from list
up_SA1035_B_D <- read.table(paste0(input_dir,'SA1035_B_treated_D/',pt_use,obs_pw,tag,'.txt'),
                           sep = '\t', header=T, check.names=F)
up_SA1035_B_D <- up_SA1035_B_D %>% dplyr::select(pathway, padj)
colnames(up_SA1035_B_D)[which(colnames(up_SA1035_B_D)=='padj')] <- tag

tag <- 'SA1035_E_A'
up_SA1035_E_A <- read.table(paste0(input_dir,'SA1035_E_treated_A/',pt_use,obs_pw,tag,'.txt'),
                            sep = '\t', header=T, check.names=F)
up_SA1035_E_A <- up_SA1035_E_A %>% dplyr::select(pathway, padj)
colnames(up_SA1035_E_A)[which(colnames(up_SA1035_E_A)=='padj')] <- tag

tag <- 'SA1035_B_A' # no kegg down-pathway
up_SA1035_B_A <- read.table(paste0(input_dir,'SA1035_B_treated_A/',pt_use,obs_pw,tag,'.txt'),
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


tag <- 'SA609_UTTTT_UT'
upSA609_UTTTT_UT <- read.table(paste0(input_dir,'SA609_UTTTT_UT/',pt_use,obs_pw,tag,'.txt'),
                              sep = '\t', header=T, check.names=F)
upSA609_UTTTT_UT <- upSA609_UTTTT_UT %>% dplyr::select(pathway, padj)
colnames(upSA609_UTTTT_UT)[which(colnames(upSA609_UTTTT_UT)=='padj')] <- tag



df_ls <- list(up_SA535_GF_D,up_SA535_E_D,up_SA535_E_C,
              up_SA1035_E_D,up_SA1035_E_A,up_SA1035_B_A,
              upSA609_R_UTTTT_E,upSA609_R_UTTT_E,upSA609_UTTTT_UT)

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

clusters <- c(rep('SA535',3), rep('SA1035',3),rep('SA609',3))
for(df in df_ls){
  pathway_stat <- pathway_stat %>% left_join(df, by = 'pathway')
}
# pathway_stat$pathway <- rownames(pathway_stat)
save_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_total/deg_analysis/summary/'
if(!file.exists(save_dir)){
  dir.create(save_dir)
}
pathway_stat$pathway <- gsub(pattern_use,'',pathway_stat$pathway)
write.csv(pathway_stat,file = paste0(save_dir,pattern_use,'_summary_SA535_SA1035_SA609_up_pathways.csv'), row.names=F, quote=F)
dim(pathway_stat)

pathway_stat <- column_to_rownames(pathway_stat, 'pathway')


hm <- ComplexHeatmap::Heatmap(as.matrix(pathway_stat), na_col = "white",
                              show_column_names=T, 
                              show_row_names = T,
                              cluster_rows=F,cluster_columns=F,
                              column_split = clusters,
                              name = 'P_adj', 
                              column_names_gp = grid::gpar(fontsize = 10), 
                              row_names_gp = grid::gpar(fontsize = 15)) 


png(paste0(save_dir,pattern_use,'_summary_SA535_SA1035_SA609_up_pathways.png'), height = 2*30*dim(pathway_stat)[1]+100, width=2*1200,res = 2*72)
print(hm)
dev.off()


library(stringr)
library(ComplexHeatmap)
library(grid)
input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_total/deg_analysis/'
tag <- 'SA535_GF_D'
up_SA535_GF_D <- read.csv(paste0(input_dir,'SA535_GF_treated_D/',pt_use,'SA535_up_pathway.csv'),
                          row.names=1, check.names=F, stringsAsFactors=F)
pt_ls <- rownames(up_SA535_GF_D)
pt_ls <- str_replace_all(tolower(pt_ls), "^hallmark_", "")
rownames(up_SA535_GF_D)  <- pt_ls
# rownames(up_SA535_GF_D)  <- paste0('SA535_',pt_ls)
dim(up_SA535_GF_D)
max(up_SA535_GF_D, na.rm = T)

top_pathway_ext <- read.table(paste0(input_dir,'SA535_GF_treated_D/',pt_use,'top_up_pathway_cls_',tag,'.txt'),
                            sep = '\t', header=T, check.names=F)
logFC_df <- read.table(paste0(input_dir,'SA535_GF_treated_D/total_markers.txt'),
                          header=T,row.names=1, sep='\t', check.names=F, stringsAsFactors=F)
write.csv(pw_df, file = paste0(input_dir,'SA535_GF_treated_D/',pt_use,'SA535_up_pathway.csv'), row.names=T, quote=F)

#SA1035 pathway 
input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_rna/deg_analysis/'
tag <- 'SA1035_E_A'
up_SA1035_E_A <- read.csv(paste0(input_dir,'SA1035_E_treated_A/',pt_use,'SA1035_up_pathway.csv'),
                          row.names=1, check.names=F, stringsAsFactors=F)
pt_ls <- rownames(up_SA1035_E_A)
pt_ls <- str_replace_all(tolower(pt_ls), "^hallmark_", "")
rownames(up_SA1035_E_A)  <- pt_ls  #paste0('SA1035_',pt_ls)
dim(up_SA1035_E_A)
max(up_SA1035_E_A, na.rm = T)
View(up_SA1035_E_A['tgf_beta_signaling',])
top_pathway_ext <- read.table(paste0(input_dir,'SA1035_E_treated_A/',pt_use,'top_up_pathway_cls_',tag,'.txt'),
                              sep = '\t', header=T, check.names=F)
logFC_df <- read.table(paste0(input_dir,'SA1035_E_treated_A/total_markers.txt'),
                       header=T,row.names=1, sep='\t', check.names=F, stringsAsFactors=F)

write.csv(pw_df, file = paste0(input_dir,'SA1035_E_treated_A/',pt_use,'SA1035_up_pathway.csv'), row.names=T, quote=F)



#SA609
input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA609_rna/deg_analysis/'
tag <- 'SA609_R_UTTTT_E'
tag1 <- 'SA609_R_E'
upSA609_R_UTTTT_E <- read.csv(paste0(input_dir,'SA609_R_UTTTT_E/',pt_use,'SA609_up_pathway.csv'),
                              row.names=1, check.names=F, stringsAsFactors=F)
pt_ls <- rownames(upSA609_R_UTTTT_E)

pt_ls <- str_replace_all(tolower(pt_ls), "^hallmark_", "")
rownames(upSA609_R_UTTTT_E)  <- pt_ls  #paste0('SA609_',pt_ls)
max(upSA609_R_UTTTT_E, na.rm = T)

top_pathway_ext <- read.table(paste0(input_dir,'SA609_R_UTTTT_E/',pt_use,'top_up_pathway_cls_',tag1,'.txt'),
                              sep = '\t', header=T, check.names=F)
logFC_df <- read.table(paste0(input_dir,'SA609_R_UTTTT_E/total_markers.txt'),
                       header=T,row.names=1, sep='\t', check.names=F, stringsAsFactors=F)

write.csv(pw_df, file = paste0(input_dir,'SA609_R_UTTTT_E/',pt_use,'SA609_up_pathway.csv'), row.names=T, quote=F)



obs_genes <- c('tnfa_signaling_via_nfkb','tgf_beta_signaling',
               'epithelial_mesenchymal_transition','uv_response_dn',
               'mitotic_spindle')

up_SA535_GF_D <- up_SA535_GF_D[rownames(up_SA535_GF_D) %in% obs_genes,]

up_SA1035_E_A <- up_SA1035_E_A[rownames(up_SA1035_E_A) %in% obs_genes,]

upSA609_R_UTTTT_E <- upSA609_R_UTTTT_E[rownames(upSA609_R_UTTTT_E) %in% obs_genes,]
max(upSA609_R_UTTTT_E, na.rm = T)
dim(up_SA535_GF_D)
rownames(up_SA535_GF_D) <- paste0('SA535_',rownames(up_SA535_GF_D))
rownames(up_SA1035_E_A) <- paste0('SA1035_',rownames(up_SA1035_E_A))
rownames(upSA609_R_UTTTT_E) <- paste0('SA609_',rownames(upSA609_R_UTTTT_E))
dim(up_SA1035_E_A)
dim(upSA609_R_UTTTT_E)
combined_genes_df <- combine_dataframe(up_SA535_GF_D, up_SA1035_E_A, upSA609_R_UTTTT_E)

pt_clusters <- c(rep(c('SA535','SA1035','SA609'),c(5,5,5)))

genes_rv <- c()
for(g in colnames(combined_genes_df)){
  if(sum(is.na(combined_genes_df[,g]))==nrow(combined_genes_df)){
    genes_rv <- c(genes_rv,g)
  }
}

combined_genes_df <- combined_genes_df[,!colnames(combined_genes_df) %in% genes_rv]
dim(combined_genes_df)
hm_up <- ComplexHeatmap::Heatmap(as.matrix(combined_genes_df), na_col = "white",
                                 show_column_names=T, 
                                 show_row_names = T,
                                 cluster_rows=F,cluster_columns=F,
                                 row_split  = pt_clusters,
                                 column_names_gp = grid::gpar(fontsize = 4.5), 
                                 row_names_gp = grid::gpar(fontsize = 9),
                                 name = "avg_logFC")   
save_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_total/deg_analysis/'

png(paste0(save_dir,'HALLMARK_marker_genes_SA535_SA1035_SA609_5_resistant_pathways.png'), height = 2*30*dim(combined_genes_df)[1]+100, width=2*1300,res = 2*72)
print(hm_up)
dev.off()



combine_dataframe <- function(df1, df2, df3){
  genes <- union(colnames(df1),colnames(df2))
  genes <- unique(union(genes,colnames(df3)))
  print(length(genes))
  genes_add <- setdiff(genes, colnames(df1))
  length(genes_add)
  df1[,genes_add] <- NA
  df1 <- df1[,genes]
  
  genes_add <- setdiff(genes,colnames(df2))
  length(genes_add)
  df2[,genes_add] <- NA
  df2 <- df2[,genes]
  
  genes_add <- setdiff(genes, colnames(df3))
  length(genes_add)
  df3[,genes_add] <- NA
  df3 <- df3[,genes]
  combine_pw <- rbind(df1, df2)
  combine_pw <- rbind(combine_pw, df3)
  print(dim(combine_pw))
  return(combine_pw)
  
}