
suppressPackageStartupMessages({
  library(tidyverse)
  library(annotables)
  library(dplyr)
  library(ggplot2)
})  


## Ki67

input_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/' #scripts/bulk_rna/
df <- data.table::fread(paste0(input_dir, 'materials/tumour_volumes/Fig2_ki67_annotation_DY.csv'))
View(head(df))
colnames(df)

df <- df %>%
  dplyr::rename(percent_Ki67=`Ki67 %`, SA_id=`SA ID`) %>%
  # dplyr::filter(!SA_id %in% c('SA575','SA609','SA1139')) %>%
  dplyr::filter(!SA_id %in% c('SA575')) %>%
  dplyr::mutate(SA_ids=
                  case_when(
                    grepl('SA535',SA_id) ~ 'SA535',
                    grepl('SA919',SA_id) ~ 'SA919',
                    grepl('SA1142',SA_id) ~ 'SA1142',
                    TRUE ~ SA_id
                  ))
unique(df$SA_ids)
dim(df)
p <- ggplot(df, aes(x=Damian_annotation, y=percent_Ki67)) + #, colour = SA_id
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(size=1.8, aes(colour = SA_ids), position=position_jitter(0.2))  + 
  theme_bw() + 
  labs(x='Tumour classification', y='% cells with Ki67')
p

p <- ggplot(df, aes(x=1, y=percent_Ki67)) + #, colour = SA_id
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(size=1, aes(colour = SA_ids), position=position_jitter(0.2))  + 
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  labs(x='Tumour classification', y='% cells with Ki67')
p


unique(df$Damian_annotation)
df_nonmet <- df %>%
  dplyr::filter(Damian_annotation=='non-met')
df_met <- df %>%
  dplyr::filter(Damian_annotation=='met')

  
ks.test(df_nonmet$percent_Ki67, df_met$percent_Ki67, alternative = "two.sided")

var.test(percent_Ki67 ~ Damian_annotation, df, alternative = "two.sided")
# or Method 2
wilcox.test(df_nonmet$percent_Ki67, df_met$percent_Ki67, alternative = "two.sided")

median(df_nonmet$percent_Ki67)
sd(df_nonmet$percent_Ki67)
median(df_met$percent_Ki67)
sd(df_met$percent_Ki67)

median(df$percent_Ki67)
sd(df$percent_Ki67)
round(100*sum(df$percent_Ki67>0)/dim(df)[1],2) #99%

png(paste0(input_dir,"figures/Ki67_met_nonmet.png"), height = 2*300, width=2*450, res = 2*72)
print(p)
dev.off()

png(paste0(input_dir,"figures/Ki67_all_tumours.png"), height = 2*300, width=2*250, res = 2*72)
print(p)
dev.off()


stat_df <- get_bootstrap_stat_sampling(df_nonmet$percent_Ki67, df_met$percent_Ki67, 
                                       sampling_fraction=0.8, nsamples=100, alternative_theory="two.sided")
stat_df

cox
median()

get_bootstrap_stat_sampling <- function(cis_genes, trans_genes, 
                                        sampling_fraction=0.7, nsamples=1000, alternative_theory="two.sided"){
  set.seed(42)
  # if(is.null(genome_genes)){
  #   # ref_dif <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/'
  #   # genome_genes_df <- read.csv(paste0(ref_dif, 'Symbol_ensembl.csv'), check.names = F, stringsAsFactors = F)  
  #   # # dim(genome_genes_df)
  #   # genome_genes <- unique(genome_genes_df$Symbol) # entire genes set
  #   
  #   ref <- annotables::grch38 %>%
  #     dplyr::select(ensgene,symbol) %>%
  #     dplyr::rename(gene_id=ensgene)
  #   ref <- ref[!duplicated(ref$gene_id),]
  #   genome_genes <- ref$gene_id
  # }
  nb_sampled_cis <- round(length(cis_genes) * sampling_fraction)
  nb_sampled_cis <- length(cis_genes)
  # cis_samples <- lapply(1:nsamples, function(i) sample(cis_genes, size=nb_sampled_cis, replace = T))
  # trans_samples <- lapply(1:nsamples, function(i) sample(trans_genes, size=nb_sampled_cis, replace = T))
  cis_samples <- list()
  trans_samples <- list()
  for(i in seq(nsamples)){
    cis_samples[[i]] <- sample(cis_genes, size=nb_sampled_cis, replace = T)
    trans_samples[[i]] <- sample(trans_genes, size=nb_sampled_cis, replace = T)
  }
  
  stat_df <- tibble::tibble()
  
  for(k in seq(nsamples)){
    out_stat <- ks.test(cis_samples[[k]], trans_samples[[k]], alternative=alternative_theory)
    
    out_vals <- tibble::tibble('idx'=k, 'stat'=out_stat$statistic, 'p_val'=out_stat$p.value)
    stat_df <- dplyr::bind_rows(stat_df, out_vals)
  }
  
  stat_df <- stat_df %>%
    dplyr::mutate(is_signf=
                    case_when(
                      p_val < 0.05 ~ 'T',
                      TRUE ~ 'F'
                    )
    )
  # summary(as.factor(stat_df$is_signf))
  summary_df <- stat_df %>%
    dplyr::group_by(is_signf) %>%
    dplyr::summarise(nb_val=n()) %>%
    dplyr::mutate(p_val=1-round(nb_val/dim(stat_df)[1],4))
  print(summary_df)
  # return(list(CI=as.numeric(r['our']),pval=as.numeric(pval)))
  return(summary_df)
}



## checking Vim marker
input_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/' #scripts/bulk_rna/
df <- data.table::fread(paste0(input_dir, 'materials/tumour_volumes/HE_score_grading_TMA_metastasis_proj.csv'))
dim(df)
df <- df %>%
  dplyr::rename(BlockID=V2, SAID=V3, mouseID=V4) %>%
  dplyr::filter(SAID=='SA919') %>%
  dplyr::select(BlockID, SAID, mouseID, `vimentin %`,`vimentin intensity`)
View(df)
colnames(df)
df$`vimentin %`  
df$`vimentin intensity` 
  
  
  