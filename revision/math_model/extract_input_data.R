
library(dplyr)
input_dir <- '/home/htran/Projects/hakwoo_project/rscript/math_model/data_met_proj/'
dir.create(input_dir)
input_fn <- paste0(input_dir,'S9_Table1_ovariancancer.csv')

clones_total <- data.table::fread(input_fn) %>% as.data.frame()
for(p in unique(clones_df$patient_id)){
  clones_df <- clones_total %>% 
    dplyr::filter(patient_id==p)
  print(dim(clones_df))
  data.table::fwrite(clones_df, paste0(input_dir,'S9_T1_OVA_patient_',p,'_raw.csv'))
  
  pt1 <- clones_df
  pt1 <- pt1 %>%
    dplyr::select(paper_id, clone_id, prevalence)
  pt1 <- pt1 %>%
    tidyr::pivot_wider(names_from = clone_id, values_from = prevalence)
  pt1 <- pt1 %>%
    dplyr::rename(tumour_id = paper_id)
  data.table::fwrite(pt1, paste0(input_dir,'S9_T1_OVA_patient_',p,'.csv'))
  
}

length(unique(clones_df$clone_id))
length(unique(clones_df$paper_id))
unique(clones_df$prevalence)
View(head(clones_df))


# %     y = matrix of clone frequencies
# %         Each row is a different tumor; each column is a different clone.
# %         Each row will be normalized so that it sums to 100%.

# %     m = row index of primary tumor
# %         Should be an integer between 1 and the number of tumors.
# %         If the primary tumor is not known or included, set m=0 instead.

pt1 <- clones_df

pt1 <- pt1 %>%
  dplyr::select(paper_id, clone_id, prevalence)

pt1 <- pt1 %>%
  tidyr::pivot_wider(names_from = clone_id, values_from = prevalence)
View(pt1)
pt1 <- pt1 %>%
  dplyr::rename(tumour_id = paper_id)

data.table::fwrite(pt1, '~/Desktop/S9_clone_proportions_metastasis_patient1.csv')

