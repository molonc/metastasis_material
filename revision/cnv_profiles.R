

input_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/'
source(paste0(input_dir,"scripts/corrupt_tree/src/cn_change/utils.R"))
save_dir <- paste0(input_dir,'revision/CN_profile/')


datatag <- 'SA535'
copynumber_fn <- paste0(input_dir,'materials/dlp_trees/',datatag,'/total_merged_filtered_states.csv.gz')
stat_tmp <- get_median_genotype_v3(copynumber_fn, 
                       datatag, save_dir,
                       cellclone_fn=NULL, library_grouping_fn=NULL)

dim(stat_tmp)
head(stat_tmp)

# res_total <- list()
series <- c('SA919','SA535')
datatag <- 'SA919'
input_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/'
results_dir <- paste0(input_dir,"materials/dlp_trees/")
calcul_distance=T
distance_type='Hamming'
for(datatag in series){
  print(datatag)
  copynumber_fn <- paste0(results_dir, datatag, '/total_merged_filtered_states.csv.gz')
  
  cellclone_fn <- paste0(results_dir, datatag, '/cell_clones.csv.gz')
  # tmp <- get_median_genotype_v2(datatag, results_dir, copynumber_fn, cellclone_fn,
  #                               calcul_distance=T, distance_type='Manhattan')
  tmp <- get_median_genotype_v2(datatag, results_dir, copynumber_fn, cellclone_fn,
                                calcul_distance=T, distance_type='Hamming')
  # res_total[[datatag]] <- tmp
}
rm(tmp)


stat <- tibble::tibble()
for(datatag in series){
  stat <- dplyr::bind_rows(stat, res_total[[datatag]]$out_mtx)
}

dim(stat)
head(stat)
unique(stat$treatment_status)
unique(stat$datatag)
data.table::fwrite(stat, paste0(save_dir,'total_stat_median_CN_distance_6series.csv.gz'))

p_manhattan <- ggplot(stat, aes(x = datatag, y=CNA_Distance, color = datatag)) +
  geom_boxplot() + #alpha = 0.6, size = 0.3
  facet_wrap(~factor(treatment_status, levels = c("UnRx","Rx")), scales = "free_x",
             strip.position = "top",
             # switch = "x",
             nrow = 2, drop = T) + 
  theme_bw()+
  theme(legend.position = 'none',
        strip.text.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=11, family='Helvetica', color='black', angle = 20),
        axis.title.x = element_text(size=11, family='Helvetica', color='black')) + #
  labs(x=NULL, y='CNA distances between clones')

p_manhattan
output_dir <- "/home/htran/Projects/farhia_project/drug_resistant_material/materials/umap_figs/main_fig2/"
saveRDS(p_manhattan, paste0(output_dir,"median_CN_distance_6series.rds"))
stat1 <- pt4_df %>%
  # select(-SourceClone,-TargetClone)%>%
  # group_by(datatag)%>%
  summarise(median_distance=round(median(CNA_Distance),2),
            sd=round(sd(CNA_Distance),2))
stat1


# TO DO: get heatmap dlp
