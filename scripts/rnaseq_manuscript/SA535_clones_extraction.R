library(dplyr)
input_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/clonealign_v3/'
sce <- readRDS(paste0(input_dir, 'combined_SA535_clones.rds'))
dim(sce)

table(sce$clone, sce$mouse_id)
unique(sce$pdxid)
sce$desc <- paste0(sce$pdxid,"_",sce$Site_origin)
table(sce$clone, sce$desc)

stat <- data.frame(desc=sce$desc, clone_id=sce$clone)
dim(stat)
head(stat)
stat <- stat %>%
  dplyr::filter(clone_id!='unassigned')

stat1 <- stat %>%
  dplyr::group_by(clone_id) %>%
  dplyr::summarise(nb_cells=n()) %>%
  dplyr::filter(nb_cells>30)
stat1

stat <- stat %>%
  dplyr::filter(clone_id %in% stat1$clone_id)

stat1 <- stat %>%
  dplyr::group_by(clone_id, desc) %>%
  dplyr::summarise(nb_cells=n())

for(d in unique(stat1$desc)){
  tmp <- stat1 %>%
    dplyr::filter(desc==d)
  tmp$clone_info <- paste0(tmp$clone_id, '(',tmp$nb_cells,')')
  print(paste0(tmp$clone_info, collapse = ' '))
}
