library(dplyr)
input_dir <- '/Users/hoatran/Documents/BCCRC_projects/hakwoo_project/code/metastasis_material/materials/'
# input_dir <- '/home/htran/storage/datasets/metastasis_results/bulk_SA919/SA919_bulk_results_v2/'
bulk <- data.table::fread(paste0(input_dir,'bulkRNAseq/SA919/library_groupings_bulk_SA919.csv')) %>% as.data.frame()
dim(bulk)
bulk$sample_id
length(unique(bulk$sample_id))
View(bulk)
data.table::fwrite(bulk, paste0(paste0(input_dir,'bulkRNAseq_dlp_SA919_for_nanopore.csv')))

# input_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_wholedata/'
clones_cnv <- data.table::fread(paste0(input_dir,'dlp_trees/SA919/cell_clones.csv.gz')) %>% as.data.frame()
dim(clones_cnv)
clones_cnv$sample_id <- get_sample_id(clones_cnv$cell_id)

clones_cnv$sample_id[1]
length(unique(clones_cnv$sample_id))
stat_cnv <- clones_cnv %>%
  dplyr::group_by(sample_id, clone_id) %>%
  dplyr::summarise(nb_cells_dlp=n()) %>%
  dplyr::filter(clone_id!='None' & nb_cells_dlp>=20)
length(unique(stat_cnv$sample_id))
dim(stat_cnv)
head(stat_cnv)
stat_cnv$clone_desc <- paste0(stat_cnv$clone_id,'(',stat_cnv$nb_cells_dlp,')')
stat_cnv <- stat_cnv %>%
  dplyr::select(sample_id, clone_desc)

stat_cnv <- stat_cnv %>% 
  dplyr::group_by(sample_id) %>% 
  # dplyr::summarise(dlp_prevalence=paste(clone_desc, sep=' '))
  dplyr::mutate(dlp_prevalence = paste0(clone_desc, collapse = " ")) 

stat_cnv <- stat_cnv[!duplicated(stat_cnv$sample_id),]
stat_cnv <- stat_cnv %>% 
  dplyr::select(-clone_desc)
dim(stat_cnv)
View(stat_cnv)
bulk <- bulk %>% left_join(stat_cnv, by='sample_id')
length(unique(stat_cnv$sample_id))

get_sample_id <- function(cell_ids) {
    labels <- sapply(strsplit(cell_ids, "-"), function(x) {
        return(x[1])
    })
    return(as.character(labels))
}
  