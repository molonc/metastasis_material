library(dplyr)
input_dir <- '/home/htran/storage/datasets/metastasis_results/bulk_SA919/SA919_bulk_results_v2/'
bulk <- data.table::fread(paste0(input_dir,'library_groupings_bulk_SA919.csv')) %>% as.data.frame()
dim(bulk)
bulk$sample_id
View(bulk)
data.table::fwrite(bulk, paste0(paste0('/home/htran/storage/bulk_dlp_SA919.csv')))

input_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_wholedata/'
clones_cnv <- data.table::fread(paste0(input_dir,'cell_clones.csv')) %>% as.data.frame()
dim(clones_cnv)



clones_cnv$sample_id <- get_sample_id(clones_cnv$cell_id)

clones_cnv$sample_id[1]

stat_cnv <- clones_cnv %>%
  dplyr::group_by(sample_id, clone_id) %>%
  dplyr::summarise(nb_cells_dlp=n()) %>%
  dplyr::filter(clone_id!='None' & nb_cells_dlp>=30)

dim(stat_cnv)
head(stat_cnv)
bulk <- bulk %>% left_join(stat_cnv, by='sample_id')
length(unique(stat_cnv$sample_id))

get_sample_id <- function(cell_ids) {
    labels <- sapply(strsplit(cell_ids, "-"), function(x) {
        return(x[1])
    })
    return(as.character(labels))
}
  