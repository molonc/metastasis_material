get_sample_id <- function(cell_ids) {
  labels <- sapply(strsplit(cell_ids, "-"), function(x) {
    return(x[1])
  })
  return(as.character(labels))
}
results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_wholedata/'
cell_clones <- read.delim(paste0(results_dir,'cell_clones.csv'), 
                          stringsAsFactors = FALSE, 
                          sep = ",", check.names = F)
dim(cell_clones)
grouping_df <- read.csv(paste0(results_dir,'library_groupings.csv'),
                        header=T,check.names=F, stringsAsFactors=F)
# colnames(grouping_df)[which(names(grouping_df) == "pdxid")] <- "timepoint"
# colnames(grouping_df)[which(names(grouping_df) == "mainsite")] <- "treatmentSt"
colnames(grouping_df)[which(names(grouping_df) == "grouping")] <- "library_id"

cell_clones$library_id <- get_library_labels_v2(cell_clones$cell_id)
cell_clones$sample_id <- get_sample_id(cell_clones$cell_id)
cell_clones <- cell_clones %>% left_join(grouping_df, by = c("library_id","sample_id"))

dim(cell_clones)
data.table::fwrite(cell_clones, paste0('/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_wholedata/prevalences/cell_clones_SA919.csv'))
df <- as.data.frame(table(cell_clones$library_id, cell_clones$sample_id, cell_clones$clone_id))
View(df)
colnames(df) <- c('library_id','sample_id','clone_id','nb_cells')
df <- df %>%
  dplyr::filter(clone_id!='None' & nb_cells>=50)
# df <- df[order(df$clone_id, df$nb_cells),]
df <- df %>% left_join(grouping_df, by = c("library_id","sample_id"))
data.table::fwrite(df, paste0('/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_wholedata/prevalences/SA919_library_info_HoaTran.csv'))

