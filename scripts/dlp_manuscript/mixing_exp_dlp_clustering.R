library(dplyr)

# Get filtered files, get list of library ids
# If lack of files --> download data from azure 
# Filter data 
# Get # cells, put into table
# Run hdbscan for combined files

script_dir <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/dlp_trees/SA919_mixing_experiment/'

meta_samples <- data.table::fread(paste0(script_dir, 'Metastasis_Hakwoo_mixing_exp_SA919_results.csv'))
dim(meta_samples)
colnames(meta_samples)
meta_samples$jira_ticket <- NULL
length(unique(meta_samples$library_id))
tks <- data.table::fread(paste0(script_dir, 'jira_tickets_library_grouping.csv'))
dim(tks)
colnames(tks)
tks <- tks %>%
  dplyr::select(library_id, jira_ticket)
length(unique(tks$library_id))
meta_samples <- meta_samples %>%
  dplyr::left_join(tks, by=c('library_id'))
data.table::fwrite(meta_samples, paste0(script_dir, 'mixing_exp_SA919_results_total.csv'))

meta_samples <- meta_samples %>%
  dplyr::select(library_id,jira_ticket)
write.table(meta_samples, "/home/htran/Projects/hakwoo_project/metastasis_material/scripts/corrupt_tree/src/downloading_scripts/jira_tickets_mixing_exp.txt",
            row.names=F, sep = ',', col.names = F, quote=F)


meta_samples$AT_ID
meta_samples <- meta_samples %>%
  dplyr::rename(sample_id=AT_ID)
save_dir = '/home/htran/storage/datasets/metastasis_results/dlp_SA919_mixing_exp/'
fns <- list.files(save_dir)
fns <- fns[grepl('_filtered_states.csv.gz',fns)]
fns <- gsub('_filtered_states.csv.gz','',fns)
meta_samples <- meta_samples %>%
  dplyr::filter(!library_id %in% fns)
dim(meta_samples)
meta_samples <- meta_samples %>%
  dplyr::select(library_id,sample_id)
data.table::fwrite(meta_samples, paste0(script_dir, 'library_groupings_extra.csv'))
