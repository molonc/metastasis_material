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
