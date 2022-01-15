input_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA919/clonealign/'

obs_sid <- 'SA919X7XB05692'

raw_segments <- readRDS(paste0(input_dir,'raw_segments.rds'))
colnames(raw_segments)
raw_segments$sample <- get_sample_id(raw_segments$cell_names, cores_use=6)

raw_segments1 <- raw_segments %>%
  dplyr::filter(sample==obs_sid)
dim(raw_segments1)
dim(raw_segments)
saveRDS(raw_segments1,file = paste0(input_dir,obs_sid,'/',obs_sid,'.rds'))


results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_wholedata/'
clones_df <- data.table::fread(paste0(results_dir,'cell_clones.csv'))
dim(clones_df)
clones_df$sample <- get_sample_id(clones_df$cell_id)
obs_sid <- 'SA919X7XB05604'
clones_df <- clones_df %>%
  dplyr::filter(sample==obs_sid)
table(clones_df$sample,clones_df$clone_id)
