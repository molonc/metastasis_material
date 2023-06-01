
input_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_filtered_v3/'
output_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_filtered/'

filter_df <- data.table::fread(paste0(input_dir,'corrupt_tree_sync-cnv_features_original.csv'), stringsAsFactors = F)
print(dim(filter_df))


# Get data 2161
s2161 <- read.csv(paste0(input_dir,'A108836B_filtered_states.csv'), stringsAsFactors = F, check.names = F)
dim(s2161)
cells_use <- colnames(s2161)
s2161_filtered <- filter_df[filter_df$cells %in% paste0('cell_',cells_use),]
dim(s2161_filtered)
length(unique(s2161_filtered$cells))

filtered <- read.csv(paste0(input_dir,'filtered.csv'))
dim(filtered)
colnames(filtered)
loci_use <- unique(filtered$loci)
length(loci_use)
s2161_filtered <- s2161_filtered[s2161_filtered$loci %in% loci_use,]
length(unique(s2161_filtered$loci))
data.table::fwrite(s2161_filtered, paste0(output_dir,'corrupt_grow/corrupt_tree_sync-cnv_features_2161.csv'), row.names = F, quote = F)


filtered <- read.table(paste0(output_dir,'/corrupt_grow/filtered_cells.txt'), header = FALSE, sep = "\t")
filter_filtered <- filter_df[filter_df$cells %in% paste0('cell_',filtered[,1]),]
filter_outliers <- filter_df[filter_df$cells %in% paste0('cell_',outliers[,1]),]
print('Divide filter loci mtx into filtered cells and outlier cells mtx')
print(dim(filter_filtered))
print(dim(filter_outliers))
data.table::fwrite(filter_filtered, output_fn, row.names = F, quote = F)
data.table::fwrite(filter_outliers, paste0(output_dir,'/corrupt_grow/corrupt_tree_sync-cnv_features_outliers.csv'), row.names = F, quote = F)
print('Save data, Completed')