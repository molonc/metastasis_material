# Detect mouse cells

# Sohrab's version, SA919
results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local/'
input_dir <- '/home/htran/storage/datasets/hakwoo_metastasis/'
cell_clones <- paste0(results_dir,'tree_cut_out/cell_clones_if_0.02_af_0.75_p0.75_e0.04.csv')
series_tag <- 'SA919'

# Tyler left, right padding version, SA919
# cell_clones <- paste0(results_dir, 'cell_clones.csv')
# input_dir <- '/home/htran/storage/datasets/hakwoo_metastasis/'
# results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local_tyler_left_padding/'
# series_tag <- 'SA919'

# Farhia, SA1035
# results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_new_encode/SA1035_whole_local/'
# cell_clones <- paste0(results_dir, 'tree_cut_out/cell_clones_if_0.1_af_0.65_p0.75_e0.04.csv')
# input_dir <- '/home/htran/storage/datasets/drug_resistance_DLP/SA1035/'
# series_tag <- 'SA1035'


source_dir <- '/home/htran/Projects/hakwoo_project/corrupt_tree/src/testing_only/'
source(paste0(source_dir,'cell_feature_utils.R'))

cellclones <- read.csv(cell_clones, check.names = F, stringsAsFactors = FALSE)

predict_state <- read.csv(paste0(input_dir,'cells_features/total_cell_predict_state.csv'), check.names = F, stringsAsFactors = FALSE)

copynumber_orig <- paste0(results_dir, 'total_merged_filtered_states_original.csv')
copy_number_orig <- read.delim(copynumber_orig, check.names = FALSE, stringsAsFactors = FALSE, sep = ",")
total_cells <- colnames(copy_number_orig)
length(total_cells)


dim(predict_state)
copynumber <- paste0(results_dir, 'total_merged_filtered_states.csv')
copy_number <- read.csv(copynumber, header=T, row.names=1, check.names = F,stringsAsFactors = FALSE)
filtered_cells <- colnames(copy_number)
length(filtered_cells)


# 1047 cells, 377 cells without profile 
# 11/ 1047 cells are mouse cells
# outlier_cells <- total_cells[!total_cells %in% cells_use]
# length(outlier_cells)
# predict_state_oc <- predict_state[predict_state$cell_id %in% outlier_cells,]
# dim(predict_state_oc)
# predict_state_oc <- predict_state_oc[(!is.na(predict_state_oc$fastqscreen_grch37) ),]  
# predict_state_oc <- predict_state_oc[(!is.na(predict_state_oc$fastqscreen_mm10) ),]
# is_mouse_cells <- predict_state_oc$fastqscreen_mm10 > predict_state_oc$fastqscreen_grch37
# sum(is_mouse_cells)
# sum(is.na(predict_state_oc$fastqscreen_mm10))


# none_cells <- filtered_cells[!filtered_cells %in% cellclones$cell_id]
# length(none_cells)
# predict_state_nc <- predict_state[predict_state$cell_id %in% none_cells,]
# predict_state_nc <- predict_state_nc[,colnames(predict_state_nc)[2:dim(predict_state_nc)[2]]]
# dim(predict_state_nc)
# sum(predict_state_nc$total_mapped_reads >= 500000)
# is_mouse_cells <- predict_state_nc$fastqscreen_mm10 > predict_state_nc$fastqscreen_grch37
# sum(is_mouse_cells)  # 37 / 609 cells
# sum(predict_state_nc$is_contaminated=="True")  # 287 / 609 cells


# predict_state <- predict_state[,colnames(predict_state)[2:dim(predict_state)[2]]]
predict_state_orig <- predict_state
# predict_state <- predict_state_orig
predict_state <- predict_state[predict_state$cell_id %in% total_cells,]
dim(predict_state)

predict_state$low_mapped_read <- predict_state$total_mapped_reads < 500000
stat <- summary(predict_state$low_mapped_read)

filter_cond <- predict_state$is_mouse_cells | predict_state$low_mapped_read
sum(filter_cond==T)
predict_state <- predict_state[!predict_state$is_mouse_cells,]

predict_state <- predict_state[!predict_state$low_mapped_read,]
copy_number_orig <- copy_number_orig[,colnames(copy_number_orig) %in% predict_state$cell_id]

write.csv(copy_number_orig, paste0(results_dir, 'total_merged_filtered_states_original_v2.csv'),  row.names = T, quote=F)

copy_number <- copy_number[,colnames(copy_number) %in% predict_state$cell_id]
sum((colnames(copy_number) %in% predict_state$cell_id))
write.csv(copy_number, paste0(results_dir, 'total_merged_filtered_states_v2.csv'),  row.names = T, quote=F)

dim(copy_number)
dim(copy_number_orig)
output_file <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Sohrab_filtered/bin_cnvs_corrupt_double_padding_original.csv'
pad_mod_cnv_mat(copy_number_orig, output_file, pad_left = TRUE)

mat <- copy_number_orig
prob = .90
jump_rank <- compute_jump_cells(mat)
# # Filter cells by jump/stuff
bad_cells <- jump_rank %>% dplyr::filter( njumps > quantile(x=njumps, probs = prob) | avgCNA > quantile(x=avgCNA, probs = prob)) %>% dplyr::select(cell_id)

# bad_cells <- jump_rank %>% dplyr::filter( njumps > quantile(x=njumps, probs = prob)) %>% dplyr::select(cell_id)

bad_cells <- bad_cells$cell_id
print("Number of bad cells: ")
print(length(bad_cells))
cc <- mat[, colnames(mat) %ni% bad_cells]
print(paste0("Filtered output: ", dim(cc)[1]," ",dim(cc)[2]))
out_fn <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Sohrab_filtered/total_merged_filtered_states.csv'
write.csv(cc, out_fn,  row.names = T, quote=F)

res_dir2 <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Sohrab_filtered/'
bin_orig <- data.table::fread(paste0(res_dir2,'bin_cnvs_corrupt_double_padding_original.csv'), stringsAsFactors = F)
# get padding, whole dataset
bin_orig <- bin_orig[bin_orig$cells %in% colnames(cc),]
length(unique(bin_orig$cells))

data.table::fwrite(bin_orig, paste0(res_dir2,'bin_cnvs_corrupt_double_padding.csv'), row.names = F, quote = F)  


xstring <- 'total_reads'
ystring <- 'total_mapped_reads'
plottype <- 'low_mapped_read'  # cell clone or sth else
save_dir <- paste0(results_dir,'cell_features/')
pmap <- plot_scatter_function_plottype(predict_state, xstring, ystring, plottype, 
                                       plottitle=paste0("Mapped read filtering"), 
                                       xlabel=paste0('Total reads',' (low=',stat[3],' high=',stat[2],')'), ylabel='Total mapped reads', 
                                       save_dir, 'mapped_reads_filter', colorcode="#2F4F4F")



predict_state <- predict_state[predict_state$total_mapped_reads >= 500000,]
sum(predict_state$total_mapped_reads >= 500000)
summary(predict_state$total_mapped_reads)

sum(is.na(predict_state))


sum(is.na(predict_state$fastqscreen_mm10))
sum(is.na(predict_state$fastqscreen_grch37))
lib_ids <- strsplit(as.character(err_cells_mouse),"-")
length(lib_ids)
nbcores <- 5
ps <- c()
lids <- parallel::mclapply(lib_ids, function(f) {
  ps <- c(ps, as.character(f[2]))
}, mc.cores = nbcores)
print(unique(lids))


predict_state$is_mouse_cells <- (predict_state$fastqscreen_mm10 > predict_state$fastqscreen_grch37)
sum(predict_state$is_mouse_cells) # 1025 cells in total of 5451 cells
stat <- summary(predict_state$is_mouse_cells)
sum(predict_state$is_contaminated=="True") #2196 cells in total of 5451 cells
xl <- c(min(predict_state$fastqscreen_grch37, predict_state$fastqscreen_mm10), max(predict_state$fastqscreen_mm10,predict_state$fastqscreen_grch37))
p <- plot_scatter_function_plottype(predict_state, 'fastqscreen_grch37', 'fastqscreen_mm10', 'is_mouse_cells', 
                                    plottitle=paste0("Detection",'(mouse=',stat[3],' human=',stat[2],')'), 
                                    xlabel='human mapped reads', ylabel='mouse mapped reads', 
                                    save_dir, 'mouse_detection', colorcode="blue", xl=xl, yl=xl) 

xstring <- 'total_reads'
ystring <- 'total_mapped_reads'
plottype <- 'is_mouse_cells'  
save_dir <- paste0(results_dir,'cell_features/')
pm <- plot_scatter_function_plottype(predict_state, xstring, ystring, plottype, 
                                       plottitle=paste0("Mapped read filtering"), 
                                       xlabel=paste0('Total reads',' (mouse=',stat[3],' human=',stat[2],')'), ylabel='Total mapped reads', 
                                       save_dir, 'mouse_detection_reads', colorcode="#2F4F4F")



sum((predict_state$is_contaminated=='True') && (predict_state$is_mouse_cells=TRUE))
xstring <- 'total_reads'
ystring <- 'total_mapped_reads'
plottype <- 'is_contaminated'  
stat <- summary(as.factor(predict_state$is_contaminated))
save_dir <- paste0(results_dir,'cell_features/')
pm <- plot_scatter_function_plottype(predict_state, xstring, ystring, plottype, 
                                     plottitle=paste0("Mapped read filtering"), 
                                     xlabel=paste0('Total reads',' (cont=',stat[2],' not_cont=',stat[1],')'), ylabel='Total mapped reads', 
                                     save_dir, 'contaminate', colorcode="#2F4F4F")

# colnames(predict_state)[1]
predict_state <- predict_state[,colnames(predict_state)[2:dim(predict_state)[2]]]
predict_state <- predict_state %>% inner_join(cellclones, by = "cell_id")
rownames(predict_state) <- predict_state$cell_id
summary(as.factor(predict_state$clone_id))
predict_state$human_subtr_mouse <- as.numeric(predict_state$fastqscreen_grch37) - as.numeric(predict_state$fastqscreen_mm10)
predict_state$human_plus_mouse <- as.numeric(predict_state$fastqscreen_grch37) + as.numeric(predict_state$fastqscreen_mm10)
predict_state$doublet_ratio <- predict_state$human_subtr_mouse / predict_state$human_plus_mouse

xl <- c(min(predict_state$human_subtr_mouse, predict_state$human_plus_mouse), max(predict_state$human_subtr_mouse,predict_state$human_plus_mouse))
p <- plot_scatter_function_plottype(predict_state, 'is_mouse_cells', 'doublet_ratio', 'is_mouse_cells', 
                                    plottitle="Mouse Detection", 
                                    xlabel='is_mouse_cells', ylabel='doublet_ratio', 
                                    paste0(results_dir,'cell_features/'), 'mouse_detection_summary', colorcode="blue") 



# check align 90% 
xstring <- 'total_reads'
ystring <- 'total_mapped_reads'
plottype <- 'is_mouse_cells'  # cell clone or sth else
save_dir <- paste0(results_dir,'cell_features/')
pmap <- plot_scatter_function_plottype(predict_state, xstring, ystring, plottype, 
                                       plottitle=paste0("Mouse Cell Detection"), 
                                       xlabel=paste0('Total reads',' (n=',dim(predict_state)[1],' cells)'), ylabel='Total mapped reads', 
                                       save_dir, 'mouse_detection_mapped_reads', colorcode="#2F4F4F")



# Debug metrics 
metrics_A96117A <- read.csv(paste0(input_dir,'A96117A/obs_cells_features/A96117A_cell_predict_state.csv'), check.names = F, stringsAsFactors = FALSE)
metrics_mA96117A <- read.csv(paste0(input_dir,'A96117A/obs_cells_features/testing/A96117A_metrics.csv'), check.names = F, stringsAsFactors = FALSE)
metrics_mA96117A <- metrics_mA96117A[metrics_mA96117A$cell_id %in% metrics_A96117A$cell_id,]
dim(metrics_A96117A)
dim(metrics_mA96117A)
metrics_A96117A
metrics_A98306A <- read.csv(paste0(input_dir,'A98306A/obs_cells_features/A98306A_cell_predict_state.csv'), check.names = F, stringsAsFactors = FALSE)

dim(metrics_A96117A)
dim(metrics_A98306A)
'total_mapped_reads' %in% colnames(metrics_mA96117A) 
sum(metrics_mA96117A$total_mapped_reads>500000)
'total_mapped_reads' %in% colnames(metrics_A98306A) 
summary(metrics_mA96117A$total_mapped_reads)
metrics_A96117A$total_mapped_reads_hmmcopy
features <- colnames(predict_state)
features[!features %in% colnames(metrics_mA96117A)]
colnames(metrics_mA96117A)[!colnames(metrics_mA96117A) %in% features]
'sample_id' %in% colnames(metrics_mA96117A)
predict_state$`biobloom_GRCh37-lite`

metrics_mA96117A$sample_id <- rep('SA919X7XB05604', length(metrics_mA96117A$cell_id))
metrics_mA96117A$`Unnamed: 0` <- rep('NA', length(metrics_mA96117A$cell_id))
metrics_mA96117A$index_sequence <- rep('NA', length(metrics_mA96117A$cell_id))
metrics_mA96117A$biobloom_noMatch <- rep('NA', length(metrics_mA96117A$cell_id))
metrics_mA96117A$biobloom_GCF_002021735.1_Okis_V1_genomic <- rep('NA', length(metrics_mA96117A$cell_id))
metrics_mA96117A$biobloom_multiMatch <- rep('NA', length(metrics_mA96117A$cell_id))
metrics_mA96117A$biobloom_mm10_build38_mouse <- rep('NA', length(metrics_mA96117A$cell_id))
metrics_mA96117A$`biobloom_GRCh37-lite` <- rep('NA', length(metrics_mA96117A$cell_id))


metrics_mA96117A
write.csv(metrics_mA96117A, paste0(input_dir,'A96117A/obs_cells_features/A96117A_cell_predict_state.csv'), quote=F, row.names=F)
dim(predict_state)
predict_state_backup <- predict_state
predict_state <- predict_state[!predict_state$cell_id %in% metrics_mA96117A$cell_id,]
metrics_mA96117A <- metrics_mA96117A[,colnames(predict_state)]
predict_state <- rbind(predict_state, metrics_mA96117A)
write.csv(predict_state, paste0(input_dir,'cells_features/total_cell_predict_state_v2.csv'), quote=F, row.names=F)
saveRDS(predict_state, file = paste0(input_dir,'cells_features/total_cell_predict_state_v2.rds'))

sum(is.na(predict_state$total_mapped_reads))

metrics_A96200A <- read.csv(paste0(input_dir,'A96200A/obs_cells_features/testing/A96200A_metrics.csv'), check.names = F, stringsAsFactors = FALSE)
metrics_A96204B <- read.csv(paste0(input_dir,'A96204B/obs_cells_features/testing/A96204B_metrics.csv'), check.names = F, stringsAsFactors = FALSE)
sum(is.na(metrics_A96200A))
dim(metrics_A96200A)
dim(metrics_A96204B)


predict_state <- read.csv(paste0(input_dir,'cells_features/total_cell_predict_state_v2.csv'), stringsAsFactors = F, check.names = F)
dim(predict_state)
sum(is.na(metrics_A96204B$fastqscreen_grch37))
sum(is.na(metrics_A96204B$fastqscreen_mm10))
summary(metrics_A96204B$fastqscreen_grch37)
predict_state <- predict_state[!predict_state$cell_id %in% metrics_A96200A$cell_id,]
predict_state <- predict_state[!predict_state$cell_id %in% metrics_A96204B$cell_id,]
c1 <- features[!features %in% colnames(metrics_A96200A)]
c2 <- features[!features %in% colnames(metrics_A96204B)]

metrics_A96200A$sample_id <- rep('SA919X7XB05402', length(metrics_A96200A$cell_id))
metrics_A96200A$`Unnamed: 0` <- rep('NA', length(metrics_A96200A$cell_id))
metrics_A96200A$index_sequence <- rep('NA', length(metrics_A96200A$cell_id))
metrics_A96200A$biobloom_noMatch <- rep('NA', length(metrics_A96200A$cell_id))
metrics_A96200A$biobloom_GCF_002021735.1_Okis_V1_genomic <- rep('NA', length(metrics_A96200A$cell_id))
metrics_A96200A$biobloom_multiMatch <- rep('NA', length(metrics_A96200A$cell_id))
metrics_A96200A$biobloom_mm10_build38_mouse <- rep('NA', length(metrics_A96200A$cell_id))
metrics_A96200A$`biobloom_GRCh37-lite` <- rep('NA', length(metrics_A96200A$cell_id))


metrics_A96204B$sample_id <- rep('SA919X7XB05691', length(metrics_A96204B$cell_id))
metrics_A96204B$`Unnamed: 0` <- rep('NA', length(metrics_A96204B$cell_id))
metrics_A96204B$index_sequence <- rep('NA', length(metrics_A96204B$cell_id))
metrics_A96204B$biobloom_noMatch <- rep('NA', length(metrics_A96204B$cell_id))
metrics_A96204B$biobloom_GCF_002021735.1_Okis_V1_genomic <- rep('NA', length(metrics_A96204B$cell_id))
metrics_A96204B$biobloom_multiMatch <- rep('NA', length(metrics_A96204B$cell_id))
metrics_A96204B$biobloom_mm10_build38_mouse <- rep('NA', length(metrics_A96204B$cell_id))
metrics_A96204B$`biobloom_GRCh37-lite` <- rep('NA', length(metrics_A96204B$cell_id))

predict_state <- rbind(predict_state, metrics_A96200A)
predict_state <- rbind(predict_state, metrics_A96204B)

