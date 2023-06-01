input_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local_v4/'
results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local_tyler_left_padding/'

source_dir <- '/home/htran/Projects/hakwoo_project/corrupt_tree/src/testing_only/'
source(paste0(source_dir,'cell_feature_utils.R'))
series_tag <- 'SA919'

# cell_clones <- paste0(results_dir, 'cell_clones.csv')
# cellclones <- read.csv(cell_clones, check.names = F, stringsAsFactors = FALSE)
# View(head(cellclones))
#
# filtered_metrics <- read.csv(paste0(results_dir,'cell_features/total_filtered_metrics.csv'), check.names = F, stringsAsFactors = FALSE)

# filtered_reads <- read.csv(paste0(results_dir,'cell_features/total_filtered_reads.csv'), check.names = F, stringsAsFactors = FALSE)
predict_state <- read.csv(paste0(results_dir,'cell_features/total_cell_predict_state.csv'), check.names = F, stringsAsFactors = FALSE)

copynumber <- paste0(results_dir, 'total_merged_filtered_states.csv')
copy_number <- read.csv(copynumber, header=T, row.names=1, check.names = F,stringsAsFactors = FALSE)
filtered_cells <- colnames(copy_number)
length(filtered_cells)

cn_outliers <- paste0(input_dir, 'corrupt_grow/total_merged_filtered_states_outlier_cells.csv')
outliers_df <- read.csv(cn_outliers, header=T, row.names=1, check.names = F,stringsAsFactors = FALSE)
oulier_cells <- colnames(outliers_df)
length(oulier_cells)
rownames(outliers_df)[1:3]

outliers_df <- outliers_df[rownames(copy_number),]

total_states <- cbind(copy_number, outliers_df)
dim(total_states)


sum(oulier_cells %in% predict_state$cell_id) 
cells_use <- c(filtered_cells, oulier_cells)
cell_quality <- c(rep('Filtered',length(filtered_cells)),rep('Outlier',length(oulier_cells),))
quality_df <- data.frame(cell_id=cells_use,cell_classified=cell_quality)
dim(quality_df)

sum(filtered_cells %in% predict_state$cell_id) 
predict_state <- predict_state[,colnames(predict_state)[2:dim(predict_state)[2]]]

predict_state <- predict_state %>% inner_join(quality_df, by = "cell_id")
dim(predict_state)






predict_state$is_contaminated <- ifelse(predict_state$is_contaminated=="",'Uncheck',predict_state$is_contaminated)
predict_state$is_s_phase_prob <- ifelse(is.na(predict_state$is_s_phase_prob),0,predict_state$is_s_phase_prob)
predict_state$total_reads <- ifelse(is.na(predict_state$total_reads),0,predict_state$total_reads)
predict_state$total_mapped_reads <- ifelse(is.na(predict_state$total_mapped_reads),0,predict_state$total_mapped_reads)
predict_state$unmapped_reads <- ifelse(is.na(predict_state$unmapped_reads),0,predict_state$unmapped_reads)

predict_state <- droplevels(predict_state)
save_dir <- paste0(results_dir,'cell_features/outlier_filtered/')
t <- summary(as.factor(predict_state$cell_classified))
class(t)
names(t[1])
tag <- paste0('F_',t[1],' O_',t[2])
plot_outlier_cells_features(predict_state, tag, series_tag='SA919', clone_id='cq', save_dir, save_data=F)

rownames(quality_df) <- quality_df$cell_id
quality_df$cell_classified
vquality <- quality_df[,'cell_classified', drop=F]
print("Plot heatmap")
library(ComplexHeatmap)
cn_colours <- structure(
  c(
    "#3182BD", "#9ECAE1", "#CCCCCC", "#FDCC8A", "#FC8D59", "#E34A33",
    "#B30000", "#980043", "#DD1C77", "#DF65B0", "#C994C7", "#D4B9DA"
  ),
  names=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11+")
)
sum(rownames(copy_number) == rownames(total_states))
sum(rownames(outliers_df) == rownames(total_states))
total_states <- format_copynumber_values(total_states)
outliers_df <- format_copynumber_values(outliers_df)
phm <- Heatmap(name="CN",
             # column_title = "Cell Quality",
             as.matrix(total_states),
             col=cn_colours,
             na_col="white",
             show_column_names=F, 
             show_row_names = F,
             cluster_rows=F,cluster_columns=F,
             column_split=vquality,
             use_raster=TRUE,
             raster_quality=5) 
save_fn <- 'heatmap_total'
png(paste0(save_dir,save_fn,".png"), height = 2*600, width=2*1000,res = 2*72)
print(phm)
dev.off()
saveRDS(phm, file = paste0(save_dir,save_fn, '_plots.rds'))

phm_outliers <- Heatmap(name="CN",
               column_title = "Outlier Cells (pjump>0.9, pmean>0.9)",
               as.matrix(outliers_df),
               col=cn_colours,
               na_col="white",
               show_column_names=F, 
               show_row_names = F,
               cluster_rows=F,cluster_columns=F,
               column_split=quality_df[colnames(outliers_df),'cell_classified', drop=F],
               use_raster=TRUE,
               raster_quality=5) 
# phm_outliers
save_fn <- 'heatmap_outlier'
png(paste0(save_dir,save_fn,".png"), height = 2*600, width=2*1000,res = 2*72)
print(phm_outliers)
dev.off()
saveRDS(phm_outliers, file = paste0(save_dir,save_fn, '_plots.rds'))



# compare filtered sohrab version and tyler version
tyler_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local_tyler/'
s_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local/'

tyler_filtered_df <- read.csv(paste0(tyler_dir, 'filtered.csv'), check.names = F, stringsAsFactors = FALSE)
tyler_average_df <- read.csv(paste0(tyler_dir, 'average.csv'), check.names = F, stringsAsFactors = FALSE)

tyler_sync_df <- read.csv(paste0(tyler_dir, 'corrupt_tree_sync-cnv_features.csv'), check.names = F, stringsAsFactors = FALSE)

s_df <- read.csv(paste0(tyler_dir, 'filtered.csv'), check.names = F, stringsAsFactors = FALSE)

dim(tyler_sync_df)
View(tyler_sync_df[1:3,1:3])
View(copy_number[1:3,1:3])
copynumber <- paste0(tyler_dir, 'bin_cnvs_corrupt_double_padding.csv')
copy_number <- read.csv(copynumber, header=T, check.names = F,stringsAsFactors = FALSE)
loci_ls <- unique(tyler_filtered_df$loci)
length(loci_ls)
sum(loci_ls %in% paste0('locus_',rownames(copy_number)))




