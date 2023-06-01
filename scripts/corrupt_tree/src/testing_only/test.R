save_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_rna/clonealign/SA1035X4XB02879/'
gene_cn_summarized <- feather::read_feather(paste0(save_dir,'SA1035X4XB02879.feather'))
colnames(gene_cn_summarized)
fn <- '/home/htran/storage/datasets/drug_resistance_DLP/SA535/A95625A/A95625A_SA535X5XB02891_filtered_states.csv'
fn1 <- '/home/htran/storage/datasets/drug_resistance_DLP/SA535/A95625A/A95625A2891_filtered_states.csv'

data_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535/SA535_total/data/'

lib_A95625A2891 <- read.csv(paste0(data_dir,'A95625A2891_filtered_states.csv'), check.names = F, stringsAsFactors = F)
colnames(lib_A95625A2891)[1:3]
dim(lib_A95625A2891)
lib_A95625A_2498 <- read.csv(paste0(data_dir,'A95625A_filtered_states.csv'), check.names = F, stringsAsFactors = F)
colnames(lib_A95625A_2498)[1:3]
dim(lib_A95625A_2498)

cells_use <- colnames(lib_A95625A)
length(cells_use)
library(stringr)
cellidx <- str_replace_all(cells_use, "A95625A", "A95625A2891")
t <- unique(get_lib_id(colnames(lib_A95625A2891)))
t

colnames(lib_A95625A) <- cellidx
dim(lib_A95625A)
write.csv(lib_A95625A, file = paste0(data_dir,'A95625A2891_filtered_states.csv'), quote=F, row.names = F)


save_dir <- '/home/htran/storage/datasets/hakwoo_metastasis/A98181A/'
cell_predict_df <- read.csv(paste0(save_dir,'cell_predict_df.csv'))
dim(cell_predict_df)


filtered_metrics <- read.csv(paste0(save_dir,'filtered_metrics.csv'))
dim(filtered_metrics)

rv_cells <- setdiff(cell_predict_df$cell_id, filtered_metrics$cell_id)
length(rv_cells) 

sum(rv_cells %in% cell_predict_df$cell_id)


get_sample_id <- function(cell_ids) {
  labels <- sapply(strsplit(cell_ids, "-"), function(x) {
    return(x[1])
  })
  return(labels)
}
get_lib_id <- function(cell_ids) {
  labels <- sapply(strsplit(cell_ids, "-"), function(x) {
    return(x[2])
  })
  return(labels)
}
metric_df <- read.csv(paste0(save_dir,'annotation/A98181A_metrics.csv'))
dim(metric_df)
rownames(metric_df) <- metric_df$cell_id
metric_rv <- metric_df[metric_df$cell_id %in% rv_cells,]
metric_filtered <- metric_df[metric_df$cell_id %in% filtered_cells,]
dim(metric_rv)
unique(metric_rv$experimental_condition)
summary(as.factor(metric_rv$cell_call))
summary(as.factor(metric_filtered$cell_call))

rv_cells[1:5]
filtered_cells[1:5]
srv <- unique(get_sample_id(rv_cells))
sf <- unique(get_sample_id(filtered_cells))

lrv <- unique(get_lib_id(rv_cells))
lf <- unique(get_lib_id(filtered_cells))
summary(metric_rv$quality)


unique(filtered_metrics$cell_call)
unique(cell_predict_df$is_filtered)

summary(as.factor(cell_predict_df$cell_call))
summary(as.factor(filtered_metrics$cell_call))

cell_predict_df$is_filtered <- 'bad_quality'
filtered_cells <- intersect(filtered_metrics$cell_id,cell_predict_df$cell_id)

rownames(cell_predict_df) <- cell_predict_df$cell_id
cell_predict_df[filtered_cells,'is_filtered'] <- 'good_quality'

predict_state <- cell_predict_df
pmap <- plot_scatter_function_plottype(predict_state, xstring, ystring, plottype, 
                                       plottitle=paste0("Mapped read filtering"), 
                                       xlabel=paste0('Total reads'), ylabel='Total mapped reads', 
                                       save_dir, 'mapped_reads_filter', colorcode="#2F4F4F")
xstring <- 'fastqscreen_grch37'
ystring <- 'fastqscreen_mm10'

pmouse <- plot_scatter_function_plottype(predict_state, xstring, ystring, plottype, 
                                       plottitle=paste0("Mapped read filtering"), 
                                       xlabel=paste0('fastqscreen_grch37'), ylabel='fastqscreen_mm10', 
                                       save_dir, 'mouse_detection', colorcode="#2F4F4F")

t <- predict_state[predict_state$cell_id %in% rv_cells,]
View(t[,c('is_s_phase')])
