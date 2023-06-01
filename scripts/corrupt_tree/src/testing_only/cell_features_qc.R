results_dir <- '/home/htran/storage/datasets/metastasis_results/SA535_new_encoding/SA535_whole_local_filter09/'
input_dir <- '/home/htran/storage/datasets/hakwoo_metastasis/SA535/'
series_tag <- 'SA1035'


source_dir <- '/home/htran/Projects/farhia_project/rscript/dlp/visualize_tree/'
source(paste0(source_dir,'cell_feature_utils.R'))


save_dir <- paste0(results_dir,'cell_features_qc/')
if(!file.exists(save_dir)){
  dir.create(save_dir)
}
cellclones <- read.csv(cell_clones, check.names = F, stringsAsFactors = FALSE)

predict_state <- read.csv(paste0(input_dir,'cells_features/total_predict_state.csv'), check.names = F, stringsAsFactors = FALSE)
dim(predict_state)
copynumber_orig <- paste0(results_dir, 'total_merged_filtered_states_original.csv')
copy_number_orig <- read.delim(copynumber_orig, check.names = FALSE, stringsAsFactors = FALSE, sep = ",")
total_cells <- colnames(copy_number_orig)
length(total_cells)


copynumber <- paste0(results_dir, 'total_merged_filtered_states.csv')
copy_number <- read.csv(copynumber, header=T, row.names=1, check.names = F,stringsAsFactors = FALSE)
filtered_cells <- colnames(copy_number)
length(filtered_cells)

predict_state <- predict_state[predict_state$cell_id %in% total_cells,]
dim(predict_state)

sum(is.na(predict_state$total_reads))
sum(is.na(predict_state$total_mapped_reads))
sum(is.na(predict_state$is_contaminated))
sum(is.na(predict_state$fastqscreen_mm10))
sum(is.na(predict_state$fastqscreen_grch37))


predict_state$low_mapped_read <- predict_state$total_mapped_reads < 500000
stat <- summary(predict_state$low_mapped_read)
xstring <- 'total_reads'
ystring <- 'total_mapped_reads'
plottype <- 'low_mapped_read'  # cell clone or sth else

# plottype <- 'is_filtered'
# p <- ggplot(predict_state, aes_string(x=xstring, y=ystring)) 
# p
pmap <- plot_scatter_function_plottype(predict_state, xstring, ystring, plottype, 
                                       plottitle=paste0("Mapped read filtering"), 
                                       xlabel=paste0('Total reads',' (low=',stat[3],' high=',stat[2],')'), ylabel='Total mapped reads', 
                                       save_dir, 'mapped_reads_filter', colorcode="#2F4F4F")





xstring <- 'total_reads'
ystring <- 'total_mapped_reads'
plottype <- 'is_contaminated'  
stat <- summary(as.factor(predict_state$is_contaminated))
pm <- plot_scatter_function_plottype(predict_state, xstring, ystring, plottype, 
                                     plottitle=paste0("Mapped read filtering"), 
                                     xlabel=paste0('Total reads',' (cont=',stat[2],' not_cont=',stat[1],')'), ylabel='Total mapped reads', 
                                     save_dir, 'contaminate', colorcode="#2F4F4F")


predict_state$is_mouse_cells <- (predict_state$fastqscreen_mm10 >= predict_state$fastqscreen_grch37)
sum(predict_state$is_mouse_cells) # 1025 cells in total of 5451 cells
stat <- summary(predict_state$is_mouse_cells)

xl <- c(min(predict_state$fastqscreen_grch37, predict_state$fastqscreen_mm10), 
        max(predict_state$fastqscreen_mm10,predict_state$fastqscreen_grch37))
p <- plot_scatter_function_plottype(predict_state, 'fastqscreen_grch37', 'fastqscreen_mm10', 'is_mouse_cells', 
                                    plottitle=paste0("Detection",'(mouse=',stat[3],' human=',stat[2],')'), 
                                    xlabel='human mapped reads', ylabel='mouse mapped reads', 
                                    save_dir, 'mouse_detection', colorcode="blue", xl=xl, yl=xl) 

xstring <- 'total_reads'
ystring <- 'total_mapped_reads'
plottype <- 'is_mouse_cells'  
pm <- plot_scatter_function_plottype(predict_state, xstring, ystring, plottype, 
                                     plottitle=paste0("Mapped read filtering"), 
                                     xlabel=paste0('Total reads',' (mouse=',stat[3],' human=',stat[2],')'), ylabel='Total mapped reads', 
                                     save_dir, 'mouse_detection_reads', colorcode="#2F4F4F")

predict_state <- predict_state[!predict_state$is_mouse_cells,]
# 
predict_state <- predict_state[!predict_state$low_mapped_read,]
# dim(predict_state)


