suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(scales)
  library(circlize)
  options(dplyr.summarise.inform = FALSE)
  options(tidyverse.quiet = TRUE)
})
# viz_kflux <- function(prop_df){
#   plotttile <- 'Seeding influx estimate kj'
#   prop_df$kflux <- ifelse(is.na(prop_df$kflux),0,prop_df$kflux)
#   prop_df$origin <- factor(prop_df$origin, levels = unique(prop_df$origin))
#   pk <- prop_df %>%
#     ggplot(aes(x=origin, y=kflux)) +
#     geom_bar(stat="identity") + 
#     coord_flip() +
#     # scale_x_reverse() +
#     # scale_x_discrete(position = "top") +
#     theme_bw()
#   pk <- pk + labs(y=plotttile, x='')
#   return(pk)
# }
# get_proportion_mathmodel <- function(meta_data, datatag, save_dir){
#   input_dir <- '/home/htran/Projects/hakwoo_project/rscript/math_model/data_met_proj_v3/'
#   dir.create(input_dir)
#   
#   meta_data <- meta_data %>%
#     dplyr::filter(origin!="Tumor_Recur")
#   
#   dim(meta_data)
#   
#   meta_data$passage <- gsub(datatag,'',meta_data$sample_id)
#   meta_data$passage <- stringr::str_sub(meta_data$passage, 1, 2)
#   unique(meta_data$passage)
#   unique(meta_data$mainsite)
#   # colnames(meta_data)+
#   # summary(as.factor(meta_data$clone_id))
#   # View(head(meta_data))
#   meta_data$origin <- paste0(meta_data$origin,'_',meta_data$passage)  # add passage to tumour origin site, for SA535: only one passage
#   meta_data$origin <- paste0(meta_data$origin,'_',meta_data$pdxid)
#   unique(meta_data$origin)
#   # meta_data$origin <- ifelse(grepl('Primary',meta_data$origin),'Primary', meta_data$origin)    
#   for(m in unique(meta_data$pdxid)){
#     # df <- meta_data %>% 
#     #   tidyr::pivot_wider(id_cols=c('origin'), 
#     #                    names_from='clone_id', values_from = 'origin',
#     #                    values_fn = list(origin = length))
#     # t <- with(meta_data, table(origin, clone_id))
#     # df <- with(meta_data, prop.table(origin, clone_id))
#     # df <- prop.table(meta_data$origin, meta_data$clone_id)
#     df <- meta_data %>%
#       dplyr::filter(pdxid==m) %>%
#       dplyr::group_by(origin,clone_id) %>%
#       dplyr::summarize(sum=n()) %>%
#       tidyr::pivot_wider(id_cols="origin",names_from="clone_id",values_from="sum")%>%
#       tibble::column_to_rownames('origin')
#     # df[is.na(df)] <- 0
#     df[is.na(df)] <- 1
#     
#     df <- df/rowSums(df)
#     # df1[df1==0] <- 0 # just formatting
#     print(dim(df))
#     # 
#     df$clone_id <- rownames(df)
#     df <- df %>%
#       dplyr::select(clone_id, everything())
#     # View(df)
#     # m <- 'total'
#     # data.table::fwrite(df, paste0(input_dir, datatag, '_', m, '_raw.csv'))
#     data.table::fwrite(df, paste0(input_dir, datatag, '_', m, '.csv'))
#     
#     # df <- data.table::fread(paste0(input_dir, datatag, '_', m, '.csv')) %>% as.data.frame()
#     # View(df)
#     # df <- df %>% tibble::column_to_rownames('clone_id')
#   }
#   
# }

viz_heatmap_frequencies <- function(prop_df, datatag, save_dir){
  prop_df1 <- prop_df %>%
    select(-kflux, -m_param) %>%
    tibble::column_to_rownames('origin')
  # View(prop_df1)
  cell_func = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.4f", prop_df1[i, j]), x, y, gp = gpar(fontsize = 10))
  }
  # library(ComplexHeatmap)
  
  col_fun = colorRamp2(c(0, 0.3,0.7, 1), c("#ADDFFF","cyan","#2B65EC", "#0000FF"))
  # col_fun(seq(-3, 3))
  row_ha = rowAnnotation(Seeding_influx_estimate_kj = anno_barplot(prop_df$kflux), width = unit(5, "cm")) #foo2 = prop_df$kflux, 
  col_ha = HeatmapAnnotation(Seeding_freq = anno_boxplot(prop_df1, height = unit(5, "cm"),col_title_rot = 90))
  p <- ComplexHeatmap::Heatmap(as.matrix(prop_df1), na_col = "white",
                               show_column_names=T,
                               show_row_names = T,
                               cluster_rows=F,cluster_columns=F,
                               name = "Clone Fraction", 
                               # row_order = sort(rownames(df)),
                               # row_split= annotation_row,
                               row_title_rot = 90,
                               # row_gap = unit(2, "mm"),
                               # column_split = annotation_col, 
                               # column_title = paste0(datatag, " frequencies"),
                               column_gap = unit(2, "mm"),
                               column_names_gp = grid::gpar(fontsize = 15),
                               row_names_gp = grid::gpar(fontsize = 12),
                               show_heatmap_legend = F,
                               top_annotation=col_ha,
                               # left_annotation = left_anno,
                               right_annotation = row_ha,
                               col = col_fun,
                               # cell_fun = cell_func,
                               row_dend_reorder=F)
  # png(paste0(save_dir,datatag,"_clone_prevalence_hm.png"), height = 2*(30*dim(df)[1]+50), width=2*(50*dim(df)[2]+50), res = 2*72)
  # print(p)
  # dev.off()
  return(p)
  
}


# datatag <- 'SA919' 'SA535'
plot_heatmap_prevalence <- function(datatag, results_dir){
  ## results_dir: 'yourdir/data_met_proj_v4/results/'
  # results_dir <- '/home/htran/storage/datasets/metastasis_results/math_model/data_met_proj_v4/results/'
  base_dir <- dirname(results_dir)
  output_dir <- paste0(base_dir, '/figs/')
  if(!file.exists(output_dir)){ 
    dir.create(output_dir)
  }
  
  fns <- list.files(results_dir)
  fns <- fns[grepl(datatag, fns)]
  print(fns)
  # input_fn <- 'SA919_total.csv'
  # input_fn <- 'SA535_total_raw_6.csv'
  prop_df <- data.table::fread(paste0(results_dir, input_fn)) %>% as.data.frame()
  # View(prop_df)
  # colnames(prop_df)
  if(datatag=='SA919'){
    label <- 'X0847' # make label shorter  
  }else{
    label <- '_X0011_' # make label shorter 
  }
  prop_df <- prop_df %>%
    rename(origin=clone_id) %>%. # the true identity column name
    mutate(origin=gsub(label,'',origin))
    
  phm <- viz_heatmap_frequencies(prop_df, datatag, output_dir)
  # png(paste0(output_dir,gsub('.csv','',input_fn), "_kflux_viz.png"), height = 2*400, width=2*700, res = 2*72)
  # print(phm)
  # dev.off()
  
  
  ggsave(paste0(output_dir,gsub('.csv','',input_fn), "_kflux_viz.png"),
         plot = phm,
         height = 6,
         width = 8,
         # useDingbats=F,
         dpi=250)
  return(phm)
  
}  

datatag <- 'SA919'
plot_heatmap_prevalence(datatag, results_dir)

datatag <- 'SA535'
plot_heatmap_prevalence(datatag, results_dir)