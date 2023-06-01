
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
}) 


viz_proportion_s_phase_per_clone <- function(datatag){
  
}

viz_rep_bk_per_clone <- function(datatag){
  out_bkmetric_fn <- paste0(results_dir, datatag,'_RT_cell_metrics_clones.csv.gz')
  bk_metric_df <- data.table::fread(out_bkmetric_fn)
  dim(bk_metric_df)
  head(bk_metric_df)
  clones_g1 <- c('J','P','T')
  clones_g2 <- c('N','S','I','K','L','M','Q','R')
  bk_metric_df <- bk_metric_df %>%
    # dplyr::filter(rep_bk>=50) %>%
    dplyr::mutate(
    clone_type = case_when(
      clone_id %in% clones_g1 ~ "increase_in_met",
      clone_id %in% clones_g2 ~ "decrease_in_met",
      TRUE ~ "other"
    ))
  bk_metric_df$rep_bk <- ifelse(bk_metric_df$rep_bk>600, 600, bk_metric_df$rep_bk)
  # bk_metric_df <- bk_metric_df %>%
  #   # dplyr::filter(rep_bk>=50) %>%
  #   dplyr::mutate(rep_bk1 = case_when(
  #                 rep_bk > 600 ~ 600,
  #                   TRUE ~ rep_bk
  #                 )
  #   )
      
  p <- ggplot(data=bk_metric_df, aes(x=clone_id, y=rep_bk, color=clone_id)) +
    geom_jitter(position=position_jitter(0.2)) +
    stat_summary(fun=median, geom="point", shape=18,
                 size=1.5, color="black")+
    facet_wrap(cell_type_status ~ clone_type , strip.position = "top", scales = "free_x", drop = F)
  p
  
  
}
viz_rep_fraction_per_clone <- function(datatag){
  out_metric_fn <- paste0(results_dir, datatag,'_RT_cell_metrics_RT_frac_clones.csv.gz')
  frac_metric_df <- data.table::fread(out_metric_fn)
  dim(frac_metric_df)
  head(frac_metric_df)
  clones_g1 <- c('J','P','T')
  clones_g2 <- c('N','S','I','K','L','M','Q','R')
  frac_metric_df <- frac_metric_df %>%
    # dplyr::filter(cell_frac_rep>=0.025 & cell_frac_rep<=0.975) %>%
    dplyr::mutate(
      clone_type = case_when(
        clone_id %in% clones_g1 ~ "slow_grow",
        clone_id %in% clones_g2 ~ "fast_grow",
        TRUE ~ "other"
      )
    )
  p <- ggplot(data=frac_metric_df, aes(x=clone_id, y=cell_frac_rep, color=clone_id)) +
    geom_jitter(position=position_jitter(0.2)) +
    stat_summary(fun=median, geom="point", shape=18,
                 size=1.5, color="black")+
    facet_wrap(cell_type_status ~ clone_type , strip.position = "top", scales = "free_x", drop = F)
  p
}

load_data <- function(datatag, results_dir){
 
  
}

results_dir <- '/home/htran/storage/datasets/metastasis_results/replication_timing/SA535X4XB05649_RT_results/'
datatag <- 'SA535X4XB05649'

results_dir <- '/home/htran/storage/datasets/metastasis_results/replication_timing/RT_results/temp/'
datatag <- ''
