
# Sohrab's version, SA919
# results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local/'
# input_dir <- '/home/htran/storage/datasets/hakwoo_metastasis/'
# cell_clones <- paste0(results_dir,'tree_cut_out/cell_clones_if_0.02_af_0.75_p0.75_e0.04.csv')
# series_tag <- 'SA919'

# Tyler left, right padding version, SA919
# cell_clones <- paste0(results_dir, 'cell_clones.csv')
# input_dir <- '/home/htran/storage/datasets/hakwoo_metastasis/'
# results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local_tyler_left_padding/'
# series_tag <- 'SA919'

# Farhia, SA1035
results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_new_encode/SA1035_whole_local/'
cell_clones <- paste0(results_dir, 'tree_cut_out/cell_clones_if_0.1_af_0.65_p0.75_e0.04.csv')
input_dir <- '/home/htran/storage/datasets/drug_resistance_DLP/SA1035/'

source_dir <- '/home/htran/Projects/hakwoo_project/corrupt_tree/src/testing_only/'
source(paste0(source_dir,'cell_feature_utils.R'))
series_tag <- 'SA1035'
cellclones <- read.csv(cell_clones, check.names = F, stringsAsFactors = FALSE)
ex_clones <- c('I','J')
cellclones <- cellclones[!cellclones$clone_id %in% ex_clones,]
unique(cellclones$clone_id)
# filtered_metrics <- read.csv(paste0(results_dir,'cell_features/total_filtered_metrics.csv'), check.names = F, stringsAsFactors = FALSE)
# 
# filtered_reads <- read.csv(paste0(results_dir,'cell_features/total_filtered_reads.csv'), check.names = F, stringsAsFactors = FALSE)
predict_state <- read.csv(paste0(input_dir,'cells_features/total_cell_predict_state.csv'), check.names = F, stringsAsFactors = FALSE)
# dim(filtered_metrics)
# features <- colnames(predict_state)
# write(features, paste0(results_dir,'cells_features/features_metrics.txt'))


# dim(filtered_reads)
dim(predict_state)
copynumber <- paste0(results_dir, 'total_merged_filtered_states.csv')
copy_number <- read.csv(copynumber, header=T, row.names=1, check.names = F,stringsAsFactors = FALSE)
cells_use <- colnames(copy_number)
length(cells_use)




predict_state$clone_id <- ifelse(is.na(predict_state$clone_id),'None',predict_state$clone_id)
predict_state <- predict_state[predict_state$cell_id %in% cells_use,]





# sum(predict_state$is_contaminated=="")
predict_state$is_contaminated <- ifelse(predict_state$is_contaminated=="",'Unkn',predict_state$is_contaminated)
predict_state$is_s_phase_prob <- ifelse(is.na(predict_state$is_s_phase_prob),0,predict_state$is_s_phase_prob)
predict_state$total_reads <- ifelse(is.na(predict_state$total_reads),0,predict_state$total_reads)
predict_state$total_mapped_reads <- ifelse(is.na(predict_state$total_mapped_reads),0,predict_state$total_mapped_reads)
predict_state$unmapped_reads <- ifelse(is.na(predict_state$unmapped_reads),0,predict_state$unmapped_reads)
dim(predict_state)
summary(as.factor(predict_state$clone_id))
# predict_state <- predict_state[predict_state$clone_id != 'None',]
obs_clones <- unique(predict_state$clone_id)
obs_clones <- sort(obs_clones)
# obs_clones <- c('A','B','C','D','None')
# obs_clones <- c('A','B','C','D','E','F','G','None')
# obs_clones <- c('A','C','G','B','D','E','F')
# plt_ls <- list()


tr <- c(min(predict_state$total_reads),max(predict_state$total_reads))
tmr <- c(min(predict_state$total_mapped_reads),max(predict_state$total_mapped_reads))
tunmr <- c(min(predict_state$unmapped_reads),max(predict_state$unmapped_reads))
tmedian <- c(min(predict_state$median_hmmcopy_reads_per_bin),max(predict_state$median_hmmcopy_reads_per_bin))
tmean <- c(min(predict_state$mean_copy),max(predict_state$mean_copy))


idx <- 0
pcMapls <- list()
pcUnMapls <- list()
pcQualityls <- list()
pcExpls <- list()
pcSphasels <- list()
pcMeanls <- list()
for(clone_id in obs_clones){
  save_dir <- paste0(results_dir,'cell_features/')
  predict_state_tmp <- predict_state[predict_state$clone_id==clone_id,]
  predict_state_tmp <- droplevels(predict_state_tmp)
  print(dim(predict_state_tmp))
  plts <- plot_cells_features(predict_state_tmp, tr, tmr, tunmr, tmedian, tmean,
                              series_tag='SA1035', clone_id=clone_id, save_dir)
  idx <- idx+1
  pcMapls[[idx]] <- plts$pmap
  pcUnMapls[[idx]] <- plts$punmap
  pcExpls[[idx]] <- plts$pexp
  pcSphasels[[idx]] <- plts$psphase
  pcMeanls[[idx]] <- plts$pquality_mean
  pcQualityls[[idx]] <- plts$pquality
  # plt_ls[[idx]] <- plts
}



# names(plt_ls) <- obs_clones
nc <- 4
pcMap <- plot_grid(plotlist = pcMapls, ncol = nc, align = 'hv')
pcUnMap <- plot_grid(plotlist = pcUnMapls, ncol = nc, align = 'hv')
pcQuality <- plot_grid(plotlist = pcQualityls, ncol = nc, align = 'hv')
pcExp <- plot_grid(plotlist = pcExpls, ncol = nc, align = 'hv')
pcSphase <- plot_grid(plotlist = pcSphasels, ncol = nc, align = 'hv')
pcMean <- plot_grid(plotlist = pcMeanls, ncol = nc, align = 'hv')

clone_ids <- 'SA1035'
plots_features(pcMap, pcUnMap, pcQuality, pcSphase, pcMean, results_dir, clone_ids, ht=380, wd=950)




# pcMap <- plot_grid(plotlist = list(plt_ls[[1]]$pmap, plt_ls[[2]]$pmap, plt_ls[[3]]$pmap), ncol = 3,align = 'hv')
# pcUnMap <- plot_grid(plotlist = list(plt_ls[[1]]$punmap, plt_ls[[2]]$punmap, plt_ls[[3]]$punmap), ncol = 3,align = 'hv')
# pcQuality <- plot_grid(plotlist = list(plt_ls[[1]]$pquality, plt_ls[[2]]$pquality, plt_ls[[3]]$pquality), ncol = 3,align = 'hv')
# pcExp <- plot_grid(plotlist = list(plt_ls[[1]]$pexp, plt_ls[[2]]$pexp, plt_ls[[3]]$pexp), ncol = 3,align = 'hv')
# pcSphase <- plot_grid(plotlist = list(plt_ls[[1]]$psphase, plt_ls[[2]]$psphase, plt_ls[[3]]$psphase), ncol = 3,align = 'hv')
# pcMean <- plot_grid(plotlist = list(plt_ls[[1]]$pquality_mean, plt_ls[[2]]$pquality_mean, plt_ls[[3]]$pquality_mean), ncol = 3,align = 'hv')
# 
# clone_ids <- 'ACG'
# plots_features(pcMap, pcUnMap, pcQuality, pcSphase, pcMean, results_dir, clone_ids, ht=250, wd=760)
# 
# 
# pcMap <- plot_grid(plotlist = list(plt_ls[[4]]$pmap, plt_ls[[5]]$pmap, 
#                                    plt_ls[[6]]$pmap, plt_ls[[7]]$pmap),
#                   ncol = 4,align = 'hv')
# pcUnMap <- plot_grid(plotlist = list(plt_ls[[4]]$punmap, plt_ls[[5]]$punmap, 
#                                      plt_ls[[6]]$punmap, plt_ls[[7]]$punmap),
#                      ncol = 4,align = 'hv')
# 
# pcQuality <- plot_grid(plotlist = list(plt_ls[[4]]$pquality, plt_ls[[5]]$pquality, 
#                                        plt_ls[[6]]$pquality, plt_ls[[7]]$pquality),
#                        ncol = 4,align = 'hv')
# 
# pcExp <- plot_grid(plotlist = list(plt_ls[[4]]$pexp, plt_ls[[5]]$pexp, 
#                                    plt_ls[[6]]$pexp, plt_ls[[7]]$pexp), 
#                   ncol = 4,align = 'hv')
# pcSphase <- plot_grid(plotlist = list(plt_ls[[4]]$psphase, psphaseD, 
#                                       psphaseE, psphaseF), 
#                       ncol = 4,align = 'hv')
# pcMean <- plot_grid(plotlist = list(plt_ls[[4]]$pquality_mean, plt_ls[[5]]$pquality_mean, 
#                                     plt_ls[[6]]$pquality_mean, plt_ls[[7]]$pquality_mean), 
#                     ncol = 4,align = 'hv')

# png(paste0(results_dir,"pcSphase",clone_ids,".png"), height = 2*ht, width=2*wd,res = 2*72)
# print(pcSphase)
# dev.off()

# Get outlier cells 
# Get filtered cells 
# add one more column denote 
