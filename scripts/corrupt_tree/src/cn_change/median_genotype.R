
library(dplyr)
library(tidyr)
library(ComplexHeatmap)

# mg_mat <- median_cnv
# results_dir <- save_dir
compute_dist_mat <- function(mg_mat, results_dir, use_hamming = FALSE) {
  # print("Testing")
  if(is.null(mg_mat)){
    # print("Testing 1")
    # return(NULL)
    stop('Check input data')
  } else if(nrow(mg_mat)==1){
    # print("Testing 2")
    # return(matrix(0, nrow = 1, ncol = 1))
    stop('Check input data')
  } else{
    # print("Testing 3")
    clone_lbs <- colnames(mg_mat)
    out_mtx <- tibble::tibble()
    n_clones <- ncol(mg_mat)
    dist_to_median <- matrix(NA, nrow = n_clones, ncol = n_clones)
    for (j in seq(n_clones)) {
      for (i in seq(n_clones)) {
        if (i<j) {  #(i + j)<=(n_clones+1)
          if (use_hamming) {
            distance_type='Hamming'
            dist_to_median[j,i] <- mean(mg_mat[ ,i] != mg_mat[ ,j])
          } else {
            distance_type='Manhattan'
            dist_to_median[j,i] <- mean(abs(mg_mat[ ,i] - mg_mat[ ,j])) #The Manhattan distance as the sum of absolute differences
          }
          distji <- c(clone_lbs[j],clone_lbs[i],round(dist_to_median[j,i],3))
          names(distji) <- c('SourceClone','TargetClone','CNA_Distance')
          out_mtx <- dplyr::bind_rows(out_mtx,distji)
        }
      }
    }
    rownames(dist_to_median) <- colnames(mg_mat)
    colnames(dist_to_median) <- colnames(mg_mat)
    data.table::fwrite(as.data.frame(dist_to_median), paste0(results_dir,'cn_distance_',distance_type,'.csv'), row.names=T)
    data.table::fwrite(out_mtx, paste0(results_dir,'cn_distance_',distance_type,'_output.csv'))
    return(list(dist_to_median=dist_to_median,out_mtx=out_mtx))
  }
  
}
get_chr_infos <- function(chr_desc) {
  genomic_regions <- strsplit(chr_desc, "_")
  chrs <- sapply(genomic_regions, function(x) {
    return(x[1])
  })
  starts <- sapply(genomic_regions, function(x) {
    return(x[2])
  })
  ends <- sapply(genomic_regions, function(x) {
    return(x[3])
  })
  return(list(chr=as.character(chrs),start=as.character(starts),end=as.character(ends)))
}
# dis_mtx is a dataframe with rownames and colnames are clone labels
# distance_type: Manhattan or Hamming
# dis_mtx <- res$dist_to_median
plot_heatmap_genotype <- function(dis_mtx, distance_type, datatag, results_dir){
  cell_func = function(j, i, x, y, width, height, fill) {
    str <- ''
    if(!is.na(dis_mtx[i, j])){
      str <- sprintf("%.2f", dis_mtx[i, j])
    }
    grid.text(str, x, y, gp = gpar(fontsize = 10))
  }
  # dis_mtx[dis_mtx==0]
  # View(dis_mtx)
  row_lbs <- as.character(rownames(dis_mtx))
  # distance between median CNA profiles of clones
  plottitle <- paste0(distance_type," copy number distance ",datatag)
  p <- ComplexHeatmap::Heatmap(as.matrix(dis_mtx), na_col = "white",
                               show_column_names=T,
                               show_row_names = T,
                               cluster_rows=F,
                               cluster_columns=F,
                               name = paste0(distance_type," dist"), 
                               # row_order = sort(rownames(test)),
                               row_split=row_lbs,
                               # row_title = "%s",
                               # column_title_side = T,
                               # row_title_rot = 90,
                               row_names_side = "left",
                               column_names_side = "bottom",
                               column_names_rot = 0,
                               row_gap = unit(1, "mm"),
                               column_split = row_lbs,
                               row_title = "Clones",
                               column_title = plottitle,
                               column_gap = unit(1, "mm"),
                               column_names_gp = grid::gpar(fontsize = 13, fontface = "bold"),
                               column_title_gp = gpar(fontsize = 13, fontface = "bold"),
                               row_names_gp = grid::gpar(fontsize = 13, fontface = "bold"),
                               show_heatmap_legend = T,
                               # top_annotation=top_anno,
                               # left_annotation = left_anno,
                               cell_fun = cell_func,
  )#row_dend_reorder=F
  # p
  
  png(paste0(results_dir, datatag, distance_type,'_distance_','.png'), height = 2*500, width=2*600, res = 2*72)
  print(p)
  dev.off()
  return(p)
}

# Get mode values
calc_mode <- function(x) {
  keys <- unique(x)
  keys[which.max(tabulate(match(x, keys)))]
}


# Input data: copy number matrix, cells names in colnames, chr labels in rownames, values are copy number states
# Ex:            c1 c2 c3
# chr1_1_50000   1  3  2
source("/home/htran/Projects/hakwoo_project/corrupt_tree/src/cn_change/utils.R")

results_dir <- '/home/htran/storage/datasets/metastasis_results/SA535_new_encoding/SA535_wholedata_v2/'
save_dir <- paste0(results_dir,'CN_profile/')

copynumber_fn <- paste0(results_dir,'total_merged_filtered_states.csv')
datatag <- 'SA535'

get_median_genotype_v3(copynumber_fn, 
                       datatag, save_dir,
                       cellclone_fn=NULL, library_grouping_fn=NULL)
  


results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_wholedata/'
save_dir <- paste0(results_dir,'CN_profile/')
dir.create(save_dir)
copynumber_fn <- paste0(results_dir,'total_merged_filtered_states.csv')
datatag <- 'SA919'
cellclone_fn <- paste0(results_dir,'cell_clones.csv')
library_grouping_fn <- paste0(results_dir,'library_groupings.csv')
get_median_genotype_v3(copynumber_fn, 
                       datatag, save_dir,
                       cellclone_fn=NULL, library_grouping_fn=NULL)



results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA501_v2/'
datatag <- 'SA501'
res_501 <- get_median_genotype_v2(datatag, results_dir, copynumber_fn=NULL, cellclone_fn=NULL,
                      calcul_distance=T, distance_type='Manhattan')

saveRDS(res_501, paste0(save_dir, datatag, '_results_CNA_distance.rds'))
min(res_501$out_mtx$CNA_Distance)
max(res_501$out_mtx$CNA_Distance)
results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA604_v2/'
datatag <- 'SA604'
res_SA604 <- get_median_genotype_v2(datatag, results_dir, copynumber_fn=NULL, cellclone_fn=NULL,
                                  calcul_distance=T, distance_type='Manhattan')
saveRDS(res_SA604, paste0(save_dir, datatag, '_results_CNA_distance.rds'))
min(res_SA604$out_mtx$CNA_Distance)
max(res_SA604$out_mtx$CNA_Distance)

results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA530_v2/'
datatag <- 'SA530'
res_SA530 <- get_median_genotype_v2(datatag, results_dir, copynumber_fn=NULL, cellclone_fn=NULL,
                                    calcul_distance=T, distance_type='Manhattan')
saveRDS(res_SA530, paste0(save_dir, datatag, '_results_CNA_distance.rds'))



results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_new_encode/SA1035_Tyler_v2/'
datatag <- 'SA1035'
res_SA1035 <- get_median_genotype_v2(datatag, results_dir, copynumber_fn=NULL, cellclone_fn=NULL,
                                    calcul_distance=T, distance_type='Manhattan')
saveRDS(res_SA1035, paste0(save_dir, datatag, '_results_CNA_distance.rds'))
min(res_SA1035$out_mtx$CNA_Distance)
max(res_SA1035$out_mtx$CNA_Distance)


# chrs <- data.table::fread('/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_new_encode/SA1035_Tyler_v2/CN_profile/median_cnv.csv')
# length(unique(chrs$chr_desc))
results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA609/'
datatag <- 'SA609'
base_dir <- '/home/htran/storage/raw_DLP/drug_resistance_DLP/SA609/reads_clones/'
df <- data.table::fread(paste0(base_dir,'SA609_cisplatin_segments_total.csv')) #filtered_reads_total.csv.gz
dim(df)
head(df)
df$chr_desc <- paste0(df$chr,'_',df$start,'_',df$end)
df$median_cn <- df$state
df$clone_id <- gsub('clone_','',df$clone_id)
res_SA609 <- get_median_genotype_v2(datatag, results_dir, copynumber_fn=NULL, cellclone_fn=NULL,
                                     calcul_distance=T, distance_type='Manhattan')
saveRDS(res_SA609, paste0(save_dir, datatag, '_results_CNA_distance.rds'))
res_SA609 <- res
min(res_SA609$out_mtx$CNA_Distance)
max(res_SA609$out_mtx$CNA_Distance)


results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535/SA535_fitness/'
datatag <- 'SA535'
res_SA535 <- get_median_genotype_v2(datatag, results_dir, copynumber_fn=NULL, cellclone_fn=NULL,
                                    calcul_distance=T, distance_type='Manhattan')
saveRDS(res_SA535, paste0(save_dir, datatag, '_results_CNA_distance.rds'))

res_501$out_mtx$datatag <- 'SA501'
res_501$out_mtx$treatment_status <- 'UnRx'
res_SA530$out_mtx$datatag <- 'SA530'
res_SA530$out_mtx$treatment_status <- 'UnRx'
res_SA604$out_mtx$datatag <- 'SA604'
res_SA604$out_mtx$treatment_status <- 'UnRx'

res_SA609$out_mtx$datatag <- 'SA609'
res_SA609$out_mtx$treatment_status <- 'Rx_RxH_UnRx'
res_SA535$out_mtx$datatag <- 'SA535'
res_SA535$out_mtx$treatment_status <- 'Rx_RxH_UnRx'
res_SA1035$out_mtx$datatag <- 'SA1035'
res_SA1035$out_mtx$treatment_status <- 'Rx_RxH_UnRx'

stat <- tibble::tibble()
stat <- dplyr::bind_rows(list(res_501$out_mtx, res_SA530$out_mtx, res_SA604$out_mtx,
                              res_SA609$out_mtx, res_SA535$out_mtx, res_SA1035$out_mtx))
dim(stat)
head(stat)
unique(stat$treatment_status)
data.table::fwrite(stat, paste0(save_dir,'total_stat_median_CN_distance_6series.csv'))
save_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA609/CN_profile/'
stat <- data.table::fread(paste0(save_dir,'total_stat_median_CN_distance_6series.csv')) %>% as.data.frame()
stat$CNA_Distance <- as.numeric(stat$CNA_Distance)
series <- c('SA501','SA530','SA604','SA609','SA535','SA1035')
stat$datatag <- factor(stat$datatag, levels = series)
p <- ggplot(stat, aes(x = datatag, y=CNA_Distance, color = datatag)) +
     geom_boxplot() + #alpha = 0.6, size = 0.3
     facet_wrap(~factor(treatment_status, levels = c('UnRx','Rx_RxH_UnRx')), scales = "free_x",
           strip.position = "bottom",
           # switch = "x",
           nrow = 1, drop = T) + 
  theme_bw()

png(paste0(save_dir,"median_CN_distance_6series.png"), height = 2*300, width=2*500, res = 2*72)
print(p)
dev.off()

# hm_ls <- res_501$hm + res_SA530$hm + res_SA604$hm + res_SA609$hm + res_SA535$hm + res_SA1035$hm
stat1 <- stat %>%
  # select(-SourceClone,-TargetClone)%>%
  group_by(datatag)%>%
  summarise(mean_distance=round(mean(CNA_Distance),2),
            sd=round(sd(CNA_Distance),2))
stat1
series <- c('SA501','SA530','SA604','SA609','SA535','SA1035')
stat1 <- as.data.frame(stat1)
stat1$datatag <- factor(stat1$datatag, levels = series)
rownames(stat1) <- stat1$datatag
stat1 <- stat1[series,]
for(d in stat1$datatag){
  print(paste0('Series: ',d,'; mean CNA distance: ',stat1[d,'mean_distance'], '; sd CNA distance: ',stat1[d,'sd']))
}

stat2 <- stat1 %>%
  filter(datatag %in% c('SA609','SA535','SA1035'))
View(stat2)
stat2$idx <- paste0(stat2$datatag, stat2$SourceClone, stat2$TargetClone)
stat3 <- stat2 %>%
  filter(idx %in% c('SA609AH','SA609HA','SA535AG','SA535GA','SA1035HE','SA1035EH'))
View(stat3)

