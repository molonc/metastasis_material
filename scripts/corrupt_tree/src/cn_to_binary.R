suppressPackageStartupMessages({
  require("dplyr")
  require("tidyr")
  require("reshape2")
  require("data.table")
  require("optparse")
})

option_list <- list(make_option(c("-i", "--input"), type="character", default=NULL, help="filtered_cnv_input", metavar="character"),
                    make_option(c("-o", "--outfile"), type="character", default=NULL, help="tip_probability_output", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

parse_bin_names <- function(bin_names) {
  bin_names <- gsub('locus_', '', bin_names)
  chr <- gsub('([0-9]+|X|Y)_[0-9]+_[0-9]+', '\\1', bin_names)
  start <- as.numeric(gsub('([0-9]+|X|Y)_([0-9]+)_[0-9]+', '\\2', bin_names))
  end <- as.numeric(gsub('([0-9]+|X|Y)_([0-9]+)_([0-9]+)', '\\3', bin_names))
  data.frame(chr = chr, start = start, end = end)
}

sort_mat_by_bins <- function(the_mat) {
  # prevent scientific notation
  options(scipen=999)
  
  chr = gsub('([0-9]+|X|Y)_[0-9]+_[0-9]+', '\\1', rownames(the_mat))
  start = as.numeric(gsub('([0-9]+|X|Y)_([0-9]+)_[0-9]+', '\\2', rownames(the_mat)))
  cnv_txt = data.frame(chr=chr, start=start, stringsAsFactors = FALSE)
  the_mat <- cbind(cnv_txt, the_mat)
  
  # Sort the matrix by their chromosome (just inside the chromosome)
  the_mat$chr[the_mat$chr == 'X'] <- '40'
  the_mat$chr[the_mat$chr == 'Y'] <- '55'
  the_mat$chr <- as.numeric(the_mat$chr)
  the_mat <- the_mat[order(the_mat$chr, the_mat$start), ]
  the_mat$chr[the_mat$chr == '40'] <- 'X'
  the_mat$chr[the_mat$chr == '55'] <- 'Y'
  
  # Remove chr, start, end
  the_mat$chr <- NULL
  the_mat$start <- NULL
  the_mat
}

get_single_feature_bins <- function(dat, mat_delta) {
  original_bins <- row.names(dat)
  delta_bin_names <- row.names(mat_delta)
  delta_bin_dat <- data.frame(delta_bin_names = delta_bin_names, left_bin = original_bins[-c(length(original_bins))], right_bin = original_bins[-c(1)], stringsAsFactors = F)
  stopifnot(all(delta_bin_dat$delta_bin_names == delta_bin_dat$left_bin))
  
  cnv_txt_left <- parse_bin_names(delta_bin_dat$left_bin)
  cnv_txt_right <- parse_bin_names(delta_bin_dat$right_bin)
  
  delta_bin_dat$dist <- cnv_txt_right$start -  cnv_txt_left$start
  # Set to inf for those that span multiple chromosomes
  delta_bin_dat$dist <- abs(ifelse(cnv_txt_left$chr == cnv_txt_right$chr, 1, Inf) * delta_bin_dat$dist)
  
  stopifnot(min(delta_bin_dat$dist, na.rm = T) == 500000)
  
  # Keep bins that only cover 500K bases
  single_feature_bins <- delta_bin_dat$delta_bin_names[delta_bin_dat$dist == 500000]
  length(single_feature_bins)
  # nrow(delta_bin_dat)
  single_feature_bins
}

pip_CNV_2_corrupt <- function(input_data, output_file) {
  
  dat <- input_data
  dat <- sort_mat_by_bins(dat)
  
  # Compute diff
  mat_delta <- abs(dat[-c(nrow(dat)), ] - dat[-c(1), ])
  mat_delta[mat_delta > 1] <- 1
  
  # keep_bins <- get_single_feature_bins(dat = dat, mat_delta = mat_delta)
  # mat_delta <- mat_delta[rownames(mat_delta) %in% keep_bins, ]
  binary_mat <- as.matrix(mat_delta)
  mat_delta$loci <- rownames(mat_delta)
  rownames(binary_mat) <- mat_delta$loci
  rownames(mat_delta) <- NULL
  
  # cells,loci,tipInclusionProbabilities
  #write.table(binary_mat, file.path(out_dir, 'merged_sorted_binary_states.csv'), sep = ",", row.names = T, quote = F)
  mat_delta <- reshape2::melt(mat_delta, id.vars = c('loci'), variable.name=c('cells'), value.name='tipInclusionProbabilities')
  mat_delta <- mat_delta[, c(2,1,3)]
  
  data.table::fwrite(mat_delta, output_file, row.names = F, quote = F)
  
}

# Hoa Tran 
# get_noisy_cells <- function(cnv_df, pb=0.9){
#   avg_ploidy <- colMeans(cnv_df)
#   thirdql <- quantile(avg_ploidy,  probs = pb)
#   # avg_ploidy[avg_ploidy>as.numeric(thirdql)]
#   is_avg_ploidy = avg_ploidy > as.numeric(thirdql)
#   
#   change_cnv <- c()
#   for(c in 1:ncol(cnv_df)){
#     count_chg <- 0
#     currv <- cnv_df[1,c]
#     for(r in 2:nrow(cnv_df)){
#       if(cnv_df[r,c]!=currv){
#         currv <- cnv_df[r,c]
#         count_chg <- count_chg + 1
#       }
#     }  
#     change_cnv <- c(change_cnv, count_chg)
#   }
#   # print(change_cnv)
#   is_change <- change_cnv > as.numeric(quantile(change_cnv,  probs = pb))
#   is_noise <- is_avg_ploidy & is_change
#   print(sum(is_noise))
#   cnv_stat <- data.frame(avg_ploidy=avg_ploidy,is_avg_ploidy=is_avg_ploidy,
#                          nb_changes=change_cnv,is_noise=is_noise, row.names = colnames(cnv_df))
#   return(cnv_stat)
# }

# cnv_stat <- get_noisy_cells(filtered_cnv_dat, pb=0.9)
# print(paste0("Before removing noise: ",dim(filtered_cnv_dat)[1]," ",dim(filtered_cnv_dat)[2]))
# cells_use <- colnames(filtered_cnv_dat) %in% rownames(cnv_stat)
# sum(cells_use)
# 
# filtered_cnv_dat <- filtered_cnv_dat[,cells_use]
# filtered_cnv_dat <- filtered_cnv_dat[,!cnv_stat[colnames(filtered_cnv_dat),'is_noise']]
# 
# print(paste0("After removing noise: ",dim(filtered_cnv_dat)[1]," ",dim(filtered_cnv_dat)[2]))
# write.table(filtered_cnv_dat, file = input_fn, quote = F, col.names = NA,)



  # this.dir <- '/home/htran/Projects/hakwoo_project/'
  # setwd(this.dir)
  # data_dir = '/home/htran/storage/datasets/metastasis_results/results/'
  # input_dir <- paste0(data_dir,'X0847_2112251/')
  # output_dir <- paste0(data_dir,'X0847_2112251/')
  # input_fn <- paste0(input_dir,"total_merged_filtered_states.csv")
  # outfile <- paste0(output_dir,"merged_corrupt_no_padding.csv")
  # 
  # filtered_cnv_dat <- read.delim(input_fn, check.names = FALSE, stringsAsFactors = FALSE, sep = ",")
  # filtered_cnv_dat <- as.matrix(filtered_cnv_dat)
  # pip_CNV_2_corrupt(filtered_cnv_dat, outfile)
  
  filtered_cnv_dat <- read.delim(opt$input, check.names = FALSE, stringsAsFactors = FALSE, sep = ",")
  filtered_cnv_dat <- as.matrix(filtered_cnv_dat)
  pip_CNV_2_corrupt(filtered_cnv_dat, opt$outfile)
  
  # filtered_cnv_dat <- read.delim(opt$input, check.names = FALSE, stringsAsFactors = FALSE, sep = ",")
  # filtered_cnv_dat <- as.matrix(filtered_cnv_dat)
  # pip_CNV_2_corrupt(filtered_cnv_dat, opt$outfile)
