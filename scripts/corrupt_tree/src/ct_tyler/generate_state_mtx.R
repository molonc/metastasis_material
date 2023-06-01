suppressPackageStartupMessages({
  require("data.table")
  require("dplyr")
  require("optparse")
})

# option_list <- list(make_option(c("-i", "--inputdir"), type="character", default=NULL, help="state_cnv_dir", metavar="character"),
#                     make_option(c("-o", "--outputdir"), type="character", default=NULL, help="output_dir", metavar="character"),
#                     make_option(c("-d", "--outdir"), type="character", default=NULL, help="output_dir", metavar="character")
# )

option_list <- list(make_option(c("-i", "--inputfile"), type="character", default=NULL, help="features_mtx_file", metavar="character"),
                    make_option(c("-o", "--outputfile"), type="character", default=NULL, help="output_file", metavar="character"),
                    make_option(c("-d", "--outdir"), type="character", default=NULL, help="output_dir", metavar="character"),
                    make_option(c("-p", "--probOutlier"), type="character", default='0.9', help="probability to be an outlier", metavar="character")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Tyler left, right padding version, SA919
# input_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local_tyler_left_padding/'
# series_tag <- 'SA919'
# output_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_filtered/'


# SA1035 Tyler
# input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_new_encode/SA1035_filtered/'
# output_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_new_encode/SA1035_Tyler/'
# series_tag <- 'SA1035_Tyler'

# copynumber <- paste0(results_dir, 'total_merged_filtered_states_chr.csv')
# copy_number_chr <- read.csv(copynumber, header=T, check.names = F,stringsAsFactors = FALSE)
# 
# unique(copy_number_chr$width)
# rownames(copy_number_chr)[1:3]
# copy_number_chr$end[1:3]
# copy_number_chr

# input_dir <- opt$inputdir
# output_dir <- opt$outputdir
# copynumber_orig <- paste0(input_dir, 'total_merged_filtered_states_original.csv')
# copy_number <- read.csv(copynumber_orig, header=T, row.names=1, check.names = F,stringsAsFactors = FALSE)
# dim(copy_number)
# View(copy_number[1:3,1:3])
# get_chr_start_end(copy_number, save_dir)

parse_bin_names <- function(bin_names, as_factors = FALSE) {
  # Remove corrupt_tree locus tag if it's there
  bin_names <- gsub('locus_', '', bin_names)
  chr <- gsub('([0-9]+|X|Y)_[0-9]+_[0-9]+', '\\1', bin_names)
  start <- as.numeric(gsub('([0-9]+|X|Y)_([0-9]+)_[0-9]+', '\\2', bin_names))
  end <- as.numeric(gsub('([0-9]+|X|Y)_([0-9]+)_([0-9]+)', '\\3', bin_names))
  data.frame(chr = chr, start = start, end = end, stringsAsFactors = as_factors)
}
### Matrix functions
################################################################################################################################################
sort_mat_by_bins <- function(the_mat) {
  # prevent scientific notation
  options(scipen=999)
  options(stringsAsFactors=FALSE) 
  cnv_txt <- parse_bin_names(rownames(the_mat), as_factors = F)
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
  the_mat$end <- NULL
  
  the_mat
}

# remove odd cells
# Does not count whole-chromosome events
compute_jump_cells <- function(mat) {
  # Cells that have too many CN up and downs may represent douplet? or s-phases
  # Exclude them after a threshold
  # dat <- load_new_cn_data(datatag)
  dat <- sort_mat_by_bins(mat)
  
  # Bins X cells
  mat_delta <- abs(dat[1:(nrow(dat)-1), ] - dat[2:nrow(dat), ])
  mat_delta[mat_delta > 1] <- 1
  njumps <- colSums(mat_delta)
  avgCNA <- colMeans(dat)
  # res <- sort(njumps, decreasing = T)
  stopifnot(colnames(mat_delta) == colnames(dat))
  data.frame(cell_id = colnames(dat), njumps = njumps, avgCNA = avgCNA, stringsAsFactors = F) %>% dplyr::arrange(desc(njumps))
}

# output_dir <- '/home/htran/storage/datasets/metastasis_results/SA535_new_encoding/SA535_wholedata_v2/corrupt_grow/'
get_outliers <- function(output_dir, prob=.90, probAVG=.90){
  cnv_path <- paste0(output_dir,'/total_merged_filtered_states_original.csv')
  mat <- read.delim(cnv_path, check.names = FALSE, stringsAsFactors = FALSE, sep = ",")
  if(prob > 0){
    jump_rank <- compute_jump_cells(mat)
    # # Filter cells by jump/stuff
    bad_cells <- jump_rank %>% dplyr::filter( njumps > quantile(x=njumps, probs = prob) | avgCNA > quantile(x=avgCNA, probs = probAVG)) %>% dplyr::select(cell_id)
    # bad_cells <- jump_rank %>% dplyr::filter( njumps > quantile(x=njumps, probs = prob)) %>% dplyr::select(cell_id)
    
    # Retrieve cells with large nb jumps but not early replicated cells - with low value of CN
    # bad_jump_cells <- jump_rank %>% dplyr::filter( njumps > quantile(x=njumps, probs = prob) & avgCNA < quantile(x=avgCNA, probs = probAVG)) %>% dplyr::select(cell_id)
    
    bad_cells <- bad_cells$cell_id
    print(paste0("Number of bad cells: ",length(bad_cells)))
    if(length(bad_cells)>0){
      cc <- mat[,!colnames(mat) %in% bad_cells]
      outlier_cc <- mat[, colnames(mat) %in% bad_cells]
      print(paste0("Outlier cells mtx: ", dim(outlier_cc)[1]," ",dim(outlier_cc)[2]))
      corrupt_grow_dir <- paste0(output_dir,'/corrupt_grow/')
      if (!file.exists(corrupt_grow_dir)){
        dir.create(corrupt_grow_dir)
      }
      # write.csv(outlier_cc, paste0(corrupt_grow_dir,'total_merged_filtered_states_outlier_cells.csv'),  row.names = T, quote=F) #get dir
      write(colnames(outlier_cc), paste0(corrupt_grow_dir,'outlier_cells.txt'), sep = "\t")
      # write(colnames(outlier_cc), paste0(output_dir,'outlier_jump_cells.txt'), sep = "\t")
    } else{
      cc <- mat
    }
    
    print(paste0("Filtered output mtx: ", dim(cc)[1]," ",dim(cc)[2]))
    write(colnames(cc), paste0(corrupt_grow_dir,'filtered_cells.txt'), sep = "\t")
    data.table::fwrite(cc, paste0(output_dir,'/total_merged_filtered_states.csv'),  row.names = T, quote=F)
  }  
}
get_filtered_mtx <- function(input_fn, output_fn, output_dir, probOutlier=.90){
  get_outliers(output_dir, prob=probOutlier, probAVG=probOutlier)
  # sync_cnv_features mtx
  filter_df <- data.table::fread(input_fn, stringsAsFactors = F)
  print(dim(filter_df))
  
  outliers <- read.table(paste0(output_dir,'/corrupt_grow/outlier_cells.txt'), header = FALSE, sep = "\t")
  filtered <- read.table(paste0(output_dir,'/corrupt_grow/filtered_cells.txt'), header = FALSE, sep = "\t")
  filter_filtered <- filter_df[filter_df$cells %in% paste0('cell_',filtered[,1]),]
  filter_outliers <- filter_df[filter_df$cells %in% paste0('cell_',outliers[,1]),]
  print('Divide filter loci mtx into filtered cells and outlier cells mtx')
  print(dim(filter_filtered))
  print(dim(filter_outliers))
  data.table::fwrite(filter_filtered, output_fn, row.names = F, quote = F)
  data.table::fwrite(filter_outliers, paste0(output_dir,'/corrupt_grow/corrupt_tree_sync-cnv_features_outliers.csv'), row.names = F, quote = F)
  print('Save data, Completed')
}




# copynumber_filtered <- paste0(save_dir, 'total_merged_filtered_states.csv')
# copy_number_fd <- read.csv(copynumber_filtered, header=T, row.names=1, check.names = F,stringsAsFactors = FALSE)
# filtered_cells <- colnames(copy_number_fd)
# length(filtered_cells)
get_filtered_outliers_state_mtx <- function(input_dir, output_dir, copy_number=NULL){
  if(is.null(copy_number)){
    copy_number <- read.csv(paste0(output_dir, 'total_merged_filtered_states_original_chr.csv'),  
                            check.names = F, stringsAsFactors = FALSE)
    
  }
  if(!file.exists(paste0(output_dir, 'corrupt_grow/'))){
    dir.create(paste0(output_dir, 'corrupt_grow/'))
  }
  outliers <- read.table(paste0(input_dir,'corrupt_grow/outlier_cells.txt'), header = FALSE, sep = "\t")
  filtered <- read.table(paste0(input_dir,'corrupt_grow/filtered_cells.txt'), header = FALSE, sep = "\t")
  
  filtered_cells <- filtered[,1]
  outlier_cells <- outliers[,1]
  cells_info <- data.frame(cell_id=c(filtered_cells, outlier_cells),
                           cell_classified=c(rep('filtered',length(filtered_cells)),
                                             rep('outlier',length(outlier_cells))))
  write.csv(cells_info, paste0(output_dir, 'corrupt_grow/cells_info.csv'),  row.names = F, quote=F)
  head(cells_info)
  
  
  
  copy_number_outliers <- copy_number[,!colnames(copy_number) %in% filtered_cells]
  dim(copy_number_outliers)
  # filtered_cells <- c('chr','start','end','width',filtered_cells)
  
  copy_number_filtered <- copy_number[,!colnames(copy_number) %in% outlier_cells]
  dim(copy_number_filtered)
  write.csv(copy_number_filtered, paste0(save_dir, 'total_merged_filtered_states_chr.csv'),  row.names = F, quote=F)
  write.csv(copy_number_outliers, paste0(save_dir, 'corrupt_grow/total_merged_filtered_states_outliers_chr.csv'),  row.names = F, quote=F)
}

# filtered_orig <- paste0(save_dir, 'filtered.csv')
# filtered_df <- read.csv(filtered_orig, header=T, check.names = F,stringsAsFactors = FALSE)
# View(head(filtered_df[1:3,1:3]))
# filtered_df1 <- filtered_df[filtered_df$cells %in% paste0('cell_',filtered_cells),]
# dim(filtered_df1)
# length(unique(filtered_df2$cells))
# write.csv(filtered_df1, paste0(save_dir, 'filtered.csv'),  row.names = F, quote=F)
# 
# filtered_df2 <- filtered_df[filtered_df$cells %in% paste0('cell_',outlier_cells),]
# write.csv(filtered_df2, paste0(save_dir, 'corrupt_grow/filtered_outliers.csv'),  row.names = F, quote=F)
# 
# length(unique(filtered_df2$loci))
# 

print(opt$inputfile)
print(opt$outputfile)
print(opt$outdir)
print(opt$probOutlier)
get_filtered_mtx(opt$inputfile, opt$outputfile, opt$outdir, as.double(opt$probOutlier))


# input_dir <- '/home/htran/storage/datasets/metastasis_results/SA535_new_encoding/SA535_wholedata_v2/'
# 
# outliers <- read.table(paste0(input_dir,'corrupt_grow/outlier_cells.txt'), header = FALSE, sep = "\t")
# filtered <- read.table(paste0(input_dir,'corrupt_grow/filtered_cells.txt'), header = FALSE, sep = "\t")
# 
# filtered_cells <- filtered[,1]
# outlier_cells <- outliers[,1]
# cells_info <- data.frame(cell_id=c(as.character(filtered_cells), as.character(outlier_cells)),
#                          cell_classified=c(rep('filtered',length(filtered_cells)),
#                                            rep('outlier',length(outlier_cells))))
# write.csv(cells_info, paste0(input_dir, 'corrupt_grow/total_filtered_cells.csv'),  row.names = F, quote=F)
# head(cells_info)
# 


# filter_df <- data.table::fread(paste0(output_dir,'corrupt_tree_sync-cnv_features_outliers.csv'), stringsAsFactors = F)
# print(dim(filter_df))
# outliers <- read.table(paste0(output_dir,'outlier_jump_cells.txt'), header = FALSE, sep = "\t")
# data.table::fwrite(filter_outliers, paste0(output_dir,'corrupt_tree_sync-cnv_features_jumps_outliers.csv'), row.names = F, quote = F)

