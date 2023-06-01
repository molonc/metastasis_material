suppressPackageStartupMessages({
  # require("tools")
  # require("phytools")
  # require("igraph")
  # require("RColorBrewer")
  require("outliers")
  require("data.table")
  require("dplyr")
  require("optparse")
  require('tidyr')
  require('tibble')
  require('ggplot2')
  require('cowplot')
  # require('ggtree')
  # library("yarrr")
  # library('ggrepel')
})

## Original script: Sohrab Salehi 
## Modified: Hoa Tran
## Citation: Please cite this manuscript if you use this script: https://github.com/UBC-Stat-ML/sitkatree


## Note: improve the encoding output, carefully selecting genomic bins
## We're missing synch events across the genome, i.e., those that happen in multiple cells, in more than one bin. Here we add the padding at the left
## to take into account the synch events.


# input : the dir+ file name for copy number file: ex: mydir/total_merged_filtered_states.csv, a matrix: cells in columns, 
# chr_start_end position in row names 
## ex:            AT128-A98207B-R65-C38,AT128-A98207B-R65-C53
# 1_2000001_2500000   1, 2
# 1_3000001_350000    3, 3

# output filename: name of binary encoding output file name 
# f: filter threshold: ex: f=0.9, remove cells with large total copy number values than 90% of cells in population, 
# remove cells with number of copy number change - jumps larger than 90% of number of CN change in population

# example run: 
# Rscript blob_cn_to_binary.R -i mydir/total_merged_filtered_states.csv -o mydir/binary_encoded.csv -f 0.9

option_list <- list(make_option(c("-i", "--input"), type="character", 
                                default=NULL, help="filtered_cnv_input", metavar="character"),
                    make_option(c("-o", "--outfile"), 
                                type="character", default=NULL, help="binary encoding output file name", metavar="character"),
                    make_option(c("-f", "--filterthrs"), type="character", 
                                default='0', help="outlier_threshold jumps, value from 0 to 1, 0 means take into account total cells", metavar="character")
)




opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


## Constants
# 1. Check if chromosome linings are included X
# 2. Remove the bins spanning multiple bins
get_single_feature_bins <- function(dat, mat_delta) {
  #original_bins <- rownames(sort_mat_by_bins(load_new_cn_data(datatag))) 
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


get_pretty_names_for_loci <- function(str_array) {
  str_array <- gsub('locus_', '', str_array)
  chr = gsub('([0-9]+|X|Y)_.*', '\\1', str_array)
  str_array <- gsub('.*_([0-9]+)_.*', '\\1', str_array)
  str_array <- substr(str_array, 1, 4)
  str_array <- paste0('chr', chr, '_', str_array, '')
  gsub('chrroot_root', 'root', str_array)
}

addrow <- function(orig, tmp) {
  if (is.null(orig)) {
    orig <- tmp
  } else {
    orig <- rbind(orig, tmp)
  }
  orig
}


format_copynumber_matrix <- function(copynumber) {
  # bin_ids <- paste(
  #   copynumber$chr, copynumber$start, copynumber$end, sep="_")
  
  # rownames(copynumber) <- bin_ids
  # copynumber <- subset(copynumber, select=-c(chr, start, end, width))
  copynumber <- as.matrix(copynumber)
  
  return(copynumber)
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

remove_pad_bins <- function(the_mat) {
  pad_index <- get_pad_bin_index_for_mat(the_mat = the_mat)
  if (length(pad_index) > 0)
    the_mat <- the_mat[-pad_index, ]
  
  the_mat
}

get_pad_bin_index_for_mat <- function(the_mat) {
  dd <- dim(the_mat)
  if (max(dd)/min(dd) > 1000) {
    bin_names <- unique(the_mat$loci)
  } else {
    bin_names <- rownames(the_mat)
  }
  
  get_pad_bin_index(bin_names)
}

get_pad_bin_index <- function(bin_names) {
  bin_dat <- parse_bin_names(bin_names)
  which(abs(bin_dat$start - bin_dat$end) == 1)
}


data_frame_to_list <- function(clustering) {
  clusts <- unique(clustering$clone_id)
  res <- list()
  for (cl in clusts) {
    res[[cl]] <- clustering$cell_id[clustering$clone_id == cl]
  }
  res
}


prepare_mat_by_chr <- function(the_mat, update_col_names=TRUE) {
  the_array = colnames(the_mat)
  chunks <- unlist(strsplit(the_array, '_'))
  chrs <- matrix(chunks, nrow=length(the_array), ncol=3, byrow = T)[, 1]
  unique(chrs)
  the_mat <- the_mat[, order_by_chr_name(chrs)]
  if (update_col_names) {
    colnames(the_mat) <- strip_chr_names(the_array = colnames(the_mat))
  } 
  
  the_mat
}

### Util_utils
################################################################################################################################################

get_ploidy_for_mat <- function(the_mat) {
  dmat <- wide_to_tidy_for_mat(the_mat)
  pa <- dmat %>% dplyr::count(cell_names, copy_number) %>% dplyr::group_by(cell_names) %>% dplyr::filter(n == max(n)) %>% dplyr::rename(ploidy = copy_number) %>% dplyr::select(cell_names, ploidy) %>% as.data.frame()
  nrow(pa)
  dim(the_mat)
  pa$cell_names <- as.character(pa$cell_names)
  pa
}

wide_to_tidy_for_mat <- function(mat, chr_filter = NULL) {
  cnv_txt <- parse_bin_names(rownames(mat))
  mat <- cbind(cnv_txt, mat)
  if (!is.null(chr_filter))
    mat <- mat[mat$chr %in% chr_filter, ]
  reshape2::melt(mat, id.vars = c('chr', 'start', 'end'), value.name='copy_number', variable.name = 'cell_names')
}


get_exp_path_for_datatag <- function(dt) {
  sa609_config <- load_main_configs(datatag = dt)
  dt_batch_path <- file.path('~/Desktop/SC-1311', dt, 'batch_runs', sa609_config$all_cells$batch_name)
  file.path(dt_batch_path, 'outputs', sa609_config$all_cells$exp_dir)
}


get_libid_from_cell_names <- function(cell_names) {
  # TODO: add the condition for SA906 and SA666
  gsub('SA([0-9]|[A-Z]|[a-z])+-(A([0-9]|[A-Z])+)-.*', '\\2', cell_names)
}



collapse_bin_names <- function(bin_dat) {
  paste0(bin_dat$chr, '_', bin_dat$start, '_', bin_dat$end)
}


parse_bin_names <- function(bin_names, as_factors = FALSE) {
  # Remove corrupt_tree locus tag if it's there
  bin_names <- gsub('locus_', '', bin_names)
  chr <- gsub('([0-9]+|X|Y)_[0-9]+_[0-9]+', '\\1', bin_names)
  start <- as.numeric(gsub('([0-9]+|X|Y)_([0-9]+)_[0-9]+', '\\2', bin_names))
  end <- as.numeric(gsub('([0-9]+|X|Y)_([0-9]+)_([0-9]+)', '\\3', bin_names))
  data.frame(chr = chr, start = start, end = end, stringsAsFactors = as_factors)
}

get_desc_names <- function(datatag, the_graph, the_node) {
  desc <- get_decendents(the_node, the_graph = the_graph)
  remove_na_from_list(node_number_to_label(the_graph, desc))
}


### Misc
################################################################################################################################################
remove_na_from_list <- function(some_list) {
  # Why doesn't it keep names?
  tmp <- lapply(names(some_list), function(x) {some_list[[x]][!is.na(some_list[[x]])] })
  names(tmp) <- names(some_list)
  tmp
}


convert_to_freq <- function(time_series_dat, total_cell_counts_breakdown) {
  if (is.null(total_cell_counts_breakdown)) {
    print('Normalising to themselves')
    for (i in 2:ncol(time_series_dat)) {
      the_sum <- sum(time_series_dat[, i])
      if (the_sum > 0)
        time_series_dat[, i] <- time_series_dat[, i]/the_sum
    }
  } else {
    print('Normalising to Total count')
    for (i in 2:ncol(time_series_dat)) { 
      time_series_dat[, i] <- time_series_dat[, i]/total_cell_counts_breakdown$Freq[i-1]
    }
  }
  
  time_series_dat
}

order_cols <- function(ts) {
  x_cols = grep('X[0-9][0-9]*', colnames(ts))
  temp = ts[, x_cols, drop =F]
  temp = temp[, order(as.numeric(gsub('X', '', colnames(temp)))), drop=F]
  temp$clone_id <- ts$clone_id
  temp = temp[, c( ncol(temp), seq(ncol(temp)-1)), drop=F]
  temp
}


order_by_chr_name <- function(the_array) {
  the_array[the_array == 'X'] <- 44
  the_array[the_array == 'Y'] <- 100
  
  order(as.numeric(the_array))
}

strip_chr_names <- function(the_array) {
  print('In strip_chr_names')
  print(the_array)
  chunks <- unlist(strsplit(the_array, '_'))
  chrs <- matrix(chunks, nrow=length(the_array), ncol=3, byrow = T)[, 1]
  uchr <- unique(chrs)
  
  for (cc in uchr) {
    cc_idx <- which(chrs == cc)
    the_array[cc_idx] <- ""
    the_array[cc_idx[1]] <- cc
  }
  the_array
}


# library(tools)





'%ni%' <- function(x,y)!('%in%'(x,y))

# fasttable <- function(...) as.data.frame(fread(..., stringsAsFactors = F))

len <- function(...) length(...)

lenu <- function(...) length(unique(...))

headm <- function(amat, n = 10, ncol= 10) head(amat[, 1:ncol], n = n)

# generate_random_str <- function(n = 1) {
#   a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
#   paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
# }



# Account for synchronized cnvs
pad_mod_cnv_mat <- function(Y, output_file=NULL, pad_left = TRUE) {
  # TODO: We're missing synch events across the genome, i.e., those that happen in multiple cells, in more than one bin. They stand out in say SA919_Hoa
  # Binary matrix B, CNV matrix Y, cell c, locus l, set of cells C
  # 1. For each locus l, after removal of multi-chr-spanning loci, 
  # 1.a. find set of cells C for which B_l,[C] == 1 & |unique(Y_l,[C])| > 1
  # Note: this will fail at chr ends, handle via chr-mod-padding
  # 2. Add a new bin such that
  # B_l',i =  1, if Y_l,i == j, where j \prop Cat(unique(Y_l,[C]); proportions of diff. CNVs in [C])
  #           0, otherwise
  # Note: consider using B after jitterfix
  # Note: this will miss whole chr events -> Add a locus to B at either end of each chr where it's 1 for the mod and 0 o.w. (consider setting the second most frequent CNV bins to 1)
  # 3. Only use bins with |max(j)| < p%
  # Only consider loci that have |B_l,. == 1| > 90%
  # 4. Set l' (the name of the new locus in B) to chr_{(l+1)_(start) - l_start}_{0}
  
  # Inputs the_mat: Y; 
  # out_dir <- file.path(dirname(cnv_path), 'results'); dir.create(out_dir, showWarnings = F, recursive = T)
  # Y <- fasttable(cnv_path)
  # rownames(Y) <- Y$V1; Y$V1 <- NULL
  
  
  Y <- sort_mat_by_bins(Y)
  bin_dat <- parse_bin_names(rownames(Y))
  # head(bin_dat)
  ########################################
  B <- NULL
  uchrs <- unique(bin_dat$chr)
  for (chr in uchrs) {
    print(sprintf('chr == %s', chr))
    sub_Y <- Y[bin_dat$chr == chr, ]
    tmp <- abs(sub_Y[-c(nrow(sub_Y)), ] - sub_Y[-c(1), ])
    B <- addrow(B, tmp)
  }
  
  B[B > 1] <- 1
  
  ## Now has the right inputs, Y and B, both sorted
  # Note: For B[l, ] == 0, chromosome padding will catch them, otherwise they will be caught on a different bin
  # TODO: jitter fix B
  # TODO: DO THIS AFTER JITTERFIX in nextflow pipeline so as to not allow jitterfix to kill the signal
  # 1. Pick loci s.t. |B_l,. == 1| > 85%
  freq <- rowMeans(B)
  #L  <- which(unname(freq > .85))
  L  <- which(unname(freq > .63))
  freq[L]
  new_B <- B
  new_bins <- c()
  for (l in L) {
    # l = L[1]
    indx <- which(B[l, ] == 1)
    # 2. Add a new bin such that everything is great
    l_y <- which(rownames(Y) == rownames(B)[l])
    counts <- as.data.frame(table(as.vector(as.matrix(Y[l_y, indx]))), stringsAsFactors = F)
    counts$Freq <- counts$Freq/sum(counts$Freq)
    counts$Var1 <- as.numeric(counts$Var1)
    counts <- counts[order(counts$Freq, decreasing = T), ]
    if (nrow(counts) > 1) {
      if (counts$Freq[2] > .05) {
        print(sprintf('%s: L_B = %d, L_Y = %d',rownames(B)[l], l, l_y))
        tmp <- B[l, ]
        tmp[] <- 0
        tmp[which(Y[l_y, ] == counts$Var1[2])] <- 1
        fake_start <- floor((bin_dat$start[l_y] + bin_dat$end[l_y])/3)
        rownames(tmp) <- sprintf('%s_%d_%d', bin_dat$chr[l_y], fake_start, fake_start + 1)
        print(rownames(tmp))
        new_B <- rbind(new_B, tmp)
        new_bins <- c(new_bins, rownames(tmp))
      }
    }
  }
  
  
  # Now add the mode chromosome padding: use the categorical here
  new_bins <- c()
  for (chr in uchrs) {
    print(sprintf('Adding padding for chr == %s', chr))
    # chr = "10"
    last_Y <- Y[bin_dat$chr == chr, ]
    last_Y <- last_Y[nrow(last_Y), ]
    last_bin_name <- rownames(last_Y)
    l <- which(rownames(Y) == last_bin_name) # index in bin_names
    counts <- as.data.frame(table(as.vector(as.matrix(last_Y))), stringsAsFactors = F) %>% 
      dplyr::mutate(freq = Freq/sum(Freq), j = as.numeric(Var1)) %>% 
      dplyr::arrange(desc(freq)) %>% 
      dplyr::select(j, freq) 
    
    if (nrow(counts) > 1) {
      print(sprintf('%s: L_Y = %d', last_bin_name, l))
      tmp <- B[1, ]
      tmp[] <- 0
      tmp[which(sub_Y[l_y, ] == counts$j[2])] <- 1
      rownames(tmp) <- sprintf('%s_%d_%d', bin_dat$chr[l], bin_dat$end[l] + 1, bin_dat$end[l] + 2) # It is always larger
      print(rownames(tmp))
      new_B <- rbind(new_B, tmp)
      new_bins <- c(new_bins, rownames(tmp))
    }
    
    
  }
  
  # Do we want to add padding on the left too? Technically they will reinforce the whole chromosome and initila chromosome events...
  # TODO: where are the bins that we're removing??
  if (pad_left) {
    new_bins <- c()
    for (chr in uchrs) {
      print(sprintf('Adding left-padding for chr  %s', chr))
      # chr = "10"
      first_Y <- Y[bin_dat$chr == chr, ][1, ]
      first_bin_name <- rownames(first_Y)
      l <- which(rownames(Y) == first_bin_name) # index in bin_names
      counts <- as.data.frame(table(as.vector(as.matrix(first_Y))), stringsAsFactors = F) %>% 
        dplyr::mutate(freq = Freq/sum(Freq), j = as.numeric(Var1)) %>% 
        dplyr::arrange(desc(freq)) %>% 
        dplyr::select(j, freq) 
      
      # Found multi cnvs?
      if (nrow(counts) > 1) {
        print(sprintf('%s: L_Y = %d', first_bin_name, l))
        tmp <- B[1, ]
        tmp[] <- 0
        tmp[which(sub_Y[l_y, ] == counts$j[2])] <- 1
        rownames(tmp) <- sprintf('%s_%d_%d', chr, 0, 1) # It is always smaller
        print(rownames(tmp))
        new_B <- rbind(new_B, tmp)
        new_bins <- c(new_bins, rownames(tmp))
      }
    }
  }
  
  
  # sort new_B by rows again to put the new loci in their place...
  new_B <- sort_mat_by_bins(the_mat = new_B)
  nrow(new_B) - nrow(B) 
  
  
  new_B$loci <- rownames(new_B)
  rownames(new_B) <- NULL
  
  # cells,loci,tipInclusionProbabilities
  new_B <- new_B %>% tidyr::gather(key = 'cells', value = 'tipInclusionProbabilities', -loci)
  
  new_B <- new_B[, c(2,1,3)]
  head(new_B); summary(new_B$tipInclusionProbabilities); boxplot(new_B$tipInclusionProbabilities); table(new_B$tipInclusionProbabilities)
  
  # data.table::fwrite(new_B, output_file, row.names = F, quote = F)  
  return(new_B)
  
  # corrupt grow jitter
  # ./nextflow run binary-infer-pipeline.nf -resume --tipInclusionProbabilities ~/Desktop/SC-1311/SA609/processed_data/bincount/sa609_cnvs_corrupt.csv
  # pipeline <- 'no_pt_fast'
  # print(sprintf("./nextflow run %s.nf -resume --tipInclusionProbabilities %s", pipeline, file_path))
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

# cn2binary_newencoding <- function(cnv_path, output_file, prob = .90) {
#   # TODO: compute average ploidy too for each cell
#   
#   # cnv_path <- paste0(data_dir,'total_merged_filtered_states_original.csv')
#   mat <- read.delim(cnv_path, check.names = FALSE, stringsAsFactors = FALSE, sep = ",")
#   if(prob > 0){
#     jump_rank <- compute_jump_cells(mat)
#     
#     
#     # # Filter cells by jump/stuff
#     bad_cells <- jump_rank %>% dplyr::filter( njumps > quantile(x=njumps, probs = prob) | avgCNA > quantile(x=avgCNA, probs = prob)) %>% dplyr::select(cell_id)
#     
#     # bad_cells <- jump_rank %>% dplyr::filter( njumps > quantile(x=njumps, probs = prob)) %>% dplyr::select(cell_id)
#     
#     bad_cells <- bad_cells$cell_id
#     print("Number of bad cells: ")
#     print(length(bad_cells))
#     cc <- mat[, colnames(mat) %ni% bad_cells]
#     print(paste0("Filtered output: ", dim(cc)[1]," ",dim(cc)[2]))
#     write.csv(cc, opt$input,  row.names = T, quote=F)
#     pad_mod_cnv_mat(cc, output_file, pad_left = TRUE)
#     
#   } else{
#     print("Convert data to binary without filtering data: ")
#     pad_mod_cnv_mat(mat, output_file, pad_left = TRUE)
#   }
#     
#   # The results would be in:
#   # bin_cnvs_corrupt_double_padding.csv
#   
# }

## Parameters: remove outlier cells that different than the fraction of cells, 
## ex: prob=.90, removing outliers with number of jumps in copy number, or large copy number values than 90% of data. 
## prob=0: including all cells in analysis

## cnv_path; path to copy number states file which contains the format: genomic bins x cells, and values are copy number states
# input : the dir+ file name for copy number file: ex: mydir/total_merged_filtered_states.csv, a matrix: cells in columns, 
# chr_start_end position in row names 
## ex:            AT128-A98207B-R65-C38,AT128-A98207B-R65-C53
# 1_2000001_2500000   1, 2
# 1_3000001_350000    3, 3


## output_file: binary encoding output that contains 3 columns: "cells","loci","tipInclusionProbabilities"


cn2binary_newencoding <- function(cnv_path, output_file, prob=.90) {
  probAVG <- prob
  # TODO: compute average ploidy too for each cell
  
  # cnv_path <- paste0(data_dir,'total_merged_filtered_states.csv')
  mat <- read.delim(cnv_path, check.names = FALSE, stringsAsFactors = FALSE, sep = ",")
  # View(mat[1:3,1:3])
  # Convert data to binary mode
  binary_cnvs <- pad_mod_cnv_mat(mat, pad_left = TRUE)
  data.table::fwrite(binary_cnvs, output_file, row.names = F, quote = F)
  if(prob > 0){
    jump_rank <- compute_jump_cells(mat)
    # # Filter cells by jump/stuff
    bad_cells <- jump_rank %>% dplyr::filter( njumps > quantile(x=njumps, probs = prob) | avgCNA > quantile(x=avgCNA, probs = probAVG)) %>% dplyr::select(cell_id)
    # bad_cells <- jump_rank %>% dplyr::filter( njumps > quantile(x=njumps, probs = prob)) %>% dplyr::select(cell_id)
    
    bad_cells <- bad_cells$cell_id
    print(paste0("Number of bad cells: ",length(bad_cells)))
    if(length(bad_cells)>0){
      cc <- mat[, colnames(mat) %ni% bad_cells]
      outlier_cc <- mat[, colnames(mat) %in% bad_cells]
      print(paste0("Outlier cells mtx: ", dim(outlier_cc)[1]," ",dim(outlier_cc)[2]))
      corrupt_grow_dir <- paste0(dirname(output_file),'/corrupt_grow/')
      if (!file.exists(corrupt_grow_dir)){
        dir.create(corrupt_grow_dir)
      }
      write.csv(outlier_cc, paste0(corrupt_grow_dir,'total_merged_filtered_states_outlier_cells.csv'),  row.names = T, quote=F) #get dir
      write(colnames(outlier_cc), paste0(corrupt_grow_dir,'outlier_cells.txt'), sep = "\t")
    } else{
      cc <- mat
    }
    
    print(paste0("Filtered output mtx: ", dim(cc)[1]," ",dim(cc)[2]))
    write(colnames(cc), paste0(corrupt_grow_dir,'filtered_cells.txt'), sep = "\t")
    data.table::fwrite(cc, paste0(dirname(output_file),'/total_merged_filtered_states.csv'),  row.names = T, quote=F)
    
    # if(length(bad_cells)>0){
    #   binary_cnvs_filtered <- binary_cnvs[binary_cnvs$cells %ni% bad_cells,]
    #   binary_cnvs_outliers <- binary_cnvs[binary_cnvs$cells %in% bad_cells,]
    # The filtered results would be in:
    # bin_cnvs_corrupt_double_padding.csv
    # data.table::fwrite(binary_cnvs_filtered, output_file, row.names = F, quote = F)
    
    # The filtered results would be in:
    # corrupt_grow/bin_cnvs_corrupt_double_padding_outlier_cells.csv
    # data.table::fwrite(binary_cnvs_outliers, paste0(corrupt_grow_dir,'bin_cnvs_corrupt_double_padding_outlier_cells.csv'), row.names = F, quote = F)
    # }
  } else{
    print("Convert data to binary without filtering data ")
    # data.table::fwrite(binary_cnvs, output_file, row.names = F, quote = F)
  }
}

# Convert copy number to binary using new encoding scheme
print(paste0("Filtered outlier threshold is: ", opt$filterthrs))
cn2binary_newencoding(opt$input, opt$outfile, as.numeric(opt$filterthrs))

