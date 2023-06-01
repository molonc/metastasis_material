suppressPackageStartupMessages({
  require("data.table")
  require("dplyr")
  require("optparse")
})

option_list <- list(make_option(c("-i", "--inputfile"), type="character", default=NULL, help="features_mtx_file", metavar="character"),
                    make_option(c("-o", "--outputfile"), type="character", default=NULL, help="output_file", metavar="character"),
                    make_option(c("-l", "--libid"), type="character", default=NULL, help="library_id", metavar="character"),
                    make_option(c("-di", "--inputdir"), type="character", default=NULL, help="input_dir", metavar="character"),
                    make_option(c("-do", "--outdir"), type="character", default=NULL, help="output_dir", metavar="character")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


get_cell_infos <- function(cells_info, obs_lib_ids){
  
  # outliers <- read.table(paste0(output_dir,'/corrupt_grow/outlier_cells.txt'), header = FALSE, sep = "\t")
  # filtered <- read.table(paste0(output_dir,'/corrupt_grow/filtered_cells.txt'), header = FALSE, sep = "\t")
  # outliers$cell_classified <- "outlier"
  # filtered$cell_classified <- "filtered"
  # cells_info <- rbind(outliers, filtered)
  # head(cells_info)
  unique(cells_info$cell_classified)
  
  cids <- strsplit(as.character(cells_info$cell_id),"-")
  # length(bin_ids)
  nbcores <- 5
  ps <- c()
  lids <- parallel::mclapply(cids, function(f) {
    ps <- c(ps, as.character(f[2]))
  }, mc.cores = nbcores)
  if(length(lids)!=length(cids)){
    print("Error, double check")
  }
  cells_info$lib_id <- as.character(lids)
  cells_info <- cells_info[cells_info$lib_id %in% obs_lib_ids,]
  filtered_cells <- cells_info[cells_info$cell_classified=="filtered",'cell_id']
  outlier_cells <- cells_info[cells_info$cell_classified=="outlier",'cell_id']
  res <- list(filtered_cells=filtered_cells, outlier_cells=outlier_cells)
}



get_filtered_mtx <- function(input_fn, output_dir, input_dir, pdx_id, filter_df=NULL){
  if (!file.exists(output_dir)){
    dir.create(output_dir)
  }
  corrupt_grow_dir <- paste(output_dir,'corrupt_grow/', sep="/")
  if (!file.exists(corrupt_grow_dir)){
    dir.create(corrupt_grow_dir)
  }
  
  
  # print(dim(filter_df))
  # cells_use <- unique(filter_df$cells)
  # length(cells_use)
  meta <- read.csv(paste0(input_dir,'/corrupt_grow/metadata.csv'), stringsAsFactors=F, check.names=F)
  # head(meta)
  # unique(meta$PDX_ID)
  meta <- meta[meta$PDX_ID==pdx_id,]
  obs_lib_ids <- meta$`Library ID_DLP`   #SA535
  # obs_lib_ids <- meta$`Library ID_DLP`   #SA919
  print(obs_lib_ids)
  
  # SA535
  outliers <- read.table(paste0(input_dir,'/corrupt_grow/outlier_cells.txt'), header = FALSE, sep = "\t")
  filtered <- read.table(paste0(input_dir,'/corrupt_grow/filtered_cells.txt'), header = FALSE, sep = "\t")
  outliers$cell_classified <- "outlier"
  filtered$cell_classified <- "filtered"
  cells_info <- rbind(outliers, filtered)
  # head(cells_info)
  colnames(cells_info)[which(names(cells_info) == "V1")] <- "cell_id"
  # unique(cells_info$cell_classified)
  
  # SA919
  # cells_info <- read.csv(paste0(input_dir,'/corrupt_grow/cells_info.csv'), stringsAsFactors=F, check.names=F)
  
  nbcells <- dim(cells_info)[2]
  cids <- strsplit(as.character(cells_info$cell_id),"-")
  # length(bin_ids)
  nbcores <- 5
  ps <- c()
  lids <- parallel::mclapply(cids, function(f) {
    ps <- c(ps, as.character(f[2]))
  }, mc.cores = nbcores)
  if(length(lids)!=length(cids)){
    print("Error, double check")
  }
  cells_info$lib_id <- as.character(lids)
  cells_info <- cells_info[cells_info$lib_id %in% obs_lib_ids,]
  filtered_cells <- cells_info[cells_info$cell_classified=="filtered",'cell_id']
  outlier_cells <- cells_info[cells_info$cell_classified=="outlier",'cell_id']
  
  data.table::fwrite(cells_info, paste0(corrupt_grow_dir,"cells_info.csv"), row.names = F, quote = F)
  
  
  print(length(filtered_cells))
  print(length(outlier_cells))
  copy_number <- read.delim(paste(input_dir,"total_merged_filtered_states_original.csv", sep="/"), 
                            check.names = FALSE, stringsAsFactors = FALSE, sep = ",", row.names = 1)
  if(dim(copy_number)[2]!=nbcells){
    print("Double check")
  }
  cn_filtered <- copy_number[,colnames(copy_number) %in% filtered_cells]
  cn_outliers <- copy_number[,colnames(copy_number) %in% outlier_cells]
  print('Divide filter loci mtx into filtered cells and outlier cells mtx')
  print(dim(cn_filtered))
  print(dim(cn_outliers))
  
  data.table::fwrite(cn_filtered, paste(output_dir,"total_merged_filtered_states.csv", sep="/"), row.names = T, quote = F)
  data.table::fwrite(cn_outliers, paste0(corrupt_grow_dir,'total_merged_filtered_states_outliers.csv'), row.names = T, quote = F)
  
  # or read cells_info.csv
  # filter_filtered <- filter_df[filter_df$cells %in% paste0('cell_', filtered_cells),]
  # filter_outliers <- filter_df[filter_df$cells %in% paste0('cell_', outlier_cells),]
  # print('Divide filter loci mtx into filtered cells and outlier cells mtx')
  # print(dim(filter_filtered))
  # print(dim(filter_outliers))
  # 
  # data.table::fwrite(filter_filtered, paste(output_dir,"corrupt_tree_sync-cnv_features.csv", sep="/"), row.names = F, quote = F)
  # data.table::fwrite(filter_outliers, paste0(corrupt_grow_dir,'corrupt_tree_sync-cnv_features_outliers.csv'), row.names = F, quote = F)
  print('Save data, Completed')
}


print(opt$inputfile)
print(opt$outputfile)
print(opt$inputdir)
print(opt$outdir)
print(opt$libid)

# input_dir <- "/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_filtered"
# input_fn <- paste(input_dir,"corrupt_tree_sync-cnv_features.csv", sep="/")
# base_dir <- "/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding"
# filter_df <- data.table::fread(input_fn, stringsAsFactors=F, check.names=F)

input_dir <- "/home/htran/storage/datasets/metastasis_results/SA535_new_encoding/Tyler_whole"
input_fn <- paste(input_dir,"corrupt_tree_sync-cnv_features_original.csv", sep="/")
base_dir <- "/home/htran/storage/datasets/metastasis_results/SA535_new_encoding"
filter_df <- data.table::fread(input_fn, stringsAsFactors=F, check.names=F)


meta <- read.csv(paste0(input_dir,'/corrupt_grow/metadata.csv'), stringsAsFactors=F, check.names=F)
# head(meta)
pdx_ids <- unique(meta$PDX_ID)
# pdx_ids <- c("X0847-2112251","X0847-2112252","X0847-2112253",
#              "X0847-2112254","X0847-216","X0847-2164")
# fns <- c("Tyler_2251","Tyler_2252","Tyler_2253","Tyler_2254","Tyler_216","Tyler_2164")
fns <- c("Tyler_2361","Tyler_2362","Tyler_2363","Tyler_2364")

for(i in rep(1:length(pdx_ids),1)){
  output_dir <- paste(base_dir,fns[i], sep="/")
  pdx_id <- pdx_ids[i]
  print(output_dir)
  print(pdx_id)
  get_filtered_mtx(input_fn, output_dir,input_dir, pdx_id, filter_df=NULL)
}

get_filtered_mtx(opt$inputfile, opt$outputfile, opt$outdir, opt$libid)