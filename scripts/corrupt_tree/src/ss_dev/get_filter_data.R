suppressPackageStartupMessages({
  require("data.table")
  require("dplyr")
  require("optparse")
})

option_list <- list(make_option(c("-i", "--input"), type="character", default=NULL, help="filtered_cnv_input", metavar="character"),
                    make_option(c("-o", "--outfile"), type="character", default=NULL, help="tip_probability_output", metavar="character"),
                    make_option(c("-d", "--outdir"), type="character", default=NULL, help="output_dir", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

get_filtered_mtx <- function(input_fn, output_fn, output_dir){
  filter_df <- read.csv(input_fn, stringsAsFactors = F)
  print(dim(filter_df))
  outliers <- read.table(paste0(output_dir,'/corrupt_grow/outlier_cells.txt'), header = FALSE, sep = "\t")
  filtered <- read.table(paste0(output_dir,'/corrupt_grow/filtered_cells.txt'), header = FALSE, sep = "\t")
  filter_df_filtered <- filter_df[filter_df$cells %in% paste0('cell_',filtered[,1]),]
  filter_df_outliers <- filter_df[filter_df$cells %in% paste0('cell_',outliers[,1]),]
  print('Divide filter loci mtx into filtered cells and outlier cells mtx')
  print(dim(filter_df_filtered))
  print(dim(filter_df_outliers))
  data.table::fwrite(filter_df_filtered, output_fn, row.names = F, quote = F)
  data.table::fwrite(filter_df_outliers, paste0(output_dir,'/corrupt_grow/filter_loci_outliers.csv'), row.names = F, quote = F)
  print('Save filtered data, Completed')
}

# get_filtered_mtx <- function(input_fn, output_fn, output_dir){
#   filter_df <- read.csv(input_fn, stringsAsFactors = F)
#   print(dim(filter_df))
#   outliers <- read.table(paste0(output_dir,'/corrupt_grow/outlier_cells.txt'), header = FALSE, sep = "\t")
#   filtered <- read.table(paste0(output_dir,'/corrupt_grow/filtered_cells.txt'), header = FALSE, sep = "\t")
#   filter_df_filtered <- filter_df[filter_df$cells %in% paste0('cell_',filtered[,1]),]
#   filter_df_outliers <- filter_df[filter_df$cells %in% paste0('cell_',outliers[,1]),]
#   print('Divide filter loci mtx into filtered cells and outlier cells mtx')
#   print(dim(filter_df_filtered))
#   print(dim(filter_df_outliers))
#   data.table::fwrite(filter_df_filtered, output_fn, row.names = F, quote = F)
#   data.table::fwrite(filter_df_outliers, paste0(output_dir,'/corrupt_grow/filter_loci_outliers.csv'), row.names = F, quote = F)
#   print('Save filtered data, Completed')
# }

# results_dir1 <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_new_encode/SA1035_whole_local/'
# results_dir2 <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_new_encode/SA1035_filtered/'
# results_dir3 <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_new_encode/filtered_data/'
# 
# 
# filter_df1 <- read.csv(paste0(results_dir1,'filtered.csv'), stringsAsFactors = F)
# filter_df2 <- read.csv(paste0(results_dir2,'filtered.csv'), stringsAsFactors = F)
# filter_df3 <- read.csv(paste0(results_dir3,'filtered.csv'), stringsAsFactors = F)
# l1 <- unique(filter_df1$loci)
# l2 <- unique(filter_df2$loci)
# l3 <- unique(filter_df3$loci)
# length(l1)
# length(l2)
# length(l3)
# 
# straight_df <- data.table::fread(paste0(results_dir2,'straighten_output.csv'), stringsAsFactors = F)
# outliers <- read.table(paste0(output_dir,'/corrupt_grow/outlier_cells.txt'), header = FALSE, sep = "\t")
# filtered <- read.table(paste0(output_dir,'/corrupt_grow/filtered_cells.txt'), header = FALSE, sep = "\t")
# straight_filtered <- straight_df[straight_df$cells %in% paste0('cell_',filtered[,1]),]
# straight_outliers <- straight_df[straight_df$cells %in% paste0('cell_',outliers[,1]),]
# data.table::fwrite(straight_filtered, paste0(results_dir3,'straight_filtered.csv'), row.names = F, quote = F)
# data.table::fwrite(straight_outliers, paste0(results_dir3,'straight_outliers.csv'), row.names = F, quote = F)
# 

# filtered_data
print(opt$input)
print(opt$outfile)
print(opt$outdir)
get_filtered_mtx(opt$input, opt$outfile, opt$outdir)


