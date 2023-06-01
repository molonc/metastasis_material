suppressPackageStartupMessages({
  require("data.table")
  require("dplyr")
  require("optparse")
})

# option_list <- list(make_option(c("-i", "--inputdir"), type="character", default=NULL, help="state_cnv_dir", metavar="character"),
#                     make_option(c("-o", "--outputdir"), type="character", default=NULL, help="output_dir", metavar="character"),
#                     make_option(c("-d", "--outdir"), type="character", default=NULL, help="output_dir", metavar="character")
# )

option_list <- list(make_option(c("-f", "--filtered_fn"), type="character", default=NULL, help="features_mtx_file", metavar="character"),
                    make_option(c("-o", "--output_fn"), type="character", default=NULL, help="output_file", metavar="character"),
                    make_option(c("-i", "--features_outliers_fn"), type="character", default=NULL, help="output_dir", metavar="character"),
                    make_option(c("-c", "--curr_cells_clone"), type="character", default=NULL, help="cells_clone file", metavar="character")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


# paste0(results_dir, 'filtered.csv')
# paste0(results_dir,'corrupt_grow/corrupt_tree_sync-cnv_features_outliers.csv')
# paste0(results_dir,'corrupt_grow/filtered_outliers.csv')

# 
get_filtered_outliers <- function(except_clones, curr_cell_clones, filtered_fn, features_outliers_fn, output_fn){
  cell_clones <- read.delim(curr_cell_clones, stringsAsFactors = FALSE, sep = ",", check.names = F)
  exception_cells <- cell_clones[cell_clones$clone_id %in% except_clones,'cell_id']
  print(paste0('Nb exception cells: ', length(exception_cells)))
  
  filtered_df <- read.csv(filtered_fn, check.names = F, stringsAsFactors = FALSE)
  print(dim(filtered_df))
  loci_ls <- unique(filtered_df$loci)
  print(paste0("Nb filtered loci: ",length(loci_ls)))
  
  filtered_outliers <- data.table::fread(features_outliers_fn, stringsAsFactors = F)
  filtered_outliers <- filtered_outliers[filtered_outliers$cells %in% paste0('cell_',exception_cells),]
  print('DEBUG')
  print(dim(filtered_outliers))
  print(length(unique(filtered_outliers$cells)))
  # data.table::fwrite(filtered_outliers, paste0(output_dir,'/corrupt_grow/corrupt_tree_sync-cnv_features_outliers.csv'), row.names = F, quote = F)
  outlier_loci <- unique(filtered_outliers$loci)
  print(length(outlier_loci))
  if(sum(outlier_loci %in% loci_ls)==length(loci_ls)){
    filtered_outliers <- filtered_outliers[filtered_outliers$loci %in% loci_ls,]
    print(length(unique(filtered_outliers$cells)))
    data.table::fwrite(filtered_outliers, output_fn, row.names = F, quote = F)
    print("Get filtered outliers mtx, done")
  } else{
    print("Nb features of outlier cells do not match with filtered cells, double check input data!!!")
  }
  
}
except_clones <- c('F','G')
get_filtered_outliers(except_clones, opt$curr_cells_clone, opt$filtered_fn, opt$features_outliers_fn, opt$output_fn)
