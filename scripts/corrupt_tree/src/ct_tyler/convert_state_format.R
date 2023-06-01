suppressPackageStartupMessages({
  require("data.table")
  require("dplyr")
  require("optparse")
})


option_list <- list(make_option(c("-i", "--inputfile"), type="character", default=NULL, help="features_mtx_file", metavar="character"),
                    make_option(c("-o", "--outputfile"), type="character", default=NULL, help="output_file", metavar="character"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)



get_chr_start_end <- function(inputfile, outputfile){
  print("Convert CN states normal format to Tyler mtx format")
  copy_number <- read.delim(inputfile, check.names = FALSE, stringsAsFactors = FALSE, sep = ",")
  # copy_number <- read.csv(inputfile, header=T, row.names=1, check.names = F,stringsAsFactors = FALSE)
  print("Read copy number matrix")
  print(dim(copy_number))
  features <- rownames(copy_number)
  print(features[1:3])
  print(length(features))
  bin_ids <- strsplit(as.character(features),"_")
  # length(bin_ids)
  nbcores <- 5
  ps <- c()
  chrs <- parallel::mclapply(bin_ids, function(f) {
    ps <- c(ps, as.character(f[1]))
  }, mc.cores = nbcores)
  if(length(chrs)!=length(bin_ids)){
    print("Error, double check")
  }
  # unique(chrs)
  
  ps <- c()
  starts <- parallel::mclapply(bin_ids, function(f) {
    ps <- c(ps, as.character(f[2]))
  }, mc.cores = nbcores)
  if(length(starts)!=length(bin_ids)){
    print("Error, double check")
  }
  starts <- unlist(starts, recursive = F, use.names = F)
  
  ps <- c()
  ends <- parallel::mclapply(bin_ids, function(f) {
    ps <- c(ps, as.character(f[3]))
  }, mc.cores = nbcores)
  if(length(ends)!=length(bin_ids)){
    print("Error, double check")
  }
  ends <- unlist(ends, recursive = F, use.names = F)
  
  copy_number$chr <- as.character(chrs)
  copy_number$start <- as.integer(starts)
  copy_number$end <- as.integer(ends)
  copy_number$width <- as.integer(rep(500000,length(features)))
  write.csv(copy_number, outputfile,  row.names = F, quote=F)
  
}

print(opt$inputfile)
print(opt$outputfile)
print(opt$outdir)

get_chr_start_end(opt$inputfile, opt$outputfile)
