# clone prevalence
suppressPackageStartupMessages({
  require("optparse")
})

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

source(paste0(script.basename, "/prevalence_utils.R"))


option_list <- list(make_option(c("-i", "--inputdir"), type="character", default=NULL, help="input_dir", metavar="character"),
                    make_option(c("-o", "--outputfile"), type="character", default=NULL, help="output_file", metavar="character"),
                    make_option(c("-d", "--datatag"), type="character", default='corrupt_tree', help="datatag", metavar="character"),
                    make_option(c("-c", "--cellclones"), type="character", default=NULL, help="cellclones_file", metavar="character"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
print(opt$inputdir)
print(opt$outputfile)
print(opt$datatag)
print(opt$cellclones)

# results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_new_encode/SA1035_Tyler_v2/'
# source_dir <- '/home/htran/Projects/farhia_project/rscript/dlp/visualize_tree/'
# source(paste0(source_dir,'prevalence_utils.R'))

# results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535/SA535_Tyler_v2'
# outputfile <- paste0(results_dir,'/prevalences/','timepoint_clones_prevalence.png')
# datatag = 'SA535_Tyler'

results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_wholedata/'
outputfile <- paste0(results_dir,'prevalences/clones_prevalence.csv')
cellclones <- paste0(results_dir,'cell_clones.csv')
datatag <- 'SA919'
outputfile <- paste0(results_dir,'prevalences/data_met_proj_v3/','timepoint_clones_prevalence.png')
# 


datatag <- 'SA535'
results_dir <- '/home/htran/storage/datasets/metastasis_results/SA535_new_encoding/SA535_wholedata_v2/'
outputfile <- paste0(results_dir,'prevalences/','timepoint_clones_prevalence.png')
# outputfile <- paste0(results_dir,'prevalences/data_met_proj_v3/','timepoint_clones_prevalence.png')
cellclones <- paste0(results_dir,'cell_clones.csv')
outputfile <- paste0(results_dir,'prevalences/same_clone_report/','timepoint_clones_prevalence.png')

get_prevalence <- function(results_dir, outputfile=NULL, cellclones=NULL, datatag = ''){
  
  print("Generating prevalence plots...")
  results_dir <- paste0(results_dir,'/')
  if(is.null(cellclones)){
    cellclones <- paste0(results_dir, 'cell_clones.csv')
  }
  save_dir <- paste0(dirname(outputfile),'/')
  if (!file.exists(save_dir)){
    dir.create(save_dir)
  }
  cell_clones <- read.delim(cellclones, stringsAsFactors = FALSE, sep = ",")
  excluded_clones <- c('Un','unassigned','None')
  cell_clones <- cell_clones[!cell_clones$clone_id %in% excluded_clones,]
  print(dim(cell_clones))
  print(summary(as.factor(cell_clones$clone_id)))
  # cell_clones_backup <- cell_clones
  library_grouping_fn <- paste0(results_dir,'library_groupings.csv')
  meta_data <- get_meta_data(cell_clones, results_dir, library_grouping_fn)
  print(dim(meta_data))
  # table(meta_data$mainsite, meta_data$clone_id)
  # table(meta_data$mainsite)
  # meta_data <- meta_data[meta_data$PDX=='SA535_CX',]
  # print('Verification...')
  # print(unique(meta_data$PDX))
  plot_prevalences_origin(cell_clones, results_dir, save_dir, datatag, meta_data, NULL)
  write.csv(meta_data, file = outputfile, quote = F, row.names = F)
  # cls_df <- get_pct_clusters(meta_data, save_dir)
  # colnames(cls_df)[which(names(cls_df) == "Freq")] <- "clone_freq"
  # color_ls <- assign_color_clone(cell_clones)
  # output_fn = paste0(save_dir,datatag,'_treatment_clones_prevalence.png')
  # plot_passage_clusters(output_fn, cls_df, color_ls, datatag, tag='treatment', save_dir,
  #                       xstring="passage",ystring="clone_freq",
  #                       plttitle="Clones Prevalence",
  #                       lgtitle="Clones", xlb = 'Treatment Status')
  # 
  # cls_tp_df <- get_pct_timepoint(meta_data, color_ls, save_dir)
  # colnames(cls_tp_df)[which(names(cls_tp_df) == "Freq")] <- "clone_freq"
  # if(is.null(outputfile)){
  #   outputfile = paste0(save_dir,'timepoint_clones_prevalence.png')
  # }
  # plot_passage_clusters(outputfile, cls_tp_df, color_ls, datatag, tag='timepoint', save_dir,
  #                       xstring="timepoint",ystring="clone_freq",
  #                       plttitle="Clones Prevalence",
  #                       lgtitle="Clones", xlb = 'Time Point')
  
}

# results_dir <- '/home/htran/storage/datasets/metastasis_results/SA535_new_encoding/SA535_wholedata_v2/'
# # cellclones <- paste0(results_dir,'cell_clones.csv')
# outputfile <- paste0(results_dir,'proportions/','mouse_clone_proportions.png')
# datatag <- 'SA535'

# results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_wholedata/'
# outputfile <- paste0(results_dir,'proportions/clones_prevalence.csv')
# cellclones <- paste0(results_dir,'cell_clones.csv')
# datatag <- 'SA919'



get_prevalence(opt$inputdir, opt$outputfile, cellclones=opt$cellclones, datatag = opt$datatag)

