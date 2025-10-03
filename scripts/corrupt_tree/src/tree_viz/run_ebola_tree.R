

# Remove all cells, and edge link to these cells
# Count nb cells in each loci 
# set the size of loci vertex proportional to nb cell counts

# Based on Sohrab Salehi's code, fitness paper
# Ebola tree reference: https://www.cell.com/action/showPdf?pii=S0092-8674%2815%2900690-X

suppressPackageStartupMessages({
  require("optparse")
})

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

source(paste0(script.basename, "/utils.R"))
source(paste0(script.basename, "/trajectory_utils_v2.R"))
source(paste0(script.basename, "/ebola_tree_utils.R"))
source(paste0(script.basename, "/collapse_tree_utils.R"))

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


print("Generating ebola tree plots...")
# Plot ebola tree at each time point, and whole dataset
# plot_ebola_tree_condition(opt$datatag, opt$inputdir, opt$outputfile, opt$cellclones)
plot_ebola_tree_condition(opt$datatag, opt$inputdir, opt$outputfile, opt$cellclones)


plot_summary_tree(opt$datatag, opt$inputdir, opt$outputfile, opt$cellclones)



script.basename <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/scripts/corrupt_tree/src/tree_viz'
source(paste0(script.basename, "/utils.R"))
source(paste0(script.basename, "/trajectory_utils_v2.R"))
source(paste0(script.basename, "/ebola_tree_utils.R"))
source(paste0(script.basename, "/collapse_tree_utils.R"))
# results_dir <- '/home/htran/storage/datasets/metastasis_results/SA535_new_encoding/SA535_wholedata_v2'
# datatag <- 'SA535'
# outputfile <- paste0(results_dir,'/tree_viz_dream_testing/trajectory_clone_wholedata.png')
# plot_ebola_tree_condition(datatag, results_dir, outputfile, cellclones=NULL)

# results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_wholedata/'
# datatag <- 'SA919'
# outputfile <- paste0(results_dir,'/tree_viz_dream/trajectory_clone_wholedata.png')
# plot_ebola_tree_condition(datatag, results_dir, outputfile, cellclones=NULL)



results_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/figures/fig_cnv/'
p_tree <- readRDS(paste0(results_dir, 'summary_tree_SA535.rds'))

p_tree










  
  
  
  
  
  
  
  
  
  