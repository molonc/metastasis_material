suppressPackageStartupMessages({
  library("optparse")##BiocManager::install("argparse")
  library("dplyr")##BiocManager::install("dplyr")
  # library("ggplot2")
  library("DropletUtils")  ##BiocManager::install("DropletUtils")
  library("SingleCellExperiment")##BiocManager::install("SingleCellExperiment")
})

option_list <- list(make_option(c("--input_dir"), 
                                type="character", 
                                default=NULL, 
                                help="inputfile", 
                                metavar="character"),
                    make_option(c("--output_file"), 
                                type="character", 
                                default=NULL, 
                                help="outputfile", 
                                metavar="character"))


opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


load_10x_alignment_scRNAseq <- function(input_dir, output_file){
  
  output_dir <- paste0(dirname(output_file),'/')
  if (!file.exists(output_dir)){
    dir.create(output_dir)
  }
  lib_id <- gsub('.rdata','',basename(output_file))
  ## R version
  
  if(dir.exists(input_dir) & !file.exists(output_file)){
    sce <- DropletUtils::read10xCounts(input_dir)
    dim(sce)
    
    ## python scanpy version
    # https://scanpy.readthedocs.io/en/stable/generated/scanpy.read_10x_mtx.html
    # https://scanpy.readthedocs.io/en/stable/generated/scanpy.read_10x_h5.html
    
    # metacells <- colData(sce)
    # print("Number of cells: ")
    # print(dim(metacells)[1])
    # print("Meta cells info: ")
    # print(head(metacells))
    if(grepl(lib_id,sce$Sample[1])){
      print('Correct label')
    }else{
      print(lib_id)
      stop('Double check observed library id: ')
    }
    sce$Sample <- lib_id
    # metagenes <- rowData(sce)
    # print("Number of genes: ")
    # print(dim(metagenes))
    # print("Meta genes info: ensembl gene id and gene symbols")
    # print(head(metagenes))
    rowData(sce)$Type <- NULL
    saveRDS(sce, output_file)
    print("Save file into folder: ")
    print(output_file)
  }
  
}

# load_10x_alignment_scRNAseq(input_dir, output_file)

base_dir <- '/home/htran/storage/rnaseq_datasets/hakwoo_metastasis_RNAseq/SA535_human_introns/'
fns <- list.files(base_dir)
fns <- fns[grepl('SCRNA10X_SA',fns)]
for(f in fns){
  input_dir <- paste0(base_dir, f, '/filtered_feature_bc_matrix/')
  output_file <- paste0(base_dir, f, '/',f,'.rdata')
  # load_10x_alignment_scRNAseq(input_dir, output_file)
  if(!file.exists(output_file)){
    print(f)
  }
}

input_dir <- '/home/htran/storage/rnaseq_datasets/hakwoo_metastasis_RNAseq/SA535_human_introns/SCRNA10X_SA_CHIP0077_004/filtered_feature_bc_matrix/'
output_file <- '/home/htran/storage/rnaseq_datasets/hakwoo_metastasis_RNAseq/SA535_human_introns/SCRNA10X_SA_CHIP0077_004/SCRNA10X_SA_CHIP0077_004.rdata'
load_10x_alignment_scRNAseq(input_dir, output_file)

## how to run script from command line: 
# Rscript /home/htran/Projects/spatial_genomics/merfish_materials/single_cell_analysis/rnaseq_scripts/pipeline/utils/load_10x.R --input_dir /home/htran/storage/images_dataset/merfish_rnaseq/SCRNA10X_SA_CHIP0267_001_AT22361/filtered_feature_bc_matrix/ --output_file /home/htran/storage/images_dataset/merfish_rnaseq/SCRNA10X_SA_CHIP0267_001_AT22361/SCRNA10X_SA_CHIP0267_001_AT22361.rdata
