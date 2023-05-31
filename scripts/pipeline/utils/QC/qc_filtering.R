suppressPackageStartupMessages({
  require("optparse")
  require("scater")
  require("argparse")
  require("SingleCellExperiment")
  require("stringr")
  # require("tidyverse")
})
# htran
option_list <- list(make_option(c("--input_file"), type="character", default=NULL, 
                                help="inputfile", metavar="character"),
                    make_option(c("--output_file"), 
                                type="character", 
                                default=NULL, 
                                help="output_file", 
                                metavar="character"),
                    make_option(c("--mouse_id"), 
                                type="character", 
                                default=NULL, 
                                help="mouse_id", 
                                metavar="character"),
                    make_option(c("--library_id"), 
                                type="character", 
                                default=NULL, 
                                help="library_id", 
                                metavar="character"),
                    make_option(c("--min_features"),
                                 type = "double",
                                 default = 1000,
                                 help = "min_features"),
                    make_option(c("--nmads_threshold"),
                                type = "double",
                                default = 3,
                                help = "nmads_threshold"),
                    make_option(c("--max_mito"),
                                type = "double",
                                default = 20,
                                help = "max_mito"),
                    make_option(c("-r", "--max_ribo"),
                                type = "double",
                                default = 60,
                                help = "max_ribo"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)



runQC <- function(input_file, min_features=1000, max_mito=20, max_ribo=60, 
                  output_file=NULL, library_id=NULL, nmads_threshold=3){
  output_dir <- paste0(dirname(output_file),'/')
  if (!file.exists(output_dir)){
    dir.create(output_dir)
  }
  if (file.exists(input_file)){
    print(paste0("Processing file:  ",input_file))
    sce <- readRDS(input_file)
    print(paste0('Original size: ',dim(sce)[1],' ',dim(sce)[2], ' library_id:', library_id))
    
    print('Detecting mito, ribo genes...')
    mito_genes <- str_detect(rowData(sce)$Symbol, "^MT\\-")
    print(sum(mito_genes==TRUE))
  
    ribo_genes <- str_detect(rowData(sce)$Symbol, "^RP(L|S)")  # or ^RP[L|S]?
    print(sum(ribo_genes==TRUE))
    
  
    # For scater version < 1.14.6
    # sce <- calculateQCMetrics(sce, exprs_values = "counts", feature_controls =
    #                             list(Mito=mito_genes, Ribo=ribo_genes))
  
    # For scater version >= 1.14.6
    print('Quality Control...')
    per.cell <- perCellQCMetrics(sce, subsets=list(Mito=mito_genes, Ribo=ribo_genes))
    colData(sce) <- cbind(colData(sce), per.cell)
    View(head(per.cell$detected))
    # print(paste0('After combination: ',dim(sce)[1],' ',dim(sce)[2]))
    
    # print(colnames(colData(sce)))
    # print(summary(sce$total_counts))  # sum
    # print(summary(summary(sce$pct_counts_Mito)))
    # print(summary(summary(as.factor(sce$pct_counts_Ribo))))
    # sm <- data.frame(cellname = rownames(colData(sce)),
    #                  subsets_Mito_percent=sce$pct_counts_mito,
    #                  subsets_Ribo_percent=sce$pct_counts_ribo,
    #                  total_features_by_counts=sce$total_features_by_counts,
    #                  sum=sce$total_counts,
    #                  library_id = sce$library_id)
    # write.table(sm, file = paste0(output_dir,library_id,'_summary.txt'),sep="\t",quote=F)
  
    # p1 <- hist(sce$pct_counts_Mito, breaks=300)
    # p1
  
    ## -----------------------------------------------------------------------------
    
    ## ----plot-pdata-pct-exprs-controls--------------------------------------------
    # plotColData(sce, x = "sum", y="subsets_Mito_percent")
    # plotColData(sce, x = "sum", y="subsets_Ribo_percent")
    
    # print(dim(colData(sce)))
    print('Raw data, percentage of ribo: ')
    print(summary(sce$subsets_Ribo_percent))
    keep_rb <- sce$subsets_Ribo_percent < max_ribo
    # print(paste0("Nb cells with high ribo percent than our threshold: ",sum(keep_rb==FALSE)))
    print('Raw data, percentage of mito: ')
    print(summary(sce$subsets_Mito_percent))
    keep_mt <- sce$subsets_Mito_percent < max_mito
    # if(max_mito>0){  # first option: using a fixed threshold
    #   keep_mt <- sce$subsets_Mito_percent <= max_mito
    #   # print("Filtering cells using a fixed mito percent as threshold")
    #   print(library_id)
    #   print(summary(sce$subsets_Mito_percent))
    # } else{          # second option: using second quantile value of subset mito as threshold
    #   max_mito_quantile = 0.5
    #   keep_mt <- sce$subsets_Mito_percent < quantile(sce$subsets_Mito_percent, max_mito_quantile)
    #   print("Filtering cells using second quantile value of mito percent as threshold")
    ## Another option to filter mito is:
    # mito.drop <- isOutlier(sce$pct_counts_feature_controls_Mt, nmads=3, type="higher")
    # pike.drop <- isOutlier(sce$pct_counts_feature_controls_ERCC, nmads=3, type="higher")
    # # }
  
  
    print(paste0("Nb cells with mito percent <=",max_mito,' mito threshold pct is: ',round(100*sum(keep_mt==T)/dim(sce)[2],2),'%'))
    # print(summary(sce$detected))
    keep_features <- sce$total_features_by_counts >= min_features
    # keep_features <- sce$detected >= min_features
    print(paste0("Nb cells with minimum nb features > ",min_features,' is: ',round(100*sum(keep_features==T)/dim(sce)[2],2),'% '))
    
    # keep_total <- sce$total_counts > min_pc_sum
    # print(sum(keep_total==FALSE))
    # 
    # keep_pc <- sce$detected > min_pc_detected
    # print(sum(keep_pc==FALSE))
    #
    # sce <- sce[, keep_rb & keep_mt & keep_features & keep_total & keep_pc]
    
    # Reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5112579/pdf/f1000research-5-10712.pdf
    # To obtain an adaptive threshold, we assume that most of the dataset consists of high-quality cells. 
    # We remove cells with log-library sizes that are more than 3 median absolute deviations (MADs) below the median log-library size. (A log-transformation improves resolution at small values, especially when the MAD of the raw values is comparable to or greater than the median.) We also remove cells 
    # where the log-transformed number of expressed genes is 3 MADs below the median.
    libsize.drop <- scater::isOutlier(sce$total_counts, nmads=nmads_threshold, type="lower", log=TRUE)
    feature.drop <- scater::isOutlier(sce$total_features_by_counts, nmads=nmads_threshold, type="lower", log=TRUE)
    
    print(paste0("Nb cells with library size outlier MAD < ",nmads_threshold,' median of all cells is: ',round(100*sum(feature.drop==F)/dim(sce)[2],2),'% '))
    print(paste0("Nb cells with feature outlier MAD < ",nmads_threshold,' median of all cells is: ',round(100*sum(libsize.drop==F)/dim(sce)[2],2),'% '))
    
    # Create a report
    sce$outlier_libsize <- libsize.drop
    sce$outlier_feature <- feature.drop
    sce$keep_ribo <- keep_rb
    sce$keep_mito<- keep_mt
    sce$keep_min_features <- keep_features
    meta_info <- as.data.frame(colData(sce))
    meta_info$library_label <- library_id
    write.csv(meta_info, file=paste0(output_dir,library_id,'_qc.csv'), quote=F, row.names=F)
    print(paste0("Raw data: ",dim(sce)[1]," ", dim(sce)[2]))
    sce <- sce[,keep_rb & keep_mt & keep_features & !(libsize.drop | feature.drop)]
    print(paste0("Filtered data: ",dim(sce)[1]," ", dim(sce)[2]))
    # sce <- sce[, keep_rb & keep_mt & keep_features]
    # print(paste0("Filtered sample: ",dim(sce)[1]," ", dim(sce)[2]))
    
    saveRDS(sce, file = output_file)
    # return(sce)
  }else{
    stop("Input file do not exist, pls double check input data!!!")
  }
      
}




# Rscript /home/htran/Projects/hakwoo_project/rscript/pipeline/qc_filtering.R 
# --inputfile /home/htran/storage/datasets/hakwoo_metastasis_RNAseq/SCRNA10X_SA_CHIP0190_001/SCRNA10X_SA_CHIP0190_001.rdata --outputfile /home/htran/storage/datasets/metastasis_results/rnaseq/filtered/SCRNA10X_SA_CHIP0190_001_filtered.rds --library_id SCRNA10X_SA_CHIP0190_001 --mouse_id SA919X4XB40503 --pdxid X0847_2164 
# --passage X4 --Site_origin Axillary --Grouping Metastasis --batch_info CHIP0190


print(opt$input_file)
print(opt$output_file)
print(paste0('Params: min_features:',opt$min_features,
             ' max_mito:',opt$max_mito,' max_ribo:', opt$max_ribo))

runQC(opt$input_file, opt$min_features, 
      opt$max_mito, opt$max_ribo, 
      opt$output_file, opt$library_id, opt$nmads_threshold)

# input_file <- '/home/htran/storage/datasets/hakwoo_metastasis_RNAseq/testing/human/SCRNA10X_SA_CHIP0146_004/SCRNA10X_SA_CHIP0146_004.rdata'
# output_file <- '/home/htran/storage/datasets/hakwoo_metastasis_RNAseq/testing/human/filtered/SCRNA10X_SA_CHIP0146_004_filtered.rds'
# library_id <- 'SCRNA10X_SA_CHIP0146_004'






# make_option(c("--passage"), 
#             type="character", 
#             default=NULL, 
#             help="passage", 
#             metavar="character"),
# make_option(c("--treatmentSt"), 
#             type="character", 
#             default=NULL, 
#             help="treatmentSt", 
#             metavar="character"),
# make_option(c("--Grouping"), 
#             type="character", 
#             default=NULL, 
#             help="Grouping", 
#             metavar="character"),
# make_option(c("--batch_info"), 
#             type="character", 
#             default=NULL, 
#             help="batch_info", 
#             metavar="character"),
# ,
# make_option(c("--pdxid"), 
#             type="character", 
#             default=NULL, 
#             help="pdxid", 
#             metavar="character"),
# make_option(c("--min_pc_detected"),
#             type = "double",
#             default = 100,
#             help = "min_pc_detected"),
# make_option(c("--min_pc_sum"),
#             type = "double",
#             default = 500,
#             help = "min_pc_sum")