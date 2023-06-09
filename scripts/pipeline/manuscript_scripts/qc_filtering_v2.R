suppressPackageStartupMessages({
  # require("optparse")
  require("scater")
  require("argparse")
  require("SingleCellExperiment")
  require("stringr")
  require("dplyr")
  require("data.table")
  # require("tidyverse")
})


option_list <- list(make_option(c("--inputfile"), 
                                type="character", 
                                default=NULL, 
                                help="inputfile", 
                                metavar="character"),
                    make_option(c("--outputfile"), 
                                type="character", 
                                default=NULL, 
                                help="outputfile", 
                                metavar="character"),
                    make_option(c("--mouse_id"), 
                                type="character", 
                                default=NULL, 
                                help="mouse_id", 
                                metavar="character"),
                    make_option(c("--pdxid"), 
                                type="character", 
                                default=NULL, 
                                help="pdxid", 
                                metavar="character"),
                    make_option(c("--library_id"), 
                                type="character", 
                                default=NULL, 
                                help="library_id", 
                                metavar="character"),
                    make_option(c("--passage"), 
                                type="character", 
                                default=NULL, 
                                help="passage", 
                                metavar="character"),
                    make_option(c("--Site_origin"), 
                                type="character", 
                                default=NULL, 
                                help="Site_origin", 
                                metavar="character"),
                    make_option(c("--Grouping"), 
                                type="character", 
                                default=NULL, 
                                help="Grouping", 
                                metavar="character"),
                    make_option(c("--batch_info"), 
                                type="character", 
                                default=NULL, 
                                help="batch_info", 
                                metavar="character"),
                    make_option(c("--min_features"),
                                 type = "double",
                                 default = 1000,
                                 help = "min_features"),
                    make_option(c("--max_mito"),
                                type = "double",
                                default = 20,
                                help = "max_mito"),
                    make_option(c("-r", "--max_ribo"),
                                type = "double",
                                default = 60,
                                help = "max_ribo"),
                    make_option(c("--min_pc_detected"),
                                type = "double",
                                default = 100,
                                help = "min_pc_detected"),
                    make_option(c("--min_pc_sum"),
                                  type = "double",
                                  default = 500,
                                  help = "min_pc_sum"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)




remove_mouse_cells <- function(sce, mouse_fn){
  print(paste0('Raw data: ',dim(sce)[1],' ',dim(sce)[2]))
  if(file.exists(mouse_fn)){
    df <- data.table::fread(mouse_fn)
    sce <- sce[, !sce$Barcode %in% df$Barcode]
    print(paste0('After excluding mouse cells: ',dim(sce)[1],' ',dim(sce)[2]))  
  }else{
    print('Do not exist mouse cells')
    print(mouse_fn)
  }
  return(sce)
}

# For scater version < 1.14.6
runQC <- function(sce, output_file=NULL, output_dir=NULL, library_id=NULL, min_features=1000, 
                  max_mito=20, max_ribo=60, min_pc_detected=200,
                  min_pc_sum=1000,  
                  nmads_threshold=3){

  mito_genes <- str_detect(toupper(rowData(sce)$Symbol), "^MT\\-")
  print(sum(mito_genes==TRUE))

  ribo_genes <- str_detect(toupper(rowData(sce)$Symbol), "^RP(L|S)")  # or ^RP[L|S]?
  sum(ribo_genes==TRUE)
  print(paste0('Original size: ',dim(sce)[1],' ',dim(sce)[2], ' library_id: ', library_id))

  # For scater version < 1.14.6
  # sce <- calculateQCMetrics(sce, exprs_values = "counts", feature_controls =
  #                             list(Mito=mito_genes, Ribo=ribo_genes))

  # For scater version >= 1.14.6
  per.cell <- perCellQCMetrics(sce, subsets=list(Mito=mito_genes, Ribo=ribo_genes))
  colData(sce) <- cbind(colData(sce), per.cell)
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
  # sce$subsets_Mito_sum
  keep_rb <- sce$subsets_Ribo_percent <= max_ribo
  print(paste0("Nb cells with high ribo percent than our max threshold ",max_ribo,': ',sum(keep_rb==FALSE)))
  # print(sum(keep_rb==T))
  # print(paste0("Nb cells with high ribo percent than our threshold: ",sum(keep_rb==FALSE)))
  
  # max_mito <- 45
  
  if(max_mito>0){  # first option: using a fixed threshold
    keep_mt <- sce$subsets_Mito_percent <= max_mito
    # print("Filtering cells using a fixed mito percent as threshold")
    print(library_id)
    print(summary(sce$subsets_Mito_percent))
  } else{          # second option: using second quantile value of subset mito as threshold
    max_mito_quantile = 0.5
    keep_mt <- sce$subsets_Mito_percent < quantile(sce$subsets_Mito_percent, max_mito_quantile)
    print("Filtering cells using second quantile value of mito percent as threshold")
  }
  

  print(paste0("Nb cells with high mito percent than our threshold ",max_mito,': ',sum(keep_mt==FALSE)))
  # summary(sce$detected)
  if(!'total_features_by_counts' %in% colnames(colData(sce))){
    sce$total_features_by_counts <- sce$detected
  }
  # keep_features <- sce$detected >= min_features
  keep_features <- sce$total_features_by_counts >= min_features
  # keep_features <- sce$detected >= min_features
  print(paste0("Nb cells with minimum nb features ",min_features,": ",sum(keep_features==T)))

  # keep_total <- sce$total_counts > min_pc_sum
  # print(sum(keep_total==FALSE))
  
  # keep_pc <- sce$detected > min_pc_detected
  # print(sum(keep_pc==FALSE))
  #
  # Reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5112579/pdf/f1000research-5-10712.pdf
  # To obtain an adaptive threshold, we assume that most of the dataset consists of high-quality cells. 
  # We remove cells with log-library sizes that are more than 3 median absolute deviations (MADs) below the median log-library size. (A log-transformation improves resolution at small values, especially when the MAD of the raw values is comparable to or greater than the median.) We also remove cells 
  # where the log-transformed number of expressed genes is 3 MADs below the median.
  # libsize.drop <- scater::isOutlier(sce$sum, nmads=nmads_threshold, type="lower", log=TRUE)
  # feature.drop <- scater::isOutlier(sce$detected, nmads=nmads_threshold, type="lower", log=TRUE)
  
  ## for old version or mouse samples 
  # libsize.drop <- scater::isOutlier(sce$total, nmads=nmads_threshold, type="lower", log=TRUE)
  # feature.drop <- scater::isOutlier(sce$detected, nmads=nmads_threshold, type="lower", log=TRUE)
  # 
  if(!'total_counts' %in% colnames(colData(sce))){
    sce$total_counts <- sce$total
  }
  
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
  # sce <- sce[,keep_rb & keep_mt & keep_features & !(libsize.drop | feature.drop)]
  sce <- sce[,keep_rb & keep_mt & keep_features & !(libsize.drop | feature.drop)]
  print(paste0("Filtered data: ",dim(sce)[1]," ", dim(sce)[2]))
  
  # sce <- sce[, keep_rb & keep_mt & keep_features]  #& keep_total & keep_pc
  # print(paste0("Filtered sample: ",dim(sce)[1]," ", dim(sce)[2]))
  saveRDS(sce, file = output_file)
  # return(sce)
      
}

input_dir <- '/home/htran/storage/rnaseq_datasets/hakwoo_metastasis_RNAseq/'
fns <- basename(list.dirs(paste0(input_dir, 'SA535_human_introns/')))
fns <- fns[grepl('SCRNA', fns)]
fns
# fns <- c('SCRNA10X_SA_CHIP0220_001','SCRNA10X_SA_CHIP0220_002')
output_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/filtered_introns/'
  
for(library_id in fns){
  # library_id <- 'SCRNA10X_SA_CHIP0176_001'
  print('')
  print('_____________________________')
  print(library_id)
  sce_human_fn <- paste0(input_dir, 'SA535_human_introns/', library_id, '/', library_id, '.rdata')
  sce_mouse_fn <- paste0(input_dir, 'mouse_cells/', library_id, '_mouse_cells.csv.gz')
  if(file.exists(sce_human_fn)){
    sce <- readRDS(sce_human_fn)
    if(dim(sce)[2]>0){
      sce <- remove_mouse_cells(sce, sce_mouse_fn)
      output_file <- paste0(output_dir, library_id,'.rds')
      runQC(sce, output_file, output_dir, library_id, min_features=1000, 
            max_mito=20, max_ribo=60, min_pc_detected=200,
            min_pc_sum=1000,  
            nmads_threshold=3)
    }  
  }
  
  
  
}


meta <- data.table::fread('/home/htran/Projects/hakwoo_project/metastasis_material/materials/10x/SA535_10x_metadata.csv')

meta <- meta %>% 
  dplyr::filter(mouse_id=='SA535X4XB05462')
dim(meta)
meta$library_id

# base_dir <- '/home/htran/storage/rnaseq_datasets/hakwoo_metastasis_RNAseq/'
output_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/filtered_introns/'
sce_list <- list()
for(lid in meta$library_id){
  sce <- readRDS(paste0(output_dir, lid, '.rds'))
  sce_list[[lid]] <- sce
}
assay_counts <- do.call(cbind, lapply(sce_list, function(y) assay(y, 'counts')))
names(assay_counts) <- 'counts'
colData_names <- c('Barcode','Sample')

colData_list <- do.call(DelayedArray::rbind, 
                        lapply(sce_list, function(y) colData(y)[, colData_names, drop = FALSE]))
sce_combine <- SingleCellExperiment::SingleCellExperiment(assays = list(counts=assay_counts), 
                                                          colData = colData_list)
dim(sce_combine)
save_dir <- '/home/htran/storage/rnaseq_datasets/testing_space/clonealign2/SA535X4XB05462/'
saveRDS(sce_combine, paste0(save_dir, 'SA535X4XB05462_introns_sce.rds'))
assayNames(sce_combine)
colnames(sce_combine) <- sce_combine$Barcode
mtx <- as.data.frame(counts(sce_combine))
dim(mtx)
rownames(mtx)[1]
cnv <- data.table::fread(paste0(save_dir, 'SA535X4XB05462_clones_cnv.csv.gz'))
genes_used <- intersect(cnv$V1, rownames(mtx))
mtx <- mtx[genes_used,]
head(cnv)
data.table::fwrite(mtx, paste0(save_dir, 'SA535X4XB05462_expr_introns.csv.gz'), row.names = T)
dim(mtx)
# /home/htran/anaconda3/envs/sisyphus/bin/python exe_TreeAlign.py --expr_fn /home/htran/storage/rnaseq_datasets/testing_space/clonealign2/SA535X4XB05462/SA535X4XB05462_expr_introns.csv.gz --cnv_fn /home/htran/storage/rnaseq_datasets/testing_space/clonealign2/SA535X4XB05462/SA535X4XB05462_clones_cnv.csv.gz --cell_clones_fn /home/htran/storage/rnaseq_datasets/testing_space/clonealign2/SA535X4XB05462/SA535X4XB05462_cell_clones.csv.gz --datatag SA535X4XB05462 --save_dir /home/htran/storage/rnaseq_datasets/testing_space/clonealign2/SA535X4XB05462/
# meta <- data.table::fread('/home/htran/storage/datasets/drug_resistance/rna_results/SA604_rna/normalized/SA604_sctransform_colData.csv.gz')
# 
# head(meta)
# table(meta$library_id, meta$Sample)
# unique(meta$material)
# meta <- meta %>%
#   dplyr::group_by(tenx,Sample) %>%
#   dplyr::summarise(nb_sample=n())
# 
# dim(meta)
# View(meta)
# 
# colnames(meta) <- c('sample_id','library_id','nb_cells')
# 
# sids <- c('SA501X2XB00096','SA530X3XB03295')
# lids <- c('SCRNA10X_SA_CHIP0008_000','SCRNA10X_SA_CHIP0164_002')
# 
# extra <- tibble(sample_id=sids, library_id=lids)
# extra$nb_cells <- NULL
# t <- dplyr::bind_rows(meta, extra)
# View(t)
# data.table::fwrite(t, '/home/htran/Projects/farhia_project/drug_resistant_material/materials/metadata_drug_resistance/untreated_Pts_SA501_SA530_SA604.csv')
# /usr/local/bin/Rscript /home/htran/Projects/farhia_project/rnaseq/pipeline/utils/identify_doublet_v2.R --sce_file /home/htran/storage/images_dataset/merfish_rnaseq/SITTA12_XP6873/analysis/filtered/SITTA12_XP6873_filtered.rds --sce_output /home/htran/Projects/farhia_project/rnaseq/pipeline/utils/identify_doublet_v2.R --sce_file /home/htran/storage/images_dataset/merfish_rnaseq/SITTA12_XP6873/analysis/filtered/SITTA12_XP6873_f_dt.rds --doublet_score 0 --output_csv /home/htran/Projects/farhia_project/rnaseq/pipeline/utils/identify_doublet_v2.R --sce_file /home/htran/storage/images_dataset/merfish_rnaseq/SITTA12_XP6873/analysis/filtered/SITTA12_XP6873_filtered.rds --sce_output /home/htran/Projects/farhia_project/rnaseq/pipeline/utils/identify_doublet_v2.R --sce_file /home/htran/storage/images_dataset/merfish_rnaseq/SITTA12_XP6873/analysis/filtered/SITTA12_XP6873__doublet_filtered.csv --summary_csv /home/htran/storage/images_dataset/merfish_rnaseq/SITTA12_XP6873/analysis/filtered/SITTA12_XP6873_qc.csv

# checking <- function(){
#   p <- plotColData(sce, x = "sum", y="subsets_Mito_sum", colour_by="subsets_Mito_percent")
#   p
#   p1 <- plotColData(sce, x = "sum", y="detected", colour_by="detected") 
#   png(paste0(output_dir,"sum_mito.png"), height = 2*300, width=2*550,res = 2*72)
#   print(p)
#   dev.off()
#   png(paste0(output_dir,"sum_detected.png"), height = 2*300, width=2*550,res = 2*72)
#   print(p1)
#   dev.off()
#   t <- counts(sce)
#   t[1:3,1:3]
#   library(scran)
#   mito_genes <- str_detect(res1$symb, "^MT\\-")
#   "^RP(L|S)"
#   mt <- grep("^MT\\-", res1$symb, value = F)
#   length(mt)
#   
#   mito_genes <- str_detect(rowData(sce)$Symbol, "^MT\\-")
#   sum(mito_genes==TRUE)
#   length(mito_genes)
#   ribo_genes <- str_detect(rowData(sce)$Symbol, "^RP(L|S)")  # or ^RP[L|S]?
#   sum(ribo_genes==TRUE)
#   
#   
#   sce1 <- sce[(!mito_genes) & (!ribo_genes), ]
#   dim(sce1)
#   
#   
#   is_high_mito <- sce1$subsets_Mito_percent > median(sce1$subsets_Mito_percent)
#   out <- findMarkers(sce1, groups=is_high_mito)
#   names(out)
#   class(out[[1]])
#   res <- as.data.frame(out[[1]])
#   res$ens <- rownames(res)
#   res$ens[1:3]
#   meta_genes <- data.frame(ens=rownames(rowData(sce)),symb=rowData(sce)$Symbol, row.names = rownames(rowData(sce)), stringsAsFactors = F)
#   dim(t)
#   
#   res <- dplyr::left_join(res, meta_genes, by=c('ens'))
#   summary(res$p.value)
#   summary(res$logFC.TRUE)
#   summary(res$Top)
#   res1 <- res[res$p.value<0.025 & abs(res$logFC.TRUE) > 0.25,]
#   dim(res1)
#   View(head(res1))
#   mito_genes <- str_detect(res1$symb, "^MT\\-")
#   "^RP(L|S)"
#   mt <- grep("^MT\\-", res1$symb, value = F)
#   length(mt)
# 
#   res2 <- res1[res1$logFC.TRUE>0.25,]
#   res2 <- res2[order(res2$logFC.TRUE,decreasing = T),]
#   View(res2[1:15,])
#   dim(res2)
#   res3 <- res1[res1$logFC.TRUE<(-0.25),]
#   sum(mito_genes==TRUE)
# }
min_features=1000
max_mito = 35
max_ribo=60
# library_id <- 'SCRNA10X_SA_CHIP0077_001'
# library_id <- 'SCRNA10X_SA_CHIP0077_003'
library_id <- 'SCRNA10X_SA_CHIP0078_001'
output_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA919/filtered/'
output_file <- paste0(output_dir, library_id,'_filtered.rds')
input_dir <- '/home/htran/storage/rnaseq_datasets/hakwoo_metastasis_RNAseq/SA919/'
inputfile <- paste0(input_dir, library_id,'/',library_id,'.rdata')
run_qc_filtering <- function(meta_info, inputfile, outputfile, 
                             min_features=1000, max_mito = 20, max_ribo=60, 
                             min_pc_detected=100, min_pc_sum=500){
  output_dir <- paste0(dirname(outputfile),'/')
  if (!file.exists(output_dir)){
    dir.create(output_dir)
  }
  
  if (file.exists(inputfile)){
    print(paste0("Processing file:  ",inputfile))
    sce <- readRDS(inputfile)
    print(paste0("SCE:  ",dim(sce)[1]," ",dim(sce)[2]))
    
    # sce$mouse_id <- meta_info$mouse_id
    # sce$pdxid <- meta_info$pdxid
    # sce$passage <- meta_info$passage
    # # sce$Site_origin <- meta_info$Site_origin
    # # sce$Grouping <- meta_info$Grouping
    # sce$library_label <- meta_info$library_id
    # sce$batch_info <- meta_info$batch_info
    
    # library_id <- meta_info$library_id
    sce <- runQC(sce, min_features, max_mito, max_ribo, min_pc_detected, 
                 min_pc_sum, outputfile, output_dir, meta_info$library_id)
    
  }  
}

# Rscript /home/htran/Projects/hakwoo_project/rscript/pipeline/qc_filtering.R 
# --inputfile /home/htran/storage/datasets/hakwoo_metastasis_RNAseq/SCRNA10X_SA_CHIP0190_001/SCRNA10X_SA_CHIP0190_001.rdata --outputfile /home/htran/storage/datasets/metastasis_results/rnaseq/filtered/SCRNA10X_SA_CHIP0190_001_filtered.rds --library_id SCRNA10X_SA_CHIP0190_001 --mouse_id SA919X4XB40503 --pdxid X0847_2164 
# --passage X4 --Site_origin Axillary --Grouping Metastasis --batch_info CHIP0190

# meta_info <- list(mouse_id="SA919X4XB09552", pdxid="X0847_2161", passage="X4",
#                   Site_origin="Primary", Grouping="Primary",
#                   library_id="SCRNA10X_SA_CHIP0191_002", batch_info="CHIP0191")
# sce_rds_input <- '/home/htran/storage/datasets/hakwoo_metastasis_RNAseq/'
# lid = "SCRNA10X_SA_CHIP0078_001"
# output_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq/'
# inputfile = paste0(sce_rds_input,'SCRNA10X_SA_CHIP0078_001/SCRNA10X_SA_CHIP0078_001.rdata')
# outputfile = paste0(output_dir, lid, '_filtered.rds')

print(opt$inputfile)
print(opt$outputfile)
print(paste0('Params: min_features:',opt$min_features,' max_mito:',opt$max_mito,' max_ribo:', opt$max_ribo,
             ' min_pc_detected',opt$min_pc_detected,' min_pc_sum:',opt$min_pc_sum))

# meta_info <- list(mouse_id=opt$mouse_id, pdxid=opt$pdxid, passage=opt$passage,
#                   Site_origin=opt$Site_origin, Grouping=opt$Grouping, 
#                   library_id=opt$library_id, batch_info=opt$batch_info)

# print(meta_info)
meta_info <- NULL
run_qc_filtering(meta_info, opt$inputfile, opt$outputfile, 
                 opt$min_features, opt$max_mito, opt$max_ribo, 
                 opt$min_pc_detected, opt$min_pc_sum)








# run_QC_Cluster <- function(sce, min_features=1000, 
#                            max_mito = 25, max_mito_quantile=0.75, 
#                            max_ribo=60, 
#                            min_pc_sum=1000, output_dir=NULL){
#   
#   
#   mito_genes <- str_detect(rowData(sce)$Symbol, "^MT-")  #"^MT\\-"
#   sum(mito_genes==TRUE)
#   
#   ribo_genes <- str_detect(rowData(sce)$Symbol, "^RP(L|S)")  # or ^RP[L|S]?
#   sum(ribo_genes==TRUE)
#   
#   dim(sce)
#   
#   # For scater version < 1.14.6
#   sce <- calculateQCMetrics(sce, exprs_values = "counts", feature_controls =
#                               list(Mito=mito_genes, Ribo=ribo_genes))
#   
#   # For scater version >= 1.14.6 
#   # per.cell <- perCellQCMetrics(sce, subsets=list(Mito=mito_genes, Ribo=ribo_genes))
#   # print(summary(sce$total_counts))  # sum
#   # print(summary(summary(sce$pct_counts_Mito)))
#   # print(summary(summary(as.factor(sce$pct_counts_Ribo))))
#   sm <- data.frame(cellname = rownames(colData(sce)),
#                    subsets_Mito_percent=sce$pct_counts_Mito,
#                    subsets_Ribo_percent=sce$pct_counts_Ribo,
#                    total_features_by_counts=sce$total_features_by_counts,
#                    sum=sce$total_counts,
#                    id = sce$id, 
#                    treatment_status = sce$treatment_status,
#                    libid = sce$libid)
#   write.table(sm, file = paste0(output_dir,unique(sce$id),'_summary.txt'),sep="\t",quote=F,col.names=NA)
#   
# } 
