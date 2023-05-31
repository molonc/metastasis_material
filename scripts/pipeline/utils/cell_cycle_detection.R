#' Cell cycle prediction with cyclone on SingleCellExperiment
suppressPackageStartupMessages({
  library(tidyverse)
  library(SingleCellExperiment)
  library(scater)
  library(data.table)
  library(methods)
  library(scran)
  library(argparse)
  library(org.Hs.eg.db)
  library(stringr)
})

parser <- ArgumentParser(description = "Run cyclone on SingleCellExperiment")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS")
parser$add_argument('--ncpus', type='character',
                    help="Number of cores to use.",default = '5')
parser$add_argument('--output_cc', type = 'character', metavar = 'FILE',default = NULL,
                    help="Output path for cell cycle assignments.")
parser$add_argument('--output_file', type = 'character', metavar = 'FILE',default = NULL,
                    help="Output path for sce file.")
args <- parser$parse_args()

sce_path <- args$sce
if(is.null(args$output_cc)){
  args$output_cc <- paste0(dirname(args$sce),'/cellcycle_scores.rds')
}
output_dir <- paste0(dirname(args$output_cc),'/')
if (!file.exists(output_dir)){
  dir.create(output_dir)
}

base_name <- basename(args$sce)
base_name <- gsub('_filtered.rds','',base_name)
print(paste0("base_name is: ", base_name))

# sce_path <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_cys/filtered/TENX048_filtered.rds'
sce <- readRDS(sce_path)
print('Raw data: ')
print(dim(sce))
hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))

# mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", 
#                                 package="scran"))


# Run cyclone
if (!str_detect(rownames(sce)[1], "^ENSG")) {
  
  print(rownames(sce)[1])
  print('Using ensemble as genes ids')
  if('ID' %in% colnames(rowData(sce))){
    print(rowData(sce)$ID[1])
    rownames(sce) <- rowData(sce)$ID  
    gene_ids <- rownames(sce)
    sce1 <- sce
  }else{
    stop('Check ensemble genes id assignment')
    # gene_ids <- mapIds(org.Hs.eg.db, rownames(sce), 'ENSEMBL', 'SYMBOL')
    # sce1 <- sce[!is.na(gene_ids),]
    # gene_ids <- unname(gene_ids[!is.na(gene_ids)][rownames(sce1)])
  }
  
} else {
  gene_ids <- rownames(sce)
  sce1 <- sce
}

assignments <- cyclone(sce1, hs.pairs, gene.names=gene_ids, 
                       min.iter = 10, verbose = TRUE, BPPARAM = MulticoreParam(as.integer(args$ncpus)))
# assignments <- cyclone(sce, hs.pairs, gene.names=gene_ids, 
#                        min.iter = 10, verbose = TRUE, BPPARAM = MulticoreParam(5))
sce$cell_cycle_phases <- assignments$phases
print('Detect cell cycles...')
print(sce_path)
print(summary(sce$cell_cycle_phases))
# p <- plot(assignments$score$G1, assignments$score$G2M,
#      xlab="G1 score", ylab="G2/M score", pch=16)
# png(paste0(output_dir,base_name,"_cell_cycle.png"), height = 2*250, width=2*300,res = 2*72)
#   print(p)
# dev.off()
# Write outputs
saveRDS(assignments, file = args$output_cc)
print(dim(sce))
if(!is.null(args$output_file)){
  saveRDS(sce, file = args$output_file)  
}

cat("Completed.\n")


# sce <- readRDS('/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/filtered/SCRNA10X_SA_CHIP0176_001_f_dt.rds')
# summary(as.factor(sce$cell_cycle_phases))
# Rscript /home/htran/Projects/farhia_project/rnaseq/pipeline/utils/cell_cycle_detection.R --sce /home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v7/SA535_v7/sce/combined_SA535_cisplatin.rds 
# Cell cycle scores

# output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v7/SA535_v7/sce/'
# assignment_cc <- readRDS(paste0(output_dir,'cellcycle_scores.rds'))
# sce <- readRDS(paste0(output_dir,'combined_SA535_cisplatin.rds'))
# dim(assignment_cc)
# assignment_cc$phases[1:5]
# assignment_cc$scores$G1[1:5]
# assignment_cc$scores$G2M[1:5]
# 
# 
# cell_cycle_info <- data.frame(S_score=assignment_cc$scores$S, G2M_score=assignment_cc$scores$G2M,
#                               G1_score=assignment_cc$scores$G1,cell_cycle_phases=assignment_cc$phases,
#                               cell_id=colnames(sce), treatmentSt=sce$treat, passage=sce$timepoint, 
#                               mouse_id=sce$sample, library_id=sce$library_id)
# dim(cell_cycle_info)
# save_dir = '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/encoder_trajectory/'
# data.table::fwrite(cell_cycle_info, paste0(output_dir,'cell_cycle_info.csv.gz'))
# data.table::fwrite(cell_cycle_info, paste0(save_dir,'cell_cycle_info.csv.gz'))
# View(head(cell_cycle_info))

