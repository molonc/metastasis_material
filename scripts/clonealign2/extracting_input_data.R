


input_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/clonealign_v3/'
sid <- 'SA535X4XB05462'

sce <- readRDS(paste0(input_dir, '/',sid, '/',sid,'_sce.rds'))
dim(sce)
sce[1:3,1:3]
dim(counts(sce))
data.table::fwrite(as.matrix(counts(sce)), paste0(input_dir, '/',sid, '/',sid,'_expr.csv.gz'))

sce <- readRDS(paste0(input_dir, '/',sid, '/',sid,'_processed_data.rds'))
dim(sce$gene_expression_data)
expr <- t(sce$gene_expression_data)
class(expr)
View(expr[1:3,1:3])
dim(expr)
expr <- as.data.frame(expr)
data.table::fwrite(expr, paste0(input_dir, '/',sid, '/',sid,'_expr.csv.gz'), row.names=T)
t <- data.table::fread(paste0(input_dir, '/',sid, '/',sid,'_expr.csv.gz'))
t[1:3,1:3]
cnv <- sce$copy_number_data
class(cnv)
dim(cnv)
cnv[1:5,1:5]

clones_df <- data.frame(cell_id=paste0('cell_',colnames(cnv)),clone_id=colnames(cnv))
colnames(cnv) <- clones_df$cell_id
data.table::fwrite(clones_df, paste0(input_dir, '/',sid, '/',sid,'_cell_clones.csv.gz'))
t <- data.table::fread(paste0(input_dir, '/',sid, '/',sid,'_cell_clones.csv.gz'))
head(t)

cnv <- as.data.frame(cnv)
data.table::fwrite(cnv, paste0(input_dir, '/',sid, '/',sid,'_clones_cnv.csv.gz'),row.names = T)
t <- data.table::fread(paste0(input_dir, '/',sid, '/',sid,'_clones_cnv.csv.gz'))
t[1:3,1:3]
input_dir <- '/home/htran/storage/datasets/metastasis_results/SA535_new_encoding/SA535_wholedata_v2/'
cnv <- data.table::fread(paste0(input_dir, 'total_merged_filtered_states.csv'))
rownames(cnv) <- cnv$V1
cnv$V1 <- NULL
dim(cnv)
clones <- data.table::fread(paste0(input_dir, 'cell_clones.csv'))
dim(clones)
sum(colnames(cnv) %in% clones$cell_id)
cells_use <- colnames(cnv)[colnames(cnv) %in% clones$cell_id]
length(cells_use)
cnv <- as.data.frame(cnv)
cnv <- cnv[,cells_use]
dim(cnv)
