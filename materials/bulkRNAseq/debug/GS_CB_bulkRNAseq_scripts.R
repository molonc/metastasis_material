library(DESeq2)
library(EnhancedVolcano)
library(textshape)

save_dir <- "~/Documents/BCCRC_projects/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/debug/"
counts <- read.csv(file = paste0(save_dir, "SA919_New_Norm_rounded_PM_CB.csv"))
# counts <- column_to_rownames(counts, loc = 4)
dim(counts)
head(counts)
counts <- column_to_rownames(counts)
head(counts)
counts <- counts %>%
  tibble::column_to_rownames('ens_gene_id')

colSums(counts)
# coldata = read.csv(file = paste0(save_dir, "9_Samp_Meta_Cond7_PM.csv"))
coldata = data.table::fread(paste0(save_dir, "9_Samp_Meta_Cond7_PM.csv"))
ds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~Condition)
colnames(ds) <- colnames(counts)
ds <- DESeq(ds)
res <- results(ds, c("Condition","Metastasis","Primary"))
dim(res)
class(res)
summary(res$log2FoldChange)
deseq2_out <- as.data.frame(res)
dim(deseq2_out)
deseq2_out <- deseq2_out %>%
  dplyr::filter(!is.na(log2FoldChange))


cnv_genes1 <- intersect(cnv$ensembl_gene_id, rownames(deseq2_out))
deseq2_out$ensembl_gene_id <- rownames(deseq2_out)

deseq2_out <- deseq2_out %>%
  dplyr::filter(ensembl_gene_id %in% cnv_genes1)
dim(deseq2_out)
sum(deseq2_out$pvalue<0.05)
summary(deseq2_out$log2FoldChange)
sum(abs(deseq2_out$log2FoldChange)>=1)
head(res)
summary(res)
norm_df <- counts(ds, normalized=TRUE)
dim(norm_df)
colSums(norm_df)
write.csv(as.matrix(res_g_ID_names), file = "P7_PM_CB_New_Norm_Meta_vs_Prim_DESeq_All.csv")
sig <- res_g_ID_names[ which(res_g_ID_names$padj < 0.05), ]
write.csv(as.matrix(sig), file = "P7_PM_CB_New_Norm_Meta_vs_Prim_DESeq_Sig.csv")


