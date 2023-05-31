source(paste0("/home/htran/Projects/farhia_project/rscript/pipeline/utils/normalize_utils.R"))
# BiocManager::install("batchelor")
library(batchelor)
library(BiocSingular)
library(dplyr)
## SEG list in ensemblGene ID
data("segList_ensemblGeneID", package = "scMerge") 

input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
output_dir <- paste0(input_dir,'rnaseq_v6/normalization_evaluation_v1/')
dir.create(output_dir)
datatag <- 'SA609'
sce_fn <- paste0(input_dir,'rnaseq_v6/',datatag,'-v6/total_sce_treated.rds')
sce <- readRDS(sce_fn)
dim(sce)


# print("Cosine normalization as logcounts input")
# cosdata <- cosineNorm(as.matrix(log2(counts(sce)+1)), BPPARAM = SerialParam())
# colnames(cosdata) <- colnames(sce)
# rownames(cosdata) <- rownames(sce)
# logcounts(sce) <- as.matrix(cosdata)

print("Simple logcounts input")
logcounts(sce) <- as.matrix(log2(counts(sce)+1))
  

metacells_fn <-  paste0(input_dir,'SA609_rna/snakemake_10x/SA609_10x.csv')
sample_df <- read.csv(metacells_fn)
# head(sample_df)
rownames(sample_df) <- sample_df$mouse_id
print(dim(sample_df))
sample_df <- sample_df %>%
              dplyr::select(mouse_id, batch_info) %>%
              dplyr::rename(batch=batch_info)
sce$batch <- 'None'
# metacells <- as.data.frame(colData(sce))
# metacells <- metacells %>% left_join(sample_df, by=c("sample"="mouse_id"))
# colData(sce) <- as.matrix(metacells)
# sapply(sample_df$mouse_id, function(s) {
#   sce[,sce$sample==s]$batch <- sample_df[s,'batch']
# })
sce$batch <- sample_df[sce$sample,'batch']
scSEG_df <- read.csv(paste0(output_dir,'segIndx_df_filtered_80.csv'),check.names = F, stringsAsFactors = F)
dim(scSEG_df)
colnames(scSEG_df)
# View(head(scSEG_df))
summary(scSEG_df$segIdx)
scSEG_df <- scSEG_df %>%
  dplyr::filter(segIdx >= median(segIdx))
scSEG_df <- scSEG_df[order(scSEG_df$segIdx,decreasing = T),]  
ntop <- 300
if(ntop>nrow(scSEG_df)){
  ntop <- nrow(scSEG_df)
}
scSEG_df <- scSEG_df[1:ntop,]
dim(scSEG_df)
# View(head(scSEG_df))
stable_genes <- intersect(scSEG_df$gene_ens, rownames(sce))
print(length(stable_genes))

summary(as.factor(sce$batch))
kmeansK_params <- c()
for(b in unique(sce$batch)){
  # sce_tmp <- sce[,sce$batch==b]
  # length(unique(sce_tmp$clone))
  kmeansK_params <- c(kmeansK_params, length(unique(sce[,sce$batch==b]$clone)))
}

t1 = Sys.time()

scMerge_res <- scMerge::scMerge(
  sce_combine = sce, 
  ctl = stable_genes,   #segList_ensemblGeneID$human$human_scSEG,
  assay_name = "scMerge_fast",
  replicate_prop = 0.4,
  cell_type = NULL, # unsupervised
  kmeansK = kmeansK_params,
  verbose=T,
  BSPARAM = IrlbaParam(), 
  svd_k = 20)

t2 = Sys.time()
print(t2-t1)


saveRDS(scMerge_res, file=paste0(output_dir,datatag,"_scMerge_correction_without_cosine.rds"))


t <- assay(scMerge_res, "scMerge_fast")
print(max(t))
print(min(t))

# Get PCA coordinates on normalized data
# normalized_data <- as.data.frame(scMerge_res@assays$data$scMerge_res) 
# normalized_data_t <- t(normalized_data)
# pca_mat <- stats::prcomp(normalized_data_t, rank = npcs, retx=TRUE, center = TRUE, scale. = FALSE)
# pca_mat <- as.data.frame(pca_mat$x)



