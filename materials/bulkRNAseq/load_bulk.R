library(scran)
library(dplyr)
library(ggplot2)
# library(scater)


# BiocManager::install('scran')
# scran::calculateSumFactors()
# scran:::computeSumFactors()

ngenes <- 10000
ncells <- 20
mu <- 2^runif(ngenes, -1, 5)
gene.counts <- matrix(rnbinom(ngenes*ncells, mu=mu, size=10), nrow=ngenes)

# library(org.Mm.eg.db)
# all.ensembl <- unique(toTable(org.Mm.egENSEMBL)$ensembl_id)
rownames(gene.counts) <- paste0('G',seq(1:dim(gene.counts)[1]))

nspikes <- 100
ncells <- 20
mu <- 2^runif(nspikes, -1, 5)
spike.counts <- matrix(rnbinom(nspikes*ncells, mu=mu, size=10), nrow=nspikes)
rownames(spike.counts) <- paste0("ERCC-", seq_len(nspikes))
all.counts <- rbind(gene.counts, spike.counts)
library(scran)
sce <- SingleCellExperiment(list(counts=all.counts))
# isSpike(sce, "MySpike") <- grep("^ERCC", rownames(sce))
calculateSumFactors(sce)
summary(sizeFactors(sce))
colnames(colData(sce))
dim(sce)
sce
scran::computeSumFactors(sce)
sce
summary(sce$sizeFactor)

raw_counts <- as.data.frame(counts(sce))
dim(raw_counts)
t <- raw_counts[1:5,1:2]
t <- raw_counts
dim(t)
s <- colSums(t)
t1 <- t / sce$sizeFactor
t1 <- t / sce$sizeFactor[1:2]
dim(t1)
t1


M1<-matrix(rpois(80,10),ncol=4)
V1<-1:4
sweep(M1,2,V1,FUN="/")

t1 <- sweep(t,2,sce$sizeFactor[1:2],FUN="/")

t1 <- sweep(t,2,sce$sizeFactor,FUN="/")
t1
s1 <- colSums(t1)
s
s1
t[1:5,1:2]
t1[1:5,1:2]
# ggplot(data = )




input_dir <- '/Users/htran/Documents/storage_tmp/metastasis_trees/SA919_10x/'

df <- data.table::fread(paste0(input_dir, 'SA919_total_raw_counts.csv.gz')) 
dim(df)
rownames(df) <- df$ens_gene_id
df$ens_gene_id <- NULL
head(df)
colSums(df)
colnames(df)
rownames(df)[1:3]
meta_genes <- data.frame(ens_gene_id=rownames(df))
sce <- SingleCellExperiment::SingleCellExperiment(list(counts=as.matrix(df)))
rowData(sce) <- meta_genes
sce
rownames(sce) <- rownames(df)
sce <- scran::computeSumFactors(sce)
sce$sizeFactor

raw_counts <- as.data.frame(counts(sce))
rownames(raw_counts)[1:3]
s <- colSums(raw_counts)
norm_counts_mtx <- sweep(raw_counts,2,sce$sizeFactor,FUN="/")
s1 <- colSums(norm_counts_mtx)

dim(norm_counts_mtx)
norm_counts_mtx <- as.data.frame(norm_counts_mtx)
norm_counts_mtx$ens_gene_id <- rownames(norm_counts_mtx)
norm_counts_mtx$ens_gene_id[1:2]

data.table::fwrite(norm_counts_mtx, paste0(input_dir, 'SA919_total_normalized_sizefactor.csv.gz'))


t <- log2(norm_counts_mtx[1:5,1:2]+1)
t
fc <- t$SA919X3XB08939 - t$SA919X4XB09563


ref <- annotables::grch38 %>%
  dplyr::select(ensgene, symbol)
get_ens_gene_ids <- function(gene_ids){
  genes <- lapply(strsplit(gene_ids,'\\.'), function(x){
    return(x[1])
  })#, mc.cores=3)
  return(as.character(genes))
}
genes_used <- get_ens_gene_ids(rownames(df))
genes_used[1:3]
colnames(t)  
dim(t)

meta_genes <- data.frame(ens_gene_id=genes_used)
dim(meta_genes)
meta_genes <- meta_genes %>%
  dplyr::left_join(ref, by=c('ens_gene_id'='ensgene'))
dim(meta_genes)

length(unique(meta_genes$symbol))

data.table::fwrite(meta_genes, paste0(input_dir, 'meta_genes.csv.gz'))




norm_counts_mtx <- data.table::fread(paste0(input_dir, 'SA919_total_normalized_sizefactor.csv.gz'))

norm_counts_mtx$ens_gene_id[1:10]
norm_counts_mtx$ens_gene_id <- get_ens_gene_ids(norm_counts_mtx$ens_gene_id)


sum(is.na(meta_genes$symbol))
meta_genes <- meta_genes %>%
  dplyr::filter(!is.na(symbol) & symbol!='')
dim(meta_genes)
norm_counts_mtx <- norm_counts_mtx %>%
  dplyr::inner_join(meta_genes, by=c('ens_gene_id'))%>%
  dplyr::filter(!is.na(symbol)) 
dim(norm_counts_mtx)
colnames(norm_counts_mtx)
norm_counts_mtx$symbol[1:10]
sum(is.na(norm_counts_mtx$symbol))
## TO DO: computing mean values of gene symbol
stat <- norm_counts_mtx %>%
  tidyr::pivot_longer(cols = starts_with("SA"), names_to = 'sample', values_to='exp') %>%
  as.data.frame()
colnames(stat)
head(stat)
sum(is.na(stat$symbol))
sum(is.na(stat$sample))
sum(is.null(stat$symbol))
## For duplicated gene symbols
stat1 <- stat %>% 
  dplyr::group_by(sample, symbol) %>%
  dplyr::summarise(exp=mean(exp))
dim(stat1)
View(head(stat1))
dim(stat)
length(unique(stat1$symbol))
dim(comp51BB)
comp51BB <- data.table::fread(paste0(input_dir, 'bulk_SA919/deg_analysis_X7_introns/SA919_X7_introns_M51_SupraSpinal_B_vs_M51_Primary_B/de_significant_genes.csv.gz'))

comp52CB <- data.table::fread(paste0(input_dir, 'deg_analysis_X7_introns/SA919_X7_introns_M52_SupraSpinal_C_vs_M52_Primary_B/de_significant_genes.csv.gz'))
meta_samples <- data.table::fread(paste0(input_dir, 'library_groupings_bulk_SA919.csv'))
# meta_samples SA919X7XB05691 SupraSpinal M51 versus SA919X7XB05402 Primary M51
dim(comp52CB)
dim(comp51BB)
length(intersect(comp51BB$gene_symb, stat1$symbol))

bulk51BB <- stat1 %>%
  dplyr::filter(sample %in% c('SA919X7XB05691','SA919X7XB05402')) %>%
  tidyr::pivot_wider(names_from = 'sample', values_from='exp') %>%
  dplyr::filter(!is.na(symbol) & symbol!='' & symbol %in% comp51BB$gene_symb) %>%
  dplyr::mutate(bulk_log2FC=log2(SA919X7XB05691+1) - log2(SA919X7XB05402+1))
  

bulk52CB <- stat1 %>%
  dplyr::filter(sample %in% c('SA919X7XB05604','SA919X7XB05378')) %>%
  tidyr::pivot_wider(names_from = 'sample', values_from='exp') %>%
  dplyr::filter(!is.na(symbol) & symbol!='' & symbol %in% comp52CB$gene_symb) %>%
  dplyr::mutate(bulk_log2FC=log2(SA919X7XB05604+1) - log2(SA919X7XB05378+1))

View(head(bulk51BB))
sum(bulk51BB$symbol=='')
dim(bulk51BB)

comp51BB$avg_log2FC
comp52CB$avg_log2FC
colnames(comp51BB)

comp51BB <- comp51BB %>%
  dplyr::select(avg_log2FC, gene_symb)

bulk51BB <- bulk51BB %>%
  inner_join(comp51BB, by=c('symbol'='gene_symb'))


comp52CB <- comp52CB %>%
  dplyr::select(avg_log2FC, gene_symb)
bulk52CB <- bulk52CB %>%
  inner_join(comp52CB, by=c('symbol'='gene_symb'))

cor_val <- cor(bulk52CB$bulk_log2FC, bulk52CB$avg_log2FC, method='pearson')
cor_val


cor_val <- cor(bulk51BB$bulk_log2FC, bulk51BB$avg_log2FC, method='pearson')
cor_val

p <- ggplot(data=bulk51BB, aes(x=avg_log2FC,y=bulk_log2FC))+
      geom_point(size=1, color='purple')
p


p <- ggplot(data=bulk52CB, aes(x=avg_log2FC,y=bulk_log2FC))+
  geom_point(size=1, color='purple') + 
  geom_hline(yintercept=0, linetype="dashed", color = "red") + 
  geom_vline(xintercept=0, linetype="dashed", color = "red")
p

bulk52CB1 <- bulk52CB %>%
  dplyr::filter((avg_log2FC>0 & bulk_log2FC>0) | (avg_log2FC<0 & bulk_log2FC<0))
round(100*dim(bulk52CB1)[1]/dim(bulk52CB)[1],1)





library(dplyr)
input_dir <- '/Users/htran/Documents/storage_tmp/metastasis_trees/SA919_10x/'

# df <- data.table::fread(paste0(input_dir, 'bulk_SA919/GS_results/Final_Sig_g_IDnames_Meta_vs_Prim_GS_Passage7.csv')) 
df <- data.table::fread(paste0(input_dir, 'bulk_SA919/GS_results/Passage7_CloneCMeta_vs_CloneBPrim_DESeq_Sig_GS.csv')) 
dim(df)
sum(df$log2FoldChange>0)
summary(as.numeric(df$log2FoldChange))
df1 <- df %>%
  dplyr::filter(as.numeric(log2FoldChange)>=0.5)
dim(df1)

summary(as.numeric(df1$log2FoldChange))



comp52CB <- data.table::fread(paste0(input_dir, 'bulk_SA919/deg_analysis_X7_introns/SA919_X7_introns_M52_SupraSpinal_C_vs_M52_Primary_B/de_significant_genes.csv.gz'))
length(intersect(comp52CB$gene_symb, df$Gene.name))

comp52CB$gene_symb[grepl('MY',comp52CB$gene_symb)]

df$Gene.name[grepl('MY',df$Gene.name)]

library(SingleCellExperiment)
sce <- readRDS(paste0(input_dir, 'SA919_10x_introns_mode/intron_alignment/SCRNA10X_SA_CHIP0077_001.rdata'))
dim(sce)
metagenes <- as.data.frame(rowData(sce))
dim(metagenes)
length(intersect(metagenes$Symbol, df$Gene.name))
dim(df)
