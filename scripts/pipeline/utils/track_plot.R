
suppressPackageStartupMessages({
  library(tidyverse)
  library(annotables)
  library(dplyr)
  library(ggplot2)
  # library(scMerge)
})
cnv_fn <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_cys/clonealign/whole_data/SA535_whole_data.csv'
obs_clones=c('A','B')
additional_genes = NULL
input_deg_fn <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_cys/deg_analysis/SA535_A_B/de_significant_genes.txt'
save_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_cys/deg_analysis/SA535_A_B/genes_cn_correlation/'
annots <- annotables::grch37

annots <- dplyr::select(annots, ensembl_gene_id = ensgene, chr, start, end) 

cnv_mat_backup <- cnv_mat
cnv_mat$ensembl_gene_id <- rownames(cnv_mat)
# View(head(annots))
cnv_mat <- inner_join(cnv_mat, annots)
View(head(cnv_mat))
intersect_genes <- intersect(cnv_mat$ensembl_gene_id, genes_map_symb$ensembl_gene_id)
length(intersect_genes)

dim(cnv_mat)
cnv_mat_chr1 <- cnv_mat[cnv_mat$chr=="1",]
View(head(cnv_mat_chr1))
cnv_mat_chr1 <- cnv_mat_chr1[order(cnv_mat_chr1$start),]
dim(cnv_mat_chr1)
deg_df_backup <- deg_df
View(head(deg_df))
deg_df_chr1 <- deg_df[,colnames(deg_df) %in% c('ensembl_gene_id','gene_symbol')]
dim(deg_df_chr1)
deg_df_chr1 <- inner_join(deg_df_chr1, annots)
deg_df_chr1 <- deg_df_chr1[deg_df_chr1$chr=='1',]
deg_df_chr1 <- deg_df_chr1[order(deg_df_chr1$start),]

genes_use <- intersect(cnv_mat_chr1$ensembl_gene_id, deg_df_chr1$ensembl_gene_id)
length(genes_use)

saveRDS(cnv_mat, file=paste0(save_dir,'cnv_mat_AB.rds'))
colnames(cnv_mat)
View(head(cnv_mat))
genes_use <- intersect(rownames(cnv_mat),rownames(z_score))
length(genes_use)
obs_cnv_mat <- cnv_mat[genes_use,]
obs_score <- z_score[genes_use,]
dim(obs_cnv_mat)
colnames(obs_cnv_mat)
obs_cnv_mat$ens_gene <- rownames(obs_cnv_mat)

library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)
obs_cnv_mat_2 <- obs_cnv_mat %>%
  pivot_longer(!ens_gene, names_to = "clone", values_to = "cn")

dim(obs_cnv_mat_2)
View(head(obs_cnv_mat_2))
obs_cnv_mat_2$copy_number <- ifelse(obs_cnv_mat_2$cn==1, -2, 
                           ifelse(obs_cnv_mat_2$cn==2, -1,
                                  ifelse(obs_cnv_mat_2$cn==3, 0, 
                                         ifelse(obs_cnv_mat_2$cn==4, 1,
                                                ifelse(obs_cnv_mat_2$cn>=5, 2,NA)))))
summary(as.factor(obs_cnv_mat_2$copy_number))


dim(obs_score)
obs_cnv_matA <- obs_cnv_mat[,colnames(obs_cnv_mat)=="A"]
names(obs_cnv_matA) <- rownames(obs_cnv_mat)
obs_cnv_matA[1:3]
obs_cnv_matB <- obs_cnv_mat[,colnames(obs_cnv_mat)=="B"]
names(obs_cnv_matB) <- rownames(obs_cnv_mat)
obs_scoreA <- obs_score[,colnames(obs_score) %in% colnames(sceA)]
obs_scoreB <- obs_score[,colnames(obs_score) %in% colnames(sceB)]

obs_scoreA_med <- apply(obs_scoreA, 1, median)
obs_scoreB_med <- apply(obs_scoreB, 1, median)
length(obs_scoreA_med)
obs_scoreA_med[1:2]
genes_stat <- data.frame(z_score=c(as.numeric(obs_scoreA_med),as.numeric(obs_scoreB_med)),
                         ens_gene=c(names(obs_scoreA_med),names(obs_scoreB_med)),
                         clone=c(rep('A', length(obs_scoreA_med)),rep('B',length(obs_scoreB_med))))
dim(genes_stat)
obs_cnv_mat_22 <- inner_join(obs_cnv_mat_2, genes_stat, by=c('ens_gene','clone'))
dim(obs_cnv_mat_22)
View(head(obs_cnv_mat_22))
colnames(obs_cnv_mat_22)
summary(obs_cnv_mat_22$z_score)
summary(as.factor(obs_cnv_mat_22$copy_number))
dim(obs_cnv_mat_22)


ctotal <- cor(as.numeric(obs_cnv_mat_22$copy_number), as.numeric(obs_cnv_mat_22$z_score), method = 'pearson')
obs_cnv_mat_22A <- obs_cnv_mat_22[obs_cnv_mat_22$clone=='A',]
cA <- cor(as.numeric(obs_cnv_mat_22A$copy_number), as.numeric(obs_cnv_mat_22A$z_score),method='pearson')
obs_cnv_mat_22B <- obs_cnv_mat_22[obs_cnv_mat_22$clone=='B',]
cB <- cor(as.numeric(obs_cnv_mat_22B$copy_number), as.numeric(obs_cnv_mat_22B$z_score),method='pearson')

cAB <- cor(as.numeric(obs_cnv_mat_22A$copy_number), as.numeric(obs_cnv_mat_22B$z_score),method='pearson') #


obs_cnv_mat_22$copy_number <- as.factor(obs_cnv_mat_22$copy_number)
corr_plot <- grouped_boxplot_correlation(obs_cnv_mat_22, 'copy_number', 'z_score', 'clone', 
                                  plttitle="Genes Copy Number Correlation",
                                  xlabel="Copy number", ylabel="Median Z-score genes")
corr_plot
png(paste(save_dir,"gene_CN_correlation_A_B_SA535.png",sep="/"), height = 2*400, width=2*600, res = 2*72)
print(corr_plot)
dev.off()

obs_cnv_mat_22$copy_number <- as.numeric()
linearMod <- lm(as.numeric(copy_number) ~ as.numeric(z_score), data=obs_cnv_mat_22B)  # build linear regression model on full data
print(linearMod)
t <- summary(linearMod)
obs_cnv_mat_22$copy_number <- as.numeric(obs_cnv_mat_22$copy_number)
p <- ggplot(obs_cnv_mat_22A, aes(copy_number, z_score)) +
  geom_point() +
  stat_smooth(method = lm)



# TO DO: add significant level
signif_cn <- c()
cns <- unique(obs_cnv_mat_22$copy_number)
for(cn in cns){
  res <- t.test(z_score ~ clone, data=obs_cnv_mat_22[obs_cnv_mat_22$copy_number==cn,])
  signif_cn <- c(signif_cn, res$p.value < 0.05)
}
names(signif_cn) <- cns



stat <- data.frame(cnvA=as.numeric(obs_cnv_matA[genes_use]),
                   cnvB=as.numeric(obs_cnv_matB[genes_use]),
                   scoreA=as.numeric(obs_scoreA_med[genes_use]),
                   scoreB=as.numeric(obs_scoreB_med[genes_use]),
                   row.names = genes_use)
dim(stat)  


corr_plotA <- ggplot(stat, aes(x = cnvA, y = scoreA, colour = cnvA)) +
  geom_jitter(position=position_jitter(0.3), size=2)
corr_plotA
stat$cnvAA <- stat$cnvA
corr_plotA1 <- ggplot(stat, aes(x = cnvA, y = scoreB, colour = cnvA)) +
  geom_jitter(position=position_jitter(0.3), size=2)
corr_plotA1
corr_plotA <- boxplot_correlation(stat1, 'cnvA', 'scoreA', 'cnvA', 
                    plttitle="Genes CN correlation in clone A",
                    xlabel="Copy number in A", ylabel="Median Z-score genes in A")

corr_plotB <- boxplot_correlation(stat1, 'cnvB', 'scoreB', 'cnvB', 
                                  plttitle="Genes CN correlation in clone B",
                                  xlabel="Copy number in B", ylabel="Median Z-score genes in B")
png(paste(save_dir,"clusterA_gene_CN_correlation.png",sep="/"), height = 2*500, width=2*500, res = 2*72)
print(corr_plotA)
dev.off()

png(paste(save_dir,"clusterB_gene_CN_correlation.png",sep="/"), height = 2*500, width=2*500, res = 2*72)
print(corr_plotB)
dev.off()


library(plyr)
convert_values <- colwise(function(x){
  return(ifelse(x==1, -2, 
                ifelse(x==2, -1,
                       ifelse(x==3, 0, 
                              ifelse(x==4, 1,
                                     ifelse(x>=5, 2,x))))))})


stat1 <- stat[,colnames(stat) %in% c('cnvA','cnvB')]
stat1 <- convert_values(stat1)
summary(as.factor(stat1$cnvA))
stat1$scoreA <- stat$scoreA
stat1$scoreB <- stat$scoreB
M1[1:3,1:3]
stat1$cnvA <- as.factor(stat1$cnvA)
stat1$cnvB <- as.factor(stat1$cnvB)

cA <- cor(as.numeric(stat1$cnvA), stat1$scoreA, method = 'pearson')
cA1 <- cor(as.numeric(stat1$cnvA), stat1$scoreB, method = 'pearson')
cB1 <- cor(as.numeric(stat1$cnvB), stat1$scoreA, method = 'pearson')
cB <- cor(as.numeric(stat1$cnvB), stat1$scoreB, method = 'pearson')


cnv_mat_chr1_ext <- cnv_mat_chr1[!(cnv_mat_chr1$ensembl_gene_id %in% genes_use), ]
dim(cnv_mat_chr1_ext)
View(cnv_mat_chr1_ext[1:50,])

input_dir <- "/home/htran/storage/datasets/drug_resistance/dlp_results/"
genes_map_symb <- read.csv(paste0(input_dir, "biodatabase/meta_genes.csv"), header=T, stringsAsFactors=F)
rownames(genes_map_symb) <- genes_map_symb$gene_ens
colnames(genes_map_symb)[which(colnames(genes_map_symb)=='gene_ens')] <- 'ensembl_gene_id'
cnv_mat_chr1_ext <- inner_join(cnv_mat_chr1_ext, genes_map_symb)
cnv_mat_chr1_ext$distance <- cnv_mat_chr1_ext$end - cnv_mat_chr1_ext$start

write.csv(cnv_mat_chr1_ext,file=paste0(save_dir,'cnv_mat_chr1_ext.csv'), row.names=F, quote=F)
deg_df_chr1_ext <- deg_df_chr1[!(deg_df_chr1$ensembl_gene_id %in% genes_use),]
dim(deg_df_chr1_ext)
View(deg_df_chr1_ext[1:50,])
deg_df_chr1_ext$distance <- deg_df_chr1_ext$end - deg_df_chr1_ext$start

write.csv(deg_df_chr1_ext,file=paste0(save_dir,'deg_df_chr1_ext.csv'), row.names=F, quote=F)


sceA <- sce[,sce$clone_id =='A',]
dim(sceA)
sceB <- sce[,sce$clone_id =='B',]
dim(sceB)

rownames(sce1)[1:3]
genes_use <- intersect(rownames(sceA),cnv_mat$ensembl_gene_id)
length(genes_use)
sceA <- sceA[rownames(sceA) %in% genes_use,]
sceA <- sceA[genes_use,]
sceB <- sceB[rownames(sceB) %in% genes_use,]
sceB <- sceB[genes_use,]
lcA <- logcounts(sceA)
lcB <- logcounts(sceB)
dim(lcB)
class(lcA)

median_genesA <- apply(lcA, 1, median)
length(median_genesA)
median_genesB <- apply(lcB, 1, median)
length(median_genesB)
summary(median_genesA)
summary(median_genesB)
median_genesB[1:3]
View(head(cnv_mat))
rownames(cnv_mat) <- cnv_mat$ensembl_gene_id
cnv_mat_use <- cnv_mat[cnv_mat$ensembl_gene_id %in% genes_use,]
dim(cnv_mat_use)
View(head(cnv_mat_use))

cnv_mat_use <- cnv_mat_use %>% distinct(ensembl_gene_id, .keep_all = TRUE)
rownames(cnv_mat_use) <- cnv_mat_use$ensembl_gene_id
cnv_mat_use <- cnv_mat_use[genes_use,]
stat <- data.frame(gene_id=names(median_genesA),
                   median_genesA=as.numeric(median_genesA),
                   median_genesB=as.numeric(median_genesB),
                   cnA=cnv_mat_use$A,
                   cnB=cnv_mat_use$B, 
                   row.names = names(median_genesA))
summary()
cor(stat$stateA, stat$median_genesA)
cor(stat$stateB, stat$median_genesB)
cor(stat$stateA, stat$median_genesB)
cor(stat$stateB, stat$median_genesA)

stat_backup <- stat
summary(stat$cnA)
cnAA <- c()
convert_df <- data.frame(cnA=c(0,1,2,3,4,5,6,7,8,9),
                      stateA=c(-2,-1,0,1,1,1,2,2,2,2))
convert_dfB <- data.frame(cnB=c(0,1,2,3,4,5,6,7,8,9),
                         stateB=c(-2,-1,0,1,1,1,2,2,2,2))
stat <- inner_join(stat, convert_df, by='cnA')
stat <- inner_join(stat, convert_dfB, by='cnB')
dim(stat)
p <- scatter.smooth(x=stat$stateA, y=stat$median_genesA, main="Copy Number ~ Median genes")
p1 <- scatter.smooth(x=stat$stateA, y=stat$median_genesB, main="Copy Number ~ Median genes")

linearMod <- lm(stateA ~ median_genesA, data=stat)  # build linear regression model on full data
print(linearMod)
linearMod2 <- lm(stateA ~ median_genesB, data=stat)  # build linear regression model on full data
print(linearMod2)

boxplot_correlation <- function(df, xstring, ystring, plottype, 
                                plttitle="Clone A - genes A",
                                xlabel=" ",ylabel=" ") {
  
  p <- ggplot(df, aes_string(x=xstring, y=ystring, fill=plottype)) +
    # stat_boxplot( aes_string(x=xstring, y=ystring), geom='errorbar', linetype=1, width=0.4)+  #whiskers
    geom_boxplot( aes_string(x=xstring, y=ystring),outlier.shape=1) +    
    stat_summary(fun=mean, geom="point", size=1) #+ 
    # stat_summary(fun.data = mean_se, geom = "errorbar")
  p <- p + labs(x=xlabel, y=ylabel, title=plttitle)
  p <- p + theme(
    legend.title = element_text(size=18), 
    legend.text = element_text(size=16), 
    plot.title = element_text(color="black", size=20, hjust = 0.5)
  )
  return(p)
}

grouped_boxplot_correlation <- function(df, xstring, ystring, plottype, 
                                plttitle="Clone A - genes A",
                                xlabel=" ",ylabel=" ") {
  
  p <- ggplot(df, aes_string(x=xstring, y=ystring, fill=plottype)) +
    stat_boxplot( aes_string(x=xstring, y=ystring), geom='errorbar', 
                  linetype=1, width=0.2, position=position_dodge(.9)) +  #whiskers
    geom_boxplot( aes_string(x=xstring, y=ystring, fill=plottype),
                  outlier.shape=1, position=position_dodge(.9)) +    
    stat_summary(fun=mean, geom="point", size=1, aes_string(group=plottype),
                 position=position_dodge(.9)) + 
    stat_summary(fun.data = mean_se, geom = "errorbar",
                 position=position_dodge(.9))
  p <- p + labs(x=xlabel, y=ylabel, title=plttitle)
  p <- p + theme(
    legend.title = element_text(size=18), 
    legend.text = element_text(size=16), 
    plot.title = element_text(color="black", size=20, hjust = 0.5)
  )
  return(p)
}

library("ggpubr")
dev.off()
p <- ggscatter(stat, x = "stateA", y = "median_genesA", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "stateA", ylab = "median_genesA")





plot_CNV(input_deg_fn, cnv_fn, obs_clones, sample_name='SA535',
                     additional_genes = NULL, n_genes_to_annotate = 60)


rm(obs_clones)
obs_clones=c('A','G')
input_deg_fn <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_cys/deg_analysis/SA535_A_G/de_significant_genes.txt'

rm(obs_clones)
obs_clones=c('DE','B')
input_deg_fn <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_cys/deg_analysis/SA535_DE_B/de_significant_genes.txt'

rm(obs_clones)
obs_clones=c('DE','G')
input_deg_fn <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_cys/deg_analysis/SA535_DE_G/de_significant_genes.txt'


plot_CNV(input_deg_fn, cnv_fn, obs_clones, sample_name='SA535',
         additional_genes = NULL, n_genes_to_annotate = 60)




plot_CNV <- function(input_deg_fn, cnv_fn, obs_clones=c('A','B'), sample_name='',
                     additional_genes = NULL, n_genes_to_annotate = 50){
  save_dir <- paste0(dirname(input_deg_fn),'/')
  if (!file.exists(save_dir)){
    dir.create(save_dir)
  }
  df_cnv <- read.csv(cnv_fn, check.names=F, stringsAsFactors=F)
  print(summary(as.factor(df_cnv$cluster)))
  
  
  
  if(sample_name=='SA535'){
    clone_names <- data.frame(cluster=c("A", "B", "C", "D", "E", "F"), 
                              clone=c("A", "B", "C", "DE", "F", "G"))
    df_cnv <- df_cnv %>% left_join(clone_names, by = "cluster")
  } else{
    colnames(df_cnv)[which(colnames(df_cnv)=='cluster')] <- 'clone'
  }
  # View(head(df_cnv))
  cnv_mat <- get_variance_genes_genome(df_cnv, obs_clones, 0.5)
  # dim(cnv_mat)
  # View(head(cnv_mat))
  track_plot2 <- plot_DE_genes_chr(save_dir, cnv_mat, input_deg_fn, obs_clones, 
                                   sample_name, NULL, n_genes_to_annotate)
  
  
  annots <- annotables::grch37
  
  annots <- dplyr::select(annots, ensembl_gene_id = ensgene, chr, start, end) 
  
  df_cnv <- inner_join(df_cnv, annots)
  
  df_cnv <- dplyr::select(df_cnv, ensembl_gene_id, chr, start) %>% 
    group_by(chr) %>% 
    dplyr::mutate(start_order = rank(start)) %>% 
    ungroup() %>% 
    inner_join(df_cnv)
  
  chr_levels <- c(as.character(1:23), "X", "Y")
  
  df_cnv$chr <- factor(df_cnv$chr, levels = chr_levels)
  # print(summary(as.factor(df_cnv$chr)))
  df_cnv <- drop_na(df_cnv)
  df_cnv$cnv <- as.character(round(df_cnv$mean_cnmode))
  df_cnv$cnv <- as.factor(df_cnv$cnv)
  
  
  # summary(as.factor(df_cnv$cnv))
  
  # MA: 11 Apr 2020, setting the copy number colors that were used in the heatmap
  # TO COME BACK
  #print("data frame")
  #cnv_colors <- c('#4880B8', '#A7C9DF','#CCCCCC','#F5CE93','#ED9364','#D2553E','#A42116','#8B1A43','#CB3576','#D06CAD','#C196C4','#D0BAD8')
  #cnv_cols <- data.frame(cn=0:11,color<-cnv_cols)
  #colnames(cnv_cols) <- c("cn","color")
  #levels(cnv_cols$color) <- cnv_colors
  
  # cnv_colors <- c('#4880B8', '#A7C9DF','#CCCCCC','#F5CE93','#ED9364',
  #                 '#D2553E','#A42116','#8B1A43','#CB3576','#D06CAD','#C196C4','#D0BAD8')
  # 
  cnv_cols <- structure(
    c(
      "#3182BD", "#9ECAE1", "#CCCCCC", "#FDCC8A", "#FC8D59", "#E34A33",
      "#B30000", "#980043", "#DD1C77", "#DF65B0", "#C994C7", "#D4B9DA"
    ),
    names=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")
  )
  # cnv_cols <- c('0'='#4880B8', '1'='#A7C9DF','2'='#CCCCCC','3'='#F5CE93','4'='#ED9364',
  #               '5'='#D2553E','6'='#A42116','7'='#8B1A43','8'='#CB3576','9'='#D06CAD','10'='#C196C4','11'='#D0BAD8')
  #levels(cnv_cols) <- 0:11
  #  cnv_cols <- c("0" = "#2166ac",
  #                "1" = "#92c5de", 
  #                "2" = "grey80", 
  #                "3" = "#f4a582", 
  #                "4" = "#d6604d",
  #                "5" = "#b2182b",
  #                "6+" = "#67001f")
  
  
  # MA: removing the 6+ restriction (clonealign still has that restriction though)
  # df_cnv$cnv[round(df_cnv$median_cnmode) >= 6] <- "6+"
  #!chr %in% c("X","Y")
  # View(head(df_cnv))
  
  levels(df_cnv$cnv) <- 0:(length(cnv_cols)-1)
  
  cnv_plot <- dplyr::filter(df_cnv, clone %in% obs_clones, chr !="Y") %>% #%in% c("X","Y")
    ggplot(aes(x = start_order, y = clone, fill = cnv)) +
    geom_raster() +
    facet_wrap(~ chr, scales = "free_x", nrow = 1, switch = "x") +
    # theme(legend.position = "top", axis.text.x = element_blank()) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_manual(values=cnv_cols, name = "CN", guide = 'legend',
                      labels = 0:(length(cnv_cols)-1),drop=FALSE)  +
    theme(legend.position = "bottom") +
    # scale_fill_manual(values = cnv_cols, name = "Copy Number", labels = 0:(length(cnv_cols)-1),drop=FALSE) +  #, labels = 0:(length(cnv_cols)-1),drop=FALSE) +
    labs(x = "chromosome") +
    theme(strip.background = element_rect(fill = 'white', colour = 'white'),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size=10), 
          strip.placement = "outside",
          legend.text = element_text(size=10),
          legend.title = element_text(size=11),
          # legend.position = "none",
          panel.spacing = unit(c(0.1), 'cm')) +   # MA: was 0.2
    theme(text = element_text(size = 18)) +
    labs(x = "Chromosome", y = "Clone")+
    guides(fill = guide_legend(nrow = 1)) +
    scale_x_continuous(expand = c(0,0))
  
  # cnv_plot
  main_plot <- cowplot::plot_grid(
    track_plot2 + theme(axis.title.x = element_blank(),
                        strip.text.x = element_blank()),
    cnv_plot,
    ncol = 1,
    rel_heights = c(3.5,1),
    axis = "b",
    align = 'v'
  )
  
  # main_plot
  base_name <- paste(obs_clones, collapse = "_vs_")
  png(paste0(save_dir,"DE_cnv_expression_",sample_name,"_",base_name,".png"), height = 2*600, width=2*1200,res = 2*72)
  print(main_plot)
  dev.off()
  
  
}


get_variance_genes_genome <- function(df_cnv, obs_clones=c('A','B'),
                                      default_threshold=0.5){
  
  cnv <- dplyr::filter(df_cnv, use_gene) %>%
    dplyr::rename(copy_number=median_cnmode) %>% 
    dplyr::select(ensembl_gene_id, clone, copy_number) %>% 
    dplyr::filter(clone %in% obs_clones) %>%
    spread(clone, copy_number)
  
  cnv_mat <- cnv %>%
    as.data.frame %>%
    column_to_rownames("ensembl_gene_id") 
  
  # third_quantile <- as.numeric(quantile(t, probs = 0.75))
  # print(third_quantile)
  # if(third_quantile > default_threshold){
  #   default_threshold = third_quantile
  # }
  
  var_genes <- apply(cnv_mat, 1, var)
  cnv_mat <- cnv_mat[var_genes >= default_threshold,]
  print(colnames(cnv_mat))
  print(dim(cnv_mat))
  return(cnv_mat)
}



