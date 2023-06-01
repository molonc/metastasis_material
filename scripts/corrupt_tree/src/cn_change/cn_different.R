# windowsFonts()

# get input median CN profile data 

# get clones at the same branch area, ex: group 1: A, B, C, D,...

# Compare clonal profiles
suppressPackageStartupMessages({
  # library(ggrastr)
  library(tidyr)
  library(dplyr)
  library(ggplot2)
})
# options(scipen=999)

# library(extrafont)
# font_import(prompt=F) # import all your fonts
# fonts() #get a list of fonts
# # fonttable()
# 
# b <- ggplot(mtcars, aes(x=wt, y=mpg)) + geom_point() +
#   ggtitle("Fuel Efficiency of 32 Cars") +
#   xlab("Weight (x1000 lb)") + ylab("Miles per Gallon") #+
#   # theme(text=element_text(size=16,  family="Comic Sans MS"))
# ggsave(
#   filename = paste0(results_dir,"testing.pdf"),
#   plot = b,
#   height = 4,
#   width = 6,
#   useDingbats=F)



source("/home/htran/Projects/hakwoo_project/corrupt_tree/src/cn_change/utils.R")
results_dir <- '/home/htran/storage/datasets/metastasis_results/SA535_new_encoding/SA535_wholedata_v2/'
save_dir <- paste0(results_dir,'CN_profile/')
dir.create(save_dir)
median_cnv <- read.csv(paste0(save_dir,'median_cnv.csv'), check.names = F, stringsAsFactors = F)
meta_cells <- read.csv(paste0(save_dir,'meta_cells.csv'), check.names = F, stringsAsFactors = F)

head(median_cnv)
dim(meta_cells)
View(meta_cells)

df_cnv <- median_cnv
# How to add number of cells here 
# cl1 <- 'A'
# cl2 <- 'B'
# t1 <- 'Primary'
# t2 <- 'Metastasis'
# desc1 <- paste0(tolower(t1)," cells (n=",meta_cells[meta_cells$clone_id==cl1 & meta_cells$mainsite==t1,'nb_cells'],")")
# desc2 <- paste0(tolower(t2)," cells (n=",meta_cells[meta_cells$clone_id==cl2 & meta_cells$mainsite==t2,'nb_cells'],")")
# print(desc1)
# print(desc2)
# viz_cn_change(median_cnv, c(paste0(cl1,'_',t1),desc1), c(paste0(cl2,'_',t2),desc2), paste0(save_dir,cl1,'_',t1,'_vs_',cl2,'_',t2,'/'))

# for SA535, it would be great if you can add comparison between primary vs metastasis for each of clone G and H.
cl1_ls <- c("A","E","N","O","K","S","G","H")
cl2_ls <- c("B","F","O","P","L","T","G","H")
t1 <- c("Primary","Metastasis","Primary","Metastasis","Metastasis","Primary","Primary","Primary")
t2 <- c("Metastasis","Metastasis","Metastasis","Metastasis","Metastasis","Metastasis","Metastasis","Metastasis")
unique(meta_cells$mainsite)

cl1_ls <- c("N",     "D",         "G",         "K")
cl2_ls <- c("O",     "E",         "H",         "L")
t1 <- c("Primary",   "Primary",   "Primary",   "Primary")
t2 <- c("Metastasis","Metastasis","Metastasis","Metastasis")
de <- data.frame(cl1=cl1_ls,cl2=cl2_ls,t1=t1,t2=t2)
detect_CN_change(de, median_cnv, meta_cells, save_dir)
  

de <- data.frame(cl1="K",cl2="L",t1="Metastasis",t2="Metastasis")
detect_CN_change(de, median_cnv, meta_cells, save_dir)

de <- data.frame(cl1="G",cl2="H",t1="Primary",t2="Metastasis")
detect_CN_change(de, median_cnv, meta_cells, save_dir)

colnames(median_cnv)

# viz_cn_change(median_cnv, c("E","Metastasis"), c("F","Metastasis"), save_dir)
# viz_cn_change(median_cnv, c("N","Primary"), c("O","Metastasis"))
# viz_cn_change(median_cnv, c("O","Metastasis"), c("P","Metastasis"))
# viz_cn_change(median_cnv, c("K","Metastasis"), c("L","Metastasis"))
# viz_cn_change(median_cnv, c("S","Primary"), c("T","Metastasis"))
# 
# rEF <- readRDS(paste0(save_dir,'F_de_E_vs_F_mediancn.rds'))
# ptEF_535 <- get_related_pathways(as.data.frame(rEF$cn_change), pathway)
# write.csv(ptEF_535, paste0(save_dir,'ptEF_535.csv'), quote = F, row.names = F)
# 
# rAB <- readRDS(paste0(save_dir,'B_de_A_vs_B_mediancn.rds'))
# ptAB_535 <- get_related_pathways(as.data.frame(rAB$cn_change), pathway)
# write.csv(ptAB_535, paste0(save_dir,'ptAB_535.csv'), quote = F, row.names = F)
# 
# rST <- readRDS(paste0(save_dir,'T_de_S_vs_T_mediancn.rds'))
# ptST_535 <- get_related_pathways(as.data.frame(rST$cn_change), pathway)
# write.csv(ptST_535, paste0(save_dir,'ptST_535.csv'), quote = F, row.names = F)
# 
# rNO <- readRDS(paste0(save_dir,'O_de_N_vs_O_mediancn.rds'))
# ptNO_535 <- get_related_pathways(as.data.frame(rNO$cn_change), pathway)
# write.csv(ptNO_535, paste0(save_dir,'ptNO_535.csv'), quote = F, row.names = F)
# 
# rPO <- readRDS(paste0(save_dir,'P_de_O_vs_P_mediancn.rds'))
# ptPO_535 <- get_related_pathways(as.data.frame(rPO$cn_change), pathway)
# write.csv(ptPO_535, paste0(save_dir,'ptPO_535.csv'), quote = F, row.names = F)
# 
# 
# rKL <- readRDS(paste0(save_dir,'L_de_K_vs_L_mediancn.rds'))
# ptKL_535 <- get_related_pathways(as.data.frame(rKL$cn_change), pathway)
# write.csv(ptKL_535, paste0(save_dir,'ptKL_535.csv'), quote = F, row.names = F)
# 
# 
# results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_wholedata/prevalences_SA919_02_Nov/'
# median_cnv <- read.csv(paste0(results_dir,'prevalences_SA535/median_cnv.csv'), check.names = F, stringsAsFactors = F)
# save_dir <- paste0(results_dir,'CN_profile/')
# dir.create(save_dir)





results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_wholedata/'
save_dir <- paste0(results_dir,'CN_profile/')
dir.create(save_dir)
copynumber_fn <- paste0(results_dir,'total_merged_filtered_states.csv')
datatag <- 'SA919'
cellclone_fn <- paste0(results_dir,'cell_clones.csv')
library_grouping_fn <- paste0(results_dir,'library_groupings.csv')
get_median_genotype_v3(library_grouping_fn, copynumber_fn, cellclone_fn, datatag, save_dir)

median_cnv <- read.csv(paste0(save_dir,'median_cnv.csv'), check.names = F, stringsAsFactors = F)
dim(median_cnv)
head(median_cnv)

meta_cells <- read.csv(paste0(save_dir,'meta_cells.csv'), check.names = F, stringsAsFactors = F)
head(meta_cells)

cl1_ls <- c("A","A","A","B","B","A","A")
cl2_ls <- c("B","B","B","C","C","C","C")
t1 <- c("Primary","Primary","Metastasis","Primary","Metastasis","Primary","Metastasis")
t2 <- c("Primary","Metastasis","Metastasis","Metastasis","Metastasis","Metastasis","Metastasis")
de <- data.frame(cl1=cl1_ls,cl2=cl2_ls,t1=t1,t2=t2)


cl1_ls <- c("A")
cl2_ls <- c("A")
t1 <- c("Primary")
t2 <- c("Metastasis")
de <- data.frame(cl1=cl1_ls,cl2=cl2_ls,t1=t1,t2=t2)


detect_CN_change(de, median_cnv, meta_cells, save_dir)


# viz_cn_change(median_cnv, c("A","Primary_Metastasis"), c("B","Primary_Metastasis"))
# viz_cn_change(median_cnv, c("B","Primary_Metastasis"), c("C","Metastasis"))
# viz_cn_change(median_cnv, c("A","Primary_Metastasis"), c("C","Metastasis"))
# 
# 
# tAB <- readRDS(paste0(save_dir,'B_de_A_vs_B_mediancn.rds'))
# ptwAB_919 <- get_related_pathways(as.data.frame(tAB$cn_change), pathway)
# write.csv(ptwAB_919, paste0(save_dir,'ptwAB_919.csv'), quote = F, row.names = F)
# 
# tBC <- readRDS(paste0(save_dir,'C_de_B_vs_C_mediancn.rds'))
# ptwBC_919 <- get_related_pathways(as.data.frame(tBC$cn_change), pathway)
# write.csv(ptwBC_919, paste0(save_dir,'ptwBC_919.csv'), quote = F, row.names = F)
# 
# tAC <- readRDS(paste0(save_dir,'C_de_A_vs_C_mediancn.rds'))
# ptwAC_919 <- get_related_pathways(as.data.frame(tAC$cn_change), pathway)
# write.csv(ptwAC_919, paste0(save_dir,'ptwAC_919.csv'), quote = F, row.names = F)












