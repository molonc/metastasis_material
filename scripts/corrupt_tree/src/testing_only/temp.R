results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_new_encode/SA1035_Tyler_v2/'
# results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535/SA535_total/'
cell_clones <- read.delim(paste0(results_dir,'cell_clones.csv'), stringsAsFactors = FALSE, sep = ",")
length(unique(cell_clones$clone_id))



# results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_cys/snakemake/cell_clones_Sohrab/'
# # results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535/SA535_total/'
# cell_clones <- read.delim(paste0(results_dir,'cell_clones_SA535.csv'), stringsAsFactors = FALSE, sep = ",")
# summary(as.factor(cell_clones$clone_id))

unique(cell_clones$clone_id)



results_dir <- '/home/htran/storage/datasets/metastasis_results/SA535_new_encoding/SA535_wholedata_v2/'
cell_clones <- read.delim(paste0(results_dir,'corrupt_grow/19_clones/cell_clones.csv'), stringsAsFactors = FALSE, sep = ",")
dim(cell_clones)
colnames(cell_clones)

clone_assignment <- read.delim(paste0(results_dir,'corrupt_grow/clone_assignment_v2.csv'), sep = ",")
View(clone_assignment)
dim(clone_assignment)
length(unique(clone_assignment$clone_label))

results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_wholedata/'
cell_clones <- read.delim(paste0(results_dir,'cell_clones_v1.csv'), stringsAsFactors = FALSE, sep = ",", check.names=F)

summary(as.factor(cell_clones$clone_id))

ex_clones <- c('I','D','G','K','A','L')
clones_ls <- list()
for(c in ex_clones){
  print(c)
  tmp <- read.delim(paste0(results_dir,'corrupt_grow/clone_',c,'/cell_clone_',c,'.csv'), stringsAsFactors = FALSE, sep = ",")
  print(dim(tmp))
  tmp <- tmp[,colnames(tmp) != 'clone_label']
  colnames(tmp)[which(colnames(tmp) == "clone_id")] <- "clone_label"
  tmp$clone_label <- ifelse(tmp$clone_label=="A",paste0(c,'1'),
                            ifelse(tmp$clone_label=="B",paste0(c,'2'),NA))
  print(summary(as.factor(tmp$clone_label)))
  clones_ls[[c]] <- tmp
}
sub_divided_clones <- do.call(rbind, clones_ls)
dim(sub_divided_clones)
colnames(sub_divided_clones)[which(names(sub_divided_clones) == "clone_label")] <- "clone_id"
sub_cell_clones <- cell_clones[!cell_clones$clone_id %in% ex_clones,]

cells_clone_v2 <- rbind(sub_cell_clones, sub_divided_clones)
dim(cells_clone_v2)
write.csv(cells_clone_v2,paste0(results_dir,'corrupt_grow/cell_clones.csv'), row.names = F, quote=F)


colnames(cells_clone_v2)
colnames(clone_assignment)
sum(clone_assignment$clone_id %in% unique(cells_clone_v2$clone_id))
head(cell_clones)
cells_clone_v2 <- cells_clone_v2 %>% inner_join(clone_assignment, by = "clone_id")

cells_clone_v3 <- cells_clone_v2
cells_clone_v3 <- cells_clone_v3 %>% inner_join(clone_assignment, by = "clone_id")
cells_clone_v3 <- cells_clone_v3[,colnames(cells_clone_v3) != 'clone_id']
colnames(cells_clone_v3)[which(names(cells_clone_v3) == "clone_label")] <- "clone_id"
dim(cells_clone_v3)
length(unique(cells_clone_v3$clone_id))
write.csv(cells_clone_v3, paste0(results_dir,'corrupt_grow/cell_clones_total.csv'), row.names = F, quote=F)

cell_clones$clone_id <- ifelse(cell_clones$clone_id=="D","C",cell_clones$clone_id)

typeof(cell_clones$clone_id)

cell_clones <- cell_clones[,colnames(cell_clones) != 'clone_label']
cell_clones <- cell_clones[,colnames(cell_clones) != 'clone_id']
cell_clones$clone_label <- cell_clones$clone_id

colnames(cell_clones)[which(names(cell_clones) == "clone_id")] <- "clone_label"



clone_id <- LETTERS[seq(length(unique(cell_clones$clone_id)))]
t = summary(as.factor(cell_clones$clone_id))
clone_names <- data.frame(clone_id=clone_id, clone_label=names(t))


clone_names <- data.frame(clone_id=c('A','K','G','I','F','H','C','D','E','B','J','M','O','N','L'), 
                       clone_label=c('B','B','C','C','C','D','D','D','E','A','F','G','G','G','G'))

clone_names <- data.frame(   clone_id=c('C','B','F','E','D','A','G'), 
                          clone_label=c('A','B','B','B','C','C','D'))

clone_names <- data.frame(      clone_id=c('A','B','C'), 
                             clone_label=c('A','A','B'))


clone_names <- data.frame(clone_label=c('M','K','D','F','I','C','E','L','J','A','H','B','G'), 
                          clone_id=c('A','A','A','A','B','C','C','D','E','F','F','G','H'))

clone_names <- data.frame(clone_label=c("A", "B", "C", "D", "E", "F", "G","H"), 
                             clone_id=c("A", "B", "C", "B", "B", "A", "None","A"))

clone_names <- data.frame(clone_label=c("A", "B", "C", "D", "E", "F","G", "H", "I", "J"), 
                             clone_id=c("A", "C", "B", "B", "B", "None", "B", "A", "A", "C"))

clone_names <- data.frame(clone_label=c("A", "B", "C", "D", "E", "F","G", "H", "I", "J",
                                        "K", "L", "M", "N", "O", "P","Q","R","S"), 
                             clone_id=c("A", "C", "B", "B", "B", "None", "B", "A", "A", "C"))

View(clone_names)
library(dplyr)
cell_clones <- cell_clones %>% inner_join(clone_names, by = "clone_label")
View(head(cell_clones))
cell_clones <- cell_clones[,colnames(cell_clones) != 'clone_label']
dim(cell_clones)
# colnames(cell_clones)[which(names(cell_clones) == "clone_label")] <- "clone_id"
write.csv(cell_clones,paste0(results_dir,'cell_clones.csv'), row.names = F, quote=F)
cell_clones <- cell_clones[cell_clones$clone_id!='None',]
dim(cell_clones)

# library(DOSE)
# data(geneList)
# de <- names(geneList)[1:100]
# x <- enrichDO(de)
# emapplot(x)
# class(x)

# xvfb-run -a python /home/htran/Projects/fitness_material/de_scripts/pathway_analysis_v2.py --csv total_markers_SA535_A_B.csv --minfc 0.25 --adjp 0.05 --minedge 1 --png results_pathways/SA535_A_B.png




