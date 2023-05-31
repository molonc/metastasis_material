
data_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/snakemake_dlp/"
sample_file <- "library_groupings.csv"
df <- read.csv(paste0(data_dir, sample_file), header=T, check.names=F, stringsAsFactors=F)
dim(df)
# 19 libraries

lib_file <- list.dirs(data_dir, recursive = F)
lib_file <- grep('A*[A|B]',lib_file, value=T)
lib_file <- gsub('/home/htran/storage/datasets/drug_resistance_DLP/SA535//','', lib_file)
dim(df)
length(unique(lib_file))
rownames(df) <- df$library_DLP
df <- df[lib_file,]

colnames(df)[which(names(df) == "library_DLP")] <- "libraryID"

length(unique(df$library_DLP))
df$id <- rownames(df)
data_dup <- df[gdata::duplicated2(df$library_DLP),]
names_dup <- data_dup$id
df$library_DLP <- df$libraryID
libs <- unique(df$library_DLP)
length(libs)
write(libs, paste0(data_dir,'library_ids.txt'), sep = "\t")
# str <- "["
# for(l in libs){
#   str = paste0(str,"'",l,"',")
# }
# str = paste0(str,"]")
# str


df1 <- df[!df$id %in% names_dup,]
dim(df1)
sample <- subset(df, select = c(mouse_id, PDX, passage,library_DLP,treatment_status))
colnames(df) 
dim(sample)
# rownames(sample) <- sample$mouse_id
# sample <- sample[which (sample$PDX=="SA1035"),]
# colnames(sample)[which(names(sample) == "Treatment Status")] <- "Treatment_Status"
write.csv(sample, file = paste0(data_dir, sample_file), quote=F, row.names = F)


lib_ids <- strsplit(as.character(sample$library_DLP),",")
length(lib_ids)

lib_ids[12][[1]]
colnames(sample)
err_lib <- c()
mouse_id <- c()
passage <- c()
treatment_status <- c()
PDX <- c()
library_DLP <- c()
for(i in rep(1:length(lib_ids),1)){
  if(length(lib_ids[i][[1]])==2){
    print(sample$library_DLP[i])
    err_lib <- c(err_lib,sample$library_DLP[i])
    mouse_id <- c(mouse_id, sample$mouse_id[i],sample$mouse_id[i])
    library_DLP <- c(library_DLP, lib_ids[i][[1]][1],lib_ids[i][[1]][2])
    passage <- c(passage, sample$passage[i],sample$passage[i])
    PDX <- c(PDX, sample$PDX[i],sample$PDX[i])
    treatment_status <- c(treatment_status,sample$treatment_status[i],sample$treatment_status[i])
  }
}

sample2 <- data.frame(mouse_id=mouse_id, PDX=PDX, passage=passage,
                      library_DLP=library_DLP, treatment_status=treatment_status, stringsAsFactors=F)

sample1 <- sample[!sample$library_DLP %in% err_lib, ]
dim(sample1)
colnames(sample1)
colnames(sample2)
meta_df <- rbind(sample1, sample2)
dim(meta_df)
# meta_df$library_DLP[meta_df$library_DLP=="A62397A "] <- "A62397A"
colnames(meta_df)[which(names(meta_df) == "mouse_id")] <- "sample"
colnames(meta_df)[which(names(meta_df) == "treatment_status")] <- "treatmentSt"
meta_df$libraryID <- meta_df$library_DLP
View(meta_df)
rownames(meta_df) <- meta_df$library_DLP
dim(meta_df)
length(unique(meta_df$library_DLP))
gdata::duplicated2(as.character(meta_df$library_DLP))
rownames(meta_df) <- NULL
View(meta_df)
write.csv(meta_df, file = paste0(data_dir, 'SA535_meta_data_v3.csv'), quote=F, row.names = F)


library(yaml)
library(gdata)
df$library_DLP <- df$libraryID
df <- df[,!colnames(df) %in% 'branch']
length(unique(df$library_DLP))
dim(df)
# df1 <- df[duplicated2(df$library_DLP),]
out <- yaml::as.yaml(list(samples=split(replace(sample, "mouse_id", NULL), sample$mouse_id)))
cat(out)

out <- yaml::as.yaml(list(library=split(replace(df, "library_DLP", NULL), df$library_DLP)))
cat(out)
save_dir <- "/home/htran/storage/datasets/drug_resistance/dlp_results/SA535/SA535_Tyler_v2/config/"
write.table(out, paste0(data_dir,"yaml_dlp_SA609.txt"),sep="\t",quote=F,col.names=NA)


