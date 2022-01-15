input_dir <- "/home/htran/storage/datasets/metastasis_results/SA535_new_encoding/Tyler_whole/"
data_dir <- "/home/htran/storage/datasets/hakwoo_metastasis/SA535/"
t <- read.table(paste0(input_dir,'log.txt'),sep = '\n')
View(t)
bin_ids <- strsplit(as.character(t$V1),"_")
# length(bin_ids)
nbcores <- 5
ps <- c()
libs <- parallel::mclapply(bin_ids, function(f) {
  ps <- c(ps, as.character(f[1]))
}, mc.cores = nbcores)
if(length(libs)!=length(bin_ids)){
  print("Error, double check")
}

libs <- as.character(libs)


meta <- read.csv(paste0(data_dir,'SA535_metadata.csv'), stringsAsFactors=F, check.names=F)
head(meta)

data_dir <- '/home/htran/storage/datasets/hakwoo_metastasis/'
meta <- read.csv(paste0(data_dir,'SA919X7_metadata_Hoa.csv'), stringsAsFactors=F, check.names=F)
results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_filtered/'

grouping_df <- read.csv(paste0(results_dir,'library_groupings.csv'),
                        header=T,check.names=F, stringsAsFactors=F)

total_libs <- meta$`Library ID_DLP`
total_libs[!total_libs %in% libs]
curr_lib <- grouping_df$library_labels

curr_lib <- c('A96201B','A95721B','A95731A','A98232A','A98261A','A96121A','A96253A',
              'A98277A','A98269A','A98277B','A98165A','A98197A','A98197B',
              'A98194B','A98194A','A98169A','A98174A','A98282B','A98294B',
              'A98296A','A95697A','A98163A','A98262A','A98303B','A98262B')
rownames(meta) <- meta$`Library ID_DLP`
meta <- meta[curr_lib,]
length(curr_lib)
dim(meta)
meta$Site_origin
meta$`Library ID_DLP`
meta1 <- meta[,'Site_origin', drop=F]
summary(as.factor(meta1$Site_origin))


input_2 <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_filtered/'
features_df <- data.table::fread(paste0(input_dir,'corrupt_tree_features.csv'), 
                                 stringsAsFactors = F, check.names = F)

features_df2 <- data.table::fread(paste0(input_2,'corrupt_tree_features.csv'), 
                                 stringsAsFactors = F, check.names = F)

state_df <- data.table::fread(paste0(input_dir,'bin_cnvs_corrupt_double_padding.csv'), 
                                 stringsAsFactors = F, check.names = F)

features_df[1:4,1:3]
features_df2[1:4,1:3]
dim(features_df)
rownames(state_df)[1:3]
state_df$end[1:3]
