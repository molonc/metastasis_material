suppressPackageStartupMessages({
  require("data.table")
  require("dplyr")
  require("scales")
})
script.basename <- '/home/htran/storage/datasets/metastasis_results/snv_analysis_SA919/ascn_SA919_Hoa/scripts/'
source(paste0(script.basename, "schnapp_utils.R"))
source(paste0(script.basename, "utils_heatmap.R"))


# Load SA919 data 
results_pseudobk_dir <- '/home/htran/storage/datasets/metastasis_results/snv_analysis_SA919/ascn_SA919_Hoa/'
results_dir <- results_pseudobk_dir
save_dir <- paste0(results_dir,'ascn_viz/')
if(!dir.exists(save_dir)){
  dir.create(save_dir)
}


ascn_fn = paste0(results_pseudobk_dir,'ascn_inference.rds')
## NOTE: allele object contains input data for all libraries, and columns: 
# [1] 37199267       19
# [1] "chr"             "start"           "end"            
# [4] "cell_id"         "state"           "copy"           
# [7] "alleleA"         "alleleB"         "totalcounts"    
# [10] "BAF"             "state_min"       "Min"            
# [13] "Maj"             "state_AS_phased" "state_AS"       
# [16] "LOH"             "phase"           "state_phase"    
# [19] "state_BAF"


## Data from copy number analysis
datatag <- 'SA919'
cellclone_fn <- paste0(results_dir,'cell_clones.csv') # cells assigned to clones
library_grouping_fn <- paste0(results_dir,'library_groupings.csv') # metadata file
copynumber_fn <- paste0(results_dir,'total_merged_filtered_states.csv.gz') # copy number states file with rownames are chr_start_end, and column names are cell_ids


## Loading data, extracting common cells between ascn object, and copy number cell_ids
## ascn cell_ids contains less cells profiles compared to copy number cell profiles. 
res <- get_median_LOH_genotype(ascn_fn, 
                               library_grouping_fn, cellclone_fn, 
                               datatag, save_dir)


## Form comparisons, and detect the LOH changes between pair of clones. 
# res <- readRDS(paste0(save_dir,'median_genotype.rds'))
detect_LOH_change(copynumber_fn, res, save_dir, datatag)





## Summary events
# Question: “”" is there a difference in LOH/CNLOH between metastasis and primary clones. 
# a matrix could be constructed to show pairwise differences between clones 
# - also, are any purely seen in metastasis ?“”"
# You may want to further break down the LOH into 3 classes
# deletion LOH - LOH and copy 1
# neutral LOH - LOH and copy 2
# amplified LOH - LOH and copy > 2

# dim(ascn2)
# head(ascn2)
# cl1 <- 'A'
# cl2 <- 'B'
# ascn_tmp <- ascn2 %>%
#           dplyr::filter(clone_id %in% c(cl1, cl2))
# ascn_tmp <- ascn_tmp %>%
#             tidyr::pivot_wider(names_from = 'chr_desc', values_from='mode_cn') %>% 
#             as.data.frame()
# 
# class(ascn_tmp)
# View(head(ascn_tmp1))
# rownames(ascn_tmp) <- ascn_tmp$clone_id
# ascn_tmp$clone_id <- NULL
# ascn_tmp1 <- ascn_tmp %>%
#             t() %>%
#             as.data.frame()%>%
#             tidyr::drop_na()%>%  # note: lots of NA in ascn values, need to drop them first.
#             as.data.frame()
# 
# dim(ascn_tmp1)
# colnames(ascn_tmp1)[which(colnames(ascn_tmp1)==cl1)] <- 'cl1'
# colnames(ascn_tmp1)[which(colnames(ascn_tmp1)==cl2)] <- 'cl2'
# # ascn_tmp1 <- ascn_tmp1 %>%
# #   dplyr::filter(cl1<3 & cl2 <3)
# 
# # ascn_tmp1 <- ascn_tmp1 %>%
# #   dplyr::filter((cl1-cl2) != 0)
# # View(ascn_tmp1)
