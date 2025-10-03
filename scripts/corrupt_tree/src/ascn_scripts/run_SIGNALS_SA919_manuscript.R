suppressPackageStartupMessages({
  require("data.table")
  require("dplyr")
  # require("scales")
})
# install.packages("devtools")
# devtools::install_github("shahcompbio/signals")
# install.packages(c('uwot', 'dbscan', 'ComplexHeatmap', 'ape', 'ggtree'))
# BiocManager::install("ComplexHeatmap")
# BiocManager::install("ggtree")

# install.packages('/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/results_CNA/signals_source/signals-master/', repos = NULL, type="source")



# script.basename <- '/home/htran/Projects/hakwoo_project/rscript/schnapp_haplo'
script.basename <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/results_CNA/ascn_SA919/scripts/'
source(paste0(script.basename, "/schnapp_utils.R"))
source(paste0(script.basename, "/utils_heatmap.R"))



# Load SA919 data 
# results_pseudobk_dir <- '/home/htran/storage/datasets/metastasis_results/snv_analysis_SA919/ascn_SA919_Hoa/'
results_pseudobk_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/results_CNA/ascn_SA919/'
results_dir <- results_pseudobk_dir
## allele object contains input data for all libraries, and columns: 
# [1] 37199267       19
# [1] "chr"             "start"           "end"            
# [4] "cell_id"         "state"           "copy"           
# [7] "alleleA"         "alleleB"         "totalcounts"    
# [10] "BAF"             "state_min"       "Min"            
# [13] "Maj"             "state_AS_phased" "state_AS"       
# [16] "LOH"             "phase"           "state_phase"    
# [19] "state_BAF"

# The main output is a dataframe with that is similar to the CNbins input file 
# with the following additional columns: 
# * A A allele copy number 
# * B B allele copy number 
# * state_AS_phased A|B 
# * state_min Minor allele copy number 
# * LOH =LOH if bin is LOH, NO otherwise 
# * state_phase Discretized haplotype specific states (see below) 
# * phase Whether the A allele or B allele is dominant 
# * alleleA Counts for the A allele 
# * alleleB Counts for the B allele 
# * totalcounts Total number of counts 
# * BAF B-allele frequency (alleleB / totalcounts)
# 
# state_phase has the following states:
# * A-Hom A allele is homozygous, ie LOH of B-Allele 
# * B-Hom B allele is homozygous, ie LOH of A-Allele 
# * A-Gained A allele is gained but B not lost (A > B) 
# * B-Gained B allele is gained but A not lost (B > A) 
# * Balanced A == B

ascn_fn = paste0(results_pseudobk_dir,'ascn_inference.rds')

# save_dir <- paste0(results_dir,'ascn_plot/')
save_dir <- paste0(results_dir,'ascn_viz/')
if(!dir.exists(save_dir)){
  dir.create(save_dir)
}
# newick_fn <- paste0(results_dir,'tree.newick')
# copynumber_fn <- paste0(results_dir,'total_merged_filtered_states.csv')
datatag <- 'SA919'
cellclone_fn <- paste0(results_dir,'cell_clones.csv.gz')
# cell_clones <- read.csv(paste0(results_dir,'ascn_plot/ascn_cell_clones.csv'), check.names = F,stringsAsFactors = FALSE)
# summary(as.factor(cell_clones$clone_id))
library_grouping_fn <- paste0(results_dir,'library_groupings.csv.gz')
copynumber_fn <- paste0(results_dir,'total_merged_filtered_states.csv.gz')

median_ascn <- get_median_LOH_genotype_v2(ascn_fn, 
                               library_grouping_fn, cellclone_fn, 
                               datatag, save_dir)


## From comparisons, and detect the LOH changes between pair of clones. 

detect_LOH_change_v2(copynumber_fn, median_ascn, save_dir, datatag)





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
