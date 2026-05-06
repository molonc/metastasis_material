

# R24. The manuscript would benefit from the inclusion of a figure showing the raw single cell data, including sequencing read counts and inferred copy number profiles. Such a figure would be important for assessing the noise levels in the data and for evaluating confidence in the identified clones, which underpin many of the key conclusions.
# Hoa: 
# Inferred copy number states based on read counts plot
# Boxplot for total reads
# Heatmap for all cells 
# Copy number profiles dot plots for representative cells - median profile in each clone
# Hamming distance between each cell versus median copy number profile showing small variance, i.e 2% of total mutation events due to noises, small background mutation

# R25. Many of the central results rely on the definition of clones derived from single cell clustering. 
# Given the inherent noise and uncertainty associated with single cell clustering, 
# it would be important to assess the robustness of these findings. 
# For example, the authors could test whether key conclusions like the identification of metastasizing clones 
# are preserved when using alternative clustering approaches, subsampling cells or genomic regions, 
# or applying other methods to evaluate clustering stability.
# Boxplot of Hamming/Mahattan distances between clones showing large change/driver of mutation events. 

# R14- Are there significant differences in copy number signatures between different clones?
  

get_heatmap_CNA_change <- function(){
  # git dir
  # load cna file
  # load cell clones, library grouping 
  input_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/'
  sa919_metasamples <- data.table::fread(paste0(input_dir, 'materials/dlp_trees/SA919/library_groupings.csv.gz'))
  head(sa919_metasamples)
  
}