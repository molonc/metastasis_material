##
## Get metadata file from github for SA919
## Selecting 1 library
## Downloading hmncopy file from azure
## Get different fields 
## Plot the corrected copy values, and copy number states
## Zoom in the chr 5, 7, 10
## Plot the corrected copy values, and copy number states
##
## Daniel:
## heatmap plot: because they want to: assessing the noise levels in the data
## heatmap is Alhena heatmaps from CNV calls
## you can also provide a distribution of the read counts: boxplot
## read counts are just a histogram bur reads per bin looks VERY noisy, so you need to include correction, but then it's a pseudocount, not real
## Do we use CG-corrected readcounts values for that?  No
## I would maybe just include total reads per cell
## Maybe I extract genomic regions where copy number changes, defined clones, and then use the boxplot for total reads
## And heatmap for all cells
## cell copy number profiles for a representative cell in each clone.
## + Each cell versus median copy number profile showing a small difference of noises i.e. 2%

