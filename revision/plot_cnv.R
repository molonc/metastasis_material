suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(ggplot2)
  library('ggExtra')
})


cnv_cols <- c('0'='#4880B8', '1'='#A7C9DF','2'='#CCCCCC','3'='#F5CE93','4'='#ED9364',
              
              '5'='#D2553E','6'='#A42116','7'='#8B1A43','8'='#CB3576','9'='#D06CAD',
              
              '10'='#C196C4','11'='#D0BAD8')
input_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/'
save_dir <- paste0(input_dir,'revision/cnv_profiles_res/')

df <- data.table::fread(paste0(save_dir, 'test_hmmcopy_cloneA.csv'))
obs_chrs <- c(1:22, "X")
df$state <- as.character(df$state)
levels(df$state) <- 0:(length(cnv_cols)-1)
df$chr <- factor(df$chr, levels = obs_chrs) 
p <- ggplot(df, aes(start, copy)) + 
  geom_point(aes(colour = state), size = 0.3) + 
  facet_grid(. ~ chr, scales = "free_x", space = "free_x", switch = "x") + 
  scale_x_continuous(expand = c(0, 0), breaks = NULL) + 
  scale_color_manual(values = cnv_cols, name = "Copy number ", breaks = names(cnv_cols)) + 
  theme(panel.spacing = unit(0.1, "lines"))
# p

p <- ggMarginal(p, type = "density", margins = "y", groupColour = TRUE, groupFill = TRUE)