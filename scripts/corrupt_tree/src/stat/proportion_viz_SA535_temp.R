
suppressPackageStartupMessages({
  # library(dplyr)
  # library(stringr)
  # library(RColorBrewer)
  # library(ggplot2)
  # library(viridis)
  # library("argparse")
  library("dplyr")
  library("ggplot2")
})


input_dir <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/dlp_trees/'
# input_dir <- '/Users/htran/Documents/storage_tmp/metastasis_trees/SA535/viz_graph_535/'
cell_clones <- data.table::fread(paste0(input_dir, 'SA535/cell_clones.csv.gz'))
grouping <- data.table::fread(paste0(input_dir, 'SA535/library_groupings.csv.gz'))
colors_df <- data.table::fread(paste0(input_dir, 'colorCode_clones/color_code_SA535.csv'))
grouping$pdxid <- gsub('X0011_236','M',grouping$pdxid)
sids <- sapply(strsplit(cell_clones$cell_id,'-'), function(x){
  return(x[1])
})
lids <- sapply(strsplit(cell_clones$cell_id,'-'), function(x){
  return(x[2])
})
cell_clones$library_id <- as.character(lids)
cell_clones$sample_id <- as.character(sids)
unique(cell_clones$sample_id)
cell_clones <- cell_clones %>%
  left_join(grouping, by=c('library_id'='grouping', 'sample_id'))
dim(cell_clones)

unique(grouping$mainsite)
grouping <- grouping %>%
  dplyr::filter(mainsite=="Primary")
dim(grouping)
input_dir <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/dlp_trees/'
t <- data.table::fread(paste0(input_dir, 'SA535/jira_tickets_library_groupings.csv.gz'))
dim(t)
typeof(t$library_id)
obs_libs = as.character(grouping$grouping)
t <- t %>%
  dplyr::filter(library_id %in% obs_libs)
paste(t$grouping, collapse = "' '")
paste(t$jira_ticket, collapse = "' '")

extracted_clones_df <- cell_clones %>%
   dplyr::filter(clone_id!='None' & 
                   mainsite=="Primary" & pdxid=="M2") 
dim(extracted_clones_df)
summary(as.factor(extracted_clones_df$clone_id))
summary(as.factor(extracted_clones_df$library_id))
table(extracted_clones_df$library_id, extracted_clones_df$sample_id)
unique(extracted_clones_df$library_id)

## for combined libraries
results_dir <- '/home/htran/storage/datasets/metastasis_results/SA535_new_encoding/SA535_wholedata_v2/'
copynumber <- data.table::fread(paste0(results_dir, 'total_merged_filtered_states.csv'))%>% as.data.frame()
rownames(copynumber) <- copynumber$V1
copynumber$V1 <- NULL
print(dim(copynumber))
cells_use <- extracted_clones_df$cell_id
copynumber <- copynumber[, cells_use]
data.table::fwrite(copynumber, paste0(results_dir, 'temp/total_merged_filtered_states_M2_Primary.csv'))
# stat_df <- cell_clones %>%
#   dplyr::filter(clone_id!='None') %>%
#   dplyr::group_by(pdxid, mainsite, clone_id) %>%
#   dplyr::summarise(nb_cells=n())

stat_df <- cell_clones %>%
  dplyr::filter(clone_id!='None') %>%
  dplyr::group_by(pdxid, mainsite, origin, clone_id) %>%
  dplyr::summarise(nb_cells=n())

unique(cell_clones$origin)
# stat1 <- stat_df %>%
#   dplyr::group_by(pdxid, mainsite) %>%
#   dplyr::summarise(nb_total_cells=sum(nb_cells))

stat1 <- stat_df %>%
  dplyr::group_by(pdxid, mainsite, origin) %>%
  dplyr::summarise(nb_total_cells=sum(nb_cells))

stat_df <- stat_df %>%
  left_join(stat1, by=c('pdxid','mainsite','origin'))

stat_df$prop <- round(100*stat_df$nb_cells/stat_df$nb_total_cells,1)

stat_df$mainsite <- factor(stat_df$mainsite, levels = c('Primary','Metastasis'))


my_font <- 'Helvetica'

## To Do: adding custom colors code
p <- ggplot(data=stat_df, aes(x=origin, y=prop, fill=clone_id)) +
  geom_bar(stat="identity", width=0.5)+
  # facet_grid(. ~ fov)+ 
  facet_wrap(mainsite ~ pdxid, ncol = 4, strip.position = "top") + 
  theme_bw() +
  theme(
    axis.text.x = element_text(size=14, angle = 90),
    legend.position = 'none',
    strip.background = element_rect(fill = 'white', colour = 'white'),
    text = element_text(color="black",size = 12, hjust = 0.5, family=my_font),
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    # axis.title.x = element_blank(),
    strip.text.x = element_text(color="black",size=14, hjust = 0, family=my_font),
    axis.text.y = element_text(color="black",size=9, hjust = 0.5, family=my_font),
    axis.title.y = element_text(color="black",size=11, hjust = 0.5, family=my_font),
    axis.line = element_line(colour = "black"),
    strip.placement = "outside",
    # legend.position = lg_pos,
    # legend.text=element_text(color="black",size=8, hjust = 0.5, family=my_font),
    # legend.title=element_text(color="black",size=8, hjust = 0.5, family=my_font),
    # legend.key.size=unit(0.3,"cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.spacing = unit(c(0.1), 'cm'),
    # legend.margin=margin(0,0,0,0),
    # legend.box.margin=margin(-2,-2,-2,-2)
  ) + 
  labs(title="SA535 clonal prevalence across different mouses", x='Metastasis sites', y='Clonal proportion')

datatag <- 'SA535'
out_fn <- paste0(input_dir,datatag,'clonal_prevalence.png')
png(out_fn, height = 2*500, width=2*600, res = 2*72)
print(p)
dev.off()





stat_df <- cell_clones %>%
  dplyr::filter(clone_id!='None') %>%
  dplyr::group_by(pdxid, mainsite, clone_id) %>% #origin, 
  dplyr::summarise(nb_cells=n())


stat1 <- stat_df %>%
  dplyr::group_by(pdxid, mainsite) %>%
  dplyr::summarise(nb_total_cells=sum(nb_cells))

stat_df <- stat_df %>%
  left_join(stat1, by=c('pdxid','mainsite'))

stat_df$prop <- round(100*stat_df$nb_cells/stat_df$nb_total_cells,1)

for(m in c('M1','M2','M3','M4')){
  m <- 'M1'
  tmp <- stat_df %>%
    dplyr::filter(pdxid==m & nb_cells>=5) 
  tmp1 <- tmp %>%
    dplyr::group_by(clone_id) %>%
    dplyr::summarise(nb_sites=n()) %>%
    dplyr::filter(nb_sites==2)
  print('-----------')
  print(m)
  print(tmp1$clone_id)
  # print(table(tmp$clone_id, tmp$mainsite))
  df <- as.matrix(table(tmp$clone_id, tmp$mainsite))
  print(df)
  
}


df <- as.data.frame(table(tmp$clone_id, tmp$mainsite))


View(tmp)


