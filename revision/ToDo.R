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
## read counts are just a histogram bur reads per bin looks VERY noisy, so you need to include correction, 
## but then it's a pseudocount, not real
## Do we use CG-corrected readcounts values for that?  No
## I would maybe just include total reads per cell
## Maybe I extract genomic regions where copy number changes, defined clones, and then use the boxplot for total reads
## And heatmap for all cells
## cell copy number profiles for a representative cell in each clone.
## + Each cell versus median copy number profile showing a small difference of noises i.e. 2%

suppressPackageStartupMessages({
  library(annotables)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  options(dplyr.summarise.inform = FALSE)
  options(tidyverse.quiet = TRUE)
})



## Get library with more cells in each clone 
get_representative_library <- function(datatag, copynumber_fn, save_dir, results_dir='', cellclone_fn=NULL,library_grouping_fn=NULL,
                                   calcul_distance=F, distance_type='Manhattan'){
  input_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/'
  results_dir <- paste0(input_dir,"materials/dlp_trees/")
  source(paste0(input_dir,"scripts/corrupt_tree/src/cn_change/utils.R"))
  save_dir <- paste0(input_dir,'revision/cnv_profiles_res/')
  datatag <- 'SA919'
  copynumber_fn <- paste0(results_dir, datatag, '/total_merged_filtered_states.csv.gz')
  cellclone_fn <- paste0(results_dir, datatag, '/cell_clones.csv.gz')
  
  
  results_dir <- paste0(dirname(copynumber_fn),'/')
  if(is.null(cellclone_fn)){
    cellclone_fn <- paste0(results_dir,'cell_clones.csv.gz')  
  }
  if(is.null(library_grouping_fn)){
    library_grouping_fn <- paste0(results_dir,'library_groupings.csv.gz')  
  }
  if(!dir.exists(save_dir)){
    dir.create(save_dir, recursive = T)
  }
  
  # copynumber <- read.csv(copynumber_fn, header=T, row.names = 1, check.names = F,stringsAsFactors = FALSE)
  copynumber <- as.data.frame(data.table::fread(copynumber_fn))
  rownames(copynumber) <- copynumber$V1
  copynumber$V1 <- NULL
  # cell_clones contain 2 columns of cell_id, and clone_id
  # ex:           cell_id                      clone_id
  # 1   SA535X4XB02498-A98163A-R09-C11          C
  cell_clones <- data.table::fread(cellclone_fn) %>% as.data.frame()
  dim(cell_clones)
  metasample <- data.table::fread(library_grouping_fn) %>% as.data.frame()
  dim(metasample)
  head(metasample)
  # View(metasample)
  cell_clones <- cell_clones %>%
    dplyr::filter(!clone_id %in% c('None','unassigned'))
  
  
  metasample <- metasample %>%
    dplyr::rename(library_id=grouping) %>%
    dplyr::select(library_id, sample_id)
  print(dim(copynumber))
  
  
  cell_clones$library_id <- get_library_id(cell_clones$cell_id)
  cell_clones$sample_id <- get_sample_id(cell_clones$cell_id)
  cell_clones <- cell_clones %>% left_join(metasample, by=c("library_id","sample_id"))
  
  meta_cells <- cell_clones %>%
    dplyr::group_by(clone_id, library_id) %>%
    dplyr::summarise(nb_cells=n()) %>%
    dplyr::ungroup()
  # View(meta_cells)
  
  lib_df <- tibble::tibble()
  for(cl in unique(meta_cells$clone_id)){
    tmp <- meta_cells %>%
      dplyr::filter(clone_id == cl) 
    max_nbcells <- max(tmp$nb_cells)
    tmp <- tmp %>%
      dplyr::filter(nb_cells==max_nbcells)
    lib_df <- dplyr::bind_rows(lib_df, tmp)
  }
  # lib_df
  data.table::fwrite(lib_df, paste0(save_dir, 'clone_libraries.csv'))  
  
}  

# https://molonc.atlassian.net/wiki/spaces/AP/pages/1879048200/DLP+GC+Normalization
# copy: The y-axis is not “copy” but instead the raw number of reads sampled from each location of the genome
viz_hmmcopy <- function(){
  # Load cell clones profiles 
  # Load reads files
  # Filtering cells 
  # Selecting a random cell
  # Visualizing CN state, and other fields? 
  # Get chr regions with CNA change from median profile
  # Get total reads per CNA change regions, and plot a boxplot for all cells 
  # 
  # 
  # 
  # install.packages("ggExtra")
  library('ggExtra')
  
  
  input_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/'
  save_dir <- paste0(input_dir,'revision/CN_profile/')
  results_dir <- paste0(input_dir,"materials/dlp_trees/")
  data_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/raw_dlp/SA919/'
  source(paste0(input_dir,"scripts/corrupt_tree/src/cn_change/utils.R"))
  
  datatag <- 'SA919'
  copynumber_fn <- paste0(results_dir, datatag, '/total_merged_filtered_states.csv.gz')
  cellclone_fn <- paste0(results_dir, datatag, '/cell_clones.csv.gz')
  
  cell_clones <- data.table::fread(cellclone_fn) %>% as.data.frame()
  dim(cell_clones)
  cell_clones <- cell_clones %>%
    dplyr::filter(!clone_id %in% c('None','unassigned'))
  
  ## Read data from 3 clones, get reads file, combined them into 1 file
  lib_df <- data.table::fread(paste0(save_dir, 'clone_libraries.csv'))  
  
  obs_chrs <- c(1:22, "X")
  total_reads <- tibble::tibble()
  for(obs_lib in lib_df$library_id){
    print(obs_lib)
    fns <- list.files(paste0(data_dir, obs_lib, '/results/'))
    fns <- fns[grepl('*_reads.csv.gz', fns)]
    for(f in fns){
      reads_df <- data.table::fread(paste0(data_dir, obs_lib, '/results/', f))
      # print(dim(reads_df)[1])
      reads_df <- reads_df %>%
        dplyr::filter(cell_id %in% cell_clones$cell_id &
                        chr %in% obs_chrs)
      total_reads <- dplyr::bind_rows(total_reads, reads_df)
      print(dim(total_reads)[1])
    }
  }
  
  save_dir_fg <- paste0(save_dir, 'viz_cell/')
  if(!dir.exists(save_dir_fg)){
    dir.create(save_dir_fg)
  }
  cell_clones_tmp <- cell_clones %>%
    dplyr::filter(cell_id %in% unique(total_reads$cell_id))
  
  ## Sampling cells in each clone
  for(cl in unique(cell_clones_tmp$clone_id)){
    cell_ids <- cell_clones_tmp %>%
      dplyr::filter(clone_id==cl) %>%
      dplyr::pull(cell_id)
    ## to do: get cells that have smallest distance to median profile
    cids <- sample(cell_ids, size = 3)
    print(cids)
    for(obs_cell in cids){
      # obs_cell <- reads_df$cell_id[1]
      df <- total_reads %>%
        dplyr::filter(cell_id == obs_cell)
      print(dim(df))
      data.table::fwrite(df, paste0(save_dir_fg, datatag, '_',obs_cell,'_clone',cl,'.csv.gz'))
      viz_obs_cell(df, obs_cell, paste0(datatag,'_clone',cl), save_dir_fg)
    }
  }
  
  plg <- plot_CN_state_legend()
  saveRDS(plg, paste0(save_dir_fg, 'CN_state_legend.rds'))
}


viz_obs_cell <- function(df, obs_cell, datatag, save_dir_fg){
  # save_dir_fg <- paste0(save_dir, 'viz_cell/')
  if(!dir.exists(save_dir_fg)){
    dir.create(save_dir_fg)
  }
  cnv_cols <- c('0'='#4880B8', '1'='#A7C9DF','2'='#CCCCCC','3'='#F5CE93','4'='#ED9364',
                
                '5'='#D2553E','6'='#A42116','7'='#8B1A43','8'='#CB3576','9'='#D06CAD',
                
                '10'='#C196C4','11'='#D0BAD8')
  obs_chrs <- c(1:22, "X")
  
  df <- df %>%
    dplyr::mutate(reads=case_when(
      reads > 400 ~ 400, 
      TRUE ~ reads
    ),
    copy=case_when(
      copy > 11 ~ 11, 
      TRUE ~ copy
    ))
  # data.table::fwrite(df, paste0(save_dir, 'test_hmmcopy_cloneA.csv'))
  df$state <- as.character(df$state)
  levels(df$state) <- 0:(length(cnv_cols)-1)
  df$state <- factor(df$state, levels = 0:(length(cnv_cols)-1))
  df$chr <- factor(df$chr, levels = obs_chrs) 
  
  ## Reads is noisy but density or boxplot maybe fine
  pr <- ggplot(df, aes(start, reads)) + 
    geom_point(aes(colour = state), size = 0.3) + 
    facet_grid(. ~ chr, scales = "free_x", space = "free_x", switch = "x") + 
    # scale_x_continuous(expand = c(0, 0), breaks = NULL) + 
    scale_color_manual(values = cnv_cols, name = "Copy number ", breaks = names(cnv_cols)) + 
    theme(panel.spacing = unit(0.1, "lines"))
  pr <- add_theme_cnv_plot(pr)
  # pr
  ## Copy is normalized read values
  # max(df$copy, na.rm = T)
  # sum(is.na(df$copy))
  # max(as.numeric(df$state), na.rm = T)
  # min(as.numeric(df$state), na.rm = T)
  pc <- ggplot(df, aes(start, copy)) + 
    geom_point(aes(colour = state), size = 0.3) + 
    facet_grid(. ~ chr, scales = "free_x", space = "free_x", switch = "x") + 
    # scale_x_continuous(expand = c(0, 0), breaks = NULL) + 
    scale_color_manual(values = cnv_cols, name = "Copy number ", breaks = names(cnv_cols)) + 
    theme(panel.spacing = unit(0.1, "lines"))
  pc <- add_theme_cnv_plot(pc)
  # pc
  # p <- ggMarginal(p, type = "density", margins = "y", groupColour = TRUE, groupFill = TRUE)
  saveRDS(pr, paste0(save_dir_fg, obs_cell, '_',datatag,'_reads_cnv_plt.rds'))
  saveRDS(pc, paste0(save_dir_fg, obs_cell, '_',datatag,'_copy_cnv_plt.rds'))
  
  
  p_total <- cowplot::plot_grid(pr, pc, ncol=1)
  # For quick demo
  png(paste0(save_dir_fg, obs_cell, '_',datatag,'_cnv_plt.png'), 
      height = 2*400, width=2*900,res = 2*72)
  print(p_total)
  dev.off()
}
viz_CNA_changes_distribution <- function(){
  ## CNA changes only
  cnv_mat <- get_CNA_change()
  dim(cnv_mat)
  # View(cnv_mat[100:150,])
  input_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/'
  save_dir <- paste0(input_dir,'revision/CN_profile/')
  results_dir <- paste0(input_dir,"materials/dlp_trees/")
  data_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/raw_dlp/SA919/'
  source(paste0(input_dir,"scripts/corrupt_tree/src/cn_change/utils.R"))
  
  datatag <- 'SA919'
  copynumber_fn <- paste0(results_dir, datatag, '/total_merged_filtered_states.csv.gz')
  cellclone_fn <- paste0(results_dir, datatag, '/cell_clones.csv.gz')
  
  cell_clones <- data.table::fread(cellclone_fn) %>% as.data.frame()
  dim(cell_clones)
  cell_clones <- cell_clones %>%
    dplyr::filter(!clone_id %in% c('None','unassigned'))
  
  total_reads <- total_reads %>%
    dplyr::mutate(chr_desc=paste0(chr,'_',start,'_', end)) %>%
    dplyr::filter(chr_desc %in% cnv_mat$chr_desc)
  print(dim(total_reads)[1])
  total_reads_backup <- total_reads

  
  counts_df <- total_reads %>%
    dplyr::group_by(chr, cell_id) %>%
    dplyr::summarise(mean_reads=round(sum(reads)/length(chr),2), 
                     mean_copy=round(sum(copy)/length(chr),2)) %>%
    dplyr::left_join(cell_clones, by='cell_id')
  
  # df <- total_reads %>%
  #   dplyr::filter(cell_id==total_reads$cell_id[1]) %>%
  #   dplyr::group_by(chr) %>%
  #   dplyr::summarise(total_reads=round(sum(reads)/length(chr),2), 
  #                    total_copy=round(sum(copy)/length(chr),2)) 
  
  data.table::fwrite(counts_df, paste0(save_dir, datatag, '_cn_changes_mean_reads_copy.csv'))
  counts_df <- data.table::fread(paste0(save_dir, datatag, '_cn_changes_mean_reads_copy.csv'))
  head(counts_df)
  p <- ggplot(counts_df, aes(clone_id, mean_reads, colour = clone_id)) + 
    geom_boxplot() + 
    facet_grid(. ~ chr, scales = "free") #, space = "free_x", switch = "x"
  p
  p <- add_theme(p)
  p
  counts_df1 <- counts_df %>%
    dplyr::filter(mean_copy<5)
  typeof(counts_df1$chr)
  counts_df1$chr <- paste0('chr ',counts_df1$chr)
  counts_df1$chr <- factor(counts_df1$chr, levels = c('chr 5','chr 7','chr 10'))
  p1 <- ggplot(counts_df1, aes(clone_id, mean_copy, colour = clone_id)) + 
    geom_boxplot(outlier.shape = NA) + 
    facet_grid(. ~ chr, scales = "free") + #, space = "free_x", switch = "x"
    theme_bw()
  p1 <- add_theme(p1)
  p1
  saveRDS(p1, paste0(save_dir, datatag, '_cn_change_mean_copy.rds'))
}
add_theme_cnv_plot <- function(cnv_plot){
  lg_pos <- "none"
  my_font <- "Helvetica"
  cnv_plot <- cnv_plot + theme(strip.background = element_rect(fill = 'white', colour = 'white'),
                               strip.text = element_text(color="black",size=10, hjust = 0.5, family=my_font, angle = 90),
                               text = element_text(color="black",size = 11, hjust = 0.5, family=my_font),
                               axis.text.x = element_blank(),
                               axis.ticks.x = element_blank(),
                               axis.text.y = element_text(color="black",size=11, hjust = 0.5, family=my_font),
                               axis.title.y = element_text(color="black",size=11, hjust = 0.5, family=my_font),
                               axis.title.x = element_text(color="black",size=11, hjust = 0.5, family=my_font),
                               axis.line = element_line(colour = "black"),
                               strip.placement = "outside",
                               legend.position = lg_pos,
                               legend.text=element_text(color="black",size=9, hjust = 0.5, family=my_font),
                               legend.title=element_text(color="black",size=9, hjust = 0.5, family=my_font),
                               legend.key.size=unit(0.3,"cm"),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               # panel.background = element_rect(fill = "#F8F8F8", colour = NA),
                               panel.spacing = unit(c(0.1), 'cm'),
                               legend.margin=margin(0,0,0,0),
                               legend.box.margin=margin(-2,-2,-2,-2)) # MA: was 0.2
  
  cnv_plot <- cnv_plot + guides(fill = guide_legend(nrow = 1, override.aes = list(size=0.1))) + #, override.aes = list(size=1.1)
    scale_x_continuous(expand = c(0,0))
  return(cnv_plot)
  
}
add_theme <- function(cnv_plot){
  lg_pos <- "none"
  my_font <- "Helvetica"
  cnv_plot <- cnv_plot + theme(strip.background = element_rect(fill = 'white', colour = 'white'),
                               strip.text = element_text(color="black",size=12, hjust = 0.5, family=my_font),
                               text = element_text(color="black",size = 11, hjust = 0.5, family=my_font),
                               axis.text.x = element_text(color="black",size=13, hjust = 0.5, family=my_font),
                               axis.ticks.x = element_blank(),
                               axis.text.y = element_text(color="black",size=11, hjust = 0.5, family=my_font),
                               axis.title.y = element_text(color="black",size=11, hjust = 0.5, family=my_font),
                               axis.title.x = element_text(color="black",size=11, hjust = 0.5, family=my_font),
                               axis.line = element_line(colour = "black"),
                               strip.placement = "outside",
                               legend.position = lg_pos,
                               legend.text=element_text(color="black",size=9, hjust = 0.5, family=my_font),
                               legend.title=element_text(color="black",size=9, hjust = 0.5, family=my_font),
                               legend.key.size=unit(0.3,"cm"),
                               # panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               # panel.background = element_rect(fill = "#F8F8F8", colour = NA),
                               panel.spacing = unit(c(0.1), 'cm'),
                               legend.margin=margin(0,0,0,0),
                               legend.box.margin=margin(-2,-2,-2,-2)) # MA: was 0.2
  
  # cnv_plot <- cnv_plot + guides(fill = guide_legend(nrow = 1, override.aes = list(size=0.1))) + #, override.aes = list(size=1.1)
  #   scale_x_continuous(expand = c(0,0))
  return(cnv_plot)
  
}
plot_CN_state_legend <- function(){
  my_font <- "Helvetica"
  cnv_cols <- c('0'='#4880B8', '1'='#A7C9DF','2'='#CCCCCC','3'='#F5CE93','4'='#ED9364',
                
                '5'='#D2553E','6'='#A42116','7'='#8B1A43','8'='#CB3576','9'='#D06CAD',
                
                '10'='#C196C4','11'='#D0BAD8')
  cls_df <- data.frame(cnv_val=names(cnv_cols))
  cls_df$Freq <- 1
  # cls_df <- cls_df[order(cls_df$clone_id),]
  cls_df$cnv_val <- factor(cls_df$cnv_val, levels = sort(as.numeric(cls_df$cnv_val))) #, ordered = T
  # cls_df$clone_id <- factor(cls_df$clone_id, levels = sort(cls_df$clone_id)) #, ordered = T
  p <- ggplot(cls_df, aes(fill=sort(cnv_val), y=Freq, x=cnv_val)) + 
    geom_bar(position="fill", stat="identity",width=0.1) + 
    scale_fill_manual(values = cnv_cols)+
    theme(legend.text=element_text(color="black",size=10, hjust = 0.5, family=my_font),
          legend.title=element_text(color="black",size=10, hjust = 0.5, family=my_font))
  # p
  lg <- cowplot::get_legend(p + guides(fill = guide_legend(title = "Copy number state ",nrow=1, 
                                                           title.position = "left", 
                                                           override.aes = list(shape = 0, size=0.001)))) #nrow=1
  plg <- cowplot::ggdraw() + cowplot::draw_plot(lg)
  # plg
  return(plg)
}

## Testing at reads, and copy level
get_CNA_change <- function(){
  
  input_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/'
  save_dir <- paste0(input_dir,'revision/CN_profile/')
  datatag <- 'SA919'
  median_cnv <- data.table::fread(paste0(save_dir, datatag, '_median_cnv.csv'))
  head(median_cnv)
  obs_chrs <- c('5','7','10')
  obs_clones <- c('A','B','C')
  # de_desc <- paste0(obs_clone1[1],'_vs_',obs_clone2[1])
  # cnv_mat <- median_cnv[,c(obs_clones, "chr_desc")] %>% as.data.frame()
  cnv_mat <- median_cnv
  chrs <- sapply(strsplit(cnv_mat$chr_desc , "_"), function(x) {
    return(x[1])
  })
  cnv_mat$chr <- as.character(chrs)
  unique(cnv_mat$chr)
  cnv_mat <- cnv_mat %>%
    dplyr::filter(chr %in% obs_chrs)
  tmp <- cnv_mat %>%
    tibble::column_to_rownames('chr_desc') %>%
    dplyr::select(-chr)
  
  var_genes <- matrixStats::rowVars(as.matrix(tmp))
  var_genes <- var_genes[var_genes>0]
  var_genes[1:3]
  chr_cna_change <- names(var_genes)
  # var_genes <- apply(cnv_mat[,!colnames(cnv_mat) %in% c("chr_desc")], 1, var)
  # cnv_mat <- cnv_mat[var_genes > 0,]
  
  print(dim(cnv_mat))
  print(length(chr_cna_change))
  cnv_mat <- cnv_mat %>%
    dplyr::filter(chr_desc %in% chr_cna_change)
  # View(cnv_mat)
  return(cnv_mat)
}

## See function at: 
# file_dir <- "/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/scripts/corrupt_tree/src/cn_change/extract_reads_rawdata_info.R"
extract_cells_features_manuscript_v2()


plot_R24 <- function(){
  # ## 2 x 2 reads, copy plots for clone A, B, C
  # [1] "SA919X4XB40509-A98217A-R22-C40" "SA919X4XB40509-A98217A-R25-C51"
  # [3] "SA919X4XB40509-A98217A-R20-C17"
  # [1] "SA919X7XB05691-A96204B-R37-C41" "SA919X7XB05691-A96204B-R29-C07"
  # [3] "SA919X7XB05691-A96204B-R42-C54"
  # [1] "SA919X7XB05588-A98299A-R12-C46" "SA919X7XB05588-A98299A-R20-C31"
  # [3] "SA919X7XB05588-A98299A-R11-C13"
  input_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/'
  save_dir <- paste0(input_dir,'revision/CN_profile/')
  save_dir_fg <- paste0(save_dir, 'viz_cell/')
  obs_cell <- "SA919X4XB40509-A98217A-R22-C40" # clone A
  # obs_cell <- "SA919X7XB05691-A96204B-R48-C58" # clone B
  # obs_cell <- "SA919X7XB05588-A98299A-R05-C21" # clone C
  cl <- 'A'
  pr <- readRDS(paste0(save_dir_fg, obs_cell, '_',datatag,'_clone',cl,'_reads_cnv_plt.rds'))
  pc <- readRDS(paste0(save_dir_fg, obs_cell, '_',datatag,'_clone',cl,'_copy_cnv_plt.rds'))
  plg <- readRDS(paste0(save_dir_fg, 'CN_state_legend.rds'))
  p_total <- cowplot::plot_grid(NULL, pr, pc, plg,NULL, ncol=1, rel_heights = c(0.05, 1,1, 0.1,0.03))
  p_total
  
  ## total mapped reads for 2 series
  p_total_reads <- readRDS(paste0(paste0(save_dir, 'cell_metrics_eval/','total_metrics_SA919_SA535_mapped_reads_plt.rds')))
  
  ## Copy values in chr 5, 7, 10
  ## Output of function: viz_CNA_changes_distribution() above
  p_cnchange_copy <- readRDS(paste0(save_dir, datatag, '_cn_change_mean_copy.rds'))
  
  
  
  ## Hamming distance between clones for 2 series 
  ## Heatmap of 2 series 
  distance_type='Hamming'
  datatag <- 'SA919'
  p_SA919 <- readRDS(paste0(save_dir, 'clone_distance/', datatag, distance_type,'_distance_clones','.rds'))
  p_cnv <- cowplot::plot_grid(p_total_reads, p_cnchange_copy, p_SA919, ncol=3)
  p_cnv
  
  datatag <- 'SA535'
  p_SA535 <- readRDS(paste0(save_dir, 'clone_distance/', datatag, distance_type,'_distance_clones','.rds'))
  
  p_all <- cowplot::plot_grid(p_total, p_cnv, p_SA535, ncol=1, rel_heights = c(1,0.6, 0.8),
                              labels = c('Raw read counts & normalized copy profiles',' ',' ')) +
    theme(plot.background = element_rect(fill = "white", colour = "white"))
  ggplot2::ggsave( 
    filename = paste0(save_dir,"R24.svg"), 
    plot = p_all,
    height = 11, 
    width = 9, 
  )
   
  ggplot2::ggsave( 
    filename = paste0(save_dir,"R24.png"), 
    plot = p_all,
    height = 11, 
    width = 8, 
  ) 
}
