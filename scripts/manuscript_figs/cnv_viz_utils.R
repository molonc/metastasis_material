suppressPackageStartupMessages({
  library(tidyverse)
  library(annotables)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  require(scales)
  require(ggrepel)
  require(stringr)
})
library(extrafont)
# font_import(prompt=F, paths ='/usr/share/fonts/truetype/myfonts/') # import Helvetica font
# fonts()
get_median_cnv <- function(df_cnv){
  df_cnv <- data.table::fread(df_cnv_fn, check.names = F, stringsAsFactors = F) %>% as.data.frame()
  
  print("Get median genotype")
  median_cnv <- df_cnv %>%
    dplyr::group_by(clone_id, chr_desc) %>%
    dplyr::summarise(median_cn=median(copy_number)) #mode_cn=calc_mode(copy_number)
  
  median_cnv <- median_cnv %>%
    dplyr::select(chr_desc, median_cn, clone_id) %>%
    group_by(clone_id) 
  
  
  median_cnv <- median_cnv %>%
    pivot_wider(names_from = clone_id, values_from = median_cn)
  # write.csv(median_cnv, paste0(results_dir,'median_cnv.csv'), quote = F, row.names = F)
  median_cnv <- as.data.frame(median_cnv)
  return(median_cnv)
}
get_median_genotype_v3 <- function(copynumber_fn, datatag, save_dir,
                                   cellclone_fn=NULL, library_grouping_fn=NULL){
  # if(is.null(cellclone_fn)){
  #   cellclone_fn <- paste0(save_dir,'cell_clones.csv.gz')  
  # }
  # if(is.null(library_grouping_fn)){
  #   library_grouping_fn <- paste0(save_dir,'library_groupings.csv.gz')  
  # }
  if(!dir.exists(save_dir)){
    dir.create(save_dir)
  }
  
  # copynumber <- read.csv(copynumber_fn, header=T, row.names = 1, check.names = F,stringsAsFactors = FALSE)
  copynumber <- as.data.frame(data.table::fread(copynumber_fn))
  rownames(copynumber) <- copynumber$V1
  copynumber$V1 <- NULL
  dim(copynumber)
  # cell_clones contain 2 columns of cell_id, and clone_id
  # ex:           cell_id                      clone_id
  # 1   SA535X4XB02498-A98163A-R09-C11          C
  cell_clones <- data.table::fread(cellclone_fn) %>% as.data.frame()
  dim(cell_clones)
  # metasample <- data.table::fread(library_grouping_fn) %>% as.data.frame()
  # dim(metasample)
  # head(metasample)
  cell_clones <- cell_clones %>%
    dplyr::filter(!clone_id %in% c('None','unassigned'))
  
  # if(datatag=='SA604'){
  #   cell_clones$clone_id <- ifelse(cell_clones$clone_id=='R1','K',
  #                                  ifelse(cell_clones$clone_id=='R2','L',cell_clones$clone_id))  
  # }
  
  # colnames(metasample)
  # metasample <- metasample %>%
  #   dplyr::rename(library_id=libraryID, sample_id=sample) %>%
  #   dplyr::select(library_id, sample_id)
  # print(dim(copynumber))
  
  
  cell_clones$library_id <- get_library_id(cell_clones$cell_id)
  cell_clones$sample_id <- get_sample_id(cell_clones$cell_id)
  # cell_clones <- cell_clones %>% left_join(metasample, by=c("library_id","sample_id"))
  # dim(cell_clones)
  
  copynumber <- copynumber[,colnames(copynumber)[colnames(copynumber) %in% cell_clones$cell_id]]
  copynumber$chr_desc <- rownames(copynumber)
  # dim(copynumber)
  cnv <- copynumber %>%
    pivot_longer(!chr_desc, names_to = "cell_id", values_to = "copy_number")
  
  print(dim(cnv))
  cnv <- cnv %>% left_join(cell_clones, by=c("cell_id"))
  clones_ls <- unique(cell_clones$clone_id)
  clones_ls <- clones_ls[clones_ls %in% unique(cnv$clone_id)]
  # get median genotype 
  print("Get median genotype")
  stat_cnv <- cnv %>%
    dplyr::group_by(clone_id, chr_desc) %>%
    dplyr::summarise(cnv=median(copy_number)) %>%
    dplyr::ungroup()
  dim(stat_cnv)
  # median_cnv <- stat_cnv %>%
  #   dplyr::select(chr_desc, median_cn, clone_id, clone_label) %>%
  #   group_by(clone_label) 
  # head(stat_cnv)
  stat_cnv <- stat_cnv %>%
    dplyr::rename(clone=clone_id)
  
  stat_cnv$clone <- factor(stat_cnv$clone, levels = clones_ls)
  # median_cnv_pivot <- stat_cnv %>%
  #   pivot_wider(names_from = clone_id, values_from = median_cn) %>%
  #   as.data.frame()
 
  regions <- sapply(stat_cnv$chr_desc, strsplit, "_")
  stat_cnv$chr <- unname(sapply(regions, '[[', 1))
  stat_cnv$start <- as.numeric(unname(sapply(regions, '[[', 2)))
  stat_cnv$end <- as.numeric(unname(sapply(regions, '[[', 3)))
  stat_cnv <- stat_cnv %>%
    dplyr::select(-chr_desc)
  data.table::fwrite(stat_cnv, paste0(save_dir, datatag,'_median_cnv.csv'))
  # data.table::fwrite(stat_cnv, paste0(save_dir,datatag, '_stat_median_cnv.csv'))
  # dim(stat_cnv)
  return(stat_cnv)  
}  

# df_cnv_fn <- paste0(save_dir, datatag,'_median_cnv.csv')
# df_cnv has column names: ensembl_gene_id
plot_CNV_profile <- function(df_cnv, clones=NULL, plttitle='',meta_genes=NULL, datatag=''){
  # df_cnv <- data.table::fread(df_cnv_fn) %>% as.data.frame()
  # df_cnv <- df_cnv %>%
  #   dplyr::rename(ensembl_gene_id=V1)
  
  # View(head(df_cnv))
  print(dim(df_cnv))
  # Select only genes that contain variance across all clones --> make plot clearer
  # if(is.null(meta_genes)){
  #   var_genes <- df_cnv %>%
  #     tibble::column_to_rownames('ensembl_gene_id') %>%
  #     dplyr::mutate(var_gene = rowVars(as.matrix(.), na.rm=T)) %>%
  #     dplyr::pull(var_gene)
  #   
  #   df_cnv <- df_cnv[var_genes>0,]
  # }else{
  #   df_cnv <- df_cnv %>%
  #     dplyr::filter(ensembl_gene_id %in% meta_genes$ensembl_gene_id)
  # }
  
  
  
  # df_cnv <- df_cnv %>% 
  #   dplyr::filter(clone %in% clones)  #, !chr %in% c("X","Y")
  print(dim(df_cnv))
  # sum(is.na(df_cnv$cnv))
  df_cnv <- df_cnv[!is.na(df_cnv$cnv),]
  
  print(dim(df_cnv))
  chr_levels <- c(as.character(1:22), "X")
  df_cnv <- df_cnv %>% dplyr::filter(chr %in% chr_levels)
 
  df_cnv <- df_cnv  %>%
    dplyr::select(cnv, clone, chr, start, end) %>%
    group_by(chr) %>%
    # dplyr::mutate(start_order = rank((start+end)/2)) %>%
    dplyr::mutate(start_order = rank(start)) %>%
    ungroup() 
  df_cnv$chr <- factor(df_cnv$chr, levels = chr_levels)
  # df_cnv <- tidyr::drop_na(df_cnv)
  
  # cnv_colors <- c('#4880B8', '#A7C9DF','#CCCCCC','#F5CE93','#ED9364','#D2553E','#A42116','#8B1A43','#CB3576','#D06CAD','#C196C4','#D0BAD8')
  
  cnv_cols <- c('0'='#4880B8', '1'='#A7C9DF','2'='#CCCCCC','3'='#F5CE93','4'='#ED9364',
                '5'='#D2553E','6'='#A42116','7'='#8B1A43','8'='#CB3576','9'='#D06CAD',
                '10'='#C196C4','11'='#D0BAD8')
  
  df_cnv$cnv <- as.character(round(df_cnv$cnv))
  levels(df_cnv$cnv) <- 0:(length(cnv_cols)-1)
  
  # # Just in case of SA609
  # df_cnv$clone <- ifelse(df_cnv$clone=='R','A',df_cnv$clone)
  # if(clones[1]=='R'){
  #   clones[1] <- 'A'
  # }
  if(is.null(clones)){
    clones <- unique(df_cnv$clone)  
  }
  df_cnv$clone <- factor(df_cnv$clone, levels = clones)
  print(summary(as.factor(df_cnv$clone)))
  
  cnv_plot <- df_cnv %>% 
    ggplot(aes(x = start_order, y = clone, fill = cnv)) + #, fill = cnv
    # geom_tile(aes(fill = cnv)) + #
    geom_raster() +
    # facet_wrap(~ chr, scales = "free_x", switch = "x", nrow = 1) + #"free_x" nrow = 1
    facet_grid(~ chr, scales = "free_x", switch = "x", space='free') +
    # theme(legend.position = "bottom", axis.text.x = element_blank()) +
    scale_y_discrete(expand = c(0, 0)) +
    # scale_fill_manual(values=cnv_colors, name = "Copy number", guide = 'legend',labels = 0:(length(cnv_colors)-1),drop=FALSE)  +
    #          theme(legend.position = "bottom") +
    scale_fill_manual(values = cnv_cols, name = "Copy number ",breaks = names(cnv_cols)) +  #, 
    labs(x = "Chromosome", y = "Clone", title=plttitle)
  # cnv_plot
  my_font <- "Helvetica"
  
  lg_pos <- "bottom"
  # thesis_theme <- ggplot2::theme(
  #   text = element_text(color="black",size = 8, hjust = 0.5, family=my_font),
  #   axis.title.x = element_text(color="black",size=8, hjust = 0.5, family=my_font),
  #   axis.title.y = element_text(color="black",size=8, hjust = 0.5, family=my_font),
  #   # axis.text.x = element_text(color="black",size=7, hjust = 0.5, family=my_font, angle = 90),
  #   axis.text.x = element_blank(),
  #   axis.ticks.x = element_blank(),
  #   strip.placement = "outside",
  #   axis.line = element_line(colour = "black"),
  #   axis.text.y = element_text(color="black",size=7, hjust = 0.5, family=my_font),
  #   plot.title = element_text(color="black",size=10, face="bold", hjust=0, family=my_font),
  #   legend.title=element_text(color="black",size=7, hjust = 0.5, family=my_font),
  #   legend.text=element_text(color="black",size=7, hjust = 0.5, family=my_font),
  #   strip.text.x = element_text(color="black",size=9, family=my_font),
  #   strip.text.y = element_text(color="black",size=9, family=my_font),
  #   # legend.spacing.x = unit(0.1, 'mm'),
  #   # legend.spacing.y = unit(0.1, 'mm'),
  #   # legend.key.height=unit(1,"line"),
  #   legend.position = lg_pos,
  #   panel.background = element_blank(), 
  #   panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  # )
  # cnv_theme <- theme(strip.background = element_rect(fill = 'white', colour = 'white'),
  #                    axis.text.x = element_blank(),
  #                    axis.ticks.x = element_blank(),
  #                    axis.text.y = element_text(size=12, colour = "black", hjust = 0.5),
  #                    axis.title.y = element_text(size=12, colour = "black", hjust = 0.5),
  #                    axis.line = element_line(colour = "black"),
  #                    strip.placement = "outside",
  #                    legend.position = "bottom",
  #                    legend.text = element_text(size=9),
  #                    legend.title = element_text(size=9),
  #                    legend.key.size=unit(0.3,"cm"),
  #                    panel.grid.major = element_blank(),
  #                    panel.grid.minor = element_blank(),
  #                    # panel.background = element_rect(fill = "#F8F8F8", colour = NA),
  #                    panel.spacing = unit(c(0.1), 'cm'))
  cnv_plot <- cnv_plot + theme(strip.background = element_rect(fill = 'white', colour = 'white'),
                               text = element_text(color="black",size = 11, hjust = 0.5, family=my_font),
                               axis.text.x = element_blank(),
                               plot.title = element_text(color="black",size=12, hjust = 0, family=my_font, face="bold"),
                               axis.ticks.x = element_blank(),
                               axis.text.y = element_text(color="black",size=11, hjust = 0.5, family=my_font),
                               axis.title.y = element_text(color="black",size=11, hjust = 0.5, family=my_font),
                               axis.line = element_line(colour = "black"),
                               strip.placement = "outside",
                               legend.position = "bottom",
                               legend.text=element_text(color="black",size=9, hjust = 0.5, family=my_font),
                               legend.title=element_text(color="black",size=11, hjust = 0.5, family=my_font),
                               legend.key.size=unit(0.3,"cm"),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               # panel.background = element_rect(fill = "#F8F8F8", colour = NA),
                               panel.spacing = unit(c(0.1), 'cm'),
                               legend.margin=margin(0,0,0,0),
                               legend.box.margin=margin(-2,-2,-2,-2))    # MA: was 0.2
  
  cnv_plot <- cnv_plot + guides(fill = guide_legend(nrow = 1, override.aes = list(size=2))) +  #, override.aes = list(size=1.1)
    scale_x_continuous(expand = c(0,0))
  
  
  # cnv_plot
  lg <- cowplot::get_legend(cnv_plot)
  plg <- cowplot::ggdraw() + cowplot::draw_plot(lg)
  # cnv_plot <- cnv_plot + theme(legend.position = "none")
  
  if(datatag=='SA919'){
    color_fn <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/dlp_trees/colorCode_clones/color_code_SA919.csv'
    df <- data.table::fread(color_fn)
    df$val <- 1
    print(datatag) 
    # clones_cols <- c('A'='#66C2A5','B'='#FC8D62','C'='#8DA0CB')
    # df <- tibble::tibble(clone_id = c('A','B','C'), color=rep(1, length(clones_cols)))  
  }else if(datatag=='SA535'){
    color_fn <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/dlp_trees/colorCode_clones/color_code_SA535.csv'
    df <- data.table::fread(color_fn)
    clones_cols <- df$color
    names(clones_cols) <- df$clone_id
    df$val <- 1
    print(datatag) 
  }else{
    print(datatag) 
  }
  
  plt <- df %>% 
    ggplot(aes(x = val, y = clone_id, fill = clone_id)) + #, fill = cnv
    # geom_tile(aes(fill = cnv)) + #
    geom_raster() + 
    scale_fill_manual(values = clones_cols) + 
    theme_bw() + 
    labs(x = "Chromosome", y = "Clone", title=plttitle) + 
    theme(legend.position = "none",
          text = element_text(color="black",size = 11, hjust = 0.5),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          # axis.text = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank())
  # plt
  
  # test_plot <- cowplot::plot_grid(
  #   plt,
  #   res_919$cnv_plot,
  #   nrow = 1,
  #   rel_widths = c(0.05,1.2),
  #   align = 'h'#,
  #   # labels = c('a')
  # )
  # test_plot
  
  results <- list(df_cnv=df_cnv, plg=plg, cnv_plot=cnv_plot, clone_plt=plt)
  return(results)
}

get_library_id <- function(cell_ids) {
  
  labels <- lapply(strsplit(cell_ids, "-"), function(x) {
    return(x[2])
  })
  return(as.character(labels))
}
get_sample_id <- function(cell_ids) {
  labels <- lapply(strsplit(cell_ids, "-"), function(x) {
    return(x[1])
  })
  return(as.character(labels))
}




# get_median_genotype_v2 <- function(datatag, results_dir, copynumber_fn=NULL, cellclone_fn=NULL,
#                                    calcul_distance=F, distance_type='Manhattan'){
#   
#   if(is.null(copynumber_fn)){
#     copynumber_fn <- paste0(results_dir,'total_merged_filtered_states.csv')
#   }
#   if(is.null(cellclone_fn)){
#     cellclone_fn <- paste0(results_dir,'cell_clones.csv')  
#   }
#   save_dir <- paste0(results_dir,'CN_profile/')
#   if(!dir.exists(save_dir)){
#     dir.create(save_dir)
#   }
#   
#   copynumber <- read.csv(copynumber_fn, header=T, 
#                          row.names = 1, check.names = F,stringsAsFactors = FALSE)
#   
#   # cell_clones contain 2 columns of cell_id, and clone_id
#   # ex:           cell_id                      clone_id
#   # 1   SA535X4XB02498-A98163A-R09-C11          C
#   cell_clones <- read.csv(cellclone_fn, check.names = F,stringsAsFactors = FALSE)
#   
#   print(dim(copynumber))
#   cells_use <- colnames(copynumber)[colnames(copynumber) %in% cell_clones$cell_id]
#   copynumber <- copynumber[,cells_use]
#   copynumber$chr_desc <- rownames(copynumber)
#   cnv <- copynumber %>%
#     pivot_longer(!chr_desc, names_to = "cell_id", values_to = "copy_number")
#   
#   print(dim(cnv))
#   cnv <- cnv %>% inner_join(cell_clones, by=c("cell_id"))
#   
#   # get median genotype 
#   print("Get median genotype")
#   stat_cnv <- cnv %>%
#     dplyr::group_by(clone_id, chr_desc) %>%
#     dplyr::summarise(median_cn=median(copy_number)) #,mode_cn=calc_mode(copy_number)
#   
#   
#   median_cnv <- stat_cnv %>%
#     dplyr::select(chr_desc, median_cn, clone_id) %>%
#     group_by(clone_id) 
#   
#   
#   median_cnv <- median_cnv %>%
#     pivot_wider(names_from = clone_id, values_from = median_cn)
#   
#   median_cnv <- as.data.frame(median_cnv)
#   write.csv(median_cnv, paste0(save_dir,'median_cnv.csv'), quote = F, row.names = F)
#   median_cnv <- median_cnv %>%
#     tibble::column_to_rownames('chr_desc')
#   
#   # median_cnv <- median_cnv %>% select(-c(chr_desc))  
#   # chr_infos <- get_chr_infos(median_cnv$chr_desc)
#   # median_cnv$chr <- chr_infos$chr
#   # median_cnv$start <- chr_infos$start
#   # median_cnv$end <- chr_infos$end
#   # head(median_cnv)
#   median_cnv <- median_cnv[unique(chrs$chr_desc),]
#   # dim(median_cnv) # SA609
#   # rownames(median_cnv)[1]
#   # ls_features <- c()
#   # for(i in seq(dim(median_cnv)[1])){
#   #   if(sum(is.na(median_cnv[i,]))==0){
#   #     ls_features <- c(ls_features, i)
#   #   }
#   # }
#   # length(ls_features)
#   # median_cnv <- median_cnv[ls_features,]
#   # dim(median_cnv)
#   res <- compute_dist_mat(median_cnv, save_dir, use_hamming = F)
#   head(res$dist_to_median)
#   p <- plot_heatmap_genotype(res$dist_to_median, distance_type, datatag, save_dir)
#   res$hm <- p
#   dim(res$out_mtx)
#   # if(calcul_distance){
#   #   if(distance_type=='Hamming'){
#   #     dis_mtx <- compute_dist_mat(median_cnv, save_dir, use_hamming = T)
#   #     p <- plot_heatmap_genotype(dis_mtx, distance_type, datatag, save_dir)
#   #   }else if(distance_type=='Manhattan'){
#   #     dis_mtx <- compute_dist_mat(median_cnv, save_dir, use_hamming = F)
#   #     head(dis_mtx)
#   #     dim(dis_mtx)
#   #     p <- plot_heatmap_genotype(dis_mtx, distance_type, datatag, save_dir)
#   #   }else{ # get results for both case
#   #     dis_mtx1 <- compute_dist_mat(median_cnv, save_dir, use_hamming = T)
#   #     p <- plot_heatmap_genotype(dis_mtx1, 'Hamming', datatag, save_dir)
#   #     dis_mtx2 <- compute_dist_mat(median_cnv, save_dir, use_hamming = F)
#   #     p <- plot_heatmap_genotype(dis_mtx2, 'Manhattan', datatag, save_dir)
#   #   }
#   #   
#   # }
#   return(res)
# }

plot_heatmap_genotype <- function(dis_mtx, distance_type, datatag, results_dir){
  cell_func = function(j, i, x, y, width, height, fill) {
    str <- ''
    if(!is.na(dis_mtx[i, j])){
      str <- sprintf("%.2f", dis_mtx[i, j])
    }
    grid.text(str, x, y, gp = gpar(fontsize = 10))
  }
  # dis_mtx[dis_mtx==0]
  # View(dis_mtx)
  row_lbs <- as.character(rownames(dis_mtx))
  # distance between median CNA profiles of clones
  plottitle <- paste0(distance_type," copy number distance ",datatag)
  p <- ComplexHeatmap::Heatmap(as.matrix(dis_mtx), na_col = "white",
                               show_column_names=T,
                               show_row_names = T,
                               cluster_rows=F,
                               cluster_columns=F,
                               name = paste0(distance_type," dist"), 
                               # row_order = sort(rownames(test)),
                               row_split=row_lbs,
                               # row_title = "%s",
                               # column_title_side = T,
                               # row_title_rot = 90,
                               row_names_side = "left",
                               column_names_side = "bottom",
                               column_names_rot = 0,
                               row_gap = unit(1, "mm"),
                               column_split = row_lbs,
                               row_title = "Clones",
                               column_title = plottitle,
                               column_gap = unit(1, "mm"),
                               column_names_gp = grid::gpar(fontsize = 13, fontface = "bold"),
                               column_title_gp = gpar(fontsize = 13, fontface = "bold"),
                               row_names_gp = grid::gpar(fontsize = 13, fontface = "bold"),
                               show_heatmap_legend = T,
                               # top_annotation=top_anno,
                               # left_annotation = left_anno,
                               cell_fun = cell_func,
  )#row_dend_reorder=F
  # p
  
  png(paste0(results_dir, datatag, distance_type,'_distance_','.png'), height = 2*500, width=2*600, res = 2*72)
  print(p)
  dev.off()
  return(p)
}

# Get mode values
calc_mode <- function(x) {
  keys <- unique(x)
  keys[which.max(tabulate(match(x, keys)))]
}

