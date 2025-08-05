suppressPackageStartupMessages({
  library(tidyverse)
  library(annotables)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  # library(scales)
  # library(ggrepel)
  library(stringr)
  # library(scran)
  library(SingleCellExperiment)
  # library(tximport)
  library(edgeR)
})

## Loading utility functions
script_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/'
source(paste0(script_dir, 'scripts/bulk_rna/bulk_utils.R'))


viz_genewise_dispersion <- function(df, obs_clones=c('B','C'), legend_pos='bottom'){
  # unique(df$gene_type)
  cutoff_values <- 2
  df <- df %>% 
    dplyr::mutate(
      gene_type=
        case_when(
          gene_type %in% c('trans','trans_whole_genome') ~ 'trans',
          TRUE ~ gene_type))  %>% 
    dplyr::rename(gene_wise_dispersion=tagwise_dispersion) %>% 
    dplyr::mutate(
      gene_wise_dispersion=
        case_when(
          gene_wise_dispersion > cutoff_values ~ cutoff_values, 
          TRUE ~ gene_wise_dispersion
        ))
  ptittle <- paste0(obs_clones[2],' vs. ', obs_clones[1])
  obs_clones <- paste0('Clone ', obs_clones)
  # color_use <- c('Clone A'='#66C2A5','Clone B'='#FC8D62','Clone C'='#8DA0CB') 
  
  cols_use <- c('cis'='#EE220C','trans'='#0076BA')
  print(cols_use)
  df$gene_type <- as.factor(df$gene_type)
  
  p.density <- ggplot(df, aes(x = gene_wise_dispersion, group = gene_type)) + geom_density()    
  p.density.data <- layer_data(p.density) %>%
    select(x, y, group) %>%
    mutate(gene_type = factor(group, labels = levels(df$gene_type), ordered = TRUE)) %>%
    select(-group)
  p.density.data <- p.density.data %>%
    rbind(p.density.data %>% 
            group_by(gene_type) %>% 
            filter(x == min(x)) %>% 
            mutate(y = 0) %>% 
            ungroup())
  my_font <- "Helvetica"
  p1 <- ggplot(df, aes(x = 10, y = gene_wise_dispersion)) + #-1.5
    # manually flipped density plot
    geom_polygon(data = p.density.data, aes(x = y, y = x, fill = gene_type), 
                 alpha=0.4, color = "black") + #fill = NA, 
    
    # box plot
    # geom_violin(aes(fill = cut), alpha=0.3) + #, group = cut #trim=T
    geom_boxplot(aes(fill = gene_type), width=1.5, alpha=0.4, outlier.alpha = 0) + 
    
    # vertical lines at Q1 / Q2 / Q3
    # stat_boxplot(geom = "hline", aes(yintercept = ..lower..)) +
    # stat_boxplot(geom = "hline", aes(yintercept = ..middle..)) +
    # stat_boxplot(geom = "hline", aes(yintercept = ..upper..)) +
    
    # facet_grid(cut ~ .) +
    scale_fill_manual(values=cols_use) +
    coord_flip() + 
    theme_bw() + 
    theme(legend.position = legend_pos,
          panel.grid = element_blank(), 
          text = element_text(size=10, color="black",hjust = 0.5, family=my_font),
          axis.text.x = element_text(color="black",size=10, hjust = 0.5, family=my_font),
          axis.title.x=element_text(color="black",size=11, hjust = 0.5, family=my_font),
          axis.text.y = element_text(color="black",size=10, hjust = 0.5, family=my_font),
          axis.title.y=element_text(color="black",size=11, hjust = 0.5, family=my_font),
          legend.key.size = unit(1,"line"),
          # legend.title = element_blank()
          ) +
    # scale_size(guide=FALSE)+
    labs(x='Density', y='Gene wise dispersion', title = ptittle) +
    guides(colour = guide_legend(override.aes = list(shape = 15)))
  # p1
  return(p1)
  
}
get_median_cnv_cis_genes <- function(obs_genes, obs_clones){
  script_dir <- "/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/"
  cnv <- data.table::fread(paste0(script_dir, 'SA919_DE_analysis_DESeq2_Hoa_09April2024/mapping_gene_cnv_SA919.csv.gz'))  
  cnv <- cnv[!duplicated(cnv$ensembl_gene_id),]
  dim(cnv)
  # head(cnv)
  ## need to generate a cnv profile for each patient
  rownames(cnv) <- cnv$ensembl_gene_id
  cnv <- as.data.frame(cnv)
  cols_use <- c('ensembl_gene_id','gene_symbol','chr') 
  cols_use <- c(cols_use, obs_clones)
  rownames(cnv) <- cnv$ensembl_gene_id
  cnv <- cnv %>%
    dplyr::select(all_of(cols_use))
  # genes_used <- cnv$ensembl_gene_id[rv>0]
  # length(genes_used)
  # cnv <- cnv %>%
  #   dplyr::filter(gene_symbol %in% obs_genes)
  cnv <- cnv %>%
    dplyr::filter(ensembl_gene_id %in% obs_genes)
  dim(cnv)
  
  cnv_pv <- cnv %>%
    pivot_longer(cols = all_of(obs_clones), names_to='clone', values_to = 'cnv')
  dim(cnv_pv)
  head(cnv_pv)
  stat <- cnv_pv %>%
    group_by(clone, chr) %>%
    dplyr::summarise(cnv=median(cnv))
  
  # stat
  # length(genes_used)
  return(stat)
}
viz_cnv_profiles_specific_chrs_horizontal <- function(df){
  df <- stat
  chrs <- unique(df$chr)
  chr_levels <- c(as.character(1:23), "X", "Y")
  chr_levels <- chr_levels[chr_levels %in% chrs]
  df$chr <- factor(df$chr, levels = chr_levels)
  
  vals <- sort(as.numeric(unique(df$cnv)), decreasing=F)
  df$cnv <- as.character(df$cnv)
  df$cnv <- factor(df$cnv, levels = as.character(vals))
  
  df <- df %>%
    dplyr::arrange(chr, clone)
  
  df <- df %>%
    dplyr::mutate(chr_clone = paste0(chr,'_', clone))
  df$chr_clone <- factor(df$chr_clone, levels=df$chr_clone)
  cnv_cols <- c('0'='#4880B8', '1'='#A7C9DF','2'='#CCCCCC','3'='#F5CE93','4'='#ED9364',
                '5'='#D2553E','6'='#A42116','7'='#8B1A43','8'='#CB3576','9'='#D06CAD',
                '10'='#C196C4','11'='#D0BAD8')
  
  my_font <- 'Helvetica'
  p <- ggplot(df, aes(x=chr_clone, y=1, fill= cnv)) + 
    geom_tile(colour = "white", size=1.2) + 
    # facet_wrap(Y~., scales = "free") + 
    theme_bw() + 
    scale_fill_manual(values = cnv_cols) + 
    theme(plot.title = element_blank(),
          axis.line=element_blank(),
          axis.title = element_text(size=11, family = my_font, colour = 'black'),
          axis.title.x = element_blank(),
          # axis.title.y = element_text(size=11), 
          # axis.text.y = element_text(size=11),
          axis.ticks.y = element_blank(),
          axis.text.y = element_text(size=13, family = my_font, colour = 'black'),
          # axis.text.x = element_text(size=13, family = my_font, colour = 'black'),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          # panel.background = element_rect(fill = "transparent", colour = NA),
          # plot.background = element_rect(fill = "white", colour = NA),
          panel.grid = element_blank(),
          panel.border = element_blank())+ 
    labs(y=' ') 
  # p
  return(p)
}

viz_cnv_profiles_specific_chrs <- function(df){
  chrs <- unique(df$chr)
  chr_levels <- c(as.character(1:23), "X", "Y")
  chr_levels <- chr_levels[chr_levels %in% chrs]
  df$chr <- factor(df$chr, levels = chr_levels)
  
  cnv_cols <- c('0'='#4880B8', '1'='#A7C9DF','2'='#CCCCCC','3'='#F5CE93','4'='#ED9364',
                '5'='#D2553E','6'='#A42116','7'='#8B1A43','8'='#CB3576','9'='#D06CAD',
                '10'='#C196C4','11'='#D0BAD8')
  vals <- sort(as.numeric(unique(df$cnv)), decreasing=F)
  df$cnv <- as.character(df$cnv)
  df$cnv <- factor(df$cnv, levels = as.character(vals))
  my_font <- 'Helvetica'
  p <- ggplot(df, aes(x=clone, y=chr, fill= cnv)) + 
    geom_tile(colour = "white", size=1.2) + #, alpha=1
    # facet_wrap(Y~., scales = "free") + 
    theme_bw() + 
    scale_fill_manual(values = cnv_cols) + 
    theme(axis.line=element_blank(),
          axis.title = element_text(size=10, family = my_font, colour = 'black'),
          # axis.title.x = element_blank(), 
          # axis.title.y = element_text(size=11), 
          # axis.text.y = element_text(size=11),
          axis.ticks.y = element_blank(),
          axis.text.y = element_text(size=10, family = my_font, colour = 'black'),
          axis.text.x = element_text(size=12, family = my_font, colour = 'black'),
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          panel.grid = element_blank(),
          panel.border = element_blank())+ 
    labs(x='Chromosome',y='Clone',title='') 
  
  # p
  return(p)
}

# df <- res$total_exp_df
viz_gene_exp_comparison <- function(df, obs_clones, legend_pos='bottom'){
  # df <- total_exp_df
  colors_code <- c('Metastasis'='#EE220C','Primary'='#0076BA')
  conds <- unique(df$obs_conds)
  colors_use1 <- ifelse(grepl('Metastasis',conds[1]),colors_code[1], colors_code[2])
  names(colors_use1) <- conds[1]
  colors_use2 <- colors_code[colors_code!=colors_use1]
  names(colors_use2) <- conds[2]
  colors_use <- c(colors_use1, colors_use2)
  print(colors_use)
  summ_df <- df %>% group_by(chr, obs_conds) %>% summarize(median=median(gene_exp),
                                                           # Lower=quantile(gene_exp, probs=0.25),
                                                           # Upper=quantile(gene_exp, probs=0.75)
  ) %>% 
    pivot_longer(cols = median, names_to = "quantile", values_to = "estimate")
  # head(summ_df)
  ptittle <- paste0(obs_clones[2],' vs. ', obs_clones[1])
  # Convert the variable dose from a numeric to a factor variable
  # df$chr <- as.factor(df$chr)
  chrs <- unique(df$chr)
  chrs_levels <- c(c(1:22),'X')
  chrs_levels <- chrs_levels[chrs_levels %in% chrs]
  df$chr <- factor(df$chr, levels = chrs_levels)
  summ_df$chr <- factor(as.character(summ_df$chr), levels = chrs_levels)
  print(chrs_levels)
  # minY <- min(df$gene_exp)
  # maxY <- max(df$gene_exp)+2 # just to have space at top area
  my_font <- "Helvetica"
  p <- ggplot(df, aes(x=chr, y=gene_exp, color=obs_conds)) +
    # geom_jitter(aes(color=obs_conds), shape=16, position=position_jitter(0.2)) +
    geom_violin(position=position_dodge(0.8), alpha=0) +
    geom_jitter(position=position_jitterdodge(0.2), size=0.1) +
    geom_point(data = summ_df, aes(x=chr, y=estimate, group=obs_conds), 
               position = position_dodge(width=0.8), color='black', fill='black',shape=23) + 
    scale_color_manual(values = colors_use) + 
    # ylim(minY, maxY) + 
    # geom_boxplot(width=1)
    theme_bw() + 
    theme(legend.position = legend_pos,
          panel.grid = element_blank(), 
          legend.text = element_text(color="black",size=10, hjust = 0.5, family=my_font),
          text = element_text(color="black",size=10,hjust = 0.5, family=my_font),
          axis.text.x = element_text(color="black",size=12, hjust = 0.5, family=my_font),
          axis.title.x=element_text(color="black",size=10, hjust = 0.5, family=my_font),
          axis.text.y = element_text(color="black",size=10, hjust = 0.5, family=my_font),
          axis.title.y=element_text(color="black",size=10, hjust = 0.5, family=my_font),
          legend.key.size = unit(0.6,"line"),
          legend.title = element_blank()) +
    labs(title = ptittle, x='Chromosome position', y= 'Log2(exp + 1)') +
    guides(colour = guide_legend(override.aes = list(shape = 15)))
  # p
  # geom_jitter(shape=16, position=position_jitter(0.2)) #+
  # geom_dotplot(aes(y=len), dotsize=0.5)# + #binaxis='y', stackdir='center', 
  # geom_boxplot(aes()width=0.1, fill="white")#+
  # scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  
  return(p)
}


viz_Fig6_partA <- function(){
  save_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/results_bulkRNAseq/SA919_full/'
  datatag <- 'SA919_full'
  
  obs_conds <- c('B_Metastasis','A_Primary')
  subtag <- 'Bmet_Apri'
  save_fig_dir <- paste0(save_dir, subtag, '/')
  obs_clones <- c('A','B')
  df <- data.table::fread(paste0(save_fig_dir, subtag, '_dispersion_cis_trans_genes.csv'))
  dim(df)
  head(df)
  p1 <- viz_genewise_dispersion(df, obs_clones, legend_pos='bottom')
  saveRDS(p1, paste0(save_fig_dir, subtag, '_gene_dispersion_with_legend.rds'))
  p1
  
  obs_conds <- c('C_Metastasis','A_Primary')
  subtag <- 'Cmet_Apri'
  save_fig_dir <- paste0(save_dir, subtag, '/')
  obs_clones <- c('A','C')
  df <- data.table::fread(paste0(save_fig_dir, subtag, '_dispersion_cis_trans_genes.csv'))
  dim(df)
  head(df)
  p2 <- viz_genewise_dispersion(df, obs_clones, legend_pos='bottom')
  saveRDS(p2, paste0(save_fig_dir, subtag, '_gene_dispersion_with_legend.rds'))
  p2
  
  obs_conds <- c('C_Metastasis','B_Primary')
  subtag <- 'Cmet_Bpri'
  save_fig_dir <- paste0(save_dir, subtag, '/')
  obs_clones <- c('B','C')
  df <- data.table::fread(paste0(save_fig_dir, subtag, '_dispersion_cis_trans_genes.csv'))
  dim(df)
  p3 <- viz_genewise_dispersion(df, obs_clones, legend_pos='bottom')
  saveRDS(p3, paste0(save_fig_dir, subtag, '_gene_dispersion_with_legend.rds'))
  p3
  
  
  p_total_partA <- cowplot::plot_grid(p1, 
                                      p2, 
                                      p3, nrow=1,
                                      rel_widths = c(1,1,1))
  saveRDS(p_total_partA, paste0(save_dir,"figs/","gene_exp_Fig6_partA.rds"))
  p_total_partA <- readRDS(paste0(save_dir,"figs/","gene_exp_Fig6_partA.rds"))
  png(paste0(save_dir,"figs/","gene_exp_Fig6_partA.png"), 
      height = 2*270, width=2*900, res = 2*72)
  print(p_total_partA)
  dev.off() 
  
  ggsave(  
    filename = paste0(save_dir,"figs/","gene_exp_Fig6_partA.svg"),  
    plot = p_total_partA,  
    height = 3.5, width = 11, dpi = 150)
  
}
viz_Fig6_partB <- function(){
  save_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/results_bulkRNAseq/SA919_full/'
  datatag <- 'SA919_full'
  
  obs_conds <- c('B_Metastasis','A_Primary')
  subtag <- 'Bmet_Apri'
  save_fig_dir <- paste0(save_dir, subtag, '/')
  obs_clones <- c('A','B')
  df <- data.table::fread(paste0(save_fig_dir, paste(obs_conds, collapse='_vs_'), '_total_exp.csv.gz'))
  dim(df)
  obs_genes <- unique(df$ens_gene_id)
  stat <- get_median_cnv_cis_genes(unique(df$ens_gene_id), obs_clones)
  stat
  p10 <- viz_cnv_profiles_specific_chrs_horizontal(stat)
  p11 <- viz_gene_exp_comparison(df, obs_clones, legend_pos='bottom')
  # saveRDS(p1, paste0(save_fig_dir, subtag, '_gene_exp_plt_with_legend.rds'))
  
  
  obs_conds <- c('C_Metastasis','A_Primary')
  subtag <- 'Cmet_Apri'
  save_fig_dir <- paste0(save_dir, subtag, '/')
  obs_clones <- c('A','C')
  df <- data.table::fread(paste0(save_fig_dir, paste(obs_conds, collapse='_vs_'), '_total_exp.csv.gz'))
  dim(df)
  stat <- get_median_cnv_cis_genes(unique(df$ens_gene_id), obs_clones)
  stat
  p20 <- viz_cnv_profiles_specific_chrs_horizontal(stat)
  p21 <- viz_gene_exp_comparison(df, obs_clones, legend_pos='bottom')
  # saveRDS(p2, paste0(save_fig_dir, subtag, '_gene_exp_plt_with_legend.rds'))
  # p2
  
  obs_conds <- c('C_Metastasis','B_Primary')
  subtag <- 'Cmet_Bpri'
  save_fig_dir <- paste0(save_dir, subtag, '/')
  obs_clones <- c('B','C')
  df <- data.table::fread(paste0(save_fig_dir, paste(obs_conds, collapse='_vs_'), '_total_exp.csv.gz'))
  dim(df)
  stat <- get_median_cnv_cis_genes(unique(df$ens_gene_id), obs_clones)
  stat
  p30 <- viz_cnv_profiles_specific_chrs_horizontal(stat)
  p31 <- viz_gene_exp_comparison(df, obs_clones, legend_pos='bottom')
  # saveRDS(p3, paste0(save_fig_dir, subtag, '_gene_exp_plt_with_legend.rds'))
  # p3
  
  # xpos <- 0.2
  # ypos <- 0.9
  # + theme(legend.position = c(xpos, ypos))
  p13 <- cowplot::plot_grid(p10 + theme(legend.position = 'top'),
                            p11,
                            ncol=1, rel_heights = c(0.7,3))
  p23 <- cowplot::plot_grid(p20 + theme(legend.position = 'top'),
                            p21,
                            ncol=1, rel_heights = c(0.7,3))
  
  p33 <- cowplot::plot_grid(p30 + theme(legend.position = 'top'),
                            p31,
                            ncol=1, rel_heights = c(0.7,3))
  
  # p13
  
  p_total_partB <- cowplot::plot_grid(p13, 
                                      p23,
                                      p33,
                                      nrow=1,
                                      rel_widths = c(1.8,3.1,1.8))
  p_total_partB
  saveRDS(p_total_partB, paste0(save_dir,"figs/","gene_exp_Fig6_partB.rds"))
  p_total_partB <- readRDS(paste0(save_dir,"figs/","gene_exp_Fig6_partB.rds"))
  png(paste0(save_dir,"figs/","gene_exp_Fig6_partB.png"), 
      height = 2*200, width=2*650, res = 2*72)
  print(p_total_partB)
  dev.off() 
  
    
}
viz_fig6 <- function(){
  
  save_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/results_bulkRNAseq/SA919_full/'
  datatag <- 'SA919_full'
  
  p_total_partA <- readRDS(paste0(save_dir,"figs/","gene_exp_Fig6_partA.rds"))
  p_total_partB <- readRDS(paste0(save_dir,"figs/","gene_exp_Fig6_partB.rds"))
  
  p_total_partAB <- cowplot::plot_grid(p_total_partA, 
                                       p_total_partB,
                                      ncol=1,
                                      rel_heights = c(1,1.2))
  p_total_partAB <- cowplot::plot_grid(p_total_partA, 
                                       p_total_partB,
                                       ncol=1,
                                       rel_heights = c(1,1.3))
  
  ggsave(paste0(save_dir,"Fig6_SA919_Pt1_trackplot_panelBC_1.svg"),
         plot = p_total_partAB,
         height = 10,
         width = 12,
         # useDingbats=F,
         dpi = 150
  )
  dev.off()
  
  
}


