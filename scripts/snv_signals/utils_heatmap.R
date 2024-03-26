
# source("/home/htran/Projects/hakwoo_project/corrupt_tree/src/cn_change/utils.R")
suppressPackageStartupMessages({
  # require("optparse")
  require("data.table")
  # require("feather")
  require("dplyr")
  require(scales)
  require("tidyr")
  require("ggplot2")
})

# library(extrafont)
# font_import(prompt=F) # import all your fonts
# fonts()


# df_cnv has column names: ensembl_gene_id
plot_CNV_corrupt_tree <- function(df_cnv, clones){
  # df_cnv <- read.csv(df_cnv_fn, check.names = F, stringsAsFactors = F)
  df_cnv <- df_cnv %>% 
    dplyr::select(c(clones, 'chr_desc')) %>% 
    tidyr::pivot_longer(!chr_desc, names_to = 'clone', values_to='cnv')%>% 
    as.data.frame()
  print(dim(df_cnv))
  df_cnv <- get_chr_infos(df_cnv)
    
  df_cnv <- df_cnv %>%
    group_by(chr) %>%
    # dplyr::mutate(start_order = rank((start+end)/2)) %>%
    dplyr::mutate(start_order = rank(start)) %>%
    ungroup() 
  
  chr_levels <- c(as.character(1:22), "X")
  df_cnv <- df_cnv %>% dplyr::filter(chr %in% chr_levels)
  df_cnv$chr <- factor(df_cnv$chr, levels = chr_levels)
  # summary(as.factor(df_cnv$chr))
  df_cnv <- tidyr::drop_na(df_cnv)
 
  
  cnv_cols <- c('0'='#4880B8', '1'='#A7C9DF','2'='#CCCCCC','3'='#F5CE93','4'='#ED9364',
                '5'='#D2553E','6'='#A42116','7'='#8B1A43','8'='#CB3576','9'='#D06CAD',
                '10'='#C196C4','11'='#D0BAD8')
 
  df_cnv$cnv <- as.character(round(df_cnv$cnv))
  levels(df_cnv$cnv) <- 0:(length(cnv_cols)-1)
  update_clones <- lapply(strsplit(as.character(df_cnv$clone), "_"), function(x) {
    return(as.character(x[1]))
  })
  if(length(unique(update_clones))!=1){
    df_cnv$clone <- update_clones  
    df_cnv$clone <- factor(df_cnv$clone, levels = unique(df_cnv$clone))
  }else{
    df_cnv$clone <- factor(df_cnv$clone, levels = c(clones[grepl('Primary',clones)],clones[grepl('Metastasis',clones)]))
  }
  
  
  # df_cnv$clone <- factor(df_cnv$clone, levels = c(clones[1],clones[2]))
  print(summary(df_cnv$clone))
  # saveRDS(df_cnv, paste0(save_fig_dir,'df_cnv_1035_test.rds'))
  cnv_plot <- df_cnv %>% 
    ggplot(aes(x = start_order, y = clone, fill = cnv)) + #, fill = cnv
    # geom_tile(aes(fill = cnv)) + #
    geom_raster() +
    # geom_tile() + 
    facet_wrap(~ chr, scales = "free_x", nrow = 1, drop = F, strip.position = "bottom") + # , switch = "x"
    # facet_grid(~ chr, scales = "free_x", switch = "x", space='free') +
    # theme(legend.position = "bottom", axis.text.x = element_blank()) +
    scale_y_discrete(expand = c(0, 0)) +
    # scale_fill_manual(values=cnv_colors, name = "Copy number", guide = 'legend',labels = 0:(length(cnv_colors)-1),drop=FALSE)  +
    #          theme(legend.position = "bottom") +
    scale_fill_manual(values = cnv_cols, name = "Copy number ", breaks = names(cnv_cols)) +  #
    labs(x = "Chromosome", y = "Clone",title ='Median CNV')

  lg_pos <- "bottom"
  my_font <- "Helvetica"
  thesis_theme <- theme(strip.background = element_rect(fill = 'white', colour = 'white'),
                        text = element_text(size = 8, hjust = 0.5, family=my_font),
                        plot.title = element_text(color="black",size=9, face="bold", hjust=0, family=my_font),
                        plot.subtitle = element_text(size=9, hjust=0, family=my_font),
                        axis.text.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.title.x = element_text(size=7, hjust = 0.5, family=my_font),
                        axis.text.y = element_text(size=7, hjust = 0.5, family=my_font),
                        axis.title.y = element_text(size=6, hjust = 0.5, family=my_font),
                        axis.line = element_line(colour = "black"),
                        strip.placement = "outside",
                        # strip.text.x = element_text(size = 7, hjust = 0.5, family=my_font),
                        # strip.text.x =element_blank(),
                        legend.position = lg_pos,
                        legend.text=element_text(size=6, hjust = 0.5, family=my_font),
                        legend.title=element_text(size=7, hjust = 0.5, family=my_font),
                        legend.key.size=unit(0.2,"cm"),
                        legend.margin=margin(0,0,0,0),
                        legend.box.margin=margin(-5,-5,-5,-5),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        # panel.background = element_blank(), 
                        panel.background = element_rect(fill = "#F8F8F8", colour = NA),
                        panel.spacing = unit(c(0.1), 'cm'))
  cnv_plot <- cnv_plot + thesis_theme   
  cnv_plot <- cnv_plot + guides(fill = guide_legend(nrow = 1, override.aes = list(size=0.01))) +  #, override.aes = list(size=1.1)
    scale_x_continuous(expand = c(0,0))
  # png(paste0(output_dir, 'cloneT_segment_profile.png'), height = 2*250, width=2*1000, res = 2*72)
  # print(cnv_plot)
  # dev.off()
  return(cnv_plot)
}

get_library_id <- function(cell_ids){
  lids <- lapply(strsplit(cell_ids,'-'), function(x){
    return(x[2])
  })
  lids <- as.character(lids)
  return(lids)
}

get_sample_id <- function(cell_ids){
  sids <- lapply(strsplit(cell_ids,'-'), function(x){
    return(x[1])
  })
  sids <- as.character(sids)
  return(sids)
}

get_chr_infos <- function(median_cnv) {
  chr_desc <- as.character(median_cnv$chr_desc)
  obs_chrs = as.character(c(paste0(1:22), "X"))
  
  chrs <- sapply(strsplit(chr_desc, "_"), function(x) {
    return(x[1])
  })
  starts <- sapply(strsplit(chr_desc, "_"), function(x) {
    return(x[2])
  })
  ends <- sapply(strsplit(chr_desc, "_"), function(x) {
    return(x[3])
  })
  median_cnv$chr <- as.character(chrs)
  median_cnv$start <- as.numeric(starts)
  median_cnv$end <- as.numeric(ends)
  median_cnv <- median_cnv %>% 
    dplyr::filter(chr %in% obs_chrs)
  # return(list(chr=as.character(chrs),start=as.character(starts),end=as.character(ends)))
  return(median_cnv)
}

get_median_LOH_genotype <- function(ascn_fn, 
                                    library_grouping_fn, cellclone_fn, 
                                    datatag, save_dir){
  if(file.exists(ascn_fn)){
    ascn_raw <- readRDS(ascn_fn)
  }else{
    stop("File ascn results does not exist!!!")
  }
  
  # length(ascn_raw)
  # summary(as.factor(ascn_raw$data$LOH))
  # summary(as.factor(ascn$LOH))
  ascn <- ascn_raw$data
  print(dim(ascn))
  print(colnames(ascn))
  
  # ascn %>%
  #   group_by(chr, start, end) %>%
  #   summarise(state = schnapps:::Mode(state),
  #             state_min = schnapps:::Mode(state_min),
  #             BAF = median(BAF),
  #             LOH = schnapps:::Mode(LOH),
  #             state_phase = schnapps:::Mode(state_phase),
  #             copy = median(copy)) %>%
  #   ungroup() 
  
  cell_clones <- read.csv(cellclone_fn, check.names = F, stringsAsFactors = FALSE)
  
  metasample <- read.csv(library_grouping_fn, check.names = F, stringsAsFactors = FALSE)
  # dim(metasample)
  # head(metasample)
  cell_clones <- cell_clones %>%
    dplyr::filter(!clone_id %in% c('None','unassigned'))
  
  
  metasample <- metasample %>%
    dplyr::rename(library_id=grouping) %>%
    dplyr::select(library_id, mainsite, sample_id)
  
  print(dim(cell_clones))
  cell_clones <- cell_clones %>% 
    dplyr::filter(cell_id %in% unique(ascn$cell_id))
  print(dim(cell_clones))
  
  cell_clones$library_id <- get_library_id(cell_clones$cell_id)
  cell_clones$sample_id <- get_sample_id(cell_clones$cell_id)
  cell_clones <- cell_clones %>% left_join(metasample, by=c("library_id","sample_id"))
  meta_cells <- cell_clones
  meta_cells <- meta_cells %>%
    dplyr::group_by(clone_id, mainsite) %>%
    dplyr::summarise(nb_cells=n()) %>%
    dplyr::ungroup()
  # View(meta_cells)
  
  meta_cells <- meta_cells %>%
    dplyr::filter(nb_cells>=20 & clone_id !='None')
  print(dim(meta_cells))
  # colnames(meta_cells)
  write.csv(meta_cells, paste0(save_dir,'meta_cells.csv'), quote = F, row.names = F)
  
  
  # dim(ascn)
  # head(ascn)
  # unique(ascn$LOH)
  # ascn <- ascn %>%
  #   dplyr::filter(LOH=='LOH')
  ascn <- ascn %>% inner_join(cell_clones, by = "cell_id")
  ascn <- ascn %>% inner_join(meta_cells, by=c("clone_id","mainsite"))
  dim(ascn)
  
  ascn$clone <- paste0(ascn$clone_id,'_',ascn$mainsite)
  # length(unique(ascn$clone))
  
  # ascn1 <- ascn
  print("Get median genotype")
  ascn$chr_desc <- paste0(ascn$chr,'_',ascn$start,'_',ascn$end)
  # ascn1 <- ascn1 %>%
  #   dplyr::group_by(clone_label, chr_desc) %>%
  #   dplyr::summarise(mode_cn=calc_mode(state))
  # 
  ascn <- ascn %>%
    dplyr::group_by(clone, chr_desc) %>%
    dplyr::summarise(cnv=median(state),
                     LOH = schnapps:::Mode(LOH)) %>%
    ungroup() 
  # dim(ascn1)
  ascn <- get_chr_infos(ascn)
  print(dim(ascn))
  # ascn1 <- ascn1 %>%
  #   dplyr::mutate(LOH_event = ifelse(cnv==2,'cnLOH',
  #                                    ifelse(cnv<2,'LOH', 'ampLOH')))
  
  ascn <- ascn %>%
    dplyr::mutate(LOH_event = case_when(
      LOH=='LOH' & cnv==2 ~ 'cnLOH',
      LOH=='LOH' & cnv<2 ~ 'LOH',
      LOH=='LOH' & cnv>2 ~ 'ampLOH',
      TRUE ~ 'others'
    ))
  res <- list(meta_cells=meta_cells, cell_clones=cell_clones, ascn=ascn)
  saveRDS(res, paste0(save_dir,'median_genotype.rds'))
  return(res)
}

detect_LOH_change <- function(copynumber_fn, res, save_dir, datatag){
  meta_cells <- res$meta_cells
  dim(meta_cells)
  
  ## Form comparisons, and detect the LOH changes between pair of clones. 
  if(datatag=='SA535'){
    cl1_ls <- c("A","E","N","O","K","S","G","H")
    cl2_ls <- c("B","F","O","P","L","T","G","H")
    t1 <- c("Primary","Metastasis","Primary","Metastasis","Metastasis","Primary","Primary","Primary")
    t2 <- c("Metastasis","Metastasis","Metastasis","Metastasis","Metastasis","Metastasis","Metastasis","Metastasis")
    # cl1_ls <- c("A")
    # cl2_ls <- c("B")
    # t1 <- c("Primary")
    # t2 <- c("Metastasis")
    de <- data.frame(cl1=cl1_ls,cl2=cl2_ls,t1=t1,t2=t2)
  }
  if(datatag=='SA919'){
    cl1_ls <- c("A","A","A","B","B","A","A","A","B")
    cl2_ls <- c("B","B","B","C","C","C","C","A","B")
    t1 <- c("Primary","Primary","Metastasis","Primary","Metastasis","Primary","Metastasis","Primary","Primary")
    t2 <- c("Primary","Metastasis","Metastasis","Metastasis","Metastasis","Metastasis","Metastasis","Metastasis","Metastasis")
    # cl1_ls <- c("A")
    # cl2_ls <- c("B")
    # t1 <- c("Primary")
    # t2 <- c("Metastasis")
    de <- data.frame(cl1=cl1_ls,cl2=cl2_ls,t1=t1,t2=t2)
  }
  # copynumber_fn <- paste0(dirname(save_dir),'/total_merged_filtered_states.csv')
  copynumber <- as.data.frame(data.table::fread(copynumber_fn))
  rownames(copynumber) <- copynumber$V1
  copynumber$V1 <- NULL
  filtered_bins <- rownames(copynumber)
  rm(copynumber)
  
  # median_cnv <- read.csv(paste0(dirname(save_dir),'/CN_profile/median_cnv.csv'), check.names = F, stringsAsFactors = F)
  
  for(i in seq(1:nrow(de))){
    cl1 <- de$cl1[i]
    cl2 <- de$cl2[i]
    t1 <-  de$t1[i]
    t2 <-  de$t2[i]
    # desc1 <- paste0(tolower(t1)," cells (n=",meta_cells[meta_cells$clone_id==cl1 & meta_cells$mainsite==t1,'nb_cells'],")")
    # desc2 <- paste0(tolower(t2)," cells (n=",meta_cells[meta_cells$clone_id==cl2 & meta_cells$mainsite==t2,'nb_cells'],")")
    desc1 <- paste0('Clone ',cl1,' (',meta_cells[meta_cells$clone_id==cl1 & meta_cells$mainsite==t1,'nb_cells'],' ',tolower(t1)," cells) ")
    desc2 <- paste0('Clone ',cl2,' (',meta_cells[meta_cells$clone_id==cl2 & meta_cells$mainsite==t2,'nb_cells'],' ',tolower(t2)," cells) ")
    
    print(desc1)
    print(desc2)
    viz_LOH_events(filtered_bins, res$ascn, c(paste0(cl1,'_',t1),desc1), c(paste0(cl2,'_',t2),desc2), 
                   paste0(save_dir,cl1,'_',t1,'_vs_',cl2,'_',t2,'/'))
    
  }
}

# df_cnv <- res$ascn
# obs_clone1 <- c(paste0(cl1,'_',t1),desc1)
# obs_clone2 <- c(paste0(cl2,'_',t2),desc2)
# save_dir <- paste0(save_dir,cl1,'_',t1,'_vs_',cl2,'_',t2,'/')
viz_LOH_events <- function(filtered_bins, df_cnv, obs_clone1, obs_clone2, save_dir){
  
  if(!dir.exists(save_dir)){
    dir.create(save_dir)
  }
  df_cnv <- df_cnv %>%
    dplyr::filter(clone %in% c(obs_clone1[1], obs_clone2[1]))
  
  # head(df_cnv)
  
  # genes1 <- df_cnv %>%
  #               dplyr::filter(clone==obs_clone1[1] & !is.na(cnv)) %>%
  #               dplyr::pull(chr_desc)
  # 
  # genes2 <- df_cnv %>%
  #   dplyr::filter(clone==obs_clone2[1] & !is.na(cnv)) %>%
  #   dplyr::pull(chr_desc)
  # 
  # common <- intersect(genes1,genes2)
  # length(common)
  # 
  # df_cnv <- df_cnv %>%
  #   dplyr::filter(chr_desc %in% common)
  # print(dim(df_cnv))
 
  
  
  df_cnv2 <- df_cnv %>%
    dplyr::filter(LOH_event !='others') %>%
    dplyr::select(chr_desc,clone,cnv) %>%
    dplyr::filter(chr_desc %in% filtered_bins) %>%
    tidyr::pivot_wider(names_from = 'clone', values_from='cnv') %>%
    tidyr::drop_na(obs_clone1[1],obs_clone2[1]) %>%
    as.data.frame()
  # dim(df_cnv2)
  total_LOH <- nrow(df_cnv2)
  print(paste0('Total LOH events: ',total_LOH))
  print(dim(df_cnv))
  df_cnv <- df_cnv %>%
    dplyr::filter(chr_desc %in% unique(df_cnv2$chr_desc))
  
  # CNV plot corrupt tree input data 
  median_cnv_plt <- plot_CNV_corrupt_tree(df_cnv2, c(obs_clone1[1], obs_clone2[1]))
  ascn1 <- df_cnv
  ascn1 <- ascn1 %>% 
    dplyr::select(clone, chr_desc, chr, cnv) %>% #LOH_event,
    tidyr::pivot_wider(names_from = 'clone', values_from='cnv')
  
  ascn2 <- ascn1 %>% 
    dplyr::select(-chr)%>% 
    as.data.frame()
  # print(head(ascn2))
  cond <- abs(ascn2[,obs_clone1[1]]-ascn2[,obs_clone2[1]])>=1
  ascn2 <- ascn2[cond,]
  rownames(ascn2) <- ascn2$chr_desc
  ascn2$chr_desc <- NULL
  rv <- c()
  for(x in rownames(ascn2)){
    if(ascn2[x,obs_clone1[1]]<=2 | ascn2[x,obs_clone2[1]]<=2){
      rv <- c(rv, x)
    }
  }
  nb_LOH_changes <- length(rv)
  print(paste0('Total LOH events: ',nb_LOH_changes))
  # ascn2 <- ascn2 %>%
  #   dplyr::filter(chr_desc %in% rv)
  # dim(ascn2)
  
  if(nb_LOH_changes>0){
    ascn1 <- ascn1 %>% 
      dplyr::filter(chr_desc %in% rv)
    
    
    ascn1$event_type <- ifelse(ascn1[,obs_clone1[1]]<ascn1[,obs_clone2[1]],'GainCN','LossCN')
    total_events <- nrow(ascn1)
    ascn1 <- ascn1 %>%
      dplyr::group_by(event_type, chr) %>%
      dplyr::summarise(nb_events=n()) %>%
      dplyr::mutate(pct_events=round(nb_events/total_events,3)*100)%>%
      ungroup()%>%
      as.data.frame()
    
    plottitle <- paste0('Proportion of LOH altered by CNA events, ',gsub('_',':',obs_clone1[1]),
                        ' vs. ',gsub('_',':',obs_clone2[1]))
    
    subtl <- paste0(total_events,' LOH altered by CNA events in total ',
                    total_LOH,' LOH events (',round(total_events/total_LOH*100,2),'%)')
    proportion_plt <- plot_stack_barplot(ascn1, xlabel=NULL, ylabel='(%) LOH events', subtl=subtl, output_dir=save_dir,  tag='A_B',  plottitle=plottitle,
                                         fa='event_type', xa='event_type', ya='pct_events', lg=T,colorcode=NULL)
    
    
    print("Number of LOH events:")
    print(dim(df_cnv))
    # print(head(df_cnv))
    saveRDS(df_cnv, paste0(save_dir,"cnv.rds"))
    plottitle <- paste0(obs_clone1[2],' vs. ',obs_clone2[2])
    cnv_plt <- plot_CNV(df_cnv, c(obs_clone1[1], obs_clone2[1]), plottitle, save_dir)
    main_plot <- cowplot::plot_grid(proportion_plt, cnv_plt, median_cnv_plt, ncol=1, rel_heights = c(1.5,1,0.85), 
                                    labels = c('a','b','c'), align = 'v')
    ggsave(paste0(save_dir,obs_clone1[1], '_vs_',obs_clone2[1], "_stat_cnchange.png"),
           plot = main_plot,
           height = 4,
           width = 7,
           # useDingbats=F,
           type = "cairo-png",
           dpi=250
    )
  }
  
  # ggsave(
  #   filename = paste0(save_dir,obs_clone1[1], '_vs_',obs_clone2[1], "_stat_cnchange.pdf"),
  #   plot = main_plot,
  #   height = 5,
  #   width = 8,
  #   useDingbats=F)

}  


plot_stack_barplot <- function(cls_df, xlabel, ylabel, subtl, output_dir,  tag='',  plottitle=NULL,
                               fa='event_type', xa='event_type', ya='pct_events', lg=T,colorcode=NULL){
  if(is.null(plottitle)){
    plottitle <- 'Proportion of LOH altered by CNA events'
  }  
  if(is.null(colorcode)){
    colorcode <- structure(
      c(
        "#003b00", "#b1b14e"
      ),
      names=c("GainCN","LossCN")
    )
  }
  
  if(!lg){
    lg_pos <- "none"
  }else{
    lg_pos <- "top"
  }
  # cls_df$pct <- round(cls_df[,ya],2)
  chr_levels <- c(as.character(1:22), "X")
  cls_df$chr <- factor(cls_df$chr, levels =chr_levels)
  cls_df$event_type <- factor(cls_df$event_type,levels=c("LossCN","GainCN"))
  max_pct <- max(cls_df[,ya])
  p <- ggplot2::ggplot(cls_df, ggplot2::aes_string(fill=fa, y=ya, x=xa)) + 
    ggplot2::geom_bar(position="dodge",stat="identity", width = 0.9) +
    # geom_col(width = 0.8, position = position_dodge2(width = 0.8)) + 
    # ggplot2::geom_text(aes_string(label = xa),size=4,  position=position_dodge(width=0.9), vjust = -0.5, angle=90) +  #position=position_dodge(width=0.9), vjust=0
    ggplot2::geom_text(ggplot2::aes_string(label=ya), position=position_dodge(width=0.9), 
                       size=2, color="black",  vjust=-0.1)+ #vjust=-0.3,
    # ggplot2::facet_grid( ~ chr) +#rows = vars(chr)
    ggplot2::facet_wrap(. ~ chr, drop = F, nrow = 1) +  #, scales="free_x", scales="free", space='free'
    # geom_text(aes(y=label_ypos, label=Freq), vjust=1.6, 
    #           color="white", size=3.5) +
    # ylim(yl[1],yl[2]) + 
    scale_y_continuous(breaks = seq(0, max_pct, 10), limits = c(0, max_pct+5)) +  #labels = scales::percent_format(scale = 1)
    ggplot2::scale_fill_manual(values = colorcode) 
  my_font <- "Arial"
  thesis_theme <- theme(strip.background = element_rect(fill = 'white', colour = 'white'),
                        text = element_text(size = 8, hjust = 0.5, family=my_font),
                        plot.title = element_text(color="black",size=9, face="bold", hjust=0, family=my_font),
                        plot.subtitle = element_text(size=8, hjust=0, family=my_font),
                        axis.text.x = element_text(size=5, hjust = 0.5, family=my_font, angle = 90),
                        axis.ticks.x = element_blank(),
                        axis.title.x = element_blank(),
                        axis.text.y = element_text(size=6, hjust = 0.5, family=my_font),
                        axis.title.y = element_text(size=6, hjust = 0.5, family=my_font),
                        axis.line = element_line(colour = "black"),
                        strip.placement = "outside",
                        strip.text.x = element_blank(),
                        legend.position = lg_pos,
                        legend.text=element_text(size=6, hjust = 0.5, family=my_font),
                        legend.title=element_text(size=7, hjust = 0.5, family=my_font),
                        legend.key.size=unit(0.25,"cm"),
                        legend.margin=margin(0,0,0,0),
                        legend.box.margin=margin(-5,-5,-5,-5),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        # panel.background = element_blank(), 
                        panel.background = element_rect(fill = "#F8F8F8", colour = NA),
                        panel.spacing = unit(c(0.1), 'cm'))
  p <- p + thesis_theme   
  
  p <- p + guides(fill = guide_legend(nrow = 1, override.aes = list(size=0.1))) #+  #, override.aes = list(size=1.1)
  p <- p + ggplot2::labs(x=NULL, y=ylabel, title=plottitle, subtitle = subtl, fill='LOH event altered by CNA ') 
  
  # scale_x_continuous(expand = c(0,0)) 
  # if(!is.null(yl)){
  #   p <- p + ylim(yl[1],yl[2])
  # }
  
  # saveRDS(p, paste0(output_dir,tag,"_stat.rds"))
  # png(paste0(output_dir,tag,"_stat_cnchange.png"), height = 2*210, width=2*800, res = 2*72)
  # print(p)
  # dev.off()
  # ggsave(paste0(output_dir,tag,"_stat_cnchange.png"),
  #        plot = p,
  #        height = 2.5,
  #        width = 8.5,
  #        # useDingbats=F,
  #        type = "cairo-png",
  #        dpi=200
  # )
  # ggsave(file=paste0(output_dir,tag,"_stat.pdf"), plot=p, height=3.8, width=6.5, useDingbats=F)  #, dpi=2*72
  return(p)
  
}



plot_CNV <- function(df_cnv, clones, plottitle, save_dir){
  # df_cnv <- read.csv(df_cnv_fn, check.names = F, stringsAsFactors = F)
  # df_cnv <- df_cnv %>% 
  #   dplyr::filter(clone %in% clones)  #, !chr %in% c("X","Y")
  # dim(df_cnv)
  
  df_cnv <- dplyr::select(df_cnv, cnv, clone, chr, start, end, LOH_event) %>%
    group_by(chr) %>%
    dplyr::mutate(start_order = rank(start)) %>%
    ungroup()
  
  chr_levels <- c(as.character(1:22), "X")
  df_cnv <- df_cnv %>% dplyr::filter(chr %in% chr_levels)
  df_cnv$chr <- factor(df_cnv$chr, levels = chr_levels)
  
  cnv_cols <- structure(
    c(
      "#003b00",'#A9A9A9', "#b1b14e"  #"#ff7f7f", "#4c4cff"
    ),
    names=c("ampLOH","cnLOH", "LOH")
  )
  labels_event = c("amplification LOH","neutral LOH", "deletion LOH")
  #levels(cnv_cols) <- 0:11
  #  cnv_cols <- c("0" = "#2166ac",
  #                "1" = "#92c5de", 
  #                "2" = "grey80", 
  #                "3" = "#f4a582", 
  #                "4" = "#d6604d",
  #                "5" = "#b2182b",
  #                "6+" = "#67001f")
  
  # df_cnv$cnv <- as.character(round(df_cnv$cnv))
  # levels(df_cnv$cnv) <- 0:(length(cnv_cols)-1)
  
  update_clones <- lapply(strsplit(as.character(df_cnv$clone), "_"), function(x) {
    return(as.character(x[1]))
  })
  
  if(length(unique(update_clones))!=1){
    df_cnv$clone <- update_clones  
    df_cnv$clone <- factor(df_cnv$clone, levels = unique(df_cnv$clone))
  }else{
    df_cnv$clone <- factor(df_cnv$clone, levels = c(clones[grepl('Primary',clones)],clones[grepl('Metastasis',clones)]))
  }
  print(summary(df_cnv$clone))
  # saveRDS(df_cnv, paste0(save_fig_dir,'df_cnv_1035_test.rds'))
  cnv_plot <- df_cnv %>% 
    ggplot(aes(x = start_order, y = clone, fill = LOH_event)) + #, fill = cnv
    # geom_tile(aes(fill = cnv)) + #
    geom_raster() +
    # geom_tile() + 
    facet_wrap(~ chr, scales = "free_x", nrow = 1,drop = F) + #, switch = "x"
    # facet_wrap(~ chr, nrow = 1) + #scales = "free_x", , switch = "x"
    # facet_grid(~ chr, scales = "free_x", switch = "x", space='free',drop = FALSE) +
    # theme(legend.position = "bottom", axis.text.x = element_blank()) +
    scale_y_discrete(expand = c(0, 0)) +
    # scale_fill_manual(values=cnv_colors, name = "Copy number", guide = 'legend',labels = 0:(length(cnv_colors)-1),drop=FALSE)  +
    #          theme(legend.position = "bottom") + 
    scale_fill_manual(values = cnv_cols, name = "LOH event ", breaks = names(cnv_cols), labels=labels_event) +  #
    labs(x = "Chromosome", y = "Clone", title = plottitle)
  
  lg_pos <- "bottom"
  my_font <- "Arial"
  thesis_theme <- theme(strip.background = element_rect(fill = 'white', colour = 'white'),
                        text = element_text(size = 8, hjust = 0.5, family=my_font),
                        plot.title = element_text(color="black",size=9, face="bold", hjust=0, family=my_font),
                        plot.subtitle = element_text(size=8, hjust=0, family=my_font),
                        axis.text.x = element_blank(),
                        # axis.text.x = element_text(size=7, hjust = 0.5, family=my_font, angle = 90),
                        axis.ticks.x = element_blank(),
                        axis.title.x = element_text(color="black", size=7, hjust = 0.5, family=my_font),
                        axis.text.y = element_text(size=7, hjust = 0.5, family=my_font),
                        axis.title.y = element_text(color="black", size=6, hjust = 0.5, family=my_font),
                        axis.line = element_line(colour = "black"),
                        strip.placement = "outside",
                        # strip.text.x = element_blank(),
                        strip.text.x = element_text(size = 8, hjust = 0.5, family=my_font),
                        legend.position = lg_pos,
                        legend.text=element_text(size=6, hjust = 0.5, family=my_font),
                        legend.title=element_text(size=7, hjust = 0.5, family=my_font),
                        legend.key.size=unit(0.2,"cm"),
                        legend.margin=margin(0,0,0,0),
                        legend.box.margin=margin(-5,-5,-5,-5),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        # panel.background = element_blank(),
                        panel.background = element_rect(fill = "#F8F8F8", colour = NA),
                        panel.spacing = unit(c(0.1), 'cm'))
  cnv_plot <- cnv_plot + thesis_theme   
  cnv_plot <- cnv_plot + guides(fill = guide_legend(nrow = 1, override.aes = list(size=0.02))) +  #, override.aes = list(size=1.1),
    scale_x_continuous(expand = c(0,0)) 
  
  # desc <- paste(clones, collapse = '_vs_')
  # saveRDS(cnv_plot, paste0(save_dir, 'LOH_profile_cnvplot.rds'))
  # if(nrow(df_cnv)<100){
  #   wd <- 600
  #   w <- 5
  # }else{
  #   wd <- 1000
  #   w <- 10
  # }
  # wd <- 1100
  # png(paste0(save_dir, desc,'_LOH_profile.png'), height = 2*250, width=2*wd, res = 2*72)
  # print(cnv_plot)
  # dev.off()
  
  # pdf(paste0(save_dir, 'LOH_profile.pdf'), height = 3.5, width=w, useDingbats = F)
  # print(cnv_plot)
  # dev.off()
  # BiocManager::install("Cairo")
  # library(Cairo)
  # ggsave(
  #   filename = paste0(save_dir, 'LOH_profile_cairo.png'),
  #   plot = cnv_plot,
  #   height = 3,
  #   width = w,
  #   # useDingbats=F
  #   type = "cairo-png")
  # results <- list(df_cnv=df_cnv, cnv_plot=cnv_plot)
  return(cnv_plot)
}




viz_boxplot_LOH <- function(ascn_stat, datatag, save_dir, 
                            xstring="clone_id", ystring="percent_LOH", 
                            group_string="clone_id",
                            lg_title="LOH per cell (%) ",
                            plttitle=paste0("Distribution of LOH in ",datatag),
                            ytitle="Percentage of LOH (%)",wd=580,ht=250){
  p <- ggplot(ascn_stat, aes_string(x = xstring, y = ystring, color=group_string)) +
    geom_boxplot() +
    scale_y_continuous(breaks = seq(round(min(ascn_stat[,ystring]))-5, round(max(ascn_stat[,ystring]))+5, by = 5)) + 
    theme(strip.background = element_rect(fill = 'white', colour = 'white'),
          axis.text.x = element_text(size=14, colour = "black"),
          # axis.ticks.x = element_line(size=8, colour = "black"),
          axis.text.y = element_text(size=14, colour = "black"),
          axis.title.y = element_text(size=14, colour = "black"),
          axis.title.x = element_text(size=14, colour = "black"),
          legend.position = "right",
          legend.title = element_text(size=14, colour = "black"),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = "#F8F8F8", colour = NA),
          panel.spacing = unit(0.01, 'cm')) +
    labs(x = "Clones", y= ytitle, title=plttitle) #+ 
  # guides(fill = guide_legend(override.aes = list(size=0.5)))
  
  
  png(paste0(save_dir,datatag,".png"), height = 2*ht, width=2*wd, res = 2*72)
  print(p)
  dev.off()
  
}  

plot_BAF_per_clone <- function(ascn, save_dir,datatag){
  plot_ls <- list()
  for(cl in unique(ascn$clone_id)){
    ascn_tmp <- ascn[ascn$clone_id==cl,]
    p <- plotBAFperstate(paste0(datatag,' BAF per state in clone ',cl), ascn_tmp, maxstate = 10)
    plot_ls[[cl]] <- p
    # png(paste0(save_dir,datatag,'_',cl,"_BAFperstate.png"), height = 2*450, width=2*900,res = 2*72)
    # print(p)
    # dev.off()
  }
  return(plot_ls)
}
viz_LOH <- function(ascn_stat, datatag, save_dir, 
                    xstring="cell_id", ystring="percent_LOH", 
                    group_string="percent_LOH",
                    lg_title="LOH per cell (%) ",
                    plttitle=paste0("Distribution of LOH in ",datatag),
                    ytitle="Percentage of LOH (%)"){
  ploh <- ggplot(ascn_stat, aes_string(x = xstring, y = ystring)) +
    geom_point(aes_string(colour = group_string), size = 1.8) +   # MA: size was 1
    facet_grid(~ clone_id, scales = "free_x", space='free', switch = "x") + #facet_wrap
    scale_color_continuous(name = lg_title,limits = c(min(ascn_stat[,ystring]), max(ascn_stat[,ystring])),
                           type = "viridis", guide=guide_colorbar(reverse=TRUE)) +
    theme(#strip.background = element_rect(fill = 'white', colour = 'white'),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size=14, colour = "black"),
      axis.title.y = element_text(size=14, colour = "black"),
      axis.title.x = element_text(size=14, colour = "black"),
      strip.placement = "outside",
      legend.position = "top",
      axis.line = element_line(colour = "black"),
      panel.spacing = unit(0.1, 'cm')) +
    labs(x = "Clones", y= ytitle, title=plttitle) + 
    guides(fill = guide_legend(col = 1, override.aes = list(size=1.1)))
  
  png(paste0(save_dir,datatag,"_LOH_plots.png"), height = 2*350, width=2*900, res = 2*72)
  print(ploh)
  dev.off()
  
}

# View(head(ascn_stat))

# dim(ascn_stat) 
# View(head(ascn_stat))
# summary(ascn_stat$percent_LOH)
# 
# ascn <- ascn %>%
#   group_by(chr, start, end) %>%
#   summarise(state_phase_stat = schnapps:::Mode(state_phase),
#             BAF_stat = mean(BAF, na.rm = T),
#             state_BAF_stat = mean(state_BAF, na.rm = T)) %>%
#   ungroup() 
#   

# print("Summary number of LOH events for whole series data")
# summary(ascn$state_phase_stat)
# 
# 
# print("Summary whole genome duplication WGD")
# summary(ascn$state_BAF_stat)
# 
# nb_WGD = sum(ascn$state_BAF_stat>4)  


# write.csv(ascn, file=paste0(results_pseudobk_dir,'ascn_events.csv'), quote=F, row.names = F)


plotBAFperstate <- function(plttitle, cn, minpts = 250, minfrac = 0.01, maxstate = 10, dens_adjust = 2) 
{
  if (schnapps::is.hscn(cn) | schnapps::is.ascn(cn)) {
    alleleCN <- cn$data
  }
  else {
    alleleCN <- cn
  }
  maj <- seq(0, max(alleleCN$state), 1)
  min <- seq(0, max(alleleCN$state), 1)
  allASstates <- expand.grid(state = maj, min = min) %>% dplyr::mutate(cBAF = min/state) %>% 
    dplyr::mutate(state_AS_phased = paste0(state - min, 
                                           "|", min)) %>% dplyr::mutate(Maj = state - min, 
                                                                        Min = min) %>% dplyr::select(-min) %>% dplyr::filter(Maj >= 
                                                                                                                               0, Min >= 0) %>% dplyr::filter(state <= maxstate)
  allASstates$cBAF[is.nan(allASstates$cBAF)] <- 0
  forplot <- alleleCN %>% dplyr::group_by(state_AS_phased) %>% 
    dplyr::mutate(n = dplyr::n()) %>% dplyr::ungroup() %>% 
    dplyr::mutate(f = n/dplyr::n()) %>% dplyr::filter(state <= 
                                                        maxstate, n > minpts, f > minfrac)
  allASstates <- allASstates %>% dplyr::filter(state_AS_phased %in% 
                                                 unique(forplot$state_AS_phased))
  text_fraction <- dplyr::distinct(forplot, state_AS_phased, 
                                   f) %>% dplyr::mutate(pct = paste0(100 * round(f, 2), 
                                                                     "%"), y = 0.95)
  g <- forplot %>% dplyr::mutate(cncol = paste0("CN", state)) %>% 
    ggplot2::ggplot(ggplot2::aes(x = forcats::fct_reorder(state_AS_phased, 
                                                          state), y = BAF)) + ggplot2::geom_violin(scale = "width", 
                                                                                                   col = "white", ggplot2::aes(fill = cncol), 
                                                                                                   adjust = dens_adjust) + 
    ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA, 
                          col = "white", ggplot2::aes(fill = cncol)) + ggplot2::scale_fill_manual(name = "Copy number \n state", 
                                                                                                  breaks = paste0("CN", seq(0, max(alleleCN$state, na.rm = TRUE), 
                                                                                                                            1)), labels = seq(0, max(alleleCN$state, na.rm = TRUE), 
                                                                                                                                              1), values = schnapps::scCN_cols(paste0("CN", seq(0, max(alleleCN$state, 
                                                                                                                                                                                                       na.rm = TRUE), 1)))) + cowplot::theme_cowplot() + 
    ggplot2::geom_crossbar(data = allASstates, ggplot2::aes(y = cBAF, 
                                                            ymin = cBAF, ymax = cBAF), alpha = 0.2, size = 0.2) + 
    ggplot2::geom_text(data = text_fraction, ggplot2::aes(x = state_AS_phased, 
                                                          y = y, label = pct)) + ggplot2::xlab("") + ggplot2::theme(legend.position = "bottom") + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, 
                                                       hjust = 1),
                   plot.title = element_text(color="black", size=10, hjust = 0.5),
                   legend.title = element_text(color="black", size=7, hjust = 0.5)) +
    labs(title=plttitle) + 
    guides(fill = guide_legend(nrow = 1, override.aes = list(size=1.1)))
  return(g)
}
