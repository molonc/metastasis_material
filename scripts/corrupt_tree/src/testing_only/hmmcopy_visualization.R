
get_library_id <- function(cells_name){
  lib_ids <- strsplit(as.character(cells_name),"-")
  # length(lib_ids)
  nbcores <- 5
  ps <- c()
  lids <- parallel::mclapply(lib_ids, function(f) {
    ps <- c(ps, as.character(f[2]))
  }, mc.cores = nbcores)
  
  if(length(lids)==length(cells_name)){
    return(as.character(lids))
  } else{
    return(NULL)
  }
  
}

get_reads_dir <- function(input_dir, lib_id){
  # lib_dirs <- c()
  fn1 <- paste0(input_dir,lib_id,'/','hmmcopy_autoploidy/',lib_id,'_multiplier0_reads.csv')
  fn2 <- paste0(input_dir,lib_id,'/','hmmcopy_autoploidy/',lib_id,'_reads.csv')
  fn3 <- paste0(input_dir,lib_id,'/','hmmcopy/',lib_id,'_reads.csv')
  if(file.exists(fn1)){
    # lib_dirs <- c(lib_dirs, fn1)
    ldir <- fn1
  } else if(file.exists(fn2)){
    # lib_dirs <- c(lib_dirs, fn2)
    ldir <- fn2
  } else if(file.exists(fn3)){
    # lib_dirs <- c(lib_dirs, fn3)
    ldir <- fn3
  } else{
    print(paste0('Double check lib:',lib_id))
    ldir <- NULL
  }
  return(ldir)
}

# t <- tmp[tmp$state==1,]
# dim(t)
# View(t)
# t <- t[,c('chr','start','end')]
# write.csv(t,file = paste0(save_dir,'chr3_CN1_position.csv'),quote = F)
plot_segment_profile <- function(obs_cell, reads, cn_colours, datatag, loci_use, save_dir,
                                 xlabel='chr'){
  
  tmp <- subset(reads, cell_id == obs_cell)
  tmp <- subset(tmp, chr == "3")
  # tmp$chr_label <- paste0(tmp$chr,'_',tmp$start,'_',tmp$end)
  # tmp <- tmp[tmp$chr_label %in% loci_use,]
  # tmp$state <- ifelse(tmp$state>11,11,tmp$state)
  tmp$state[tmp$state>11] <- 11
  the_mat$chr[the_mat$chr == 'X'] <- '40'
  lim <- max(quantile(tmp$copy, na.rm = TRUE, 0.99), 4)
  # tmp <- tmp %>% group_by(chr, start, end)
  # tmp$chr_val <- ifelse(tmp$chr=="X",23,
  #                       ifelse(tmp$chr=="Y",24,
  #                              ifelse(!tmp$chr %in% c("X","Y"),tmp$chr,-1)))
  # tmp$chr_val <- as.numeric(tmp$chr_val)
  # tmp <- dplyr::arrange(tmp, chr_val)
  # tmp$start <- factor(tmp$start)
  
  # order_chr <- data.frame(chr=c(as.character(rep(1:22,1)),"X","Y"),
  #                         chr_val=rep(1:24,1))
  # tmp$chr <- factor(tmp$chr, levels = order_chr$chr[order(order_chr$chr_val)])
  
  # p <- ggplot(tmp, aes(start, copy * 2, col = as.factor(state))) + geom_point(size = 0.5) + 
  #   facet_grid(cell_id ~ chr, space = "free_x", switch = "x") + 
  #   scale_x_continuous(expand = c(0, 0), breaks = NULL) +  ylim(0, lim*2) + #+ theme(limit = c(0, lim))
  #   scale_color_manual(values = cn_colours)
  
  p <- ggplot(tmp, aes(start, copy, col = as.factor(state))) + geom_point(size = 3, shape=16) + 
    facet_grid(cell_id ~ chr, scales = "free", space = "free", switch = "x") + 
    scale_x_continuous(expand = c(0, 0), breaks = NULL) +  ylim(0, lim) + #+ theme(limit = c(0, lim))
    scale_color_manual(values = cn_colours)
  
  p <- p + labs(x=xlabel)
  p <- p + theme(panel.spacing = unit(0.1, "lines"),
                 legend.title = element_text(color="black", size=9,hjust = 0.5),
                 axis.text.x = element_text(color="black", size=9, hjust = 0.5),
                 axis.text.y = element_text(color="black", size=9, hjust = 0.5))
  
  
  
  # png(paste0(save_dir,datatag,'_',obs_cell,".png"), height = 2*200, width=2*550, res = 2*72)
  # print(p)
  # dev.off()
  return(p)
}  


library(ggplot2)
library(plyr)
# library(RColorBrewer)

results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_filtered/'
input_dir <- '/home/htran/storage/datasets/hakwoo_metastasis/'
# cell_clones <- paste0(results_dir,'tree_cut_out/cell_clones_if_0.02_af_0.75_p0.75_e0.04.csv')
cell_clones <- paste0(results_dir,'cell_clones.csv')
series_tag <- 'SA919'
cellclones <- read.csv(cell_clones, check.names = F, stringsAsFactors = FALSE)
save_dir <- paste0(results_dir,'cell_features/')
if (!file.exists(save_dir)){
  dir.create(save_dir)
}

# copynumber <- paste0(results_dir, 'total_merged_filtered_states.csv')
# copy_number <- read.csv(copynumber, header=T, row.names=1, check.names = F,stringsAsFactors = FALSE)
# loci_use <- rownames(copy_number)



cn_colours <- structure(
  c(
    "#3182BD", "#9ECAE1", "#CCCCCC", "#FDCC8A", "#FC8D59", "#E34A33",
    "#B30000", "#980043", "#DD1C77", "#DF65B0", "#C994C7", "#D4B9DA"
  ),
  names = 0:11
  # names=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11+")
)


obs_cells <- c()
clones_ls <- c()
nbsample <- 4
for(c in unique(cellclones$clone_id)){
  cn <- cellclones[cellclones$clone_id==c,'cell_id']
  scn <- cn[sample(1:length(cn), nbsample, replace=FALSE)]
  obs_cells <- c(obs_cells, scn)
  clones_ls <- c(clones_ls, rep(c, nbsample))
}
meta_reads <- data.frame(cell_id=obs_cells, 
                         library_id = get_library_id(obs_cells), 
                         clone_id = clones_ls,
                         row.names = obs_cells)
# View(meta_reads)
write.csv(meta_reads, file = paste0(save_dir,'meta_reads.csv'), quote = F, row.names = T)

plots <- list()
for(lib_id in unique(meta_reads$library_id)){
  
  reads_fn <- get_reads_dir(input_dir, lib_id)
  reads <- data.table::fread(reads_fn, check.names = F, stringsAsFactors = FALSE)
  
  for(obs_cell in meta_reads[meta_reads$library_id==lib_id,'cell_id']){
    print(obs_cell)
    datatag <- meta_reads[obs_cell,'clone_id']
    p <- plot_segment_profile(obs_cell, reads, cn_colours, datatag, loci_use, save_dir, 'chr')
    plots[[obs_cell]] <- p
  }
  
}

for(cid in unique(meta_reads$clone_id)){
  cells_use <- meta_reads[meta_reads$clone_id==cid,'cell_id']
  ps <- list()
  for(c in cells_use){
    ps[[c]] <- plots[[c]]
  }
  pl <- cowplot::plot_grid(plotlist = ps, align = "v", ncol = 4, rel_heights = c(1, 1.3))
  png(paste0(save_dir,cid,'_diagnostic.png'), height = 2*500, width=2*1200, res = 2*72)
  print(pl)
  dev.off()
}  

saveRDS(plots,file = paste0(save_dir,'plots.rds'))



