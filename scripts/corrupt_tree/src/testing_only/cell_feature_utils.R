library(ggExtra)
library(ggplot2)

# heatmap_if_0.02_af_0.75_p0.75_e0.04.pdf
plot_scatter_function <- function(meta_data, xstring, ystring, plottype, 
                          plottitle="Metrics", 
                          xlabel='read', ylabel='mapped_read', 
                          save_dir, save_fn, colorcode="blue",xl=NULL, yl=NULL) {
  
  p <- ggplot(meta_data, aes_string(x=xstring, y=ystring)) +  #, group=plottype
    # geom_line(aes_string(color=plottype)) #+   #,color=colorcode
    geom_point(color=colorcode) #+   #aes_string(color=plottype)
  # scale_colour_manual(values = colorcode, guide = FALSE) 
  p <- p + labs(x=xlabel,y=ylabel,title=plottitle)
  if(!is.null(xl)){
    p <- p + xlim(xl[1],xl[2])
  }
  if(!is.null(yl)){
    p <- p + ylim(yl[1],yl[2])
  }
  
  p <- p + theme(legend.title = element_blank(), 
                 plot.title = element_text(color="black", size=8, hjust = 0.5),
                 legend.position = "none", 
                 # axis.line = element_blank(), 
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                 panel.border = element_blank(),
                 axis.text.y = element_text(color="black",size=5,angle = 30, hjust = 1),
                 axis.text.x = element_text(color="black",size=5, hjust = 1),
                 # axis.text = element_blank(),
                 #axis.title = element_blank(),
                 # axis.ticks = element_blank()
  )
  p1 <- ggMarginal(p, type="histogram",fill = "white")  #, size=10
  # png(paste0(save_dir,save_fn,".png"), height = 2*500, width=2*500,res = 2*72)
  # print(p1)
  # dev.off()
  return(p1)
  
}  

plot_scatter_function_plottype <- function(meta_data, xstring, ystring, plottype, 
                                  plottitle="Metrics", 
                                  xlabel='read', ylabel='mapped_read', 
                                  save_dir, save_fn, colorcode="blue",
                                  xl=NULL, yl=NULL) {
  
  p <- ggplot(meta_data, aes_string(x=xstring, y=ystring)) +  #, group=plottype
    # geom_line(aes_string(linetype=plottype, color=plottype))+   #,color=colorcode
    geom_point(aes_string(color=plottype),size=0.8) #+   #aes_string(color=plottype)
  # scale_colour_manual(values = colorcode, guide = FALSE) 
  p <- p + labs(x=xlabel,y=ylabel,title=plottitle)
  if(!is.null(xl)){
    p <- p + xlim(xl[1],xl[2])
  }
  if(!is.null(yl)){
    p <- p + ylim(yl[1],yl[2])
  }
  p <- p + theme(legend.title = element_text(color="black", size=8, hjust = 0.5),
                 plot.title = element_text(color="black", size=8, hjust = 0.5),
                 # legend.position = "none", 
                 # axis.line = element_blank(), 
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                 panel.border = element_blank(),
                 axis.text.y = element_text(color="black",size=5,angle = 30, hjust = 1),
                 axis.text.x = element_text(color="black",size=5, hjust = 1),
                 axis.text = element_text(color="black",size=7, hjust = 1),
                 legend.text = element_text(color="black",size=7, hjust = 1)
                 #axis.title = element_blank(),
                 # axis.ticks = element_blank()
  )
  p1 <- ggMarginal(p, type="histogram",fill = "white")  #, size=10, ,margins = 'x'
  png(paste0(save_dir,save_fn,".png"), height = 2*400, width=2*550,res = 2*72)
  print(p1)
  dev.off()
  return(p1)
  
}  

plot_scatter_function_contaminate <- function(meta_data, xstring, ystring, plottype, 
                                           plottitle="Metrics", 
                                           xlabel='read', ylabel='mapped_read', 
                                           save_dir, save_fn, colorcode="blue",
                                           xl=NULL, yl=NULL) {
  
  p <- ggplot(meta_data, aes_string(x=xstring, y=ystring)) +  #, group=plottype
    # geom_line(aes_string(linetype=plottype, color=plottype))+   #,color=colorcode
    geom_point(aes_string(color=plottype),size=0.8) #+   #aes_string(color=plottype)
    # scale_colour_manual(breaks = c("False", "Uncheck"), c("False", "True", "Uncheck")
    #                     values=c("green","blue"))  c("green","red", "blue")
  p <- p + labs(x=xlabel,y=ylabel,title=plottitle)
  if(!is.null(xl)){
    p <- p + xlim(xl[1],xl[2])
  }
  if(!is.null(yl)){
    p <- p + ylim(yl[1],yl[2])
  }
  p <- p + theme(legend.title = element_text(color="black", size=6, hjust = 0.5),
                 plot.title = element_text(color="black", size=8, hjust = 0.5),
                 # legend.position = "none", 
                 # axis.line = element_blank(), 
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                 panel.border = element_blank(),
                 axis.text.y = element_text(color="black",size=5,angle = 30, hjust = 1),
                 axis.text.x = element_text(color="black",size=5, hjust = 1),
                 axis.text = element_text(color="black",size=5, hjust = 1),
                 legend.text = element_text(color="black",size=5, hjust = 1)
                 #axis.title = element_blank(),
                 # axis.ticks = element_blank()
  )
  p1 <- ggMarginal(p, type="histogram",fill = "white",margins = 'y')  #, size=10
  # png(paste0(save_dir,save_fn,".png"), height = 2*500, width=2*500,res = 2*72)
  # print(p1)
  # dev.off()
  return(p1)
  
}  
plot_scatter_function_plottype_v2 <- function(meta_data, xstring, ystring, plottype,plottype2, 
                                           plottitle="Metrics", 
                                           xlabel='read', ylabel='mapped_read', 
                                           save_dir, save_fn, colorcode="blue",
                                           xl=NULL, yl=NULL) {
  
  p <- ggplot(meta_data, aes_string(x=xstring, y=ystring)) +  #, group=plottype
    # geom_line(aes_string(linetype=plottype, color=plottype))+   #,color=colorcode
    geom_point(aes_string(color=plottype, shape=plottype2),size=0.8) #+   #aes_string(color=plottype)
  # scale_colour_manual(values = colorcode, guide = FALSE) 
  p <- p + labs(x=xlabel,y=ylabel,title=plottitle)
  if(!is.null(xl)){
    p <- p + xlim(xl[1],xl[2])
  }
  if(!is.null(yl)){
    p <- p + ylim(yl[1],yl[2])
  }
  p <- p + theme(legend.title = element_text(color="black", size=8, hjust = 0.5),
                 plot.title = element_text(color="black", size=8, hjust = 0.5),
                 # legend.position = "none", 
                 # axis.line = element_blank(), 
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                 panel.border = element_blank(),
                 axis.text.y = element_text(color="black",size=5,angle = 30, hjust = 1),
                 axis.text.x = element_text(color="black",size=5, hjust = 1),
                 axis.text = element_text(color="black",size=7, hjust = 1),
                 legend.text = element_text(color="black",size=7, hjust = 1)
                 #axis.title = element_blank(),
                 # axis.ticks = element_blank()
  )
  p1 <- ggMarginal(p, type="histogram",fill = "white",margins = 'x')  #, size=10
  # png(paste0(save_dir,save_fn,".png"), height = 2*500, width=2*500,res = 2*72)
  # print(p1)
  # dev.off()
  return(p1)
  
}  

plot_line_function <- function(meta_data, xstring, ystring, plottype, 
                                  plottitle="Metrics", 
                                  xlabel='read', ylabel='mapped_read', 
                                  save_dir, save_fn, colorcode="blue") {
  
  p <- ggplot(meta_data, aes_string(x=xstring, y=ystring)) +  #, group=plottype
    geom_line(aes_string(linetype=plottype, color=plottype)) #+   #,color=colorcode
    # geom_point(color=colorcode) #+   #aes_string(color=plottype)
  # scale_colour_manual(values = colorcode, guide = FALSE) 
  p <- p + labs(x=xlabel,y=ylabel,title=plottitle)
  p <- p + theme(legend.title = element_text(color="black", size=8, hjust = 0.5), 
                 plot.title = element_text(color="black", size=8, hjust = 0.5),
                 # legend.position = "none", 
                 # axis.line = element_blank(), 
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                 panel.border = element_blank(),
                 axis.text.y = element_text(color="black",size=5,angle = 30, hjust = 1),
                 axis.text.x = element_text(color="black",size=5, hjust = 1),
                 # axis.text = element_blank(),
                 #axis.title = element_blank(),
                 # axis.ticks = element_blank()
  )
  # p1 <- ggMarginal(p, type="histogram",fill = "white")  #, size=10
  # png(paste0(save_dir,save_fn,".png"), height = 2*500, width=2*500,res = 2*72)
  # print(p)
  # dev.off()
  return(p)
  
}  

plot_boxplot_function <- function(meta_data, xstring, ystring, plottype, 
                               plottitle="Metrics", 
                               xlabel='read', ylabel='mapped_read', 
                               save_dir, save_fn, colorcode="blue",
                               xl=NULL, yl=NULL) {
  
  p <- ggplot(meta_data, aes_string(x=xstring, y=ystring)) +  #, group=plottype
    geom_violin(aes_string(color=plottype)) +   #,color=colorcode
    geom_point(aes_string(color=plottype)) #+   #aes_string(color=plottype)
  # scale_colour_manual(values = colorcode, guide = FALSE) 
  p <- p + labs(x=xlabel,y=ylabel,title=plottitle)
  if(!is.null(xl)){
    p <- p + xlim(xl[1],xl[2])
  }
  if(!is.null(yl)){
    p <- p + ylim(yl[1],yl[2])
  }
  p <- p + theme(legend.title = element_text(color="black", size=7, hjust = 0.5), 
                 plot.title = element_text(color="black", size=8, hjust = 0.5),
                 # legend.position = "none", 
                 # axis.line = element_blank(), 
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                 panel.border = element_blank(),
                 axis.text.y = element_text(color="black",size=5,angle = 30, hjust = 1),
                 axis.text.x = element_text(color="black",size=5, hjust = 1),
                 # axis.text = element_blank(),
                 #axis.title = element_blank(),
                 # axis.ticks = element_blank()
  )
  # p1 <- ggMarginal(p, type="histogram",fill = "white", margins = c("x"))  #, size=10
  # png(paste0(save_dir,save_fn,".png"), height = 2*500, width=2*500,res = 2*72)
  # print(p)
  # dev.off()
  return(p)
  
}  

# plts <- plot_outlier_cloneA(predict_state_tmp, tr, tmr, tunmr, tmedian, tmean,
#                             series_tag='SA1035', clone_id=clone_id, save_dir)

plot_outlier_cloneA <- function(predict_state,  tr, tmr, tunmr, tmedian, tmean,
                                series_tag='SA919', 
                                clone_id='A', save_dir='',save_data=F){
  if (!file.exists(save_dir)){
    dir.create(save_dir)
  }
  xstring <- 'total_reads'
  ystring <- 'total_mapped_reads'
  plottype <- 'cell_classified'  # cell clone or sth else
  save_fn <- paste0(clone_id,'_mapped_reads')
  pmap <- plot_scatter_function_plottype(predict_state, xstring, ystring, plottype, 
                              plottitle=paste0("Metrics ",series_tag,' Cl_ ',clone_id), 
                              xlabel=paste0('Total reads',' (n=',dim(predict_state)[1],')'), ylabel='Total mapped reads', 
                              save_dir, save_fn, colorcode="#2F4F4F", xl=tr, yl=tmr)
  
  
  
  xstring <- 'total_reads'
  ystring <- 'unmapped_reads'
  plottype <- 'cell_classified' 
  save_fn <- paste0(clone_id,'_unmapped_reads')
  punmap <- plot_scatter_function_plottype(predict_state, xstring, ystring, plottype, 
                              plottitle=paste0("Metrics ",series_tag,' Cl_ ',clone_id), 
                              xlabel=paste0('Total reads',' (n=',dim(predict_state)[1],')'), ylabel='Total Unmapped Reads', 
                              save_dir, save_fn, colorcode="#330000", xl=tr, yl=tunmr)
  
  
  xstring <- 'median_hmmcopy_reads_per_bin'
  ystring <- 'quality'
  save_fn <- paste0(clone_id,'_quality_median_read')
  predict_state1 <- predict_state
  plottype <- 'cell_classified'  # cell clone or sth else
  colnames(predict_state1)[which(names(predict_state1) == "experimental_condition")] <- "exp"
  pquality <- plot_scatter_function_plottype(predict_state1, xstring, ystring, plottype, 
                           plottitle=paste0("Metrics ",series_tag,' Cl_ ',clone_id), 
                           xlabel=paste0('Median cn reads per bin',' (n=',dim(predict_state)[1],')'), ylabel='Quality', 
                           save_dir, save_fn, colorcode="#330000", xl=tmedian, yl=c(0.75,1))
  
  
  xstring <- 'exp'
  ystring <- 'quality'
  save_fn <- paste0(clone_id,'_quality_exp')
  plottype <- 'exp' 
  pexp <- plot_boxplot_function(predict_state1, xstring, ystring, plottype, 
                              plottitle=paste0("Metrics ",series_tag,' Cl_ ',clone_id), 
                              xlabel=paste0('Experimental Condition',' (n=',dim(predict_state)[1],')'), 
                              ylabel='Quality', 
                              save_dir, save_fn, colorcode="#330000",xl=NULL, yl=c(0.75,1))
  
  
  # xstring <- 'median_hmmcopy_reads_per_bin'
  # ystring <- 'is_s_phase_prob'
  # plottype <- 'is_contaminated' 
  # save_fn <- paste0(clone_id,'_sphase_contaminate')
  # p5 <- plot_scatter_function(predict_state, xstring, ystring, plottype, 
  #                             plottitle=paste0("Predict Cell State ",series_tag), 
  #                             xlabel=paste0('Median hmmcopy reads per bin ',' (n=',dim(predict_state)[1],')'), ylabel='Is s-phase prob', 
  #                             save_dir, save_fn, colorcode="#53868B")
  # 
  
  xstring <- 'is_cont'
  ystring <- 'is_s_phase_prob'
  plottype <- 'cell_classified'  
  save_fn <- paste0(clone_id,'_sphase_contaminate')
  
  predict_state1$is_contaminated <- as.factor(predict_state1$is_contaminated)
  cont <- unique(predict_state1$is_contaminated)
  colnames(predict_state1)[which(names(predict_state1) == "is_contaminated")] <- "is_cont"
  
  # color_df <- data.frame(bk=c("False", "True", "Uncheck"),val=c("green","red", "blue"), 
  #                        row.names = c("False", "True", "Uncheck"))
  # color_df <- color_df[cont,]
  # color_df <- droplevels(color_df)
  psphase <- plot_scatter_function_contaminate(predict_state1, xstring, ystring, plottype, 
                                        plottitle=paste0("Predicted",series_tag,' Cl_ ',clone_id), 
                                        xlabel=paste0('Is-contaminated',' (n=',dim(predict_state)[1],')'), ylabel='Is s-phase prob', 
                                        save_dir, save_fn, colorcode="#53868B",xl=NULL, yl=c(0,0.5))
  
  xstring <- 'mean_copy'  #mean_hmmcopy_reads_per_bin
  ystring <- 'quality'
  save_fn <- paste0(clone_id,'_quality_mean_read')
  plottype <- 'cell_classified' 
  # xm <- c(1.5,2.45)
  if(clone_id!='None'){
    # xm <- c(1.5,2.45)
    xm <- tmean
  } else{
    xm <- NULL
  }

  pquality_mean <- plot_scatter_function_plottype(predict_state1, xstring, ystring, plottype, 
                                                     plottitle=paste0("Metrics ",series_tag,' Cl_ ',clone_id), 
                                                     xlabel=paste0('Mean CN',' (n=',dim(predict_state)[1],')'), ylabel='Quality', 
                                                     save_dir, save_fn, colorcode="#330000",xl=xm, yl=c(0.75,1))
  
  plts <- list(pmap=pmap, punmap=punmap, 
               pquality=pquality, pexp=pexp, 
               psphase=psphase,pquality_mean=pquality_mean)
  if(save_data){
    saveRDS(plts, file = paste0(save_dir,series_tag,'_clone',clone_id,'_plots_feature.rds'))
  }
  return(plts)
}




plot_outlier_cells_features <- function(predict_state, tag='', series_tag='SA919', 
                                        clone_id='cq', save_dir='', save_data=F, plotOutput=T){
  if (!file.exists(save_dir)){
    dir.create(save_dir)
  }
  xstring <- 'mean_copy'
  ystring <- 'is_s_phase_prob'
  plottype <- 'is_cont'  
  save_fn <- paste0(clone_id,'_sphase_contaminate')
  predict_state1 <- predict_state
  predict_state1$is_contaminated <- as.factor(predict_state1$is_contaminated)
  # cont <- unique(predict_state1$is_contaminated)
  colnames(predict_state1)[which(names(predict_state1) == "is_contaminated")] <- "is_cont"
  
  # color_df <- data.frame(bk=c("False", "True", "Uncheck"),val=c("green","red", "blue"), 
  #                        row.names = c("False", "True", "Uncheck"))
  # color_df <- color_df[cont,]
  # color_df <- droplevels(color_df)
  save_fn <- 'contaminated_outlier'
  psphase <- plot_scatter_function_contaminate(predict_state1, xstring, ystring, plottype, 
                                               plottitle=paste0("Predicted",series_tag,' Cl_ ',clone_id), 
                                               xlabel=paste0('Mean copy (',tag,')'), ylabel='Is s-phase prob', 
                                               save_dir, save_fn, colorcode="#53868B")
  psphase
  
  xstring <- 'total_reads'
  ystring <- 'total_mapped_reads'
  plottype <- 'cell_classified'  # cell clone or sth else
  save_fn <- paste0(clone_id,'_mapped_reads')
  pmap <- plot_scatter_function_plottype(predict_state, xstring, ystring, plottype, 
                                plottitle=paste0("Cells Metrics ",series_tag), 
                                xlabel=paste0('Total reads (',tag,')'), ylabel='Total Mapped Reads', 
                                save_dir, save_fn, colorcode="#2F4F4F", xl=c(0,7000000), yl=c(0,5500000))
  
  # pmap
  
  xstring <- 'total_reads'
  ystring <- 'unmapped_reads'
  plottype <- 'cell_classified' 
  save_fn <- paste0(clone_id,'_unmapped_reads')
  punmap <- plot_scatter_function_plottype(predict_state, xstring, ystring, plottype, 
                                         plottitle=paste0("Cells Metrics ",series_tag), 
                                         xlabel=paste0('Total reads (',tag,')'), ylabel='Total Unmapped Reads', 
                                         save_dir, save_fn, colorcode="#330000", xl=c(0,7000000), yl=c(0,3200000))
  
  # punmap
  
  xstring <- 'median_hmmcopy_reads_per_bin'
  ystring <- 'quality'
  save_fn <- paste0(clone_id,'_quality_median_read')
  predict_state1 <- predict_state
  plottype <- 'cell_classified' 
  plottype2 <- 'exp_cond' 
  colnames(predict_state1)[which(names(predict_state1) == "experimental_condition")] <- "exp_cond"
  pquality_median <- plot_scatter_function_plottype_v2(predict_state1, xstring, ystring, plottype, plottype2,
                                             plottitle=paste0("Cells Metrics ",series_tag), 
                                             xlabel=paste0('Median cn reads per bin (',tag,')'), ylabel='Quality', 
                                             save_dir, save_fn, colorcode="#330000")
  
  
  xstring <- 'mean_copy'  #mean_hmmcopy_reads_per_bin
  ystring <- 'quality'
  save_fn <- paste0(clone_id,'_quality_mean_read')
  predict_state1 <- predict_state
  plottype <- 'cell_classified' 
  plottype2 <- 'exp_cond' 
  colnames(predict_state1)[which(names(predict_state1) == "experimental_condition")] <- "exp_cond"
  pquality_mean <- plot_scatter_function_plottype_v2(predict_state1, xstring, ystring, plottype, plottype2,
                                                plottitle=paste0("Cells Metrics ",series_tag), 
                                                xlabel=paste0('Mean CN (',tag,')'), ylabel='Quality', 
                                                save_dir, save_fn, colorcode="#330000")
  
  # pquality_mean
  
  xstring <- 'mean_hmmcopy_reads_per_bin'  #
  ystring <- 'quality'
  save_fn <- paste0(clone_id,'_quality_mean_read_perbin')
  predict_state1 <- predict_state
  plottype <- 'cell_classified' 
  plottype2 <- 'exp_cond' 
  colnames(predict_state1)[which(names(predict_state1) == "experimental_condition")] <- "exp_cond"
  pquality_mean_perbin <- plot_scatter_function_plottype_v2(predict_state1, xstring, ystring, plottype, plottype2,
                                                     plottitle=paste0("Cell Metrics ",series_tag), 
                                                     xlabel=paste0('Mean cn reads per bin (',tag,')'), ylabel='Quality', 
                                                     save_dir, save_fn, colorcode="#330000")
  
  # pquality_mean_perbin
  # xstring <- 'exp_cond'
  # ystring <- 'quality'
  # save_fn <- paste0(clone_id,'_quality_exp')
  # plottype <- 'cell_classified' 
  # pexp <- plot_boxplot_function(predict_state1, xstring, ystring, plottype, 
  #                               plottitle=paste0("Filtered Metrics ",series_tag), 
  #                               xlabel=paste0('Experimental Condition (',tag,')'), 
  #                               ylabel='Quality', 
  #                               save_dir, save_fn, colorcode="#330000")
  # 
  # pexp
  
  
  # xstring <- 'median_hmmcopy_reads_per_bin'
  # ystring <- 'is_s_phase_prob'
  # plottype <- 'is_contaminated' 
  # save_fn <- paste0(clone_id,'_sphase_contaminate')
  # p5 <- plot_scatter_function(predict_state, xstring, ystring, plottype, 
  #                             plottitle=paste0("Predict Cell State ",series_tag), 
  #                             xlabel=paste0('Median hmmcopy reads per bin ',' (n=',dim(predict_state)[1],')'), ylabel='Is s-phase prob', 
  #                             save_dir, save_fn, colorcode="#53868B")
  # 
  
  xstring <- 'median_hmmcopy_reads_per_bin'
  ystring <- 'is_s_phase_prob'
  plottype <- 'cell_classified'  
  plottype2 <- 'is_contaminated'
  save_fn <- paste0(clone_id,'_sphase_contaminate')
  
  predict_state$is_contaminated <- as.factor(predict_state$is_contaminated)
  # print("Debug")
  # print(clone_id)
  # print(summary(predict_state$is_contaminated))
  psphase <- plot_scatter_function_plottype_v2(predict_state, xstring, ystring, plottype, plottype2,
                                            plottitle=paste0("Predicted Cell State ",series_tag), 
                                            xlabel=paste0('Median hmmcopy reads per bin (',tag,')'), ylabel='Is s-phase probability', 
                                            save_dir, save_fn, colorcode="#53868B")
  # psphase
  
  plts <- list(pmap=pmap, punmap=punmap, 
               pquality_median=pquality_median, 
               pquality_mean= pquality_mean,
               pquality_mean_perbin = pquality_mean_perbin,
               psphase=psphase)
  
  if(plotOutput){
    pcm <- plot_grid(plotlist = list(pmap, punmap, psphase), ncol = 3,align = 'hv')
    pcq <- plot_grid(plotlist = list(pquality_mean, pquality_mean_perbin, pquality_median), ncol = 3,align = 'hv')
    png(paste0(save_dir,"cellquality_pcm.png"), height = 2*300, width=2*900,res = 2*72)
    print(pcm)
    dev.off()
    
    png(paste0(save_dir,"cellquality_pcq.png"), height = 2*300, width=2*900,res = 2*72)
    print(pcq)
    dev.off()
    
  }
  if(save_data){
    saveRDS(plts, file = paste0(save_dir,series_tag, clone_id,'_plots_feature.rds'))
  }
  return(plts)
}

plots_features <- function(pcMap, pcUnMap, pcQuality, pcSphase, pcMean,
                           results_dir, clone_ids, ht=300, wd=900){
  save_dir <- paste0(results_dir,'cell_features/')
  
  png(paste0(save_dir,"pcMap",clone_ids,".png"), height = 2*ht, width=2*wd,res = 2*72)
  print(pcMap)
  dev.off()
  
  png(paste0(save_dir,"pcUnMap",clone_ids,".png"), height = 2*ht, width=2*wd,res = 2*72)
  print(pcUnMap)
  dev.off()
  
  png(paste0(save_dir,"pcQuality",clone_ids,".png"), height = 2*ht, width=2*wd,res = 2*72)
  print(pcQuality)
  dev.off()
  
  png(paste0(save_dir,"pcExp",clone_ids,".png"), height = 2*ht, width=2*wd,res = 2*72)
  print(pcExp)
  dev.off()
  
  png(paste0(save_dir,"pcSphase",clone_ids,".png"), height = 2*ht, width=2*wd,res = 2*72)
  print(pcSphase)
  dev.off()
  
  png(paste0(save_dir,"pcMean",clone_ids,".png"), height = 2*ht, width=2*wd,res = 2*72)
  print(pcMean)
  dev.off()
  
}

format_copynumber_values <- function(copynumber) {
  copynumber[copynumber > 11] <- 11
  for(col in colnames(copynumber)) {
    values <- as.character(copynumber[, col])
    values[values == "11"] <- "11+"
    copynumber[, col] <- values
  }
  return(copynumber)
}