

get_allele_counts_v2 <- function(lib_ticket_fn, cells_use, library_ids,
                              download_dir, save_dir, prefix='total', return_data=F){
  lib_ticket_df <- read.csv(lib_ticket_fn, stringsAsFactors=F, check.names=F)
  # View(lib_ticket_df)
  rownames(lib_ticket_df) <- lib_ticket_df$library_id
  fn <- 'allele_counts.tsv'
  allele_counts_ls <- list()
  for(library_id in library_ids){
    searchpath = paste0(download_dir,lib_ticket_df[library_id,'jira_ticket'],'/results/count_haps/')
    filenames <- list.files(searchpath, pattern=fn, recursive = TRUE)
    if(!is.null(filenames)){
      f <- filenames[1]  
      print(f)
      data <- as.data.frame(fread(paste0(searchpath, f)))
      sample_id <- gsub('sample_','',dirname(f))
      data$library_id <- rep(library_id, dim(data)[1])
      data$sample_id <- rep(sample_id, dim(data)[1])
      print(dim(data))
      allele_counts_ls[[library_id]] <- data
    } else{
      print(paste0("**** ERROR, check allele output for lib id: ",library_id))
    }
  }
  
  total_allele_counts <- dplyr::bind_rows(allele_counts_ls)
  print(dim(total_allele_counts))
  colnames(total_allele_counts)[which(colnames(total_allele_counts) == "chromosome")] <- "chr"
  # total_allele_counts <- total_allele_counts[total_allele_counts$cell_id %in% cells_use,]
  # print(dim(total_allele_counts))
  saveRDS(total_allele_counts, file = paste0(save_dir,prefix,'_allele_data.rds'))
  print("Completed!")
  if(return_data){
    return(total_allele_counts)  
  }
}  
# Read allele data and combine libs into a total csv file
get_allele_counts <- function(cells_use, library_ids,
                              download_dir, save_dir, prefix='total', return_data=F){
  fn <- 'allele_counts.tsv'
  options(warn=2)
  allele_counts_ls <- list()
  for(library_id in library_ids){
    searchpath = paste0(download_dir,library_id,'/')
    filenames <- list.files(searchpath, pattern=fn, recursive = TRUE)
    if(!is.null(filenames)){
      f <- filenames[1]  
      print(f)
      data <- as.data.frame(fread(paste0(searchpath, f)))
      sample_id <- gsub('sample_','',dirname(f))
      data$library_id <- rep(library_id, dim(data)[1])
      data$sample_id <- rep(sample_id, dim(data)[1])
      print(dim(data))
      allele_counts_ls[[library_id]] <- data
    } else{
      print(paste0("**** ERROR, check allele output for lib id: ",library_id))
    }
  }
  
  total_allele_counts <- dplyr::bind_rows(allele_counts_ls)
  print(dim(total_allele_counts))
  colnames(total_allele_counts)[which(colnames(total_allele_counts) == "chromosome")] <- "chr"
  total_allele_counts <- total_allele_counts[total_allele_counts$cell_id %in% cells_use,]
  print(dim(total_allele_counts))
  saveRDS(total_allele_counts, file = paste0(save_dir,prefix,'_allele_data.rds'))
  print("Completed!")
  if(return_data){
    return(total_allele_counts)  
  }
  
}

get_hmmcopy_reads <- function(cells_use, genes_use, library_ids,
                              download_dir, save_dir, prefix='total', return_data=F){
  bincnv_ls <- list()
  features_use <- c('chr','start','end','cell_id','state','copy','reads','gc','width')
  ls_fn <- c()
  for(library_id in library_ids){
    f1 = paste0(download_dir,library_id,'/','hmmcopy','/',library_id,'_reads.csv')
    f2 = paste0(download_dir,library_id,'/','hmmcopy_autoploidy','/',library_id,'_reads.csv')
    f3 = paste0(download_dir,library_id,'/','hmmcopy_autoploidy','/',library_id,'_multiplier0_reads.csv')
    if(file.exists(f1)){
      f <- f1
    } else if(file.exists(f2)){
      f <- f2
    } else if(file.exists(f3)){
      f <- f3
    } else{
      print(paste0("**** ERROR, check allele output for lib id: ",library_id))
      return(NULL)
    }
    if(!is.null(f)){
      print(paste0("Read lib id: ",library_id))
      print(f)
      ls_fn <- c(ls_fn,f)
    } 
  }
  if(length(ls_fn)==length(library_ids)){
    for(f in ls_fn){
      hmmread_data <- as.data.frame(fread(f))
      hmmread_data <- hmmread_data[,colnames(hmmread_data) %in% features_use]
      hmmread_data$library_id <- rep(library_id, dim(hmmread_data)[1])
      hmmread_data$chr_lb <- paste0(hmmread_data$chr,'_',hmmread_data$start,'_',hmmread_data$end)
      hmmread_data <- hmmread_data[hmmread_data$chr_lb %in% genes_use,]
      bincnv_ls[[library_id]] <- hmmread_data
    }
    total_bincnv <- dplyr::bind_rows(bincnv_ls)
    print(dim(total_bincnv))
    total_bincnv <- total_bincnv[total_bincnv$cell_id %in% cells_use,]
    print(dim(total_bincnv))
    saveRDS(total_bincnv, file = paste0(save_dir,prefix,'_CNbins_data.rds'))
    print("Completed!")
  }  
  
  if(return_data){
    return(total_bincnv)  
  }
}    

# allele_counts_fn <- paste0(save_dir,prefix,'_allele_data.rds')
# total_bincnv_fn <- paste0(save_dir,prefix,'_CNbins_data.rds')
# tag='total_filtered'
# output_fn <- paste0(save_dir,tag,'_allele_data.rds')
# results_dir <- '/home/htran/storage/datasets/metastasis_results/SA535_new_encoding/Tyler_whole'
# filtered_cells_fn <- paste0(results_dir,'/corrupt_grow/filtered_cells.txt')

filter_data_v2 <- function(allele_counts_fn, total_bincnv_fn, filtered_cells_fn, output_fn,
                           return_data=F, tag='total_filtered'){
  # Read allele data
  save_dir <- paste0(dirname(output_fn),'/')
  if (!file.exists(save_dir)){
    dir.create(save_dir)
  }
  allele_counts <- readRDS(allele_counts_fn)
  print(paste0("Allele counts: ",dim(allele_counts)[1],' ',dim(allele_counts)[2]))
  total_bincnv <- readRDS(total_bincnv_fn)
  print(paste0("total_bincnv: ",dim(total_bincnv)[1],' ',dim(total_bincnv)[2]))
  total_filtered_cn <- read.table(filtered_cells_fn, header=F)
  colnames(total_filtered_cn) <- 'cell_id'
  print(paste0("Nb filtered cells:  ",dim(total_filtered_cn)[1]))
  
  cells_use <- intersect(unique(allele_counts$cell_id),total_filtered_cn$cell_id)
  print(paste0('DEBUG 1: Nb observed cells after allele counts data intersection: ',length(cells_use)))
  cells_use <- intersect(cells_use, unique(total_bincnv$cell_id))
  print(paste0('DEBUG 1: Nb observed cells after bincn data intersection: ',length(cells_use)))
  allele_counts <- allele_counts[allele_counts$cell_id %in% cells_use,]
  total_bincnv <- total_bincnv[total_bincnv$cell_id %in% cells_use,]
  print("After filtering: ")
  print(paste0("Allele counts: ",dim(allele_counts)[1],' ',dim(allele_counts)[2]))
  print(paste0("total_bincnv: ",dim(total_bincnv)[1],' ',dim(total_bincnv)[2]))
  # if(dim(allele_counts)[1]!=dim(allele_counts))
  # Save data
  saveRDS(allele_counts, file = output_fn)

  # data.table::fwrite(total_bincnv, paste0(save_dir,'total_merged_filtered_states_cnbins.csv'),
  #                    row.names = F, quote = F)
  saveRDS(total_bincnv, file = paste0(save_dir,tag,'_CNbins_data.rds'))
  if(return_data){
    return(list('allele_counts'=allele_counts,
                'total_bincnv'=total_bincnv))
  }
}


filter_data <- function(save_dir,prefix, 
                        prefix_filtered='total_filtered',return_data=F){
  # Read allele data
  allele_counts <- readRDS(paste0(save_dir,prefix,'_allele_data.rds'))
  total_bincnv <- readRDS(paste0(save_dir,prefix,'_CNbins_data.rds'))
  total_filtered_cn <- read.csv(paste0(results_dir,'total_merged_filtered_states.csv'),
                                header=T, row.names=1, check.names = F,stringsAsFactors = FALSE)
  
  cells_use <- intersect(unique(allele_counts$cell_id),colnames(total_filtered_cn))
  cells_use <- intersect(cells_use, unique(total_bincnv$cell_id))
  print(paste0('Nb observed cells: ',length(cells_use)))
  allele_counts <- allele_counts[allele_counts$cell_id %in% cells_use,]
  total_bincnv <- total_bincnv[total_bincnv$cell_id %in% cells_use,]
  print(dim(allele_counts))
  print(dim(total_bincnv))
  prefix_filtered <- 'total_filtered'
  # Save data
  saveRDS(allele_counts, file = paste0(save_dir,prefix_filtered,'_allele_data.rds'))
  
  # data.table::fwrite(total_bincnv, paste0(save_dir,'total_merged_filtered_states_cnbins.csv'), 
  #                    row.names = F, quote = F)
  saveRDS(total_bincnv, file = paste0(save_dir,prefix_filtered,'_CNbins_data.rds'))
  if(return_data){
    return(list('allele_counts'=allele_counts,
                'total_bincnv'=total_bincnv))
  }
}

combine_BAF <- function(save_dir,prefix='total',
                        prefix_filtered='total_filtered',
                        return_data=T,filtering_data=F){
  # Filtering data first
  if(filtering_data){
    filter_data(save_dir,prefix)
  }
  print("Read data allele_counts")
  # Read filtered output
  allele_counts <- readRDS(paste0(save_dir,prefix_filtered,'_allele_data.rds'))
  print("Read data total_bincnv")
  total_bincnv <- readRDS(paste0(save_dir,prefix_filtered,'_CNbins_data.rds'))
  
  
  # First we need to do some data wrangling to merge the allele_data table 
  # and the CNbins data table so that we get a table that contains BAF and 
  # total copy number state values per bin.
  
  print("Combining...")
  CNBAF <- combineBAFCN(allele_counts, total_bincnv, binsize = 0.5e6)
  rm(allele_counts)
  rm(total_bincnv)
  CNBAF$chr_lb <- paste0(CNBAF$chr,'_',CNBAF$start,'_',CNBAF$end)
  print("Save data as rds...")
  fn <- paste0(save_dir,prefix,'_CNBAF.rds')
  saveRDS(CNBAF, file = fn)
  
  length(unique(CNBAF$chr_lb)) #?, filtered chr total: 4375
  if(file.exists(fn)){
    print(paste0('Save output of combined CNBAF as: ',fn))
  }
  
  if(return_data){
    return(CNBAF)
  }
  # length(rownames(total_filtered_cn))  #4375
  # d <- setdiff(unique(CNBAF$chr_lb),rownames(total_filtered_cn))
  # d <- rownames(total_filtered_cn)[rownames(total_filtered_cn) %ni% unique(CNBAF$chr_lb)]
  # To ensure that there is signal in the data to perform the inference 
  # it is useful to plot the BAF vs CN state. This can be done as follows. 
  # If there is signal in the data we should see clusters of points 
  # representing different allele-specific states.
  
  # p1 <- plotCNBAF(CNBAF)
  # png(paste0(save_dir,paste0(prefix,"_BAF_CNstate.png")), height = 2*460, width=2*900,res = 2*72)
  # print(p1)
  # dev.off()
  # # 
}

run_alleleSC_v2 <- function(allele_counts_fn, total_bincnv_fn, output_fn, prefix='ascn'){
  
  save_dir <- paste0(dirname(output_fn),'/')
  if (!file.exists(save_dir)){
    dir.create(save_dir)
  }
  if(file.exists(allele_counts_fn) & file.exists(total_bincnv_fn)){
    print("Loading...")
    allele_data <- readRDS(allele_counts_fn)
    print(dim(allele_data))
    CNbins <- readRDS(total_bincnv_fn)
    print(dim(CNbins))
  } else{
    stop("Files dont exist, double check input data!!!")
    return(FALSE)
  }
  
  allele_data <- schnapps::format_haplotypes_dlp(allele_data, CNbins)
  # colnames(CNBAF)
  #HMM algo perform allele-specific copy number inference 
  print("Execute allele specific...")
  ascn <- schnapps::callAlleleSpecificCN(CNbins, allele_data)
  saveRDS(ascn, file = output_fn)
  
  # # After performing this inference, to QC the results it is useful 
  # # to plot a different representation of the BAF. Here we plot the BAF 
  # # as a function of the inferred states. 
  # # The black lines indicate where we should see the mean BAF based on the state.
  
  p2 <- schnapps::plotBAFperstate(ascn, maxstate = 10)
  png(paste0(save_dir,prefix,"_BAFperstate.png"), height = 2*480, width=2*900,res = 2*72)
  print(p2)
  dev.off()
  return(TRUE)
  
}


# older version
run_alleleSC <- function(CNBAF, save_dir, prefix='total'){
  
  if(is.null(CNBAF)){
    baf_fn <- paste0(save_dir,prefix,'_CNBAF.rds')
    if(file.exists(baf_fn)){
      print("Loading...")
      CNBAF <- readRDS(baf_fn)
    } else{
      print("combined baf file dont exist!!!")
      return(FALSE)
    }
  }
  
  
  # colnames(CNBAF)
  #HMM algo perform allele-specific copy number inference 
  print("Execute allele specific...")
  CNbins <- callAlleleSpecificCNHMM(CNBAF=CNBAF,ncores = 5)
  fn <- paste0(save_dir,prefix,'_alleleCNbins_out.rds')
  saveRDS(CNbins, file = fn)
  if(file.exists(fn)){
    print(paste0('Save output of ashmm as: ',fn))
  }
  # # After performing this inference, to QC the results it is useful 
  # # to plot a different representation of the BAF. Here we plot the BAF 
  # # as a function of the inferred states. 
  # # The black lines indicate where we should see the mean BAF based on the state.
  
  p2 <- plotBAFperstate(CNbins)
  png(paste0(save_dir,prefix,"_BAFperstate.png"), height = 2*460, width=2*900,res = 2*72)
  print(p2)
  dev.off()
  return(TRUE)
  
}


visualize_hm_clustering <- function(CNbins,save_dir,prefix='total'){
  # 
  # # After having ensured the results make sense, you can plot a heatmap 
  # # of the states across all cells with the following.
  # # This will cluster the cell using umap and hdbscan.
  # 
  if(is.null(CNbins)){
    fn <- paste0(save_dir,prefix,'_alleleCNbins_out.rds')
    CNbins <- readRDS(fn)
  }
  
  p31 <- plotHeatmap(CNbins)
  png(paste0(save_dir,prefix,"_BAF_hm_dbscan_clustering.png"), height = 2*900, width=2*1400,res = 2*72)
  print(p31)
  dev.off()
  saveRDS(p31, file = paste0(save_dir,prefix,'_plots_p31_dbscan_clustering.rds'))
  # 
  # saveRDS(clusters, file = paste0(save_dir,prefix,'_clusters_hdbscan.rds'))
  # saveRDS(clustering_results, file = paste0(save_dir,prefix,'_clustering_results_hdbscan.rds'))
  
}

# CNbins <- ascn$data
# summary(CNbins$copy)
# # class(CNbins)
# clones_labels_df <- clones_df
visualize_hm_each_clone <- function(CNbins, clones_labels_df, save_dir, prefix='total'){
  if(is.null(CNbins)){
    fn <- paste0(save_dir, prefix, '_alleleCNbins_out.rds')
    CNbins <- readRDS(fn)
  }
  CNbins <- CNbins %>% inner_join(clones_labels_df, by = "cell_id")
  
  
  CNbins <- CNbins %>%
    group_by(chr, start, end, clone_id) %>%
    summarise(state = schnapps:::Mode(state),
              state_phase = schnapps:::Mode(state_phase),
              state_min = schnapps:::Mode(state_min),
              BAF = median(BAF),
              copy = median(copy,na.rm=T)) %>%
    ungroup() %>%
    mutate(cell_id = paste0("Clone ", clone_id))
  
  plots <- list()
  c <- 0
  for (cl in unique(CNbins$cell_id)){
    cat('###', cl,' \n')
    c <- c + 1
    plots[[c]] <- schnapps::plotCNprofileBAF(CNbins, cellid = cl)
    cat(' \n \n')
  }
  
  length(plots)
  saveRDS(plots, file = paste0(save_dir,prefix,'_plots_hm_each_clones.rds'))
  saveRDS(CNbins, file = paste0(save_dir,prefix,'CNbins_each_clone.rds'))
  return(plots)
  
}  


visualize_hm_each_clone_plot <- function(plots, save_dir, prefix='total'){
  plots.combinedAB <- plot_grid(plotlist = plots[1:2], nrow = 2,align = 'hv')
  plots.combinedCD <- plot_grid(plotlist = plots[3:4], nrow = 2,align = 'hv')
  # plots.combinedEF <- plot_grid(plotlist = plots[c(5,7)], nrow = 2,align = 'hv')
  # plots.combinedGH <- plot_grid(plotlist = plots[7:length(plots)], nrow = 2,align = 'hv')
  # plots.combinedGH <- plot_grid(plotlist = plots[7], nrow = 1,align = 'hv')
  png(paste0(save_dir,prefix,"_clonesAB.png"), height = 2*600, width=2*760,res = 2*72)
  print(plots.combinedAB)
  dev.off()
  
  png(paste0(save_dir,prefix,"_clonesCD.png"), height = 2*600, width=2*760,res = 2*72)
  print(plots.combinedCD)
  dev.off()
  
  # png(paste0(save_dir,prefix,"_clonesEF.png"), height = 2*600, width=2*750,res = 2*72)
  # print(plots.combinedEF)
  # dev.off()

  # png(paste0(save_dir,prefix,"_clonesG.png"), height = 2*300, width=2*750,res = 2*72)
  # print(plots.combinedGH)
  # dev.off()
}
# newick <- paste0(results_dir, 'tree.newick')
# ascn_fn <- '/home/htran/storage/datasets/metastasis_results/pseudobk_SA535/ascn_SA535/ascn_inference.rds'
# save_dir <- '/home/htran/storage/datasets/metastasis_results/pseudobk_SA535/ascn_SA535/'
# cell_clones_fn <- '/home/htran/storage/datasets/metastasis_results/SA535_new_encoding/Tyler_whole/cell_clones.csv'
# treenewick_fn <- '/home/htran/storage/datasets/metastasis_results/SA535_new_encoding/Tyler_whole/tree.newick'
# output_fn <- '/home/htran/storage/datasets/metastasis_results/pseudobk_SA535/ascn_SA535/hm_tree.png'

# results_pseudobk_dir <-  '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_pseudobulk/'
# cells_clone_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_cys/snakemake/cell_clones_Sohrab/'
# ascn_fn <- paste0(results_pseudobk_dir, 'ascn_SA535/ascn_inference.rds')
# save_dir <- paste0(results_pseudobk_dir, 'ascn_SA535/')
# cell_clones_fn <- paste0(cells_clone_dir, 'cell_clones_SA535_Sohrab.csv')
# output_fn <- paste0(results_pseudobk_dir, 'ascn_SA535/hm_SA535_cys.png')
# prefix='SA535_cys'

# results_pseudobk_dir <-  '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_pseudobulk/'
# cells_clone_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_rna/snakemake_10x_SA1035/'
# ascn_fn <- paste0(results_pseudobk_dir, 'ascn_SA1035/ascn_inference.rds')
# save_dir <- paste0(results_pseudobk_dir, 'ascn_SA1035/')
# cell_clones_fn <- paste0(cells_clone_dir, 'cell_clones_SA1035.csv')
# output_fn <- paste0(results_pseudobk_dir, 'ascn_SA1035/hm_SA1035.png')
# prefix='SA1035'

# save_dir <- '/home/htran/storage/datasets/metastasis_results/pseudobk_SA535/ascn_SA535/'
# cell_clones_fn <- '/home/htran/storage/datasets/metastasis_results/SA535_new_encoding/SA535_wholedata_v2/cell_clones_v2.csv'
# treenewick_fn <- '/home/htran/storage/datasets/metastasis_results/SA535_new_encoding/Tyler_whole/tree.newick'
# output_fn <- '/home/htran/storage/datasets/metastasis_results/pseudobk_SA535/ascn_SA535/hm_tree.png'
# ascn_fn <- paste0(save_dir, 'ascn_inference.rds')
# prefix='SA535'
# summary(as.factor(clones_labels_df$clone_id))
# save_dir <- '/home/htran/storage/datasets/metastasis_results/snv_analysis_SA919/ascn_SA919/'
# cell_clones_fn <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_wholedata/cell_clones.csv'
# # # treenewick_fn <- '/home/htran/storage/datasets/metastasis_results/SA535_new_encoding/Tyler_whole/tree.newick'
# output_fn <- '/home/htran/storage/datasets/metastasis_results/snv_analysis_SA919/ascn_SA919/ascn_SA919_02_Nov/copy_number_ascn_not_output_corrupt_tree_hm.png'
# ascn_fn <- paste0(save_dir, 'ascn_inference.rds')
# prefix='SA919'

# sps <- get_sample_id(clones_df$cell_id)
# length(unique(sps))

visualize_hm_tree_v2 <- function(ascn_fn, cell_clones_fn, output_fn,
                                 treenewick_fn, prefix='SA535'){
  
  library(ggtree)
  library(dplyr)
  save_dir <- paste0(dirname(output_fn),'/')
  if (!file.exists(save_dir)){
    dir.create(save_dir)
  }
  if(file.exists(ascn_fn)){
    ascn <- readRDS(ascn_fn)
  }else{
    stop("File ascn results does not exist!!!")
  }
  
  length(ascn)
  # clones_labels_df <- read.csv(cell_clones_fn, check.names = F, stringsAsFactors = FALSE)
  # print(dim(clones_labels_df))
  # View(head(clones_labels_df))
  # clones_labels_df <- clones_labels_df[order(clones_labels_df$clone_id, decreasing = T),]
  
  # summary(as.factor(clones_labels_df$clone_id))
  # length(unique(clones_df$clone_id))
  # length(unique(clones_labels_df$clone_id))
  # total_filtered_cn <- read.table(filtered_cells_fn, header=F)
  # colnames(total_filtered_cn) <- 'cell_id'
  # print(paste0("Nb filtered cells:  ",dim(total_filtered_cn)[1]))
  
  # clones_df <- data.frame(cell_id=unique(ascn$data$cell_id))
  # clones_df <- clones_df %>% left_join(clones_labels_df, by = "cell_id")
  # print(dim(clones_df))
  # rownames(clones_df) <- clones_df$cell_id
  # clones_df$clone_id <- ifelse(is.na(clones_df$clone_id),'None',as.character(clones_df$clone_id))
  # clones_df <- clones_df[order(clones_df$clone_id, decreasing = T),]
  # print(summary(as.factor(clones_df$clone_id)))
  # View(head(clones_df))
  
  # ascn1 <- ascn[ascn$data$cell_id]
  # View(clones_df)
  clones_df <- clones_labels_df
  sum(unique(ascn$cell_id) %in% clones_labels_df$cell_id)
  ascn1 <- ascn
  ascn <- ascn1$data
  cells_use <- intersect(unique(ascn$cell_id),clones_labels_df$cell_id)
  clones_df <- clones_df[clones_df$cell_id %in% cells_use,]
  ascn <- ascn[ascn$cell_id %in% cells_use,]
  dim(ascn)
  set.seed(1)
  
  # p1 <- schnapps::plotHeatmap(ascn, plotcol = "state", plottree = F, 
  #                             spacer_cols = 15, clusters = clones_df)
  # p2 <- schnapps::plotHeatmap(ascn, plotcol = "state_min", plottree = F,
  #                             spacer_cols = 15, clusters = clones_df)
  # p3 <- schnapps::plotHeatmap(ascn, plotcol = "state_BAF", plottree = F, 
  #                             spacer_cols = 15, clusters = clones_df)
  # p4 <- schnapps::plotHeatmap(ascn, plotcol = 'state_phase', plottree = F, 
  #                             clusters = clones_df)
  
   
  p2 <- plotHeatmap(ascn, plotcol = "state_min", plottree = F,
                    spacer_cols = 15, clusters = clones_df)
  
  p3 <- plotHeatmap(ascn, plotcol = "state_BAF", plottree = F,
                    spacer_cols = 15, clusters = clones_df)

  p1 <- plotHeatmap(ascn, plotcol = "state", plottree = F,
                              spacer_cols = 15, clusters = clones_df)

 

  p4 <- plotHeatmap(ascn, plotcol = 'state_phase', plottree = F,
                              clusters = clones_df)

  # save_dir <- paste0(save_dir,'ascn_SA919_02_Nov/')
  # png(paste0(save_dir,prefix,"_state_min_ascn.png"), height = 2*700, width=2*1450,res = 2*72)
  # print(p2)
  # dev.off()
  
  png(paste0(save_dir,prefix,"_state_phase_ascn.png"), height = 2*700, width=2*1450,res = 2*72)
  print(p4)
  dev.off()
  
  png(output_fn, height = 2*700, width=2*1450,res = 2*72)
  print(p1)
  dev.off()
  
  png(paste0(save_dir,prefix,"_state_BAF_ascn.png"), height = 2*700, width=2*1450,res = 2*72)
  print(p3)
  dev.off()
  # The tree
  # tree <- read.tree(treenewick_fn)
  # all_cells <- grep('cell', tree$tip.label, value = T)
  # all_cells <- grep('cell', tree_cg$tip.label, value = T)
  # length(all_cells)
  # sum()
  # none_cells <- all_cells[!all_cells %in% paste0('cell_',clones_df$cell_id)]
  # length(none_cells)
  # tree_cg <- ape::drop.tip(tree, none_cells, trim.internal =F, collapse.singles = F)
  # length(tree_cg$tip.label)
  # p5 <- schnapps::plotHeatmap(ascn, plotcol = 'state_phase', plottree = T, 
  #                             clusters = clones_df, tree = tree_cg, spacer_cols = 15)
  # png(paste0(save_dir,prefix,"_state_phase_treenewick_ascn.png"), height = 2*800, width=2*1400,res = 2*72)
  # print(p5)
  # dev.off()
  
  # png(output_fn, height = 2*900, width=2*1400,res = 2*72)
  # print(p1)
  # dev.off()
  # saveRDS(list(p32,p33), file = paste0(save_dir,prefix,'_BAF_hm_tree_plot.rds'))
}



visualize_hm_tree <- function(CNbins, clones_labels_df, 
                              results_dir, save_dir, prefix='total'){
  
  library(ggtree)
  if(is.null(CNbins)){
    fn <- paste0(save_dir,prefix,'_alleleCNbins_out.rds')
    CNbins <- readRDS(fn)
  }
  # The tree
  newick <- paste0(results_dir, 'tree.newick')
  treenk <- read.tree(newick)
  
  total_filtered_cn <- read.csv(paste0(results_dir,'total_merged_filtered_states.csv'),
                                header=T, row.names=1, check.names = F,stringsAsFactors = FALSE)
  clones_df <- data.frame(cell_id=colnames(total_filtered_cn))
  clones_df <- clones_df %>% left_join(clones_labels_df, by = "cell_id")
  sum(is.na(clones_df$clone_id))
  # unique(clones_df$clone_id)
  rownames(clones_df) <- clones_df$cell_id
  clones_df$clone_id <- ifelse(is.na(clones_df$clone_id),'None',as.character(clones_df$clone_id))
  clones_df <- clones_df[order(clones_df$clone_id),]
  # p32 <- h
  p32 <- plotHeatmap(CNbins=CNbins, clusters=clones_df, 
                     tree=treenk, plotcol = 'state_phase', plottree = F)

  png(paste0(save_dir,prefix,"_BAF_hm_tree_ascn.png"), height = 2*900, width=2*1400,res = 2*72)
  print(p32)
  dev.off()
  
  # p33 <- plotHeatmap(CNbins=CNbins, clusters=clones_df, 
  #                    tree=treenk, plotcol = 'state')
  # 
  # png(paste0(save_dir,prefix,"_BAF_hm_cnstate.png"), height = 2*900, width=2*1400,res = 2*72)
  # print(p33)
  # dev.off()
  # saveRDS(list(p32,p33), file = paste0(save_dir,prefix,'_BAF_hm_tree_plot.rds'))
}


umap_clustering <- function(CNbins, n_neighbors = 30, min_dist = 0, minPts = 30, 
          seed = 1, field = "state") 
{
  message("Creating CN matrix...")
  cnmatrix <- schnapps:::createCNmatrix(CNbins, na.rm = TRUE, field = field)
  dim(cnmatrix)
  cnmatrix <- subset(cnmatrix, select = -c(chr, start, end, 
                                           width))
  cnmatrix <- t(cnmatrix)
  set.seed(seed)
  message("Calculating UMAP dimensionality reduction...")
  umapresults <- uwot::umap(cnmatrix, n_neighbors = n_neighbors, 
                            n_components = 2, min_dist = min_dist, ret_model = TRUE, 
                            ret_nn = TRUE)
  dfumap <- data.frame(umap1 = umapresults$embedding[, 1], 
                       umap2 = umapresults$embedding[, 2], cell_id = row.names(cnmatrix))
  rownames(dfumap) <- row.names(cnmatrix)
  message("Clustering cells using hdbscan...")
  hdbscanresults <- dbscan::hdbscan(dfumap[, 1:2], minPts = minPts, 
                                    gen_simplified_tree = TRUE)
  clusterids <- hdbscanresults$cluster
  clusterids[clusterids == 0] <- 702
  LETTERS702 <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, 
                                                              LETTERS)))
  dfumap$clone_id <- LETTERS702[clusterids]
  dfumap <- dfumap %>% dplyr::mutate(clone_id = ifelse(clone_id == 
                                                         "ZZ", "0", clone_id))
  tree <- ape::as.phylo(hdbscanresults$hc, use.labels = TRUE)
  tree$tip.label <- row.names(cnmatrix)[as.numeric(tree$tip.label)]
  message(paste0("Identified ", length(unique(dfumap$clone_id)), 
                 " clusters"))
  message("Distribution of clusters:")
  print(table(dfumap$clone_id))
  return(list(clustering = dfumap, hdbscanresults = hdbscanresults, 
              umapresults = umapresults, tree = tree))
}
