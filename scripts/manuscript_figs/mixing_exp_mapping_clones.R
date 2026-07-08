# Mixing exp, mapping clones back to main clones exp
library(dplyr)


## Median clonal profiles in main experiment
get_cn_profiles <- function(){
  source("/home/htran/Projects/hakwoo_project/metastasis_material/scripts/corrupt_tree/src/cn_change/cn_profile_utils.R")
  results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_wholedata/'
  save_dir <- paste0(results_dir,'CN_profile/')
  if(!dir.exists(save_dir)){
    dir.create(save_dir) 
  }
  copynumber_fn <- paste0(results_dir,'total_merged_filtered_states.csv')
  datatag <- 'SA919'
  cellclone_fn <- paste0(results_dir,'cell_clones.csv')
  library_grouping_fn <- paste0(results_dir,'library_groupings.csv')
  ref_median_clones <- get_median_genotype_v3(copynumber_fn, 
                                              datatag, save_dir,
                                              cellclone_fn=cellclone_fn, library_grouping_fn=library_grouping_fn)
  
  
  return(ref_median_clones)
}
# obs_library_id <- 'A138967B'
map_clones_to_ref_labels <- function(obs_library_id, curr_cell_clones_fn, ref_cell_clones_fn){
  
  source("/home/htran/Projects/hakwoo_project/metastasis_material/scripts/corrupt_tree/src/cn_change/cn_profile_utils.R")
  ## Median clonal profiles for each library in mixing experiment
  results_dir <- '/home/htran/storage/datasets/metastasis_results/dlp_SA919_mixing_exp/'
  curr_save_dir <- paste0(results_dir,'CN_profile/')
  if(!dir.exists(curr_save_dir)){
    dir.create(curr_save_dir) 
  }
  
  curr_cell_clones_fn1 <- paste0(results_dir, 'graph_cut/',obs_library_id,'_cell_clones.csv')
  curr_cell_clones_fn2 <- paste0(results_dir, 'graph_cut/',obs_library_id,'_cell_clones.csv.gz')
  if(file.exists(curr_cell_clones_fn1)){
    curr_cell_clones_fn <- curr_cell_clones_fn1
  }else if(file.exists(curr_cell_clones_fn2)){
    curr_cell_clones_fn <- curr_cell_clones_fn2
  }else{
    stop('Do not exist cell clonal assignment, re-run cell clustering algo first')
  }
  print(curr_cell_clones_fn)
  curr_copynumber_fn <- paste0(results_dir, obs_library_id,'_filtered_states.csv.gz')
  curr_cell_clones <- data.table::fread(curr_cell_clones_fn)
  print(summary(as.factor(curr_cell_clones$clone_id)))
  
  ## Note: Manual change clone label for this library id
  # obs_library_id <- 'A130854A'
  ## Note: clone None has similar profile as clone A
  # curr_cell_clones$clone_id <- ifelse(curr_cell_clones$clone_id=='None','A',curr_cell_clones$clone_id)

  # obs_library_id <- 'A138912A'
  ## Note: clone None has similar profile as clone A
  # curr_cell_clones$clone_id <- ifelse(curr_cell_clones$clone_id=='None','A',curr_cell_clones$clone_id)
  
  # obs_library_id <- 'A130854B'
  ## Note: clone A contains all noisy cells
  # curr_cell_clones$clone_id <- ifelse(curr_cell_clones$clone_id=='A','None',curr_cell_clones$clone_id)
  # (copynumber_fn=curr_copynumber_fn, 
  #   datatag=obs_library_id, save_dir=curr_save_dir,
  #   cell_clones=curr_cell_clones, library_grouping_fn=NULL)
  curr_median_clones <- get_median_genotype_v3(curr_copynumber_fn, 
                                            datatag=obs_library_id, curr_save_dir,
                                            cell_clones=curr_cell_clones)
  dim(curr_median_clones)
  head(curr_median_clones)
  colnames(curr_median_clones) <- ifelse(colnames(curr_median_clones)!='chr_desc',paste0(obs_library_id,'_',colnames(curr_median_clones)),'chr_desc')
  
  
  ref_median_clones_fn <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_wholedata/CN_profile/SA919_median_cnv.csv.gz'
  ref_median_clones <- data.table::fread(ref_median_clones_fn)
  # Find the closest ref clone label first, and get distance to closest one
    # ref_median_clones <- median_cnv_pivot
  # head(ref_median_clones)
  dim(ref_median_clones)
  ## Merge ref clones, and curr clones mtx into one mtx
  ref_median_clones <- ref_median_clones %>%
    tibble::column_to_rownames('chr_desc')
  colnames(ref_median_clones) <- paste0('main_',colnames(ref_median_clones))
  # curr_median_clones_backup <-curr_median_clones
  # curr_median_clones <- curr_median_clones_backup
  curr_median_clones <- curr_median_clones %>%
    tibble::column_to_rownames('chr_desc')
  dim(curr_median_clones)
  rownames(ref_median_clones)[1]
  curr_median_clones <- curr_median_clones[rownames(ref_median_clones),,drop=F]
  combined_median_cnv <- dplyr::bind_cols(ref_median_clones, curr_median_clones)
  dim(combined_median_cnv)
  res <- compute_dist_mat(obs_library_id, combined_median_cnv, curr_save_dir, use_hamming = TRUE)
  # dim(res$out_mtx)  
  dis_cnv <- res$out_mtx
  distance_Hamming_thrs <- 0.005 # less than 0.5% of difference between 2 clones, mean 20 regions
  obs_clones <- colnames(curr_median_clones)
  assigned <- tibble::tibble()
  for(c in obs_clones){
    tmp <- dis_cnv %>%
      dplyr::filter(SourceClone==c & !TargetClone %in% obs_clones)
    tmp <- tmp[which.min(tmp$CNA_Distance),]
    assigned <- dplyr::bind_rows(assigned, tmp)
  }
  print(assigned)
  data.table::fwrite(assigned, paste0(curr_save_dir,obs_library_id,'_assigned_clones.csv'))
  
  ## Updated labels
  assigned$SourceClone <- gsub(paste0(obs_library_id,'_'),'',assigned$SourceClone)
  assigned$TargetClone <- gsub('main_','',assigned$TargetClone)
  assigned$TargetClone <- ifelse(assigned$CNA_Distance>distance_Hamming_thrs,paste0('extra_',assigned$TargetClone),assigned$TargetClone)
  
  ## Note: manual assign ancestor or successor for this library
  # obs_library_id <- 'A138967B'
  # assigned$TargetClone <- ifelse(assigned$CNA_Distance>distance_Hamming_thrs,paste0(assigned$TargetClone,'0'),assigned$TargetClone)
  
  
  assigned$CNA_Distance <- NULL
  assigned <- dplyr::bind_rows(assigned, tibble::tibble(SourceClone='None',TargetClone='None'))
  curr_cell_clones <- curr_cell_clones %>%
    dplyr::left_join(assigned, by=c('clone_id'='SourceClone'))
  
  colnames(curr_cell_clones)
  print(summary(as.factor(curr_cell_clones$TargetClone)))
  # summary(as.factor(curr_cell_clones$clone_id))
  ## Note: manual label for ancestor, successor clones
  # obs_library_id <- 'A98166B'
  # curr_cell_clones$TargetClone <- ifelse(curr_cell_clones$TargetClone=='extra_A','A0',curr_cell_clones$TargetClone)

  
  # obs_library_id <- 'A130841A'
  # curr_cell_clones$TargetClone <- ifelse(curr_cell_clones$TargetClone=='extra_C','C1',curr_cell_clones$TargetClone)

  # obs_library_id <- 'A130841A'
  # curr_cell_clones$TargetClone <- ifelse(curr_cell_clones$TargetClone=='extra_A','A0',curr_cell_clones$TargetClone)
  
  # obs_library_id <- 'A138912B'
  # curr_cell_clones$TargetClone <- ifelse(curr_cell_clones$TargetClone=='extra_C','C1',curr_cell_clones$TargetClone)
  
  curr_cell_clones <- curr_cell_clones %>%
    dplyr::select(cell_id, TargetClone, clone_id) %>%
    dplyr::rename(previous_clone_id=clone_id, clone_id=TargetClone)
  print(summary(as.factor(curr_cell_clones$clone_id)))
  data.table::fwrite(curr_cell_clones, paste0(results_dir, 'graph_cut/',obs_library_id,'_cell_clones_v2.csv.gz'))
  
  
}


summary_results <- function(){
  script_dir <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/dlp_trees/SA919_mixing_experiment/'
  grouping_file <- paste0(script_dir, 'Metastasis_Hakwoo_mixing_exp_SA919_results.csv')
  lib_df <- data.table::fread(grouping_file)
  obs_libs <- unique(lib_df$library_id)
  print(obs_libs)
  colnames(lib_df)
  dim(lib_df)
  lib_df$jira_ticket <- NULL
  
  jira_tickets <- data.table::fread(paste0(script_dir, 'jira_tickets_library_grouping_SA919.csv'))
  dim(jira_tickets)
  colnames(jira_tickets)
  jira_tickets <- jira_tickets %>%
    dplyr::select(library_id, jira_ticket)
  lib_df <- lib_df %>% 
    dplyr::left_join(jira_tickets, by='library_id')
  data.table::fwrite(lib_df, grouping_file)
  
  # View(lib_df)
  qc <- data.table::fread('/home/htran/storage/datasets/metastasis_results/dlp_SA919_mixing_exp/SA919_mixing_QC_cells_2023Sep22_002931.csv.gz')
  qc$V1 <- NULL
  colnames(qc)
  dim(lib_df)
  lib_df <- lib_df %>% 
    dplyr::left_join(qc, by='library_id')
  lib_df$nb_filtered_cells
  data.table::fwrite(lib_df, grouping_file)
  
  metasample_fn <- paste0(results_dir, datatag,'_total_raw_cells.csv.gz')
  metasample_df <- data.table::fread(metasample_fn)
  lib_df <- lib_df %>% 
    dplyr::left_join(metasample_df, by='library_id')
  data.table::fwrite(lib_df, grouping_file)
  
  
  metasample_fn <- paste0(results_dir, datatag,'_prop_cell_clones_final.csv.gz')
  metasample_df <- data.table::fread(metasample_fn)
  lib_df <- lib_df %>% 
    dplyr::left_join(metasample_df, by='library_id')
  data.table::fwrite(lib_df, grouping_file)
  dim(lib_df)
  colnames(lib_df)
  t <- lib_df %>% 
    dplyr::select(transplanted_mouse, Expected, `Achieved Results`,prop_cell_clones)
  data.table::fwrite(t, paste0(script_dir, 'Metastasis_Hakwoo_mixing_exp_SA919_results_short_version.csv'))
  
  
  results_dir <- '/home/htran/storage/datasets/metastasis_results/dlp_SA919_mixing_exp/'
  datatag <- 'SA919_mixing'
  metasample_fn <- paste0(results_dir, datatag,'_prop_cell_clones_final.csv.gz')
  metasample_df <- data.table::fread(metasample_fn)
  dim(metasample_df)  
  View(metasample_df)
}

load_raw_nb_cells <- function(obs_libs, download_dir, save_data=T){
  datatag <- 'SA919_mixing'
  download_dir <- '/home/htran/storage/raw_DLP/metastasis_DLP/SA919/'
  results_dir <- '/home/htran/storage/datasets/metastasis_results/dlp_SA919_mixing_exp/'
  script_dir <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/dlp_trees/SA919_mixing_experiment/'
  grouping_file <- paste0(script_dir, 'Metastasis_Hakwoo_mixing_exp_SA919_results.csv')
  lib_df <- data.table::fread(grouping_file)
  colnames(lib_df)
  # View(lib_df)
  obs_libs <- unique(lib_df$library_id)
  
  print(paste(obs_libs, collapse = ' '))
  obs_metrics <- list()
  
  for(l in obs_libs){
    
    lib_fn0 <- paste0(download_dir,l,'/annotation/',l,'_metrics.csv') # csv.gz or csv
    lib_fn1 <- paste0(download_dir,l,'/annotation/','metrics.csv.gz') # csv.gz or csv
    # lib_fn <- paste0(download_dir,l,'/',l,'/hmmcopy/','reads.csv.gz') # csv.gz or csv
    lib_fn2 <- paste0(download_dir,l,'/annotation/',l,'_metrics.csv.gz') # csv.gz or csv
    lib_fn <- ''
    lib_fn <- ifelse(file.exists(lib_fn0),lib_fn0,
                     ifelse(file.exists(lib_fn1),lib_fn1,
                            ifelse(file.exists(lib_fn2),lib_fn2,'')))
    # if(file.exists(lib_fn0)){
    #   lib_fn <- lib_fn0
    # } else{
    #   if(file.exists(lib_fn2)){
    #     lib_fn <- lib_fn2
    #   } else{
    #     print(paste0('***ERROR: library: ', l))
    #     stop('Do not exist hmmcopy reads for library, double check input data')
    #   } 
    # } 
    if(lib_fn==''){
      stop('Do not exist hmmcopy reads for library, double check input data')
    }
    if(lib_fn==''){
      stop('Do not exist hmmcopy reads for library, double check input data')
    }
    if(lib_fn!=''){
      tmp <- data.table::fread(lib_fn)
      print(dim(tmp))
      if(dim(tmp)[1]>0){
        m <- tibble::tibble(library_id=l,nb_total_cells_raw=dim(tmp)[1])
        obs_metrics[[l]] <- m
      }else{
        print(paste0('***Warning:  Do not have any output for clone: ',c))
      }
      
    }else{
      print(paste0('***ERROR: library: ', l))
      # print('Do not exist hmmcopy reads for library, double check input data')
      stop('Do not exist hmmcopy reads for library, double check input data')
    }
    
  }
  metasample_df <- as.data.frame(dplyr::bind_rows(obs_metrics))
  if(save_data){
    metasample_fn <- paste0(results_dir, datatag,'_total_raw_cells.csv.gz')
    data.table::fwrite(metasample_df, metasample_fn)  
  }
  
  print(dim(metasample_df))
}

process_prevalence_main_exp <- function(){
  input_dir <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/dlp_trees/SA919_mixing_experiment/'
  stat_main_exp <- data.table::fread(paste0(input_dir, 'prevalence_main_exp.csv.gz'))
  dim(stat_main_exp)
  colnames(stat_main_exp)
  stat_main_exp <- stat_main_exp %>%
    dplyr::select(clone_id, pct_cells, mainsite, transplanted_mouse)
  unique(stat_main_exp$transplanted_mouse)
  obs_samples1 <- c("M2164_X4_Axillary","M2112252_X7_SupraSpinal",
                   "M2164_X4_Primary","M2112251_X7_Primary")
  obs_samples2 <- c("M2112253_X7_Primary")
  stat_main_exp1 <- stat_main_exp %>%
    dplyr::filter(transplanted_mouse %in% obs_samples1)
  
  
  ## Mixing with ratio 0.25 for this sample
  stat_main_exp2 <- stat_main_exp %>%
    dplyr::filter(transplanted_mouse %in% obs_samples2)
  stat_main_exp2$pct_cells <- stat_main_exp2$pct_cells * 0.25
  none_clone <- tibble::tibble(clone_id='None', pct_cells=100-sum(stat_main_exp2$pct_cells),
                 mainsite='Primary',transplanted_mouse='M2112253_X7_Primary')
  stat_main_exp2 <- dplyr::bind_rows(stat_main_exp2, none_clone)
  
  
  stat_main <- dplyr::bind_rows(stat_main_exp1, stat_main_exp2)
  data.table::fwrite(stat_main, paste0(input_dir, 'prevalence_main_exp_4_mixing.csv.gz'))
}

load_cell_clones_prop <- function(obs_libs, download_dir, save_data=T){
  datatag <- 'SA919_mixing'
  download_dir <- '/home/htran/storage/raw_DLP/metastasis_DLP/SA919/'
  results_dir <- '/home/htran/storage/datasets/metastasis_results/dlp_SA919_mixing_exp/'
  script_dir <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/dlp_trees/SA919_mixing_experiment/'
  grouping_file <- paste0(script_dir, 'Metastasis_Hakwoo_mixing_exp_SA919_results.csv')
  lib_df <- data.table::fread(grouping_file)
  colnames(lib_df)
  dim(lib_df)
  obs_libs <- unique(lib_df$library_id)
  
  print(paste(obs_libs, collapse = ' '))
  obs_metrics <- list()
  
  for(l in obs_libs){
    cell_clones_fn <- paste0(results_dir, 'graph_cut/',l,'_cell_clones_v2.csv.gz')
    if(!file.exists(cell_clones_fn)){
      stop('Do not exist cell_clones_fn for library, double check input data')
    }
    cell_clones_tmp <- data.table::fread(cell_clones_fn)
    print(dim(cell_clones_tmp))
    cell_clones_tmp <- cell_clones_tmp %>%
      dplyr::filter(!clone_id %in% c('None','unassigned'))
    print(summary(as.factor(cell_clones_tmp$clone_id)))
    ## summary proportion here
    nb_cells <- dim(cell_clones_tmp)[1]
    # cell_clones_tmp <- cell_clones_tmp %>%
    #   dplyr::group_by(clone_id) %>%
    #   dplyr::summarise(pct_cells=round(100*n()/nb_cells,2)) %>%
    #   dplyr::mutate(desc=paste0(clone_id,'-',pct_cells))
    # m <- tibble::tibble(library_id=l,prop_cell_clones=paste(cell_clones_tmp$desc, collapse = '_'))
    # print(m)
    
    cell_clones_tmp <- cell_clones_tmp %>%
      dplyr::group_by(clone_id) %>%
      dplyr::summarise(pct_cells=100*n()/nb_cells) #%>%
      # dplyr::mutate(desc=paste0(clone_id,'-',pct_cells))
    cell_clones_tmp$library_id <- l
    obs_metrics[[l]] <- cell_clones_tmp
     
    
  }
  metasample_df <- as.data.frame(dplyr::bind_rows(obs_metrics))
  dim(metasample_df)
  
  t <- lib_df %>% 
    dplyr::select(mainsite, transplanted_mouse, library_id)
  
  metasample_df <- metasample_df %>%
    left_join(t, by='library_id')
  if(save_data){
    # metasample_fn <- paste0(results_dir, datatag,'_prop_cell_clones_final.csv.gz')
    metasample_fn <- paste0(results_dir, datatag,'_prop_cell_clones_final_long_version.csv.gz')
    data.table::fwrite(metasample_df, metasample_fn)  
  }
  
  
  
}

plot_mixing_exp_Fig3 <- function(){
  datatag <- 'SA919_mixing'
  download_dir <- '/home/htran/storage/raw_DLP/metastasis_DLP/SA919/'
  results_dir <- '/home/htran/storage/datasets/metastasis_results/dlp_SA919_mixing_exp/'
  script_dir <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/dlp_trees/SA919_mixing_experiment/'
  script_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/dlp_trees/'
  
  metasample_main_df <- data.table::fread(paste0(script_dir,'SA919_mixing_experiment/','prevalence_main_exp_4_mixing.csv.gz'))
  metasample_df <- data.table::fread(paste0(script_dir, 'SA919_mixing_experiment/', datatag,'_prop_cell_clones_final_long_version.csv'))
  print(dim(metasample_df))
  print(dim(metasample_main_df))
  View(metasample_main_df)
  metasample_df$library_id <- NULL
  metasample_df$transplanted_mouse <- stringr::str_sub(metasample_df$transplanted_mouse, 4, length(metasample_df$transplanted_mouse))
  metasample_df$transplanted_mouse <- gsub(' ','',metasample_df$transplanted_mouse)
  
  colnames(metasample_df)
  colnames(metasample_main_df)
  ## Combining main and mixing prevalence for plotting
  metasample_df <- dplyr::bind_rows(metasample_main_df, metasample_df)
  # View(metasample_df)
  # input_dir <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/dlp_trees/'
  colors_df <- data.table::fread(paste0(script_dir, 'colorCode_clones/color_code_SA919_mixing_experiment.csv'))
  cols_use <- colors_df$color
  names(cols_use) <- colors_df$clone_id
  cols_use <- c(cols_use, 'black')
  names(cols_use) <- c(colors_df$clone_id, 'None')
  
  ## Using piechart for primary, and donut chart for met
  viz_prop_chart()
  
  
  library(ggplot2)
  my_font <- 'Helvetica'
  p <- ggplot(data=metasample_df, aes(x=transplanted_mouse, y=pct_cells, fill=clone_id)) +
    geom_bar(stat="identity", width=0.3)+
    # facet_grid(. ~ fov)+ 
    facet_wrap(~ mainsite, strip.position = "top", nrow = 2,scales = "free_y") + #ncol = 4, 
    theme_bw() +
    scale_fill_manual(values=cols_use) + 
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
    labs(title="SA919 clonal prevalence mixing exp", x='Metastasis sites', y='Clonal proportion')
  # p  
  
  ggsave(  
    filename = paste0(script_dir,"mixing_exp_clonal_prop.svg"),  
    plot = p,  
    height = 6.5,  
    width = 8.5
    #useDingbats=F
  )
  ggsave(  
    filename = paste0(script_dir,"mixing_exp_clonal_prop.png"),  
    plot = p,  
    height = 6.5,  
    width = 8.5
    #useDingbats=F
  )
  
  ## To Do: need to add SA919 main exp here, so there is same scale
  
}

viz_prop_chart <- function(){
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(tidyverse)
  # View(metasample_df)
  
  ## To Do: move other clones to main clones
  save_fig_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/main_figures/fig4_revised/'
  dir.create(save_fig_dir)
  colors_df <- data.table::fread(paste0(script_dir, 'colorCode_clones/color_code_SA919_mixing_experiment.csv'))
  cols_use <- colors_df$color
  names(cols_use) <- colors_df$clone_id
  cols_use <- c(cols_use, 'black')
  names(cols_use) <- c(colors_df$clone_id, 'None')
  
  metasample_df <- metasample_df %>%
    dplyr::rename(sample_id=transplanted_mouse, original_clone=clone_id)
  
  mapping_clones <- c("A","B","C","None","C","A","A","C","B","C")
  names(mapping_clones) <- c("A","B","C","None","C10","A01","A02","C11","B0","C0")
  mapping_df <- tibble::tibble(clone_id=mapping_clones, original_clone=names(mapping_clones))
  metasample_df <- metasample_df %>%
    dplyr::left_join(mapping_df, by='original_clone')
  unique(metasample_df$clone_id)
  metasample_df <- metasample_df %>%
    dplyr::group_by(sample_id, clone_id, mainsite, original_clone) %>%
    dplyr::summarise(pct_cells=sum(pct_cells))
  data.table::fwrite(metasample_df, paste0(save_fig_dir, 'mixing_exp_prop.csv.gz'))
  View(metasample_df)
  ## Main experiment
  main_df <- data.table::fread(paste0(script_dir,'clonal_prevalence_per_sample/','SA919_cells_proportions.csv'))
  dim(main_df)
  head(main_df)
  View(main_df)
  mts <- sapply(strsplit(main_df$sample_id, '_'), function(x){
    return(x[3])
  })
  main_df$mainsite <- as.character(mts)
  # mids <- sapply(strsplit(main_df$sample_id, '_'), function(x){
  #   return(paste0(x[2],'_',x[3]))
  # })
  mids <- main_df$sample_id
  mids <- gsub('SA919','',mids)
  mids <- gsub('mouse','',mids)
  mids <- gsub('_Primary','P',mids)
  mids <- gsub('_Metastasis','M',mids)
  main_df$sample_id <- mids
  
  total_df <- dplyr::bind_rows(main_df, metasample_df)
  
  # View(total_df)
  
  met_plt_list <- list()
  for(s in unique(total_df$sample_id)){
    df <- total_df %>%
      dplyr::filter(sample_id==s)
    t <- sum(df$pct_cells)
    
    if(t>100){ # issue with round function, slightly higher than 100%
      print(s)
      print(t)
      m <- max(df$pct_cells)
      df <- df %>%
        dplyr::mutate(pct_cells= case_when(
          pct_cells==m ~ (pct_cells- t + 100),
          TRUE ~ pct_cells
        ))
      print(sum(df$pct_cells))
    }
    
    cols_use1 <- cols_use[unique(df$clone_id)]
    # print(cols_use1)
    df <- as.data.frame(df)
    p <- NULL
    if(df$mainsite[1]=='Primary'){
      p <- clone_prevalence_piechart(df, s, save_fig_dir, cols_use1)
    }else{
      p <- clone_prevalence_donutchart(df, s, save_fig_dir, cols_use1)
    }
    met_plt_list[[s]] <- p
  }
  saveRDS(met_plt_list, paste0(save_fig_dir, 'SA919_met_plt_list_main_mixing.rds'))
  # plot one with label and one without label because label change the scale of plot
  p_total <- cowplot::plot_grid(plotlist = met_plt_list, ncol=6, labels = names(met_plt_list))
  p_total
  p_total_nolabel <- cowplot::plot_grid(plotlist = met_plt_list, ncol=6)
  p_total_nolabel
  
  p_donut_lg <- clone_prevalence_donutchart_legend(save_fig_dir)
  p_pie_lg <- clone_prevalence_piechart_legend(save_fig_dir)
  
  p_total_lg <- cowplot::plot_grid(p_pie_lg, p_donut_lg, ncol=1)
  
  ggsave(  
    filename = paste0(save_fig_dir,"clonal_prop_SA919_total.svg"),  
    plot = p_total_nolabel,  
    height = 6.5,  
    width = 6
    #useDingbats=F
  )
  ggsave(  
    filename = paste0(save_fig_dir,"clonal_prop_SA919_total.png"),  
    plot = p_total,  
    height = 10,  
    width = 16
    #useDingbats=F
  )
  ggsave(  
    filename = paste0(save_fig_dir,"prop_mainsite_legend.svg"),  
    plot = p_total_lg,  
    height = 1,  
    width = 0.5
    #useDingbats=F
  )
}
# pct_cells, clone_id, cols_use
clone_prevalence_donutchart <- function(df, datatag, save_dir, cols_use){
  # df <- data.frame(value = c(10, 30, 32, 28),
  #                  clone_id = paste0("G", 1:4))
  # Hole size
  hsize <- 1.5
  
  df <- df %>% 
    mutate(x = hsize)
  df2 <- df %>% 
    mutate(csum = rev(cumsum(rev(pct_cells))), 
           pos = pct_cells/2 + lead(csum, 1),
           pos = if_else(is.na(pos), pct_cells/2, pos))
  
  p <- ggplot(df, aes(x = hsize, y = pct_cells, fill = clone_id)) +
    geom_col(color = "black") +
    # geom_text(aes(label = pct_cells),
    #           position = position_stack(vjust = 0.5)) +
    coord_polar(theta = "y") +
    # scale_fill_brewer(palette = "GnBu") +
    scale_fill_manual(values = cols_use) + 
    xlim(c(0.2, hsize + 0.5)) +
    theme_void() + 
    theme(legend.position = 'none')
  # panel.background = element_rect(fill = "transparent", colour = NA),
  # panel.grid = element_blank(),
  # axis.title = element_blank(),
  # axis.ticks = element_blank(),
  # axis.text = element_blank()
  # p
  
  # saveRDS(p, paste0(save_dir, datatag, '_donutchart.rds'))
  return(p)
}
# pct_cells, clone_id, cols_use
clone_prevalence_piechart <- function(df, datatag, save_dir, cols_use){
  
  
  # df <- data.frame(pct_cells = c(15, 25, 32, 28),
  #                  clone_id = paste0("G", 1:4))
  # cols_use <- c('red','blue','cyan','pink')
  # names(cols_use) <- df$clone_id
  # Get the positions
  df2 <- df %>% 
    mutate(csum = rev(cumsum(rev(pct_cells))), 
           pos = pct_cells/2 + lead(csum, 1),
           pos = if_else(is.na(pos), pct_cells/2, pos))
  
  p <- ggplot(df, aes(x = "" , y = pct_cells, fill = fct_inorder(clone_id))) +
    geom_col(width = 1, color = 1) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = cols_use) + 
    # scale_fill_brewer(palette = "Pastel1") +
    # geom_text_repel(data = df2,
    #                 aes(y = pos, label = paste0(clone_id,' ',pct_cells, "%")),
    #                 size = 4.5, nudge_x = 1, show.legend = FALSE, segment.color = NA) +
    # guides(fill = guide_legend(title = "clone_id")) +
    theme_void() + 
    theme(legend.position = 'none')
  # p
  # cols_use <- rep('white', 4)
  # names(cols_use) <- df$clone_id
  # p_legend <- ggplot(df, aes(x = "" , y = pct_cells, fill= fct_inorder(clone_id))) + #, fill = fct_inorder(clone_id))
  #   geom_col(width = 1, color = 1) +
  #   coord_polar(theta = "y") +
  #   scale_fill_manual(values = cols_use)+
  #   theme_void() + 
  #   theme(legend.position = 'none')
  # p_legend
  # saveRDS(p, paste0(save_dir, datatag, '_piechart_legend_shape.rds'))
  # saveRDS(p, paste0(save_dir, datatag, '_piechart.rds'))
  # res <- list(p=p, p_legend=p_legend)
  return(p)
}

clone_prevalence_piechart_legend <- function(save_dir, cols_use=NULL){
  
  df <- data.frame(pct_cells = c(100),
                   clone_id = paste0("G", 1))
  # df <- data.frame(pct_cells = c(15, 25, 32, 28),
  #                  clone_id = paste0("G", 1:4))
  # cols_use <- c('red','blue','cyan','pink')
  # names(cols_use) <- df$clone_id
  # Get the positions
  df2 <- df %>% 
    mutate(csum = rev(cumsum(rev(pct_cells))), 
           pos = pct_cells/2 + lead(csum, 1),
           pos = if_else(is.na(pos), pct_cells/2, pos))
  
  
  cols_use <- rep('white', length(unique(df$clone_id)))
  names(cols_use) <- df$clone_id
  p_legend <- ggplot(df, aes(x = "" , y = pct_cells, fill= fct_inorder(clone_id))) + #, fill = fct_inorder(clone_id))
    geom_col(width = 1, color = 1) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = cols_use)+
    theme_void() + 
    theme(legend.position = 'none')
  # p_legend
  saveRDS(p_legend, paste0(save_dir, 'piechart_legend_shape.rds'))
  
  return(p_legend)
}
clone_prevalence_donutchart_legend <- function(save_dir, cols_use=NULL){
  # Hole size
  hsize <- 1.5
  df <- data.frame(pct_cells = c(100),
                   clone_id = paste0("G", 1))
  # df <- data.frame(pct_cells = c(15, 25, 32, 28),
  #                  clone_id = paste0("G", 1:4))
  # cols_use <- c('red','blue','cyan','pink')
  # names(cols_use) <- df$clone_id
  # Get the positions
  df2 <- df %>% 
    mutate(csum = rev(cumsum(rev(pct_cells))), 
           pos = pct_cells/2 + lead(csum, 1),
           pos = if_else(is.na(pos), pct_cells/2, pos))
  
  
  cols_use <- rep('white', length(unique(df$clone_id)))
  names(cols_use) <- df$clone_id
  # p_legend <- ggplot(df, aes(x = "" , y = pct_cells, fill= fct_inorder(clone_id))) + #, fill = fct_inorder(clone_id))
  #   geom_col(width = 1, color = 1) +
  #   coord_polar(theta = "y") +
  #   scale_fill_manual(values = cols_use)+
  #   theme_void() + 
  #   theme(legend.position = 'none')
  # p_legend
  p_legend <- ggplot(df, aes(x = hsize, y = pct_cells, fill = clone_id)) +
    geom_col(color = "black") +
    # geom_text(aes(label = pct_cells),
    #           position = position_stack(vjust = 0.5)) +
    coord_polar(theta = "y") +
    # scale_fill_brewer(palette = "GnBu") +
    scale_fill_manual(values = cols_use) + 
    xlim(c(0.2, hsize + 0.5)) +
    theme(panel.background = element_rect(fill = "transparent", colour = NA),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.position = 'none')
  # p_legend
  saveRDS(p_legend, paste0(save_dir, 'donutchart_legend_shape.rds'))  
  
  return(p_legend)
}

