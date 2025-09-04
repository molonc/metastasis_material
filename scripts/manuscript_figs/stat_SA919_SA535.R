suppressPackageStartupMessages({
  library("dplyr")
  library("ggplot2")
})

## Loading utility functions
script_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/scripts/'
source(paste0(script_dir, 'corrupt_tree/src/stat/prevalence_utils.R'))
source(paste0(script_dir, 'manuscript_figs/cnv_viz_utils.R'))

stat_SA535_Pt2 <- function(){
  input_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/dlp_trees/'
  datatag <- 'SA535'
  save_fig_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/figures/fig_cnv/'
  cellclones_fn <- paste0(input_dir, datatag, '/', 'cell_clones.csv.gz')
  cell_clones <- data.table::fread(cellclones_fn)
  excluded_clones <- c('Un','unassigned','None')
  cell_clones <- cell_clones[!cell_clones$clone_id %in% excluded_clones,]
  print(dim(cell_clones))
  print(summary(as.factor(cell_clones$clone_id)))
  
  library_grouping_fn <- paste0(input_dir, datatag, '/' ,'library_groupings.csv.gz')
  meta_data <- get_meta_data(cell_clones, results_dir, library_grouping_fn)
  print(dim(meta_data))
  colnames(meta_data)
  View(head(meta_data))
  
  # We classified 3861 high quality single cell genomes, with sample size range from XX to XX cells, 
  # with median of XX single cells per sample (Supplementary Table 4, Figure 3A) 
  # from 4 primary and 6 metastatic tumors (origin from axillary, inguinal, left and right axillary)
  unique(meta_data$origin)
  prevalence_df <- meta_data %>%
    dplyr::filter(origin!='Tumor_Recur') #%>%
    # dplyr::group_by(clone_id, pdxid, mainsite) %>%
    # dplyr::summarise(nb_cells=n())
  unique(prevalence_df$sample_id)
  t <- prevalence_df %>%
    dplyr::filter(origin!='Tumor_Recur') %>%
  dplyr::group_by(sample_id, mainsite) %>%
  dplyr::summarise(nb_cells=n())
  
  #Phylogenetic analysis revealed the existence of 
  # NN-YY CNA clones in each primary passage, 
  # with no single clone occupying >75% cells
  meta_data$pdxid
  t <- meta_data %>%
    dplyr::filter(origin!='Tumor_Recur' & mainsite=='Primary') %>%
    dplyr::group_by(sample_id, clone_id, pdxid) %>%
    dplyr::summarise(nb_cells=n()) %>%
    dplyr::group_by(sample_id, pdxid) %>%
    dplyr::summarise(nb_clones=n())
  t
  
  unique(meta_data$pdxid)
  meta_data$pdxid <- gsub('_','-',meta_data$pdxid)
  
  prevalence_df <- meta_data %>%
    dplyr::filter(origin!='Tumor_Recur' & mainsite=='Primary') %>%
    dplyr::group_by(sample_id, clone_id, pdxid) %>%
    dplyr::summarise(nb_cells=n())
  
  prevalence_df <- meta_data %>%
    dplyr::filter(mainsite %in% c('Metastasis','Primary') & origin!='Tumor_Recur') %>%
    dplyr::group_by(sample_id, clone_id, pdxid) %>%
    dplyr::summarise(nb_cells=n())
  
  # for table 5
  prevalence_df <- meta_data %>%
    dplyr::mutate(mouse_id = stringr::str_sub(pdxid, nchar(pdxid), nchar(pdxid))) %>%
    dplyr::filter(mainsite %in% c('Metastasis','Primary') & origin!='Tumor_Recur') %>%
    dplyr::mutate(sample_id=paste0(sample_id, '_mouseM',mouse_id, '_',mainsite)) %>%
    dplyr::group_by(sample_id, clone_id) %>%
    dplyr::summarise(nb_cells=n())
  dim(prevalence_df)
  length(unique(prevalence_df$sample_id))
  total_df <- prevalence_df %>%
    dplyr::group_by(sample_id) %>% #, mainsite
    dplyr::summarise(nb_cells_per_sample=sum(nb_cells)) 
  
  prevalence_df <- prevalence_df %>%
    dplyr::left_join(total_df, by='sample_id') %>%
    dplyr::mutate(pct_cells=round(100*nb_cells/nb_cells_per_sample,1))
  prevalence_df$nb_cells_per_sample <- NULL
  
  prevalence_df <- prevalence_df %>%
    dplyr::mutate(series=stringr::str_sub(sample_id, 1, 5))
  prevalence_df$patient_ID <- 'Pt2'
  View(head(prevalence_df))
  max(prevalence_df$pct_cells)
  unique(prevalence_df$series)
  dim(prevalence_df)
  t <- prevalence_df %>%
    dplyr::filter(pct_cells>=50)
  View(t)
  t <- prevalence_df %>%
    dplyr::filter(clone_id %in% c('E','H') & pdxid=='X0011_2361')
  # Metastatic clones exhibited a lower degree 
  # of polyclonality with 1-2 major clones per metastasis
  
  save_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/dlp_trees/dlp_met_proj_Sean/'
  save_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/supp_tables/'
  data.table::fwrite(prevalence_df, paste0(save_dir, 'SuppTable5_',datatag, '_cells_proportions.csv'))
  data.table::fwrite(prevalence_df, paste0(save_dir, datatag, '_cells_proportions.csv'))
  
}


stat_SA919_Pt1 <- function(){
  input_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/dlp_trees/'
  datatag <- 'SA919'
  save_fig_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/figures/fig_cnv/'
  cellclones_fn <- paste0(input_dir, datatag, '/', 'cell_clones.csv.gz')
  cell_clones <- data.table::fread(cellclones_fn)
  excluded_clones <- c('Un','unassigned','None')
  cell_clones <- cell_clones[!cell_clones$clone_id %in% excluded_clones,]
  print(dim(cell_clones))
  print(summary(as.factor(cell_clones$clone_id)))
  
  library_grouping_fn <- paste0(input_dir, datatag, '/' ,'library_groupings.csv.gz')
  meta_data <- get_meta_data(cell_clones, results_dir, library_grouping_fn)
  print(dim(meta_data))
  colnames(meta_data)
  # View(head(meta_data))
  unique(meta_data$sample_id)
  summary(as.factor(meta_data$mainsite))
  meta_data <- meta_data %>%
    dplyr::mutate(passage=stringr::str_sub(sample_id,6,7)) 
  
  
  
  prevalence_df <- meta_data %>%
    dplyr::mutate(passage=stringr::str_sub(sample_id,6,7)) %>%
    dplyr::group_by(clone_id, passage) %>% #, mainsite
    dplyr::summarise(nb_cells=n()) #%>%
    # dplyr::ungroup() %>%
  total_df <- prevalence_df %>%
    dplyr::group_by(passage) %>% #, mainsite
    dplyr::summarise(nb_cells_per_passage=sum(nb_cells)) 
  prevalence_df <- prevalence_df %>%
    dplyr::left_join(total_df, by='passage') %>%
    dplyr::mutate(pct_cells=round(100*nb_cells/nb_cells_per_passage,1))
  View(prevalence_df)
  
  # Dominant clone in each sample 
  prevalence_df <- meta_data %>%
    dplyr::filter(origin!='Tumor_Recur') %>%
    dplyr::mutate(passage=stringr::str_sub(sample_id,6,7)) %>%
    dplyr::mutate(mouse_id = stringr::str_sub(pdxid, nchar(pdxid)-1, nchar(pdxid))) %>%
    dplyr::mutate(sample_id=paste0(sample_id, '_mouseM',mouse_id, '_',mainsite,'_',passage)) %>%
    dplyr::group_by(sample_id, clone_id) %>% #, mainsite
    dplyr::summarise(nb_cells=n()) #%>%
  
  total_df <- prevalence_df %>%
    dplyr::group_by(sample_id) %>% #, mainsite
    dplyr::summarise(nb_cells_per_sample=sum(nb_cells)) 
  
  prevalence_df <- prevalence_df %>%
    dplyr::left_join(total_df, by='sample_id') %>%
    dplyr::mutate(pct_cells=round(100*nb_cells/nb_cells_per_sample,1))
  dim(prevalence_df)
  View(prevalence_df)
  length(unique(prevalence_df$sample_id))
  
  prevalence_df <- prevalence_df %>%
    dplyr::mutate(series=stringr::str_sub(sample_id, 1, 5))
  unique(prevalence_df$series)
  prevalence_df$patient_ID <- 'Pt1'
  prevalence_df$nb_cells_per_sample <- NULL
  dim(prevalence_df)
  View(head(prevalence_df))
  save_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/supp_tables/'
  data.table::fwrite(prevalence_df, paste0(save_dir, 'SuppTable5_',datatag, '_cells_proportions.csv'))
  
  
  # prevalence_df <- meta_data %>%
  #   dplyr::group_by(clone_id, pdxid) %>% #, mainsite
  #   dplyr::summarise(nb_cells=n()) %>%
  #   dplyr::ungroup() %>%
  #   dplyr::group_by(pdxid) %>% #, mainsite
  #   dplyr::summarise(nb_cells_per_mouse=n()) %>%
  #   dplyr::mutate()
  # View(prevalence_df)
  prevalence_df <- meta_data %>%
    dplyr::group_by(pdxid, sample_id) %>% #, mainsite
    dplyr::summarise(nb_cells_per_sample=n())
  
  median(prevalence_df$nb_cells_per_sample)
  round(sd(prevalence_df$nb_cells_per_sample), 1)
  
  # In the X3 primary passage, clone A dominated at the primary site (X% A, Y%B, Z% C) 
  # and no metastases were found in any of the 3 primary site replicate transplants
  unique(meta_data$mainsite)
  prevalence_df <- meta_data %>%
    dplyr::filter(mainsite=="Primary" & passage=="X3")%>%
    dplyr::group_by(clone_id) %>% #, mainsite
    dplyr::summarise(nb_cells=n()) #%>%
  
  # In primary passage X4 (X% A, Y%B, Z% C), 
  # metastatic sites were observed in 50% of transplants(2 of 4 replicate transplants)
  prevalence_df <- meta_data %>%
    dplyr::filter(mainsite=="Primary" & passage=="X4")%>%
    dplyr::group_by(clone_id) %>% #, mainsite
    dplyr::summarise(nb_cells=n()) #%>%
  total_cells <- sum(prevalence_df$nb_cells)
  prevalence_df <- prevalence_df %>%
    dplyr::mutate(pct_cells=round(100*nb_cells/total_cells,1))
  prevalence_df
  
  save_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/dlp_trees/dlp_met_proj_Sean/'
  data.table::fwrite(prevalence_df, paste0(save_dir, datatag, '_cells_proportions.csv'))
  # We noted the composition of primary site clones had evolved by X7, 
  # to contain a majority of clone B in all 4 primary sites, 
  # with a minor contribution of clone C (percentages, Figure 3C)
  prevalence_df <- meta_data %>%
    dplyr::filter(mainsite=="Primary" & passage=="X7")%>%
    dplyr::group_by(clone_id, sample_id) %>% #, mainsite
    dplyr::summarise(nb_cells=n()) #%>%
  total_df <- prevalence_df %>%
    dplyr::group_by(sample_id) %>% #, mainsite
    dplyr::summarise(nb_cells_per_sample=sum(nb_cells)) 
  
  prevalence_df <- prevalence_df %>%
    dplyr::left_join(total_df, by='sample_id') %>%
    dplyr::mutate(pct_cells=round(100*nb_cells/nb_cells_per_sample,1))
  View(prevalence_df)
  stat1 <- prevalence_df %>%
    dplyr::filter(clone_id=='B')
  round(mean(stat1$pct_cells),1)
  round(sd(stat1$pct_cells),1)
  
  stat2 <- prevalence_df %>%
    dplyr::filter(clone_id=='C')
  round(mean(stat2$pct_cells),1)
  round(sd(stat2$pct_cells),1)
  
  stat3 <- prevalence_df %>%
    dplyr::filter(clone_id=='A')
  round(mean(stat3$pct_cells),1)
  round(sd(stat3$pct_cells),1)
  
    
  # Strikingly, in two primary transplants (mouse M2, mouse M3) 
  # with minor proportions of clone C (1.1%, 4.1%), 
  # the metastases were dominated by clone C (Z% C) with only a minor contribution of clone B (x%) 
  # and no clone A at supra-spinal and ventral-spinal metastatic sites.
  unique(meta_data$pdxid)
  prevalence_df <- meta_data %>%
    dplyr::filter(mainsite=="Metastasis" & passage=="X7" &
                  origin %in% c("VentralSpinal","SupraSpinal") &
                  pdxid %in% c("X08472112253","X08472112252"))%>%
    dplyr::group_by(clone_id, sample_id) %>% #, mainsite
    dplyr::summarise(nb_cells=n())
  total_df <- prevalence_df %>%
    dplyr::group_by(sample_id) %>% #, mainsite
    dplyr::summarise(nb_cells_per_sample=sum(nb_cells)) 
  
  prevalence_df <- prevalence_df %>%
    dplyr::left_join(total_df, by='sample_id') %>%
    dplyr::mutate(pct_cells=round(100*nb_cells/nb_cells_per_sample,1))
  View(prevalence_df)
  
  
  ## Figure 3, SA919, panel A
  base_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/'
  datatag <- 'SA919'
  save_dir <- paste0(base_dir, 'figures/fig_cnv/')
  copynumber_fn <- paste0(base_dir,'materials/dlp_trees/',datatag,'/total_merged_filtered_states.csv.gz')
  cellclone_fn <- paste0(base_dir,'materials/dlp_trees/', datatag, '/','cell_clones.csv.gz')
  library_grouping_fn <- paste0(base_dir,'materials/dlp_trees/', datatag, '/','library_groupings.csv.gz')
  # color_codes_df
  df_cnv <- get_median_genotype_v3(copynumber_fn, datatag, save_dir, cellclone_fn, library_grouping_fn) 
  dim(df_cnv)
  head(df_cnv)  
  df_cnv <- df_cnv %>%
    dplyr::mutate(chr_desc=paste0(chr, '_', start, '_', end))
  
  obs_chrs <- c('7')
  obs_clones <- c('A','B')
  statBA <- get_CNA_change_regions_stat(df_cnv, obs_clones, obs_chrs)
  
  obs_chrs <- c('7')
  obs_clones <- c('A','C')
  statCA <- get_CNA_change_regions_stat(df_cnv, obs_clones, obs_chrs)
  
  obs_chrs <- c('5')
  obs_clones <- c('A','C')
  statCA <- get_CNA_change_regions_stat(df_cnv, obs_clones, obs_chrs)
  
  obs_chrs <- c('10')
  obs_clones <- c('A','C')
  statCA <- get_CNA_change_regions_stat(df_cnv, obs_clones, obs_chrs)
  
  
  
} 
create_SuppTable5 <- function(){
  save_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/supp_tables/'
  df1 <- data.table::fread(paste0(save_dir, 'SuppTable5_SA919_cells_proportions.csv'))
  df2 <- data.table::fread(paste0(save_dir, 'SuppTable5_SA535_cells_proportions.csv'))
  colnames(df2)
  df <- dplyr::bind_rows(df1, df2)
  data.table::fwrite(df, paste0(save_dir, 'SuppTable5_cells_proportions_3_Sept_2025.csv'))
  
  df1 <- data.table::fread(paste0(save_dir, 'SuppTable6_mixing_exp_SA919_results_Hakwoo_annotation.csv'))
  paste(colnames(df1), collapse = ', ')
  View(df1)
}
# obs_chrs <- c('7')
# obs_clones <- c('A','B')
get_CNA_change_regions_stat <- function(df_cnv, obs_clones, obs_chrs){
  
  cnv_vals <- df_cnv %>%
    dplyr::filter(clone %in% obs_clones) %>%
    # dplyr::select(clone, cnv, chr_desc) %>%
    tidyr::pivot_wider(names_from='clone', values_from = 'cnv') %>%
    dplyr::mutate(var_region=
                    case_when(
                      !!sym(obs_clones[2])-!!sym(obs_clones[1]) > 0 ~ 'amp',
                      !!sym(obs_clones[2])-!!sym(obs_clones[1]) < 0 ~ 'del',
                      TRUE ~ 'no_variance'
                    ))
  
  # dim(cnv_vals)
  # colnames(cnv_vals)
  # head(cnv_vals)
  # summary(as.factor(cnv_vals$var_region))
  
  # rv <- DelayedArray::rowVars(as.matrix(cnv_vals))  # median copy number profile of obs clones
  # # rv <- rowVars(as.matrix(cnv[,c(5, 6)]))  # median copy number profile of clone B, C 
  # genes_used <- cnv$gene_symbol[rv>0]
  # sum(rv>0)
  
  total_tmp <- tibble::tibble()
  for(obs_chr in obs_chrs){
    print(obs_chr)
    tmp <- cnv_vals %>%
      dplyr::filter(chr %in% obs_chr)
    print(paste0(obs_clones[1],', median value is: ',median(tmp[[obs_clones[1]]])))
    print(paste0(obs_clones[2],', median value is: ',median(tmp[[obs_clones[2]]])))
    
    obs_chr_df <- tmp %>%
      # dplyr::filter(chr %in% obs_chrs) %>%
      dplyr::group_by(var_region) %>%
      dplyr::summarise(nb_regions=n(), 
                       median_cnv_clone1=median(!!sym(obs_clones[1])),
                       median_cnv_clone2=median(!!sym(obs_clones[2])))
    
    total_regions <- sum(obs_chr_df$nb_regions)
    obs_chr_df$pct_regions <- round(100*obs_chr_df$nb_regions/total_regions,1)
    print(obs_chr_df)
    obs_chr_df$chr <- obs_chr
    total_tmp <- dplyr::bind_rows(total_tmp, obs_chr_df)
  }
  return(total_tmp)
}
