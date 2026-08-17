suppressPackageStartupMessages({
  library("dplyr")
  library("ggplot2")
  library("DESeq2")
  library("SingleCellExperiment")
})

script_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/'
source(paste0(script_dir, 'scripts/bulk_rna/bulk_utils.R'))
source(paste0(script_dir, 'scripts/bulk_rna/DE_analysis_utils.R'))
source(paste0(script_dir, 'scripts/bulk_rna/bulk_pathway_utils.R'))

# B met vs. A pri main-main
# B pri vs. A pri main-main

# B met vs. B pri main-main
# C met vs. B pri main-main

# A met vs. A pri main-main
# C met vs. A pri main-main

# C met vs. C pri mixing
# C met vs. A pri mixing
# C pri vs. A pri mixing

## DE analysis, and pathway analysis
## Noted the list of samples for comparisons


## get DE analysis, pathway analysis, cis trans analysis in one function
process_data <- function(dds, meta_genes, obs_clones, filter_conds, save_dir, exp_conds=c('main_exp','main_exp')){
  # obs_clones <- c('A','B')
  # filter_conds <- c('pri','pri')
  subtag <- paste0(obs_clones[2],'_',filter_conds[2],'_',exp_conds[2],'__',obs_clones[1],'_',filter_conds[1],'_',exp_conds[1])
  print(subtag)
  save_fig_dir <- paste0(save_dir, subtag, '/')
  if(!dir.exists(save_fig_dir)){
    dir.create(save_fig_dir)  
  }
  
  
  res <- get_DE_results(subtag, save_fig_dir, dds, meta_genes, obs_clones, 
                        exps=exp_conds, filter_conds=filter_conds)
  obs_genes_symb <- unique(res$symbol)
  print(length(obs_genes_symb))
  pw_res <- get_gprofiler_pathways_obsgenes_v2(obs_genes_symb, save_fig_dir, subtag, correction_method='gSCS',
                                               custom_id="gp__RYib_WEFZ_gEw", pathway_fn=NULL, save_data=T)
  # sum(pw_res$stat$p_value<0.05)
  
  
  
  # # get cis genes, check degree in primary, and in met samples. 
  # obs_clones <- c('A','C')
  # cnv <- get_cis_trans_gene_type(obs_clones)
  # dim(cnv)
  
  
}
double_check <- function(){
  save_dir <- "/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/results_bulkRNAseq/SA919_full/"
  bb <- data.table::fread(paste0(save_dir, 'Bmet_Bpri/Bmet_Bpri_DE_genes.csv.gz'))
  cb1 <- data.table::fread(paste0(save_dir, 'Cmet_Bpri/Cmet_Bpri_DE_genes.csv.gz'))
  dim(cb1)
  head(cb1)
  dim(bb)
  bb <- bb %>%
    filter(pvalue<0.05 & abs(log2FoldChange)>=1)
  obs_genes_symb <- unique(cb1$symbol)
  pw_res <- get_gprofiler_pathways_obsgenes_v2(obs_genes_symb, save_fig_dir, subtag, correction_method='gSCS',
                                               custom_id="gp__RYib_WEFZ_gEw", pathway_fn=NULL, save_data=F)
  sum(pw_res$stat$p_value<0.05)
  
}

extract_DE_genes_DESeq2_revision_v2 <- function(){
  save_dir <- "/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/results_bulkRNAseq/SA919_full/"
  datatag <- 'SA919_full'
  sce <- readRDS(paste0(save_dir, datatag, '_sce.rds'))
  # print(head(assay(sce), 3))
  
  dds <- create_DESeq2_object(sce)
  meta_genes <- rowData(sce) %>% as.data.frame()
  dim(meta_genes)
  head(meta_genes)
  script_dir <- "/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/"
  cnv <- data.table::fread(paste0(script_dir, 'SA919_DE_analysis_DESeq2_Hoa_09April2024/mapping_gene_cnv_SA919.csv.gz'))  
  dim(cnv)
  cnv <- cnv[!duplicated(cnv$ensembl_gene_id),]
  cols_use <- colnames(cnv)[!colnames(cnv) %in% c('gene_symbol','chr')]
  # head(cnv)
  cnv <- cnv %>%
    dplyr::select(all_of(cols_use)) %>%
    dplyr::rename(ens_gene_id=ensembl_gene_id)
  sum(meta_genes$ens_gene_id %in% cnv$ens_gene_id)
  meta_genes <- meta_genes %>%
    dplyr::left_join(cnv, by=c('ens_gene_id'))
  data.table::fwrite(meta_genes, paste0(script_dir, 'SA919_DE_analysis_DESeq2_Hoa_09April2024/meta_genes_SA919.csv.gz'))  
    
  metasamples <- tibble::tibble()
  ## B pri vs. A pri, main experiment
  obs_clones <- c('A','B')
  filter_conds <- c('Primary','Primary')
  exp_conds=c('main_exp','main_exp')
  process_data(dds, meta_genes, obs_clones, filter_conds, save_dir, exp_conds)
  ms <- tibble::tibble(obs_clone1=obs_clones[1], obs_clone2=obs_clones[2],
                       filter_cond1=filter_conds[1],filter_cond2=filter_conds[2],
                       exp_cond1=exp_conds[1],exp_cond2=exp_conds[2])
  metasamples <- dplyr::bind_rows(metasamples, ms)
  metasamples <- metasamples %>%
    dplyr::mutate(desc=paste0(obs_clone2,'_',filter_cond2,'_',exp_cond2,'__',
                              obs_clone1,'_',filter_cond1,'_',exp_cond1))
  
  # B met vs. A pri main-main
  obs_clones <- c('A','B')
  filter_conds <- c('Primary','Metastasis')
  exp_conds=c('main_exp','main_exp')
  process_data(dds, meta_genes, obs_clones, filter_conds, save_dir, exp_conds)
 
  
  # B met vs. B pri main-main
  obs_clones <- c('B','B')
  filter_conds <- c('Primary','Metastasis')
  exp_conds=c('main_exp','main_exp')
  process_data(dds, meta_genes, obs_clones, filter_conds, save_dir, exp_conds)
 
  
  # C met vs. B pri main-main
  ## note: only one sample from C at main exp, met condition
  obs_clones <- c('B','C')
  filter_conds <- c('Primary','Metastasis')
  exp_conds=c('main_exp','main_exp')
  process_data(dds, meta_genes, obs_clones, filter_conds, save_dir, exp_conds)
  
  # C met vs. B pri all main, mixing
  ## note: only one sample from C at main exp, met condition
  obs_clones <- c('B','C')
  filter_conds <- c('Primary','Metastasis')
  exp_conds=c('both','both')
  process_data(dds, meta_genes, obs_clones, filter_conds, save_dir, exp_conds)
  
  
  # A met vs. A pri main-main
  obs_clones <- c('A','A')
  filter_conds <- c('Primary','Metastasis')
  exp_conds=c('main_exp','main_exp')
  process_data(dds, meta_genes, obs_clones, filter_conds, save_dir, exp_conds)
  
  # C met vs. A pri main-main
  ## note: only one sample from C at main exp, met condition
  obs_clones <- c('A','C')
  filter_conds <- c('Primary','Metastasis')
  exp_conds=c('main_exp','main_exp')
  process_data(dds, meta_genes, obs_clones, filter_conds, save_dir, exp_conds)
  
  
  # C met vs. A pri all main, mixing
  ## note: only one sample from C at main exp, met condition
  obs_clones <- c('A','C')
  filter_conds <- c('Primary','Metastasis')
  exp_conds=c('both','both')
  process_data(dds, meta_genes, obs_clones, filter_conds, save_dir, exp_conds)
  
  # C met vs. C pri mixing
  ## note: only one sample from C at mixing exp, pri condition
  obs_clones <- c('C','C')
  filter_conds <- c('Primary','Metastasis')
  exp_conds=c('mixing_exp','mixing_exp')
  process_data(dds, meta_genes, obs_clones, filter_conds, save_dir, exp_conds)
  
  # C met vs. A pri mixing
  obs_clones <- c('A','C')
  filter_conds <- c('Primary','Metastasis')
  exp_conds=c('mixing_exp','mixing_exp')
  process_data(dds, meta_genes, obs_clones, filter_conds, save_dir, exp_conds)
  
  
  # C pri vs. A pri mixing
  obs_clones <- c('A','C')
  filter_conds <- c('Primary','Primary')
  exp_conds=c('mixing_exp','mixing_exp')
  process_data(dds, meta_genes, obs_clones, filter_conds, save_dir, exp_conds)
  
  # View(metasamples)
  data.table::fwrite(metasamples, paste0(save_dir, 'list_DE_comparisons.csv'))
}  

get_summary_exp <- function(){
  save_dir <- "/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/results_bulkRNAseq/SA919_full/"
  # comp <- c('A_Metastasis_main_exp__A_Primary_main_exp','B_Metastasis_main_exp__A_Primary_main_exp',
  #   'B_Metastasis_main_exp__B_Primary_main_exp','B_Primary_main_exp__A_Primary_main_exp',
  #   'C_Metastasis_main_exp__A_Primary_main_exp','C_Metastasis_main_exp__B_Primary_main_exp',
  #   'C_Metastasis_mixing_exp__A_Primary_mixing_exp','C_Metastasis_mixing_exp__C_Primary_mixing_exp',
  #   'C_Metastasis_both__A_Primary_both','C_Metastasis_both__B_Primary_both')
  
  
  de_comps <- data.table::fread(paste0(save_dir,'list_DE_comparisons.csv'))
  de_comps <- de_comps %>%
    dplyr::filter(exp_cond1=='main_exp' & exp_cond2=='main_exp')
  de_comps
  View(de_comps)
  dim(de_comps)
  
  stat_total_df <- tibble::tibble()
  for(de in de_comps$desc){
    print(de)
    de_fn <- paste0(save_dir, de,'/',de,'_DE_genes.csv.gz')
    if(file.exists(de_fn)){
      comp <- de_comps %>%
        dplyr::filter(desc==de)
      df <- data.table::fread(de_fn)
      print(dim(df))
      if(comp$obs_clone1==comp$obs_clone2){
        df <- df  %>%
          dplyr::mutate(gene_type='trans_gene')
      }else{
        df <- df  %>%
          dplyr::mutate(
            gene_type = case_when(
              !!sym(comp$obs_clone1)-!!sym(comp$obs_clone2) != 0 ~ 'cis_gene',
              TRUE ~ 'trans_gene'
            )
          )  
      }
      
      print(summary(as.factor(df$gene_type)))
      data.table::fwrite(df, paste0(save_dir, de,'/',de,'_DE_genes_gene_type.csv.gz'))
      df_stat <- df  %>%
        dplyr::group_by(gene_type)  %>%
        dplyr::summarise(nb_genes=n(),
                         median_log2FC=round(median(abs(log2FoldChange)),2),
                         mean_log2FC=round(mean(abs(log2FoldChange)),2),
                         sd_log2FC=round(sd(abs(log2FoldChange)),2))
      df_stat$desc <- de  
      stat_total_df <- dplyr::bind_rows(stat_total_df, df_stat)  
    }
    
  }
  data.table::fwrite(stat_total_df, paste0(save_dir,'stat_total_revision.csv'))
  View(stat_total_df)  
  
  # for(de in de_comps$desc){
  #   print(de)
  #   de_fn <- paste0(save_dir, de,'/',de,'_DE_genes.csv.gz')
  #   if(file.exists(de_fn)){
  #     print('yes')
  #   }else{
  #     print('no')
  #   }
  # }  
  # dim(de_comps)
  # de_comps <- de_comps %>%
  #   dplyr::filter(desc %in% de_comps$desc[1:10])
  de_comps$obs_clone1
  obs_gene_type <- 'cis_gene'
  de_comps_cis <- de_comps %>%
      dplyr::filter(obs_clone1!=obs_clone2)
  ncomp <- dim(de_comps_cis)[1]
  common_mtx <- matrix(NA, nrow = ncomp, ncol = ncomp)
  for (j in seq(ncomp)) {
    for (i in seq(ncomp)) {
      if (i<j) {  
        de_fni <- paste0(save_dir, de_comps_cis$desc[i],'/',de_comps_cis$desc[i],'_DE_genes_gene_type.csv.gz')
        de_fnj <- paste0(save_dir, de_comps_cis$desc[j],'/',de_comps_cis$desc[j],'_DE_genes_gene_type.csv.gz')
        if(file.exists(de_fni) & file.exists(de_fnj)){
          dfi <- data.table::fread(de_fni)
          dfj <- data.table::fread(de_fnj)
          dfi <- dfi %>%
            dplyr::filter(gene_type==obs_gene_type)
          dfj <- dfj %>%
            dplyr::filter(gene_type==obs_gene_type)
          common_genes <- intersect(dfi$ensembl_gene_id, dfj$ensembl_gene_id)
          common_mtx[j,i] <- length(common_genes)
        }  
      }
    }
  }  
  rownames(common_mtx) <- de_comps_cis$desc
  colnames(common_mtx) <- de_comps_cis$desc
  View(common_mtx)
  
  data.table::fwrite(common_mtx, paste0(save_dir,'stat_common_cis_revision.csv'), row.names = T)
  
  
  obs_gene_type <- 'trans_gene'
  ncomp <- dim(de_comps)[1]
  common_mtx <- matrix(NA, nrow = ncomp, ncol = ncomp)
  for (j in seq(ncomp)) {
    for (i in seq(ncomp)) {
      if (i<j) {  
        de_fni <- paste0(save_dir, de_comps$desc[i],'/',de_comps$desc[i],'_DE_genes_gene_type.csv.gz')
        de_fnj <- paste0(save_dir, de_comps$desc[j],'/',de_comps$desc[j],'_DE_genes_gene_type.csv.gz')
        if(file.exists(de_fni) & file.exists(de_fnj)){
          dfi <- data.table::fread(de_fni)
          dfj <- data.table::fread(de_fnj)
          dfi <- dfi %>%
            dplyr::filter(gene_type==obs_gene_type)
          dfj <- dfj %>%
            dplyr::filter(gene_type==obs_gene_type)
          common_genes <- intersect(dfi$ensembl_gene_id, dfj$ensembl_gene_id)
          common_mtx[j,i] <- length(common_genes)
        }  
      }
    }
  }  
  rownames(common_mtx) <- de_comps$desc
  colnames(common_mtx) <- de_comps$desc
  # View(common_mtx)
  data.table::fwrite(common_mtx, paste0(save_dir,'stat_common_trans_revision.csv'), row.names = T)
}

viz_upsetR_annotated_genotypes <- function(genes_df, save_fig_dir, tag=''){
  save_dir <- "/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/results_bulkRNAseq/SA919_full/"
  save_fig_dir <- paste0(save_dir,'revision/')
  dir.create(save_fig_dir)
  
  ## Visualizing cis genes upSetR
  total_df <- tibble::tibble()
  obs_gene_type <- 'cis_gene'
  for(de in de_comps$desc){
    print(de)
    de_fn <- paste0(save_dir, de,'/',de,'_DE_genes_gene_type.csv.gz')
    if(file.exists(de_fn)){
      df <- data.table::fread(de_fn)
      df <- df %>%
        dplyr::filter(gene_type==obs_gene_type)
      print(dim(df))
      df$desc <- de
      total_df <- dplyr::bind_rows(total_df, df)  
    }
  }
  
  dim(total_df)
  genes_df <- total_df %>%
    dplyr::mutate(desc=gsub('_main_exp','',desc))
  genes_df <- genes_df %>%
    dplyr::mutate(desc=gsub('__',' vs ',desc))
  library(UpSetR)
  # Get all unique elements across vectors
  print(dim(genes_df))
  all_elements <- unique(genes_df$ensembl_gene_id)
  print(length(all_elements))
  # Build a binary membership data frame (1 = present, 0 = absent)
  df <- data.frame(element = all_elements)
  lbs <- c()
  for(gt in unique(genes_df$desc)){
    # lb <- strsplit(gt,'_')[[1]][1]
    lb <- gt
    lbs <- c(lbs, lb)
    tmp <- genes_df %>%
      dplyr::filter(desc==gt)
    df[,lb] <- as.integer(all_elements %in% tmp$ensembl_gene_id)
  }
  print(lbs)
  # print(head(df))
  dim(df)
  ordered_lbs <- c("B_Primary vs A_Primary","B_Metastasis vs A_Primary",
                   "C_Metastasis vs A_Primary","C_Metastasis vs B_Primary")
  # df <- df %>%
  #   dplyr::arrange()
  df <- df[,c("element",ordered_lbs)]
 
  png(paste0(save_fig_dir, obs_gene_type,"_upSetR.png"), height = 2*450, width=2*750, res = 2*72)
  p <- UpSetR::upset(df,
                     sets = ordered_lbs,
                     order.by = "freq",
                     keep.order      = TRUE,           # preserve the order given in gene_list
                     main.bar.color = "steelblue",
                     sets.bar.color = "darkgreen",
                     text.scale = 1.5,
                     point.size = 3,
                     line.size = 1.2,
                     mainbar.y.label = "#intersect genes",
                     sets.x.label = "#genes")
  print(p)
  dev.off()
  saveRDS(p, paste0(save_fig_dir, obs_gene_type,"_upSetR.rds"))
  
  
  
  ## Visualizing trans genes upSetR
  total_df <- tibble::tibble()
  obs_gene_type <- 'trans_gene'
  for(de in de_comps$desc){
    print(de)
    de_fn <- paste0(save_dir, de,'/',de,'_DE_genes_gene_type.csv.gz')
    if(file.exists(de_fn)){
      df <- data.table::fread(de_fn)
      df <- df %>%
        dplyr::filter(gene_type==obs_gene_type)
      print(dim(df))
      df$desc <- de
      total_df <- dplyr::bind_rows(total_df, df)  
    }
  }
  
  dim(total_df)
  
  genes_df <- total_df %>%
    dplyr::mutate(desc=gsub('_main_exp','',desc))
  genes_df <- genes_df %>%
    dplyr::mutate(desc=gsub('__',' vs ',desc))
  library(UpSetR)
  # Get all unique elements across vectors
  print(dim(genes_df))
  all_elements <- unique(genes_df$ensembl_gene_id)
  print(length(all_elements))
  # Build a binary membership data frame (1 = present, 0 = absent)
  df <- data.frame(element = all_elements)
  lbs <- c()
  for(gt in unique(genes_df$desc)){
    # lb <- strsplit(gt,'_')[[1]][1]
    lb <- gt
    lbs <- c(lbs, lb)
    tmp <- genes_df %>%
      dplyr::filter(desc==gt)
    df[,lb] <- as.integer(all_elements %in% tmp$ensembl_gene_id)
  }
  print(lbs)
  print(head(df))
  dim(df)
  ordered_lbs <- c("B_Primary vs A_Primary","B_Metastasis vs A_Primary",
                   "A_Metastasis vs A_Primary","C_Metastasis vs A_Primary",
                   "B_Metastasis vs B_Primary","C_Metastasis vs B_Primary")
  
  # df <- df %>%
  #   dplyr::arrange()
  df <- df[,c("element",ordered_lbs)]
  
  png(paste0(save_fig_dir, obs_gene_type,"_upSetR.png"), height = 2*450, width=2*750, res = 2*72)
  p <- UpSetR::upset(df,
                     sets = ordered_lbs,
                     order.by = "freq",
                     keep.order      = TRUE,           # preserve the order given in gene_list
                     main.bar.color = "steelblue",
                     sets.bar.color = "darkgreen",
                     text.scale = 1.5,
                     point.size = 3,
                     line.size = 1.2,
                     mainbar.y.label = "#intersect genes",
                     sets.x.label = "#genes")
  print(p)
  dev.off()
  saveRDS(p, paste0(save_fig_dir, obs_gene_type,"_upSetR.rds"))
  
  
  ## Visualizing intersection genes
  obs_gene_type <- 'cis_gene'
  cis_plt <- readRDS(paste0(save_fig_dir, obs_gene_type,"_upSetR.rds"))

  svg(paste0(save_fig_dir, obs_gene_type,"_upSetR.svg"), height = 5, width = 10)
  cis_plt
  dev.off()
  
  ggsave(
    filename = paste0(save_fig_dir, obs_gene_type,"_upSetR.svg"),
    plot     = grid::grid.grabExpr(print(cis_plt), wrap.grobs = TRUE),
    device   = "svg",     # requires the svglite package
    width    = 10,
    height   = 5,
    units    = "in"
  )
  
  
  
  obs_gene_type <- 'trans_gene'
  trans_plt <- readRDS(paste0(save_fig_dir, obs_gene_type,"_upSetR.rds"))

  ggsave(
    filename = paste0(save_fig_dir, obs_gene_type,"_upSetR.svg"),
    plot     = grid::grid.grabExpr(print(trans_plt), wrap.grobs = TRUE),
    device   = "svg",     # requires the svglite package
    width    = 10,
    height   = 6.5,
    units    = "in"
  )
  
  
  p_total_upSet <- cowplot::plot_grid(grid::grid.grabExpr(print(cis_plt), wrap.grobs = TRUE), 
                                      grid::grid.grabExpr(print(trans_plt), wrap.grobs = TRUE), 
                                      rel_heights = c(4.5, 6), ncol=1)
  #, labels = c('Cis genes','Trans genes')
  ggsave(  
    filename = paste0(save_fig_dir,"intersection_gene_type_upSetR.svg"),  
    plot = p_total_upSet,  
    device   = "svg",     # requires the svglite package
    width    = 10,
    height   = 10,
    units    = "in", dpi = 150)  
  
  
  
}


extract_DE_genes_DESeq2_revision <- function(){
  save_dir <- "/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/results_bulkRNAseq/SA919_full/"
  datatag <- 'SA919_full'
  sce <- readRDS(paste0(save_dir, datatag, '_sce.rds'))
  dds <- create_DESeq2_object(sce)
  meta_genes <- rowData(sce) %>% as.data.frame()
  dim(meta_genes)
  
  
  ## B pri vs. A pri, main experiment
  obs_clones <- c('A','B')
  subtag <- paste0(obs_clones[2],'pri','_',obs_clones[1],'pri')
  print(subtag)
  save_fig_dir <- paste0(save_dir, subtag, '/')
  dir.create(save_fig_dir)
  res <- get_DE_results(subtag, save_fig_dir, dds, meta_genes, obs_clones, exps=c('main_exp'), filter_cond='pri_pri')
  obs_genes_symb <- unique(res$symbol)
  length(obs_genes_symb)
  pw_res <- get_gprofiler_pathways_obsgenes_v2(obs_genes_symb, save_fig_dir, subtag, correction_method='gSCS',
                                               custom_id="gp__RYib_WEFZ_gEw", pathway_fn=NULL, save_data=T)
  sum(pw_res$stat$p_value<0.05)
  
  ## B met vs. A pri , main experiment
  ## See above function extract_DE_genes_DESeq2()
  
  
  ## C pri vs. A pri, mixing experiment
  obs_clones <- c('A','C')
  subtag1 <- paste0(obs_clones[2],'met','_',obs_clones[1],'pri')
  print(subtag1)
  save_fig_dir <- paste0(save_dir, subtag1, '_mixing/')
  dir.create(save_fig_dir)
  res <- get_DE_results(subtag, save_fig_dir, dds, meta_genes, obs_clones, 
                        exps=c('mixing_exp'), filter_cond='met_pri')
  df1 <- data.table::fread(paste0(save_fig_dir, '/', subtag1, '_DE_genes.csv.gz'))
  dim(df1)
  
  obs_clones <- c('A','C')
  subtag2 <- paste0(obs_clones[2],'pri','_',obs_clones[1],'pri')
  print(subtag2)
  save_fig_dir <- paste0(save_dir, subtag2, '_mixing/')
  dir.create(save_fig_dir)
  res <- get_DE_results(subtag2, save_fig_dir, dds, meta_genes, obs_clones, 
                        exps=c('mixing_exp'), filter_cond='pri_pri')
  df2 <- data.table::fread(paste0(save_fig_dir, '/', subtag2, '_DE_genes.csv.gz'))
  dim(df2)
  
  
  
  # get cis genes, check degree in primary, and in met samples. 
  obs_clones <- c('A','C')
  cnv <- get_cis_trans_gene_type(obs_clones)
  dim(cnv)
  # head(cnv)
  ## To Do: need to add cis, trans type into DE_genes table. 
  
  
  ## stat values
  ## Total # genes, # cis genes, # trans genes
  df2_cis_total <- df2 %>%
    dplyr::filter(ensembl_gene_id %in% cnv$ensembl_gene_id)
  df2_trans_total <- df2 %>%
    dplyr::filter(!ensembl_gene_id %in% cnv$ensembl_gene_id)
  
  
  median(abs(df1_cis_total$log2FoldChange))
  mean(abs(df1_cis_total$log2FoldChange))
  sd(abs(df1_cis_total$log2FoldChange))
  
  median(abs(df1_trans_total$log2FoldChange))
}  
cis_trans_logFC <- function(){
  
  ## Get common genes between 2 comparisons
  ## Comparing variance, and gene wise dispersion between metpri, vs pripri
  obs_clones <- c('A','B')
  subtag1 <- paste0(obs_clones[2],'met','_',obs_clones[1],'pri')
  print(subtag1)
  df1 <- data.table::fread(paste0(save_dir, subtag1, '/', subtag1, '_DE_genes.csv.gz'))
  dim(df1)
  
  obs_clones <- c('A','B')
  subtag2 <- paste0(obs_clones[2],'pri','_',obs_clones[1],'pri')
  print(subtag2)
  df2 <- data.table::fread(paste0(save_dir, subtag2, '/', subtag2, '_DE_genes.csv.gz'))
  dim(df2)
  
  
  # get cis genes, check degree in primary, and in met samples. 
  obs_clones <- c('A','B')
  cnv <- get_cis_trans_gene_type(obs_clones)
  dim(cnv)
  head(cnv)
  
  common_genes <- intersect(df1$ensembl_gene_id, df2$ensembl_gene_id)
  cis_common_genes <- intersect(common_genes, cnv$ensembl_gene_id)
  trans_common_genes <- common_genes[!common_genes %in% cis_common_genes]
  
  df2_cis <- df2 %>%
    dplyr::filter(ensembl_gene_id %in% cis_common_genes)
  df2_cis_total <- df2 %>%
    dplyr::filter(ensembl_gene_id %in% cnv$ensembl_gene_id)
  df2_trans <- df2 %>%
    dplyr::filter(ensembl_gene_id %in% trans_common_genes)
  df2_trans_total <- df2 %>%
    dplyr::filter(!ensembl_gene_id %in% cnv$ensembl_gene_id)
  
  df1_cis <- df1 %>%
    dplyr::filter(ensembl_gene_id %in% cis_common_genes)
  df1_cis_total <- df1 %>%
    dplyr::filter(ensembl_gene_id %in% cnv$ensembl_gene_id)
  df1_trans <- df1 %>%
    dplyr::filter(ensembl_gene_id %in% trans_common_genes)
  df1_trans_total <- df1 %>%
    dplyr::filter(!ensembl_gene_id %in% cnv$ensembl_gene_id)
  
  t <- cnv %>%
    dplyr::filter(ensembl_gene_id %in% cis_common_genes)
  paste(t$gene_symbol, collapse = ', ')
  
  
  round(median(abs(df2_cis$log2FoldChange)),2)
  
  median(abs(df1_cis_total$log2FoldChange))
  median(abs(df2_cis_total$log2FoldChange))
  dim(df1_cis_total)
  dim(df2_cis_total)
  median(abs(df1_trans$log2FoldChange))
  median(abs(df2_trans$log2FoldChange))
  
  
  median(abs(df1_trans_total$log2FoldChange))
  median(abs(df2_trans_total$log2FoldChange))
  ## To Do: need gene wise dispersion values here 
  
  dim(df1)
  dim(df2)
  
  sum(df2_extracted$ensembl_gene_id %in% cnv$ensembl_gene_id)
  sum(df2_extracted$symbol %in% cnv$gene_symbol)
  
  
  dim(df2)[1]
  
  
  
  
  
  c1 <- c(subtag1, dim(df1)[1], dim(df1_cis)[1], round(median(abs(df1_cis$log2FoldChange)),2), 
          dim(df1_cis_total)[1], round(median(abs(df1_cis_total$log2FoldChange)),2), 
          dim(df1_trans)[1], round(median(abs(df1_trans$log2FoldChange)),2), 
          dim(df1_trans_total)[1], round(median(abs(df1_trans_total$log2FoldChange)),2))
  
  c2 <- c(subtag2, dim(df2)[1], dim(df2_cis)[1], round(median(abs(df2_cis$log2FoldChange)),2), 
          dim(df2_cis_total)[1], round(median(abs(df2_cis_total$log2FoldChange)),2), 
          dim(df2_trans)[1], round(median(abs(df2_trans$log2FoldChange)),2), 
          dim(df2_trans_total)[1], round(median(abs(df2_trans_total$log2FoldChange)),2))
  col_use <- c('comp','#DEgenes','#common_cis','medianlogFC_cis',
               '#total_cis','medianlogFC_cis',
               '#common_trans','medianlogFC_trans',
               '#total_trans','medianlogFC_trans')
  stat_de1 <- tibble::tibble(c1)
  stat_de2 <- tibble::tibble(c2)
  
  stat_de <- dplyr::bind_cols(stat_de1, stat_de2)
  stat_de <- t(stat_de)
  colnames(stat_de) <- col_use
  # stat_de
  data.table::fwrite(stat_de, paste0(save_dir, subtag2, '_mixing/', 'summary_comparisons.csv'))
  ## For common genes, get cis variance, and trans variance in 2 comparisons
  
  
  
  
}  
