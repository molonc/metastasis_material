## All DLP+ plotting
# library(ggpubr)
# BiocManager::install('ggpubr', ask = F)
base_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/'
source(paste0(base_dir, 'scripts/dlp_manuscript/cnv_viz_utils.R'))
save_dir <- paste0(base_dir, 'figures/fig_cnv/')

## Figure 3, panel A
datatag <- 'SA919'
copynumber_fn <- paste0(base_dir,'materials/dlp_trees/',datatag,'/total_merged_filtered_states.csv.gz')
cellclone_fn <- paste0(base_dir,'materials/dlp_trees/', datatag, '/','cell_clones.csv.gz')
library_grouping_fn <- paste0(base_dir,'materials/dlp_trees/', datatag, '/','library_groupings.csv.gz')
# color_codes_df
df_cnv <- get_median_genotype_v3(copynumber_fn, datatag, save_dir, cellclone_fn, library_grouping_fn) 
dim(df_cnv)
plt_title <- 'Pt1: DLP+ median copy number profiles'
plt_title <- ''
res_919 <- plot_CNV_profile(df_cnv, clones= levels(df_cnv$clone),plttitle=plt_title, meta_genes = NULL, datatag)
res_919$cnv_plot
res_919$plg
# lg <- cowplot::get_legend(res_919$cnv_plot)
# plg <- cowplot::ggdraw() + cowplot::draw_plot(lg)
dev.off()
# leg <- ggpubr::get_legend(res_919$cnv_plot)

# plg <- cowplot::ggdraw() + cowplot::draw_plot(res_919$lg) # do not work at this R version, need ggpubr package func above
main_plot <- cowplot::plot_grid(
  res_919$clone_plt,
  res_919$cnv_plot + theme(legend.position = 'none'),
  nrow = 1,
  rel_widths = c(0.05,1.2),
  align = 'h'#,
  # labels = c('a')
)
main_plot

# main_plot <- cowplot::plot_grid(
#   res_919$cnv_plot,
#   ggpubr::as_ggplot(leg),
#   ncol = 1,
#   rel_heights = c(1.6,0.3),
#   align = 'v'#,
#   # labels = c('a')
# )
saveRDS(res_919$cnv_plot, paste0(save_dir, datatag,'_cnv_profiles.rds'))
ggsave(paste0(save_dir,"Fig3_SA919_Pt1_cnv_profile.svg"),
       plot = main_plot,
       height = 4,
       width = 12,
       # useDingbats=F,
       dpi = 150
)
dev.off()


## Figure 4, panel A
datatag <- 'SA535'
copynumber_fn <- paste0(base_dir,'materials/dlp_trees/',datatag,'/total_merged_filtered_states.csv.gz')
cellclone_fn <- paste0(base_dir,'materials/dlp_trees/', datatag, '/','cell_clones.csv.gz')
library_grouping_fn <- paste0(base_dir,'materials/dlp_trees/', datatag, '/','library_groupings.csv.gz')
# color_codes_df
df_cnv <- get_median_genotype_v3(copynumber_fn, datatag, save_dir,
                                 cellclone_fn, library_grouping_fn) 
dim(df_cnv)

res_535 <- plot_CNV_profile(df_cnv, clones=unique(df_cnv$clone), plttitle=' ', meta_genes = NULL, datatag)
res_535$cnv_plot
res_535$plg
dev.off()
leg <- ggpubr::get_legend(res_535$cnv_plot)
ggpubr::as_ggplot(leg)

main_plot <- cowplot::plot_grid(
  res_535$clone_plt,
  res_535$cnv_plot + theme(legend.position = 'none'),
  nrow = 1,
  rel_widths = c(0.05,1.2),
  align = 'h'#,
  # labels = c('a')
)
main_plot
ggsave(paste0(save_dir,"Fig4_SA535_Pt2_cnv_profile_panelA.svg"),
       plot = main_plot,
       height = 4,
       width = 12,
       # useDingbats=F,
       dpi = 150
)
dev.off()




base_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/'
source(paste0(base_dir, 'scripts/manuscript_figs/plot_utils.R'))
save_dir <- paste0(base_dir, 'figures/fig_cnv/')

clones <- c('B','C')
datatag <- 'SA919'

## What kind of input data here? 
# df_cnv_fn <- paste0(save_dir, datatag, '_cnv_mat.csv.gz')
script_dir <- "/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/"
df_cnv_fn <- paste0(script_dir, 'SA919_DE_analysis_DESeq2_Hoa_09April2024/mapping_gene_cnv_SA919.csv.gz')
# df_cnv <- data.table::fread(df_cnv_fn)
# dim(df_cnv)
# cnv <- data.table::fread(paste0(script_dir, 'SA919_DE_analysis_DESeq2_Hoa_09April2024/mapping_gene_cnv_SA919.csv.gz'))  
# cnv <- cnv[!duplicated(cnv$ensembl_gene_id),]
# dim(cnv)

# a list of filtered genes, if plot all genes --> very dense area, can not see small regions with CN change
input_cnv_dir <- "/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/results_bulkRNAseq/"
datatag <- 'SA919_full'
# sce <- readRDS(paste0(input_cnv_dir, datatag,'/',datatag,'_sce.rds'))
# dim(sce)
# sce
# meta_genes <- data.table::fread(paste0(save_dir,datatag, '_meta_filtered_genes.csv')) %>% as.data.frame()
# dim(meta_genes)
# meta_genes <- rownames(sce)
meta_genes <- get_obs_genes(obs_clones)
length(meta_genes)
pcnv_919 <- plot_CNV(df_cnv_fn, clones, meta_genes)
pcnv_919$cnv_plot
pcnv_919$plg
saveRDS(pcnv_919$cnv_plot, paste0(save_dir,datatag, '_pcnv.rds'))


input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/signif_genes/'
save_fig_dir <- paste0(base_dir, 'figures/fig_cnv/')
# dir.create(save_fig_dir)
# plottitle <- 'SA609 cisplatin: logFC X7-Rx clone A vs X7-UnRx clone H'
plottitle <- NULL
# desc <- 'log2FC X7-Rx clone A vs X7-UnRx clone H'
# desc <- 'SA609: Rx X7:cloneA vs. UnRx X7:cloneH'
desc <- 'Pt1: Metastasis:cloneC vs. Primary:cloneB'
datatag <- 'Pt1'
tag <- 'Pt1'

# deg_fn <- paste0(base_dir,'SA609_rna/deg_analysis/SA609-v6/SA609_UTTTT_R_UUUUU_H/signif_genes.csv')
# deg_fn <- paste0(input_dir,'scrande_SA609_4/signif_genes.csv')
data_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/results_bulkRNAseq/SA919_full/'
obs_clones <- c('B','C')
subtag <- paste0(obs_clones[2],'met_',obs_clones[1],'pri')
deg_fn <- paste0(data_dir, subtag, '/', subtag, '_DE_genes.csv.gz')

cis_genes <- get_cis_genes(meta_genes, obs_clones)
  
  
deg_df <- read.csv(deg_fn, check.names=F, stringsAsFactors=F)
deg_df <- deg_df %>%
  dplyr::filter(ens_gene_id %in% cis_genes)
deg_df$symbol
t <- deg_df %>%
  dplyr::filter(symbol=='VIM')
t
dim(deg_df)
# exp_df <- data.table::fread(paste0(data_dir, subtag, '/', subtag, '_DE_genes.csv.gz'))
data_dir_drivernet <- "/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/preprocessed_09April2024/drivernet_demo/"
additional_genes_fn <- paste0(data_dir_drivernet, subtag,'/', 'significant_genes_DriverNet.csv')
annotated_genes_df <- data.table::fread(additional_genes_fn)
dim(annotated_genes_df)
head(annotated_genes_df)
annotated_genes_df$`p-value`
annotated_genes_df <- annotated_genes_df %>%
  dplyr::filter(`p-value`<0.05)

# additional_genes_symb <- unique(annotated_genes_df$gene_symbol)
additional_genes_symb <- c("RBM17","CEP72","EXOC3","MCM10")
# additional_genes_symb <- annotated_genes_df$gene_symbol
## get drivernet significant genes output here
# additional_genes_fn <- paste0(save_fig_dir, 'SA609_R_vs_H_annotated_genes.csv')
# subtitle <- paste0(datatag,': in cis DE genes, Rx X7:cloneA vs. UnRx X7:cloneH')
# subtitle <- paste0(tag,': Rx X7:cloneA vs. UnRx X7:cloneH, in cis DE genes')
subtitle <- desc
track_incis <- plot_gex_cnv_v2(deg_df=deg_df,
                               df_cnv = pcnv_919$df_cnv,
                               clones=NULL,
                               save_dir=save_fig_dir,
                               sample_name=NULL,
                               clone_str=desc,
                               additional_genes_symb = additional_genes_symb,
                               n_genes_to_annotate = 10, 
                               plttitle=subtitle,
                               gene_type='in-cis')
track_incis$plg
main_plot <- cowplot::plot_grid(
  track_incis$track_plot,
  pcnv_919$cnv_plot,
  ncol = 1,
  rel_heights = c(1.6,1.1),
  align = 'v'#,
  # labels = c('a')
)

main_plot
save_fig_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/results_bulkRNAseq/SA919_full/'
saveRDS(main_plot, paste0(save_fig_dir,"figs/","gene_exp_Fig6_partA.rds"))
# check vim gene
save_fn <- paste0(save_fig_dir, 'clones_CB_cis_genes_track_plot.png')
png(save_fn, height = 2*500, width=2*1300, res = 2*72)
print(main_plot)
dev.off()


ggsave(paste0(save_dir,"Fig6_SA919_Pt1_trackplot_panelA.svg"),
       plot = main_plot,
       height = 4,
       width = 12,
       # useDingbats=F,
       dpi = 150
)
dev.off()



