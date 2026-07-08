suppressPackageStartupMessages({
  library(annotables)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  options(dplyr.summarise.inform = FALSE)
  options(tidyverse.quiet = TRUE)
  library(sigminer)
})
get_chr_infos <- function(cnv) {
  chr_desc <- as.character(cnv$chr_desc)
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
  cnv$chromosome <- as.character(chrs)
  cnv$start <- as.character(starts)
  cnv$end <- as.character(ends)
  cnv <- cnv %>% 
    dplyr::filter(chromosome %in% obs_chrs)
  # return(list(chr=as.character(chrs),start=as.character(starts),end=as.character(ends)))
  return(cnv)
}
load_cnv_data <- function(){
  input_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/'
  save_dir <- paste0(input_dir,'revision/cnv_signatures/')
  results_dir <- paste0(input_dir,"materials/dlp_trees/")
  source(paste0(input_dir,"scripts/corrupt_tree/src/cn_change/utils.R"))
  
  datatag <- 'SA919'
  copynumber_fn <- paste0(results_dir, datatag, '/total_merged_filtered_states.csv.gz')
  cellclone_fn <- paste0(results_dir, datatag, '/cell_clones.csv.gz')
  
  if(!dir.exists(save_dir)){
    dir.create(save_dir, recursive = T)
  }
  
  cell_clones <- data.table::fread(cellclone_fn) %>% as.data.frame()
  dim(cell_clones)
  cell_clones <- cell_clones %>%
    dplyr::filter(!clone_id %in% c('None','unassigned'))
  
  copynumber <- as.data.frame(data.table::fread(copynumber_fn))
  rownames(copynumber) <- copynumber$V1
  copynumber$V1 <- NULL
  
  sampling_fraction <- 0.3
  # sampling_cells <- sample(cell_clones$cell_id, dim(cell_clones)[1]*sampling_fraction)
  ## Subsampling data per clone
  set.seed(42)
  sampling_cells <- c()
  for(cl in unique(cell_clones$clone_id)){
    cell_ids <- cell_clones %>%
      dplyr::filter(clone_id==cl) %>%
      dplyr::pull(cell_id)
    sampled_per_clone <- sample(cell_ids, size = length(cell_ids) * sampling_fraction)
    sampling_cells <- c(sampling_cells, sampled_per_clone)
  }
  print('Number of randomized cells is: ')
  print(length(sampling_cells))
  
  # sampling_chr_pos <- sample(rownames(copynumber), dim(copynumber)[1]*0.5)
  # length(sampling_chr_pos)
  # copynumber <- copynumber[sampling_chr_pos, sampling_cells]
  copynumber <- copynumber[, sampling_cells]
  copynumber$chr_desc <- rownames(copynumber)
  cnv <- copynumber %>%
    pivot_longer(!chr_desc, names_to = "cell_id", values_to = "segVal")
  dim(cnv)
  head(cnv)
  print(dim(cnv))
  dim(copynumber)
  cnv_chrs <- data.frame(chr_desc=copynumber$chr_desc)
  cnv_chrs <- get_chr_infos(cnv_chrs)
  head(cnv_chrs)
  cnv <- cnv %>%
    dplyr::left_join(cnv_chrs, by='chr_desc') %>%
    dplyr::rename(sample=cell_id) %>%
    dplyr::select(sample, chromosome, start, end, segVal)
  dim(cnv)
  data.table::fwrite(cnv, paste0(save_dir,'sampling_cnv_',sampling_fraction,'.csv.gz'))
  
  
}
# cnv
get_signatures <- function(){
  sampling_fraction <- 0.3
  input_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/'
  save_dir <- paste0(input_dir,'revision/cnv_signatures/')
  datatag <- 'SA919'
  results_dir <- paste0(input_dir,"materials/dlp_trees/")
  save_result_dir <- paste0(save_dir,'signatures/')
  if(!dir.exists(save_result_dir)){
    dir.create(save_result_dir, recursive = T)
  }
  cellclone_fn <- paste0(results_dir, datatag, '/cell_clones.csv.gz')
  cell_clones <- data.table::fread(cellclone_fn) %>%
    as.data.frame() %>%
    dplyr::filter(!clone_id %in% c('None','unassigned'))
  cnv <- data.table::fread(paste0(save_dir,'sampling_cnv_',sampling_fraction,'.csv.gz'))
  library(sigminer)
  ## Without ascn
  cn <- read_copynumber(
    cnv,
    seg_cols       = c("chromosome", "start", "end", "segVal"),
    genome_build   = "hg19",
    genome_measure = "wg",
    complement     = TRUE,
    add_loh        = FALSE   # <-- set FALSE; no minor_cn column needed
  )
  # Inspect the CopyNumber object
  cn
  show(cn)
  saveRDS(cn, paste0(save_result_dir,datatag,'_cn.rds'))
  # ============================================================
  # STEP 3: Tally Copy Number Features (build the sample × component matrix)
  # ============================================================
  # Method "W" = Wang et al. (2021) — 8 CN features, 48 fixed components (CN_48)
  # Method "M" = Macintyre et al. (2018) — 6 CN features, mixture modeling
  
  # Method "W" = Wang et al. (2021) — total copy-number features.
  # Method "S" (Steele 2019 / COSMIC CN48) requires allele-specific copy
  # number (a minor_cn / LOH column); this dataset has total copy number
  # only, so method "S" yields an all-zero matrix. Method "W" is the
  # appropriate choice for total-CN data.
  tally_w <- sig_tally(cn, method = "W")
  str(tally_w$nmf_matrix, max.level = 1)
  saveRDS(tally_w, paste0(save_result_dir,datatag,'_tally_w.rds'))
  # The matrix for NMF is tally_w$nmf_matrix (samples × components). Drop any
  # all-zero components (NMF cannot factorize null feature rows).
  nmf_mat <- tally_w$nmf_matrix
  nmf_mat <- nmf_mat[, colSums(nmf_mat) > 0, drop = FALSE]
  dim(nmf_mat)
  
  # ============================================================
  # STEP 4: Extract Copy Number Signatures (NMF de novo)
  # ============================================================
  # Estimate the best number of signatures (test K = 2 to 5)
  # NOTE: use nrun >= 30 in real analyses; nrun=5 here for speed
  sig_est <- sig_estimate(
    nmf_mat,
    range   = 2:5,
    nrun    = 5,
    verbose = TRUE
  )
  saveRDS(sig_est, paste0(save_result_dir,datatag,'_sig_est.rds'))
  # Plot the survey to choose optimal K
  show_sig_number_survey(sig_est)
  
  # Extract signatures with the chosen K (e.g. K=5)
  sig <- sig_extract(
    nmf_mat,
    n_sig  = 5,
    nrun   = 10,
    cores  = 1
  )
  saveRDS(sig, paste0(save_result_dir,datatag,'_sig.rds'))
  # ============================================================
  # STEP 5: Visualise Signatures
  # ============================================================
  # Plot the CN signature profiles
  show_sig_profile(sig, mode = "copynumber", style = "cosmic", x_label_angle = 90)
  
  # Plot per-sample signature exposures (contributions)
  show_sig_exposure(sig)
  
  # ============================================================
  # STEP 6: (Alternative) Re-fit Samples Against Known Signatures
  # ============================================================
  # Instead of (or after) de novo extraction you can quantify how much
  # each known signature contributes to every sample. Here we refit the
  # samples against the de novo signatures extracted in Step 4.
  # (The built-in reference DBs such as sig_db = "CNS_TCGA" use the COSMIC
  #  CN48 component scheme, which the toy method-"W" tally does not match.)
  
  act_refit <- sig_fit(
    t(nmf_mat),   # transpose: components × samples
    sig = sig,               # refit against the de novo signatures from Step 4
    return_class = "matrix"
  )
  # act_refit is a signatures × samples exposure matrix
  
  # View the exposure matrix
  # head(act_refit)
  
  # Plot the refitted exposures
  show_sig_exposure(act_refit)
  saveRDS(act_refit, paste0(save_result_dir,datatag,'_act_refit.rds'))
  # ============================================================
  # STEP 7: Compare Signatures Between Clones / Groups
  # ============================================================
  # Transpose to samples × signatures for per-sample comparison
  expo <- t(act_refit)
  
  # Assuming you have a data.frame mapping samples → clone labels:
  cell_clones_map <- cell_clones %>%
    dplyr::filter(cell_id %in% unique(cnv$sample)) %>%
    dplyr::rename(clone=clone_id, sample=cell_id)
  # Merge exposures with clone labels
  expo_df <- as.data.frame(expo)
  expo_df$sample <- rownames(expo_df)
  expo_df <- merge(expo_df, cell_clones_map, by = "sample")
  data.table::fwrite(expo_df, paste0(save_result_dir,datatag,'_expo_df.csv.gz'))
  # Compare signature exposure between clones (Wilcoxon test)
  library(dplyr)
  library(tidyr)
  
  expo_long <- expo_df |>
    pivot_longer(cols = -c(sample, clone), names_to = "signature", values_to = "exposure")
  
  # Summary statistics per clone per signature
  expo_long |>
    group_by(clone, signature) |>
    summarise(mean_exp = mean(exposure), sd_exp = sd(exposure), .groups = "drop")
  data.table::fwrite(expo_long, paste0(save_result_dir,datatag,'_expo_long.csv.gz'))

}
# ============================================================
# STEP 8: Pairwise statistical comparison of signature exposures between clones
# ============================================================
# For each signature, run a two-sample Kolmogorov-Smirnov test on the
# exposure distributions of every pair of clones.
stat_signatures <- function(){
  sampling_fraction <- 0.3
  input_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/'
  save_dir <- paste0(input_dir,'revision/cnv_signatures/')
  datatag <- 'SA919'
  save_result_dir <- paste0(save_dir,'signatures/')

  expo_long <- data.table::fread(paste0(save_result_dir,datatag,'_expo_long.csv.gz')) %>%
    as.data.frame()

  signatures <- sort(unique(expo_long$signature))
  clones     <- sort(unique(expo_long$clone))
  clone_pairs <- utils::combn(clones, 2, simplify = FALSE)

  results <- list()
  for(sig in signatures){
    sig_df <- expo_long %>% dplyr::filter(signature == sig)
    for(pair in clone_pairs){
      x <- sig_df %>% dplyr::filter(clone == pair[1]) %>% dplyr::pull(exposure)
      y <- sig_df %>% dplyr::filter(clone == pair[2]) %>% dplyr::pull(exposure)
      kt <- suppressWarnings(stats::ks.test(x, y))
      results[[length(results) + 1]] <- data.frame(
        signature = sig,
        clone1    = pair[1],
        clone2    = pair[2],
        n1        = length(x),
        n2        = length(y),
        D         = unname(kt$statistic),
        p_value   = kt$p.value,
        stringsAsFactors = FALSE
      )
    }
  }
  ks_df <- do.call(rbind, results)
  # Adjust p-values across all tests (Benjamini-Hochberg)
  ks_df$p_adj <- stats::p.adjust(ks_df$p_value, method = "BH")
  print(ks_df)
  data.table::fwrite(ks_df, paste0(save_result_dir,datatag,'_signature_ks_pairwise.csv.gz'))
  return(ks_df)
}

# ============================================================
# STEP 9: Visualise signature exposures per clone (boxplot)
# ============================================================
viz_signatures <- function(){
  input_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/'
  save_dir <- paste0(input_dir,'revision/cnv_signatures/')
  datatag <- 'SA919'
  save_result_dir <- paste0(save_dir,'signatures/')

  expo_long <- data.table::fread(paste0(save_result_dir,datatag,'_expo_long.csv.gz')) %>%
    as.data.frame()

  # Simple box plot comparison
  library(ggplot2)
  clone_colors <- c(A = "#66C2A5", B = "#FC8D62", C = "#8DA0CB")
  p <- ggplot(expo_long, aes(x = clone, y = exposure, fill = clone)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = clone_colors) +
    facet_wrap(~signature, nrow = 1, scales = "free_y") +
    theme_bw() +
    theme(legend.position = 'bottom',
          panel.grid.minor = element_blank()) + 
    labs(title = "CN Signature Exposure by Clone",
         x = "Clone", y = "Signature Exposure")
  # p
  png(paste0(save_result_dir,datatag,'_signatures.png'),
      height = 2*250, width=2*800,res = 2*72)
  print(p)
  dev.off()

  # Also save as vector graphics (.svg)
  ggsave(paste0(save_result_dir,datatag,'_signatures.svg'),
         plot = p, width = 800/72, height = 250/72, units = "in")
}

get_signatures()
stat_signatures()
viz_signatures()