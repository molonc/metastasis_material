# ============================================================
# Custom Permutation Test: CNA / SNV in-cis Gene Expression Association
# ============================================================
# Purpose: test whether a gene's OWN copy number (or SNV/mutation status)
# is associated with its OWN expression, across your clones (or samples),
# more than expected by chance -- via permutation rather than parametric
# eQTL-style tools (which assume large population sample sizes).
#
# This is appropriate for small-n designs (e.g. 3 clones: A, B, C) where
# FastQTL/Matrix eQTL's asymptotic assumptions don't hold.

library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(123)

# ============================================================
# STEPS 1-5 WRAPPED INTO ONE FUNCTION
# ============================================================
# permutation_test_correlation_CNA_gene_exp()
#
# Runs the full workflow: observed CN-expression correlation per gene,
# permutation null distribution (shuffling expression labels), empirical
# p-values, BH/FDR correction, and a diagnostic plot for one example gene.
#
# ARGUMENTS:
#   copy_number : genes (rows) x samples (columns) matrix, numeric/integer CN
#   expression  : genes (rows) x samples (columns) matrix, normalized expression
#                 (must have identical dimnames/order to copy_number)
#   n_perm      : number of permutations per gene (default 10000)
#   method      : correlation method, "spearman" (default) or "pearson"
#   plot_gene   : which gene to plot the null distribution for; defaults to
#                 the top hit (lowest p-value) if left NULL
#
# RETURNS: a list with
#   $results  -- data.frame(gene, observed_r, p_value, fdr), sorted by p_value
#   $plot     -- ggplot object of the null distribution for plot_gene
#   $null_distributions -- named list of each gene's full null vector
#                          (useful for custom plotting/inspection)

permutation_test_correlation_CNA_gene_exp <- function(copy_number,
                                                        expression,
                                                        n_perm = 10000,
                                                        method = "spearman",
                                                        plot_gene = NULL) {

  stopifnot(identical(rownames(copy_number), rownames(expression)))
  stopifnot(identical(colnames(copy_number), colnames(expression)))

  genes <- rownames(copy_number)

  # --- Step 2: observed statistic per gene ---
  obs_cor <- sapply(genes, function(g) {
    cor(copy_number[g, ], expression[g, ], method = method)
  })

  # --- Step 3: permutation null distribution (per gene) ---
  permutation_test_gene <- function(cn_vec, expr_vec) {
    obs <- cor(cn_vec, expr_vec, method = method)
    null_dist <- replicate(n_perm, {
      shuffled_expr <- sample(expr_vec)          # permute sample labels
      cor(cn_vec, shuffled_expr, method = method)
    })
    p_value <- mean(abs(null_dist) >= abs(obs))  # two-sided empirical p
    list(observed = obs, null_dist = null_dist, p_value = p_value)
  }

  gene_results <- lapply(genes, function(g) {
    permutation_test_gene(copy_number[g, ], expression[g, ])
  })
  names(gene_results) <- genes

  results_df <- data.frame(
    gene       = genes,
    observed_r = sapply(gene_results, function(x) x$observed),
    p_value    = sapply(gene_results, function(x) x$p_value)
  )

  # --- Step 4: multiple testing correction (BH / FDR) ---
  results_df$fdr <- p.adjust(results_df$p_value, method = "BH")
  results_df <- results_df[order(results_df$p_value), ]

  # --- Step 5: diagnostic plot for one gene's null distribution ---
  if (is.null(plot_gene)) plot_gene <- results_df$gene[1]
  res_example <- gene_results[[plot_gene]]

  p <- ggplot(data.frame(null_r = res_example$null_dist), aes(x = null_r)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = res_example$observed, color = "red", linewidth = 1) +
    labs(
      title    = paste("Permutation null distribution:", plot_gene),
      subtitle = paste0("Observed r = ", round(res_example$observed, 3),
                         ", empirical p = ", signif(res_example$p_value, 3)),
      x = paste0("Permuted ", method, " correlation (CN vs expression)"),
      y = "Frequency"
    ) +
    theme_bw()
  print(p)
  list(
    results             = results_df,
    plot                = p,
    null_distributions  = lapply(gene_results, function(x) x$null_dist)
  )
}

# ============================================================
# EXAMPLE USAGE
# ============================================================
genes   <- c("HOXA9","LAMB1","SEMA3A","RHEB","HSPB1","VIM","RAC1","TWIST1")
samples <- c("CloneA_1","CloneA_2","CloneB_1","CloneB_2","CloneC_1","CloneC_2")

genes   <- c("HOXA9","LAMB1","SEMA3A","RHEB","HSPB1","VIM","RAC1","TWIST1")
samples <- c("CloneA_1","CloneB_1","CloneC_1")

set.seed(42)
copy_number <- matrix(
  sample(1:6, length(genes) * length(samples), replace = TRUE),
  nrow = length(genes), dimnames = list(genes, samples)
)

# Simulate expression correlated with CN for some genes (for demo purposes)
expression <- copy_number * 1.5 + matrix(rnorm(length(genes) * length(samples), sd = 2),
                                          nrow = length(genes))
dimnames(expression) <- list(genes, samples)

out <- permutation_test_correlation_CNA_gene_exp(
  copy_number = copy_number,
  expression  = expression,
  n_perm      = 10000,
  method      = "spearman"
)

print(out$results)
print(out$plot)

# Plot a specific gene instead of the top hit:
out_hoxa9 <- permutation_test_correlation_CNA_gene_exp(
  copy_number, expression, n_perm = 10000, plot_gene = "HOXA9"
)
print(out_hoxa9$plot)

# ============================================================
# STEP 6 (ALTERNATIVE DESIGN): SNV / mutation status vs expression
# ============================================================
# If instead of continuous CN you have binary SNV/mutation status
# (mutant vs wild-type) per clone, use a mean-difference permutation
# test instead of correlation.

# Example: binary mutation status per gene per sample
mutation_status <- matrix(
  sample(c(0, 1), length(genes) * length(samples), replace = TRUE),
  nrow = length(genes), dimnames = list(genes, samples)
)

permutation_test_snv <- function(mut_vec, expr_vec, n_perm = 10000) {
  obs_diff <- mean(expr_vec[mut_vec == 1]) - mean(expr_vec[mut_vec == 0])

  null_dist <- replicate(n_perm, {
    shuffled_mut <- sample(mut_vec)
    mean(expr_vec[shuffled_mut == 1]) - mean(expr_vec[shuffled_mut == 0])
  })

  p_value <- mean(abs(null_dist) >= abs(obs_diff), na.rm = TRUE)
  list(observed_diff = obs_diff, null_dist = null_dist, p_value = p_value)
}

snv_results <- lapply(genes, function(g) {
  res <- permutation_test_snv(mutation_status[g, ], expression[g, ], n_perm = 10000)
  data.frame(gene = g, observed_diff = res$observed_diff, p_value = res$p_value)
})

snv_results_df <- do.call(rbind, snv_results)
snv_results_df$fdr <- p.adjust(snv_results_df$p_value, method = "BH")
snv_results_df <- snv_results_df[order(snv_results_df$p_value), ]

print(snv_results_df)

# ============================================================
# STEP 7: Export results
# ============================================================
write.csv(out$results, "cna_cis_permutation_results.csv", row.names = FALSE)
write.csv(snv_results_df, "snv_cis_permutation_results.csv", row.names = FALSE)

# ============================================================
# NOTES / CAVEATS FOR YOUR ACTUAL DESIGN
# ============================================================
# 1. With very few samples per clone (e.g. n=1 median CN per clone x 3 clones),
#    the number of DISTINCT permutations is limited (e.g. only 3! = 6 unique
#    orderings for 3 samples) -- the empirical p-value resolution is then
#    capped at 1/6 ~ 0.17 minimum. Consider:
#      a) Using replicate/single-cell level data (not just clone medians) to
#         increase n and permutation resolution, OR
#      b) Pooling across genes for a GLOBAL null (STEP 8 below) rather than
#         a per-gene permutation with too few unique orderings.
#
# 2. If testing MANY genes, always report FDR-adjusted p-values, not raw p.
#
# 3. Oncodrive-CIS (Tamborero et al. 2013, PLOS ONE) is the published,
#    peer-reviewed method built specifically for this CNA-in-cis-expression
#    question and additionally benchmarks against diploid/normal tissue
#    expression -- worth citing/using as the primary method, with this
#    custom permutation as either a validation step or for cases where
#    Oncodrive-CIS's TCGA-scale assumptions don't fit your PDX clone design.

# ============================================================
# STEP 8 (OPTION FOR SMALL n): Gene-pooled global permutation null
# ============================================================
# Instead of permuting each gene separately (limited by n! for small n),
# pool ALL gene-sample CN/expression pairs together and permute genes'
# expression vectors across THEMSELVES (i.e. shuffle which gene's expression
# profile is paired with which gene's CN profile). This tests whether the
# observed per-gene CN-expression correlations are more extreme, ON AVERAGE
# or IN AGGREGATE, than a random re-pairing of genes -- useful when you want
# genome-wide/gene-set-level significance rather than per-gene significance.

global_permutation_test <- function(copy_number, expression, n_perm = 10000, method = "spearman") {
  obs_cors <- sapply(seq_len(nrow(copy_number)), function(i) {
    cor(copy_number[i, ], expression[i, ], method = method)
  })
  obs_mean_abs_r <- mean(abs(obs_cors), na.rm = TRUE)

  null_means <- replicate(n_perm, {
    shuffled_idx <- sample(nrow(expression))          # shuffle which gene's
    shuffled_expression <- expression[shuffled_idx, ] # expression pairs with
                                                        # which gene's CN
    perm_cors <- sapply(seq_len(nrow(copy_number)), function(i) {
      cor(copy_number[i, ], shuffled_expression[i, ], method = method)
    })
    mean(abs(perm_cors), na.rm = TRUE)
  })

  p_value <- mean(null_means >= obs_mean_abs_r)
  list(observed_mean_abs_r = obs_mean_abs_r, null_dist = null_means, p_value = p_value)
}

global_res <- global_permutation_test(copy_number, expression, n_perm = 1000)
cat("Global CN-expression pairing test: observed mean|r| =",
    round(global_res$observed_mean_abs_r, 3),
    ", p =", signif(global_res$p_value, 3), "\n")
