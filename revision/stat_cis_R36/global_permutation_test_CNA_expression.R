# ============================================================
# Global Permutation Test: CNA-Expression Correlation ACROSS Genes
# ============================================================
# Purpose: test whether copy number is associated with expression, IN CIS,
# across your WHOLE gene set (not gene-by-gene) -- i.e. a single p-value
# summarising whether CNA-expression coupling is stronger than expected
# by chance across all genes considered together.
#
# Null model: randomly re-pair each gene's expression profile with a
# DIFFERENT gene's copy number profile (permuting gene identity), then
# see how extreme the observed aggregate correlation is relative to that
# shuffled-pairing null distribution.

library(ggplot2)

# ============================================================
# FUNCTION: permutation_test_correlation_CNA_gene_exp
# ============================================================
# ARGUMENTS:
#   copy_number : genes (rows) x samples (columns) matrix, numeric/integer CN
#   expression  : genes (rows) x samples (columns) matrix, normalized expression
#                 (must have identical row/column names & order as copy_number)
#   n_perm      : number of permutations (default 10000)
#   method      : correlation method, "spearman" (default) or "pearson"
#   summary_stat: how to aggregate per-gene correlations into one global
#                 statistic -- "mean_abs_r" (default), "mean_r", or "median_abs_r"
#
# WHAT IT DOES:
#   1. Computes each gene's own CN-vs-expression correlation (per-gene r)
#   2. Aggregates those r values into ONE global statistic
#   3. Builds a null distribution by shuffling WHICH gene's expression
#      profile is paired with WHICH gene's CN profile, recomputing the
#      global statistic each time
#   4. Returns a single empirical p-value for the whole gene set
#
# RETURNS: a list with
#   $observed_stat   -- the observed global statistic
#   $p_value         -- empirical p-value (one number, for the whole gene set)
#   $null_dist       -- vector of null statistics (length n_perm)
#   $per_gene_r      -- data.frame of each gene's own correlation (for reference/plotting only,
#                        NOT individually tested -- no per-gene p-values)
#   $plot            -- ggplot histogram of the null distribution with observed value marked

permutation_test_correlation_CNA_gene_exp <- function(copy_number,
                                                        expression,
                                                        n_perm = 10000,
                                                        method = "spearman",
                                                        summary_stat = "mean_abs_r") {

  stopifnot(identical(rownames(copy_number), rownames(expression)))
  stopifnot(identical(colnames(copy_number), colnames(expression)))

  genes <- rownames(copy_number)
  n_genes <- length(genes)

  # --- per-gene correlations (building block for the global statistic) ---
  per_gene_r <- sapply(seq_len(n_genes), function(i) {
    cor(copy_number[i, ], expression[i, ], method = method)
  })
  names(per_gene_r) <- genes

  aggregate_fun <- switch(
    summary_stat,
    mean_abs_r   = function(r) mean(abs(r), na.rm = TRUE),
    mean_r       = function(r) mean(r, na.rm = TRUE),
    median_abs_r = function(r) median(abs(r), na.rm = TRUE),
    stop("summary_stat must be one of: 'mean_abs_r', 'mean_r', 'median_abs_r'")
  )

  # --- observed global statistic ---
  observed_stat <- aggregate_fun(per_gene_r)

  # --- null distribution: shuffle which gene's expression pairs with
  #     which gene's CN, recompute the global statistic each time ---
  null_dist <- replicate(n_perm, {
    shuffled_idx <- sample(n_genes)
    shuffled_expression <- expression[shuffled_idx, , drop = FALSE]

    perm_r <- sapply(seq_len(n_genes), function(i) {
      cor(copy_number[i, ], shuffled_expression[i, ], method = method)
    })

    aggregate_fun(perm_r)
    # perm_r <- cor(copy_number, shuffled_expression, method = method)
  })

  # --- empirical p-value ---
  # one-sided: is the observed coupling MORE extreme (larger) than random pairing?
  p_value <- mean(null_dist >= observed_stat)

  # --- diagnostic plot ---
  p <- ggplot(data.frame(null_stat = null_dist), aes(x = null_stat)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = observed_stat, color = "red", linewidth = 1) +
    labs(
      title    = "Global permutation test: CNA-expression correlation across gene set",
      subtitle = paste0("Observed ", summary_stat, " = ", round(observed_stat, 3),
                         ", empirical p = ", signif(p_value, 3),
                         " (n_genes = ", n_genes, ", n_perm = ", n_perm, ")"),
      x = paste0("Permuted ", summary_stat, " (", method, " correlation, shuffled gene pairing)"),
      y = "Frequency"
    ) +
    theme_bw()

  list(
    observed_stat = observed_stat,
    p_value       = p_value,
    null_dist     = null_dist,
    per_gene_r    = data.frame(gene = genes, observed_r = per_gene_r, row.names = NULL),
    plot          = p
  )
}

# ============================================================
# EXAMPLE USAGE
# ============================================================
genes   <- c("HOXA9","LAMB1","SEMA3A","RHEB","HSPB1","VIM","RAC1","TWIST1")
samples <- c("CloneA_1","CloneA_2","CloneB_1","CloneB_2","CloneC_1","CloneC_2")

set.seed(42)
copy_number <- matrix(
  sample(1:6, length(genes) * length(samples), replace = TRUE),
  nrow = length(genes), dimnames = list(genes, samples)
)

# Simulate expression correlated with CN for demo purposes
expression <- copy_number * 1.5 + matrix(rnorm(length(genes) * length(samples), sd = 2),
                                          nrow = length(genes))
dimnames(expression) <- list(genes, samples)

out <- permutation_test_correlation_CNA_gene_exp(
  copy_number  = copy_number,
  expression   = expression,
  n_perm       = 1000,
  method       = "spearman",
  summary_stat = "mean_abs_r"
)

cat("Observed global statistic:", round(out$observed_stat, 3), "\n")
cat("Global empirical p-value :", signif(out$p_value, 3), "\n\n")

print(out$per_gene_r)   # per-gene correlations, for reference/reporting only
print(out$plot)

# ============================================================
# NOTES
# ============================================================
# - This gives ONE p-value describing whether your gene set, AS A WHOLE,
#   shows more CN-expression coupling than expected from randomly
#   re-pairing genes -- appropriate when you have few samples per gene
#   (e.g. 3 clones) and per-gene permutation would lack resolution.
# - $per_gene_r is included only as descriptive/reference output (e.g. to
#   rank or plot individual genes) -- it is NOT independently significance-
#   tested. Report only the single global p-value as your formal result.
# - "mean_abs_r" treats over- and under-expression relative to CN symmetrically
#   (i.e. both amplification->up and deletion->down count as "coupling").
#   Use "mean_r" instead if you specifically expect a consistent DIRECTION
#   (e.g. dosage effect: more copies -> more expression) across the gene set.
write.csv(out$per_gene_r, "per_gene_correlations_reference.csv", row.names = FALSE)
