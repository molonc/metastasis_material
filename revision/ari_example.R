#!/usr/bin/env Rscript
# ============================================================
# Adjusted Rand Index (ARI) — from-scratch implementation
# with the same 6-point example used in discussion.
# ============================================================

#' Compute the Adjusted Rand Index between two label vectors.
#'
#' @param true_labels  Integer or character vector of ground-truth cluster labels.
#' @param pred_labels  Integer or character vector of predicted cluster labels
#'                     (same length as true_labels).
#' @return A single numeric value: the ARI.
compute_ari <- function(true_labels, pred_labels) {
  stopifnot(length(true_labels) == length(pred_labels))

  # --- Step 1: build the contingency table -----------------------
  cont <- table(true_labels, pred_labels)

  # --- Step 2: helper — C(n, 2) = n*(n-1)/2 ---------------------
  choose2 <- function(n) n * (n - 1) / 2

  # --- Step 3: compute the sums of C(n_ij, 2) -------------------
  sum_nij <- sum(choose2(cont))         # over all cells
  sum_ai  <- sum(choose2(rowSums(cont))) # over row sums
  sum_bj  <- sum(choose2(colSums(cont))) # over column sums
  n_total <- choose2(sum(cont))          # C(n, 2)

  # --- Step 4: ARI formula --------------------------------------
  expected <- (sum_ai * sum_bj) / n_total
  max_idx  <- 0.5 * (sum_ai + sum_bj)

  ari <- (sum_nij - expected) / (max_idx - expected)
  return(ari)
}


# ================================================================
# Example: the same data points 1-6 from our discussion
# ================================================================
# True labels (U):  X = {1,2,3}   Y = {4,5,6}
# Predicted   (V):  A = {1,2}     B = {3,4,5,6}

true_labels <- c("X", "X", "X", "Y", "Y", "Y")
pred_labels <- c("A", "A", "B", "B", "B", "B")

ari <- compute_ari(true_labels, pred_labels)

cat("True labels:     ", true_labels, "\n")
cat("Predicted labels:", pred_labels, "\n")
cat("ARI =", round(ari, 4), "\n")

# ================================================================
# Cross-check with the mclust package (if installed)
# ================================================================
if (requireNamespace("mclust", quietly = TRUE)) {
  ari_pkg <- mclust::adjustedRandIndex(true_labels, pred_labels)
  cat("ARI (mclust):   ", round(ari_pkg, 4), "\n")
} else {
  cat("\nTip: install.packages('mclust') to cross-check with adjustedRandIndex()\n")
}
