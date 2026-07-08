# ============================================================
# sigminer: Copy Number Signature Analysis - Simple Example
# ============================================================
# Install (run once):
#   install.packages("sigminer")
#   remotes::install_github("ShixiangWang/copynumber")  # optional but recommended

library(sigminer)

# ============================================================
# STEP 1: Load & Prepare Copy Number Segment Data
# ============================================================
# Required columns: sample, chromosome, start, end, segVal (integer absolute CN)

# Load the built-in toy dataset (absolute copy number segments)
load(system.file("extdata", "toy_segTab.RData",
  package = "sigminer", mustWork = TRUE
))

# Peek at the data structure
head(segTabs)
# Expected columns: sample, chromosome, start, end, segVal

# Optional: add minor allele copy number (for LOH analysis)
set.seed(1234)
segTabs$minor_cn <- sample(c(0, 1), size = nrow(segTabs), replace = TRUE)
summary(segTabs$minor_cn)
length(segTabs$minor_cn)
segTabs$minor_cn[1:50]
# ============================================================
# STEP 2: Read Copy Number Data into sigminer CopyNumber Object
# ============================================================
## With ascn
cn <- read_copynumber(
  segTabs,
  seg_cols        = c("chromosome", "start", "end", "segVal"),
  genome_build    = "hg19",      # or "hg38"
  genome_measure  = "wg",        # "wg" = whole genome
  complement      = TRUE,        # fill uncalled chromosomes with normal CN=2
  add_loh         = TRUE,        # compute loss-of-heterozygosity (requires minor_cn)
  verbose         = TRUE
)

## Without ascn
cn <- read_copynumber(
  segTabs,
  seg_cols       = c("chromosome", "start", "end", "segVal"),
  genome_build   = "hg19",
  genome_measure = "wg",
  complement     = TRUE,
  add_loh        = FALSE   # <-- set FALSE; no minor_cn column needed
)
# Inspect the CopyNumber object
cn
show(cn)

# ============================================================
# STEP 3: Tally Copy Number Features (build the sample × component matrix)
# ============================================================
# Method "W" = Wang et al. (2021) — 8 CN features, 48 fixed components (CN_48)
# Method "M" = Macintyre et al. (2018) — 6 CN features, mixture modeling

tally_w <- sig_tally(cn, method = "W")   # recommended, reproducible
# tally_m <- sig_tally(cn, method = "M")   # alternative (method "M" not supported in this sigminer version)

# The key matrix for NMF is in:
#   tally_w$nmf_matrix   (samples × 48 components)
dim(tally_w$nmf_matrix)

# ============================================================
# STEP 4: Extract Copy Number Signatures (NMF de novo)
# ============================================================
# Estimate the best number of signatures (test K = 2 to 5)
# NOTE: use nrun >= 30 in real analyses; nrun=5 here for speed
sig_est <- sig_estimate(
  tally_w$nmf_matrix,
  range   = 2:5,
  nrun    = 5,
  verbose = TRUE
)

# Plot the survey to choose optimal K
show_sig_number_survey(sig_est)

# Extract signatures with the chosen K (e.g. K=3)
sig <- sig_extract(
  tally_w$nmf_matrix,
  n_sig  = 3,
  nrun   = 10,
  cores  = 1
)

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
  t(tally_w$nmf_matrix),   # transpose: components × samples
  sig = sig,               # refit against the de novo signatures from Step 4
  return_class = "matrix"
)
# act_refit is a signatures × samples exposure matrix

# View the exposure matrix
head(act_refit)

# Plot the refitted exposures
show_sig_exposure(act_refit)

# ============================================================
# STEP 7: Compare Signatures Between Clones / Groups
# ============================================================
# Transpose to samples × signatures for per-sample comparison
expo <- t(act_refit)

# Assuming you have a data.frame mapping samples → clone labels:
clone_map <- data.frame(
  sample = rownames(expo),
  clone  = sample(c("CloneA", "CloneB"), nrow(expo), replace = TRUE)
)

# Merge exposures with clone labels
expo_df <- as.data.frame(expo)
expo_df$sample <- rownames(expo_df)
expo_df <- merge(expo_df, clone_map, by = "sample")

# Compare signature exposure between clones (Wilcoxon test)
library(dplyr)
library(tidyr)

expo_long <- expo_df |>
  pivot_longer(cols = -c(sample, clone), names_to = "signature", values_to = "exposure")

# Summary statistics per clone per signature
expo_long |>
  group_by(clone, signature) |>
  summarise(mean_exp = mean(exposure), sd_exp = sd(exposure), .groups = "drop")

# Simple box plot comparison
library(ggplot2)
ggplot(expo_long, aes(x = clone, y = exposure, fill = clone)) +
  geom_boxplot() +
  facet_wrap(~signature, scales = "free_y") +
  theme_bw() +
  labs(title = "CN Signature Exposure by Clone",
       x = "Clone", y = "Signature Exposure")
