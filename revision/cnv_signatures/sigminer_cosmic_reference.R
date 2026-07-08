# ============================================================
# sigminer: Copy Number Signature Fitting with Reference (COSMIC-style) Database
# ============================================================
# Install (run once):
#   install.packages("sigminer")
#   remotes::install_github("ShixiangWang/copynumber")  # optional, recommended

library(sigminer)

# ============================================================
# STEP 1: Load Copy Number Segment Data
# ============================================================
# Required columns: sample, chromosome, start, end, segVal (integer absolute CN)
load(system.file("extdata", "toy_segTab.RData",
  package = "sigminer", mustWork = TRUE
))

head(segTabs)

# ============================================================
# STEP 2: Read into CopyNumber Object
# ============================================================
cn <- read_copynumber(
  segTabs,
  seg_cols       = c("chromosome", "start", "end", "segVal"),
  genome_build   = "hg19",       # or "hg38"
  genome_measure = "wg",
  complement     = TRUE,
  add_loh        = FALSE,        # set TRUE only if you have minor_cn
  verbose        = TRUE
)

# ============================================================
# STEP 3: Tally Components -> sample x component matrix
# ============================================================
# tally_w <- sig_tally(cn, method = "W")  # 8 features -> CN_48 components

mat <- tally_s$all_matrices$CN_48
head(mat)
length(mat)
dim(mat)
mat[1:5,1:5]
colnames(mat)
# ============================================================
# STEP 4: Available Reference (COSMIC-style) CN Signature Databases
# ============================================================
# "CNS_TCGA"      -> 48 categories, TCGA-derived  (matches CN_48 matrix above)
# "CNS_USARC"     -> 40 categories, USARC cohort
# "CNS_TCGA176"   -> 176 categories, COSMIC/PCAWG+TCGA reference (Sanger COSMIC CN sigs)
# "CNS_PCAWG176"  -> 176 categories, PCAWG reference
#
# NOTE: sigminer does not use a literal name "CNS_COSMIC".
# The COSMIC-equivalent CN signature set is "CNS_TCGA176" (or "CNS_PCAWG176"),
# both published on the Sanger COSMIC CN signatures page:
# https://cancer.sanger.ac.uk/signatures/cn/

# Inspect the reference database
db_info <- get_sig_db("CNS_TCGA176")
dim(db_info$db)          # 176 components x N reference signatures
colnames(db_info$db)     # signature names, e.g. "CN1", "CN2", ...

# ============================================================
# STEP 5: Fit Samples Against COSMIC-style Reference Signatures
# ============================================================
# IMPORTANT: For 176-category databases, you must tally with the
# matching component scheme. CN_48 (method "W") is NOT directly compatible
# with CNS_TCGA176 (176-category scheme uses method "X"/absolute+allele-based
# components). If you only have CN_48, use "CNS_TCGA" (48-cat) instead.

# --- Option A: Fit using CN_48 against CNS_TCGA (48 categories, matched) ---
act_refit_48 <- sig_fit(
  t(mat),                  # transpose: components x samples
  sig_index = "ALL",
  sig_db    = "CNS_TCGA"   # 48-category COSMIC/TCGA-derived reference
)

head(t(act_refit_48))

# --- Option B: Fit using 176-category matrix (if you tallied with method "X") ---
# tally_x <- sig_tally(cn, method = "X")  # generates allele-specific 176 components
# mat176  <- tally_x$all_matrices$CN_176
# act_refit_176 <- sig_fit(
#   t(mat176),
#   sig_index = "ALL",
#   sig_db    = "CNS_TCGA176"   # COSMIC-style 176-category reference
# )

# ============================================================
# STEP 6: Visualise Reference Signature Profiles
# ============================================================
# Plot the COSMIC/TCGA reference signature profiles themselves
# show_sig_profile(
#   get_sig_db("CNS_TCGA")$db,
#   style       = "cosmic",
#   mode        = "copynumber",
#   method      = "W",
#   check_sig_names = FALSE
# )

show_sig_profile(
  get_sig_db("CNS_TCGA")$db,
  style       = "cosmic",
  mode        = "copynumber",
  method      = "S",
  check_sig_names = FALSE
)

# Plot fitted exposures for your samples
show_sig_exposure(act_refit_48)

# ============================================================
# STEP 7: Compare Fitted Exposures Between Clones
# ============================================================
expo <- t(act_refit_48)  # samples x signatures

clone_map <- data.frame(
  sample = rownames(expo),
  clone  = sample(c("CloneA", "CloneB"), nrow(expo), replace = TRUE)
)

expo_df <- as.data.frame(expo)
expo_df$sample <- rownames(expo_df)
expo_df <- merge(expo_df, clone_map, by = "sample")

library(dplyr)
library(tidyr)
library(ggplot2)

expo_long <- expo_df |>
  pivot_longer(cols = -c(sample, clone), names_to = "signature", values_to = "exposure")

# Summary stats per clone per signature
expo_long |>
  group_by(clone, signature) |>
  summarise(mean_exp = mean(exposure), sd_exp = sd(exposure), .groups = "drop")

# Boxplot comparison
ggplot(expo_long, aes(x = clone, y = exposure, fill = clone)) +
  geom_boxplot() +
  facet_wrap(~signature, scales = "free_y") +
  theme_bw() +
  labs(title = "COSMIC/TCGA CN Signature Exposure by Clone",
       x = "Clone", y = "Exposure")

# ============================================================
# STEP 8 (Optional): Match De Novo Signatures to COSMIC References
# ============================================================
# If you extracted de novo signatures (sig_extract), compare to references:
#
# sig_denovo <- sig_extract(mat, n_sig = 3, nrun = 10)
# get_sig_similarity(sig_denovo, sig_db = "CNS_TCGA")
# -> reports cosine similarity + COSMIC aetiology links
