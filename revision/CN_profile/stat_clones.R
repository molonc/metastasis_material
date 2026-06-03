library(ggplot2)

cn_changes <- read.csv("SA919_cn_changes.csv", stringsAsFactors = FALSE)
clone_CNA <- read.csv("clone_CNA.csv", stringsAsFactors = FALSE)
clone_CNA <- clone_CNA[, c("chr", "clone_id", "CNA")]
cn_changes <- merge(cn_changes, clone_CNA, by = c("chr", "clone_id"))

ks_test_clone_pairs_by_chr <- function(df) {
  results <- list()
  for (chr in sort(unique(df$chr))) {
    sub <- df[df$chr == chr, ]
    clones <- sort(unique(sub$clone_id))
    if (length(clones) < 2) next
    pairs <- combn(clones, 2, simplify = FALSE)
    for (pair in pairs) {
      x <- sub$total_reads[sub$clone_id == pair[1]]
      y <- sub$total_reads[sub$clone_id == pair[2]]
      if (length(x) < 1 || length(y) < 1) next
      test <- suppressWarnings(ks.test(x, y))
      results[[length(results) + 1]] <- data.frame(
        chr = chr,
        clone_a = pair[1],
        clone_b = pair[2],
        n_a = length(x),
        n_b = length(y),
        mean_a = round(mean(x), 2),
        median_a = round(median(x), 2),
        # median_a_per_cell = round(median(x) / length(x), 2),
        sd_a = round(sd(x), 2),
        mean_b = round(mean(y), 2),
        median_b = round(median(y), 2),
        # median_b_per_cell = round(median(y) / length(y), 2),
        sd_b = round(sd(y), 2),
        D = round(unname(test$statistic), 2),
        p_value = test$p.value,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, results)
}

plot_median_vs_cna <- function(df, out_file = "SA919_median_vs_cna.png") {
  agg <- aggregate(total_reads ~ chr + clone_id, data = df,
                   FUN = function(v) median(v) / length(v))
  names(agg)[names(agg) == "total_reads"] <- "median_per_cell"
  chr_levels <- unique(agg$chr)
  chr_levels <- chr_levels[order(suppressWarnings(as.numeric(chr_levels)), chr_levels)]
  agg$chr <- factor(agg$chr, levels = chr_levels)
  p <- ggplot(agg, aes(x = chr, y = median_per_cell, fill = clone_id)) +
    geom_col(position = "dodge") +
    labs(x = "Chromosome", y = "Median total_reads per cell", fill = "Clone") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(out_file, p, width = 8, height = 4, dpi = 150)
  p
}

ks_results <- ks_test_clone_pairs_by_chr(cn_changes)
ks_results$result <- ifelse(ks_results$p_value < 0.05, "significant", "not significant")
write.csv(ks_results, "SA919_ks_results.csv", row.names = FALSE)
plot_median_vs_cna(cn_changes, "SA919_median_vs_cna.png")

set.seed(42)
cells_C <- unique(cn_changes$cell_id[cn_changes$clone_id == "C"])
cells_A <- unique(cn_changes$cell_id[cn_changes$clone_id == "A"])
cells_B <- unique(cn_changes$cell_id[cn_changes$clone_id == "B"])

n_sample <- round(0.8 * length(cells_C))
sampled_C <- sample(cells_C, n_sample)
sampled_A <- sample(cells_A, n_sample)
sampled_B <- sample(cells_B, n_sample)

sampled_cells <- c(sampled_A, sampled_B, sampled_C)
cn_changes_sampled <- cn_changes[cn_changes$cell_id %in% sampled_cells, ]

ks_results_sampled <- ks_test_clone_pairs_by_chr(cn_changes_sampled)
ks_results_sampled$result <- ifelse(ks_results_sampled$p_value < 0.05, "significant", "not significant")
write.csv(ks_results_sampled, "SA919_ks_results_sampled.csv", row.names = FALSE)
