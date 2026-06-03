#!/usr/bin/env Rscript
# ============================================================
# Read a Newick tree and visualize it in R, colored by cluster
# ============================================================
# Required packages:
#   install.packages("ape")
#   install.packages("BiocManager")
#   BiocManager::install("ggtree")
# ============================================================

library(ape)

# ---- 0. Paths and data -------------------------------------
input_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/revision/clustering_eval/DICE_SA919/SA919_0.5_batch/sitka/'
newick_fn <- paste0(input_dir, 'trimmed_tree.newick')
labels_fn <- paste0(input_dir, 'sitka_labels.tsv')

tree <- read.tree(newick_fn)
cat("Tree summary:\n")
print(tree)
cat("\nNumber of tips:", length(tree$tip.label), "\n\n")

# ---- 0b. Build tip-to-cluster mapping ----------------------
# Collapsed leaves: "cluster_0", "cluster_0_1", etc. -> extract cluster id
# Singleton cells: look up from kmeans_labels.tsv
labels_df <- read.delim(labels_fn, stringsAsFactors = FALSE)
cell_to_cluster <- setNames(labels_df$kmeans_cluster, labels_df$cell)

get_cluster_id <- function(tip_label) {
  # Collapsed representative: cluster_<id> or cluster_<id>_<idx>
  m <- regmatches(tip_label, regexec("^cluster_(\\d+)", tip_label))
  if (length(m[[1]]) > 1) return(as.integer(m[[1]][2]))
  # Singleton cell: look up in labels
  cid <- cell_to_cluster[tip_label]
  if (!is.na(cid)) return(as.integer(cid))
  return(NA)
}

tip_clusters <- sapply(tree$tip.label, get_cluster_id)
cluster_ids <- sort(unique(tip_clusters[!is.na(tip_clusters)]))
n_clusters <- length(cluster_ids)

cat("Cluster assignments:\n")
for (cid in cluster_ids) {
  n <- sum(tip_clusters == cid, na.rm = TRUE)
  cat(sprintf("  Cluster %d: %d tips\n", cid, n))
}
n_na <- sum(is.na(tip_clusters))
if (n_na > 0) cat(sprintf("  Unassigned: %d tips\n", n_na))

# ---- 0c. Define cluster colors ------------------------------
cluster_palette <- c(
  "#E41A1C",  # red       - cluster 0
  "#377EB8",  # blue      - cluster 1
  "#4DAF4A",  # green     - cluster 2
  "#FF7F00",  # orange    - cluster 3
  "#984EA3",  # purple    - cluster 4
  "#A65628",  # brown     - cluster 5
  "#F781BF",  # pink      - cluster 6
  "#999999",  # grey      - cluster 7
  "#66C2A5",  # teal      - cluster 8
  "#E7298A"   # magenta   - cluster 9
)
# Extend palette if needed
if (n_clusters > length(cluster_palette)) {
  cluster_palette <- rep_len(cluster_palette, n_clusters)
}

get_tip_color <- function(cid) {
  if (is.na(cid)) return("black")
  idx <- match(cid, cluster_ids)
  return(cluster_palette[idx])
}

tip_colors <- sapply(tip_clusters, get_tip_color)

# Helper: get descendant tip indices for an internal node
get_descendant_tips <- function(tree, node) {
  n_tips <- length(tree$tip.label)
  desc <- node
  tips_found <- c()
  stack <- node
  while (length(stack) > 0) {
    current <- stack[1]
    stack <- stack[-1]
    children <- tree$edge[tree$edge[, 1] == current, 2]
    for (ch in children) {
      if (ch <= n_tips) {
        tips_found <- c(tips_found, ch)
      } else {
        stack <- c(stack, ch)
      }
    }
  }
  return(tips_found)
}

# Map edge colors: color an edge by the cluster of tips below it
# If all tips below an edge belong to one cluster, color the edge
get_edge_colors <- function(tree, tip_clusters, cluster_ids) {
  n_tips <- length(tree$tip.label)
  edge_colors <- rep("grey30", nrow(tree$edge))

  for (i in seq_len(nrow(tree$edge))) {
    child <- tree$edge[i, 2]
    if (child <= n_tips) {
      edge_colors[i] <- tip_colors[child]
    } else {
      desc_tips <- get_descendant_tips(tree, child)
      desc_clusters <- unique(tip_clusters[desc_tips])
      desc_clusters <- desc_clusters[!is.na(desc_clusters)]
      if (length(desc_clusters) == 1) {
        edge_colors[i] <- get_tip_color(desc_clusters)
      }
    }
  }
  return(edge_colors)
}

edge_colors <- get_edge_colors(tree, tip_clusters, cluster_ids)

# ---- 1. Colored ape plot ------------------------------------
plot_tree_ape_colored <- function(tree, tip_colors, edge_colors,
                                  outfile = "tree_colored.png",
                                  width = 16, height = 10, res = 150) {
  if (!is.null(outfile)) png(outfile, width = width, height = height,
                             units = "in", res = res)

  par(mfrow = c(1, 2), mar = c(2, 1, 2, 4))

  # Cladogram (ignore branch lengths for better layout)
  show_labels <- length(tree$tip.label) <= 200
  plot(tree, main = "Cladogram — colored by cluster",
       show.tip.label = show_labels, cex = 0.2, edge.width = 0.3,
       tip.color = tip_colors, edge.color = edge_colors,
       use.edge.length = FALSE, label.offset = 0.5)

  # Fan layout (cladogram)
  plot(tree, type = "fan", main = "Fan layout — colored by cluster",
       show.tip.label = FALSE, cex = 0.2, edge.width = 0.3,
       tip.color = tip_colors, edge.color = edge_colors,
       use.edge.length = FALSE, label.offset = 2)

  # Add legend
  legend("bottomright",
         legend = paste("Cluster", cluster_ids),
         col = cluster_palette[seq_along(cluster_ids)],
         pch = 15, pt.cex = 1.5, cex = 0.8,
         bty = "n", title = "K-means cluster")

  if (!is.null(outfile)) {
    dev.off()
    cat("Saved:", outfile, "\n")
  }
}

ape_output <- paste0(input_dir, 'SA919_sitka_trimmed_tree_colored.png')
plot_tree_ape_colored(tree, tip_colors, edge_colors, outfile = ape_output)

# ---- 2. Colored ggtree plot --------------------------------
plot_tree_ggtree_colored <- function(tree, tip_clusters, cluster_ids,
                                      outfile = "tree_ggtree_colored.png",
                                      width = 18, height = 10, res = 150) {
  if (!requireNamespace("ggtree", quietly = TRUE)) {
    cat("\nTip: install ggtree for richer visualisations:\n")
    cat('  install.packages("BiocManager")\n')
    cat('  BiocManager::install("ggtree")\n')
    return(invisible(NULL))
  }

  library(ggtree)
  library(ggplot2)

  # Build tip annotation data frame
  tip_df <- data.frame(
    label = tree$tip.label,
    cluster = factor(ifelse(is.na(tip_clusters), "NA",
                            paste0("Cluster ", tip_clusters)),
                     levels = c(paste0("Cluster ", cluster_ids), "NA")),
    stringsAsFactors = FALSE
  )

  color_vals <- cluster_palette[seq_along(cluster_ids)]
  names(color_vals) <- paste0("Cluster ", cluster_ids)
  if (any(is.na(tip_clusters))) {
    color_vals["NA"] <- "black"
  }

  # Rectangular cladogram with colored tip points (ignore branch lengths)
  n_tips <- length(tree$tip.label)
  p1 <- ggtree(tree, linewidth = 0.3, branch.length = "none") %<+% tip_df +
    geom_tippoint(aes(color = cluster), size = ifelse(n_tips > 200, 0.3, 1.5)) +
    scale_color_manual(values = color_vals, name = "Cluster") +
    theme_tree2() +
    theme(legend.position = "right") +
    ggtitle("Sitka tree (cladogram) — colored by clone")
  if (n_tips <= 200) {
    p1 <- p1 + geom_tiplab(aes(color = cluster), size = 1.8, offset = 0.5)
  }

  # # Circular tree with colored tip points
  # p2 <- ggtree(tree, layout = "circular", linewidth = 0.3) %<+% tip_df +
  #   geom_tippoint(aes(color = cluster), size = 1.5) +
  #   scale_color_manual(values = color_vals, name = "Cluster") +
  #   theme(legend.position = "right") +
  #   ggtitle("Circular tree — colored by cluster")

  if (!is.null(outfile)) {
    ggsave(outfile, plot = p1, width = width, height = height, dpi = res)
    cat("Saved:", outfile, "\n")
    # Also save circular
    # circ_outfile <- sub("\\.png$", "_circular.png", outfile)
    # ggsave(circ_outfile, plot = p2, width = height + 2, height = height, dpi = res)
    # cat("Saved:", circ_outfile, "\n")
  } else {
    print(p1)
    # print(p2)
  }
}

ggtree_output <- paste0(input_dir, 'SA919_sitka_trimmed_tree_colored_ggtree.png')
plot_tree_ggtree_colored(tree, tip_clusters, cluster_ids, outfile = ggtree_output)

# ---- 3. Colored ggplot2 tree --------------------------------
plot_tree_ggplot2_colored <- function(tree, tip_clusters, cluster_ids,
                                      outfile = "tree_ggplot2_colored.png",
                                      width = 14, height = 10) {
  library(ggplot2)

  tree_cs <- ape::collapse.singles(tree)

  # Recompute tip_clusters for collapsed tree (order may shift)
  tc_cs <- sapply(tree_cs$tip.label, get_cluster_id)

  coords <- ape::plotPhyloCoor(tree_cs, direction = "rightwards",
                               use.edge.length = FALSE)
  n_tips  <- length(tree_cs$tip.label)
  n_nodes <- tree_cs$Nnode + n_tips

  edges <- data.frame(
    parent = tree_cs$edge[, 1],
    child  = tree_cs$edge[, 2]
  )
  edges$x_start <- coords[edges$parent, 1]
  edges$y_start <- coords[edges$parent, 2]
  edges$x_end   <- coords[edges$child,  1]
  edges$y_end   <- coords[edges$child,  2]

  # Compute edge cluster: if child is a tip, use its cluster;
  # if child is internal, check if all descendant tips share one cluster
  tc_colors_cs <- sapply(tc_cs, get_tip_color)
  edge_cluster_color <- sapply(seq_len(nrow(edges)), function(i) {
    child <- edges$child[i]
    if (child <= n_tips) {
      return(tc_colors_cs[child])
    }
    desc <- get_descendant_tips(tree_cs, child)
    desc_cl <- unique(tc_cs[desc])
    desc_cl <- desc_cl[!is.na(desc_cl)]
    if (length(desc_cl) == 1) return(get_tip_color(desc_cl))
    return("grey30")
  })
  edges$edge_color <- edge_cluster_color

  # Tips data frame
  tips_df <- data.frame(
    label   = tree_cs$tip.label,
    x       = coords[1:n_tips, 1],
    y       = coords[1:n_tips, 2],
    cluster = factor(ifelse(is.na(tc_cs), "NA",
                            paste0("Cluster ", tc_cs)),
                     levels = c(paste0("Cluster ", cluster_ids), "NA"))
  )

  color_vals <- cluster_palette[seq_along(cluster_ids)]
  names(color_vals) <- paste0("Cluster ", cluster_ids)
  if (any(is.na(tc_cs))) {
    color_vals["NA"] <- "black"
  }

  p <- ggplot() +
    # Vertical connectors (colored by edge cluster)
    geom_segment(
      data = edges,
      aes(x = x_start, y = y_start, xend = x_start, yend = y_end),
      color = edges$edge_color, linewidth = 0.4
    ) +
    # Horizontal branches (colored by edge cluster)
    geom_segment(
      data = edges,
      aes(x = x_start, y = y_end, xend = x_end, yend = y_end),
      color = edges$edge_color, linewidth = 0.5
    ) +
    # Colored tip points
    geom_point(
      data = tips_df,
      aes(x = x, y = y, color = cluster),
      size = 2
    ) +
    scale_color_manual(values = color_vals, name = "Cluster") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.15))) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "right"
    ) +
    labs(x = "",
         title = "SA919 Sitka tree (cladogram) — colored by clone")

  # Add tip labels only for small trees
  if (n_tips <= 200) {
    p <- p + geom_text(
      data = tips_df,
      aes(x = x, y = y, label = label, color = cluster),
      hjust = -0.1, size = 1.8
    )
  }

  if (!is.null(outfile)) {
    ggsave(outfile, plot = p, width = width, height = height)
    cat("Saved:", outfile, "\n")
  } else {
    print(p)
  }
}

ggplot2_output <- paste0(input_dir, 'SA919_sitka_trimmed_tree_colored_ggplot2.png')
plot_tree_ggplot2_colored(tree, tip_clusters, cluster_ids, outfile = ggplot2_output)
