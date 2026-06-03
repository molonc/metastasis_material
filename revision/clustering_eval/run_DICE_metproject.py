#!/usr/bin/env python3
"""
Run DICE phylogenetic reconstruction and k-means clustering on
MetProject SA919 total copy number data.

Uses utility functions from generate_and_cluster.py for reading
data and running k-means with silhouette selection.

Usage:
    python run_DICE_metproject.py -i <input.tsv> [-o <output_dir>] [-m <method>] [-k <max_k>]

Dependencies:
    pip install numpy pandas scikit-learn biopython
"""

import argparse
import os

import numpy as np
import pandas as pd
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade
from io import StringIO

from generate_and_cluster import (
    read_total_cn_matrix,
    run_kmeans_silhouette,
)


def run_dice(input_tsv, output_dir, method="balME", dist_type="root",
             fastme_path="fastme"):
    """Run DICE on the input TSV file.

    Returns:
        tree_path: path to the output Newick tree file, or None if failed.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Find the output tree file
    for f in os.listdir(output_dir):
        if f.endswith("_tree.nwk"):
            tree_path = os.path.join(output_dir, f)
            size = os.path.getsize(tree_path)
            if size > 0:
                print(f"  Tree file: {tree_path} ({size} bytes)")
                return tree_path
            else:
                print(f"  Warning: tree file is empty: {tree_path}")

    print("  Warning: no tree file found in output")
    return None


def save_kmeans_labels(cell_names, best_k, best_labels, output_dir):
    """Save k-means cluster assignments to a TSV file."""
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, "kmeans_labels.tsv")

    df = pd.DataFrame({
        "cell": cell_names,
        "kmeans_cluster": best_labels,
    })
    df.to_csv(out_path, sep="\t", index=False)
    print(f"  Saved k-means labels (k={best_k}) to {out_path}\n")
    return out_path


def load_tree(nwk_path):
    """Load a Newick tree file and return a Bio.Phylo tree object."""
    tree = Phylo.read(nwk_path, "newick")
    tips = tree.get_terminals()
    internals = tree.get_nonterminals()
    print(f"  Loaded tree: {len(tips)} tips, {len(internals)} internal nodes")
    return tree


def save_tree(tree, out_path):
    """Save a tree to Newick format."""
    Phylo.write(tree, out_path, "newick")
    size = os.path.getsize(out_path)
    tips = len(tree.get_terminals())
    internals = len(tree.get_nonterminals())
    print(f"  Saved pruned tree: {out_path} ({size} bytes, {tips} tips, {internals} internal)")
    return out_path


def load_kmeans_labels(labels_path):
    """Load k-means cluster labels from TSV.

    Returns:
        df: DataFrame with columns cell, kmeans_cluster
        cell_to_cluster: dict mapping cell name -> cluster id
    """
    df = pd.read_csv(labels_path, sep="\t")
    cell_to_cluster = dict(zip(df["cell"], df["kmeans_cluster"]))
    return df, cell_to_cluster


def prune_short_branches(tree, threshold):
    """Collapse internal branches shorter than threshold.

    When an internal branch is shorter than the threshold, its children
    are moved up to its parent, effectively removing that branching level
    and reducing sub-structure.

    Args:
        tree: Bio.Phylo tree object (modified in place).
        threshold: minimum branch length to keep.

    Returns:
        n_collapsed: number of internal nodes collapsed.
    """
    n_collapsed = 0
    changed = True
    while changed:
        changed = False
        for node in tree.find_clades(order="level"):
            if node.is_terminal():
                continue
            children_to_remove = []
            children_to_add = []
            for child in node.clades:
                if (child.is_terminal() or child.branch_length is None):
                    continue
                if child.branch_length < threshold:
                    children_to_remove.append(child)
                    for grandchild in child.clades:
                        if grandchild.branch_length is not None and child.branch_length is not None:
                            grandchild.branch_length += child.branch_length
                        children_to_add.append(grandchild)
                    n_collapsed += 1
                    changed = True
            for c in children_to_remove:
                node.clades.remove(c)
            node.clades.extend(children_to_add)

    return n_collapsed


def _find_pure_subclades(node, cell_set):
    """Find maximal pure sub-clades under node where ALL tips belong to cell_set.

    A pure sub-clade is an internal node whose every descendant tip is in
    cell_set. 'Maximal' means we return the highest such node — we don't
    recurse into its children.

    Returns:
        list of (node, n_tips) tuples for each maximal pure sub-clade
        with >= 2 tips.
    """
    if node.is_terminal():
        return []

    tips_under = node.get_terminals()
    tip_names = {t.name for t in tips_under}

    # If all tips under this node belong to cell_set, it's pure
    if tip_names and tip_names <= cell_set:
        if len(tip_names) >= 2:
            return [(node, len(tip_names))]
        return []

    # Otherwise recurse into children
    result = []
    for child in node.clades:
        result.extend(_find_pure_subclades(child, cell_set))
    return result


def collapse_cluster_subtrees(tree, cell_to_cluster, min_clade_size=2):
    """Collapse pure sub-clades within each cluster into representative leaves.

    For each k-means cluster, finds all maximal pure sub-clades (internal
    nodes where every tip belongs to that cluster). Each such sub-clade is
    replaced by a single representative leaf named 'cluster_<id>_<index>'.
    Isolated (singleton) tips from the cluster remain as-is.

    This preserves the overall tree topology and inter-cluster branching
    while reducing within-cluster sub-structure.

    Args:
        tree: Bio.Phylo tree object (modified in place).
        cell_to_cluster: dict mapping cell name -> cluster id.
        min_clade_size: minimum tips in a pure sub-clade to collapse (default 2).

    Returns:
        collapse_stats: dict mapping cluster_id -> dict with keys:
            n_cells, n_pure_clades, n_collapsed_tips, n_singletons.
    """
    clusters = {}
    for cell, cid in cell_to_cluster.items():
        clusters.setdefault(cid, []).append(cell)

    collapse_stats = {}
    total_collapsed = 0

    for cid in sorted(clusters.keys()):
        cells = clusters[cid]
        cell_set = set(cells)

        pure_clades = _find_pure_subclades(tree.root, cell_set)
        # Filter by minimum size
        pure_clades = [(node, n) for node, n in pure_clades if n >= min_clade_size]

        n_collapsed_tips = 0
        for idx, (node, n_tips) in enumerate(pure_clades):
            # Compute mean branch length from this node to its tips
            depths = tree.depths()
            node_depth = depths.get(node, 0)
            tip_depths = [depths.get(t, 0) for t in node.get_terminals()]
            mean_bl = np.mean(tip_depths) - node_depth if tip_depths else 0

            # Replace subtree with single representative leaf
            node.clades = []
            rep_name = f"cluster_{cid}" if len(pure_clades) == 1 else f"cluster_{cid}_{idx}"
            rep_leaf = Clade(name=rep_name, branch_length=max(mean_bl, 0.001))
            node.clades.append(rep_leaf)
            n_collapsed_tips += n_tips

        # Count remaining singleton tips
        remaining_tips = [t for t in tree.get_terminals() if t.name in cell_set]
        n_singletons = len(remaining_tips)

        collapse_stats[cid] = {
            "n_cells": len(cells),
            "n_pure_clades": len(pure_clades),
            "n_collapsed_tips": n_collapsed_tips,
            "n_singletons": n_singletons,
        }
        total_collapsed += n_collapsed_tips

    return collapse_stats


def remove_outlier_cells(tree, cell_to_cluster):
    """Remove tree tips that are not assigned to any cluster.

    Identifies outlier cells in two ways:
      1. Tips with no entry in cell_to_cluster (unassigned cells).
      2. Singleton tips that are isolated from their cluster — cells
         whose nearest-neighbour tip in the tree belongs to a different
         cluster, indicating they are misplaced outliers.

    A cell is considered an outlier if its sibling tips (children of the
    same parent node) are ALL from different clusters. This catches cells
    that sit alone in a foreign part of the tree.

    Args:
        tree: Bio.Phylo tree object (modified in place).
        cell_to_cluster: dict mapping cell name -> cluster id.

    Returns:
        removed: list of (cell_name, reason) tuples for each removed tip.
    """
    removed = []

    # Pass 1: remove tips with no cluster assignment
    unassigned = [
        tip for tip in tree.get_terminals()
        if tip.name not in cell_to_cluster
    ]
    for tip in unassigned:
        tree.prune(tip)
        removed.append((tip.name, "unassigned"))

    if unassigned:
        print(f"  Removed {len(unassigned)} unassigned tips (no cluster label)")

    # Pass 2: remove isolated singleton tips whose siblings are all foreign
    changed = True
    pass_count = 0
    while changed:
        changed = False
        pass_count += 1
        tips_to_remove = []

        for tip in tree.get_terminals():
            if tip.name not in cell_to_cluster:
                continue
            my_cluster = cell_to_cluster[tip.name]

            # Find parent of this tip
            path = tree.get_path(tip)
            if len(path) < 2:
                continue
            parent = path[-2]

            # Get sibling tips (other tips under the same parent)
            sibling_tips = [
                t for t in parent.get_terminals()
                if t.name != tip.name and t.name in cell_to_cluster
            ]

            if not sibling_tips:
                continue

            # Check if ALL siblings belong to different clusters
            sibling_clusters = {cell_to_cluster[t.name] for t in sibling_tips}
            if my_cluster not in sibling_clusters:
                tips_to_remove.append(tip)

        for tip in tips_to_remove:
            try:
                tree.prune(tip)
                removed.append((tip.name, f"isolated_pass{pass_count}"))
                changed = True
            except ValueError:
                pass

        if tips_to_remove:
            print(f"  Pass {pass_count}: removed {len(tips_to_remove)} isolated outlier tips")

    # Summary
    if removed:
        cluster_counts = {}
        for name, reason in removed:
            cid = cell_to_cluster.get(name, -1)
            cluster_counts[cid] = cluster_counts.get(cid, 0) + 1
        print(f"  Total outliers removed: {len(removed)}")
        for cid in sorted(cluster_counts.keys()):
            label = f"cluster {cid}" if cid >= 0 else "unassigned"
            print(f"    {label}: {cluster_counts[cid]} cells")
    else:
        print("  No outlier cells found")

    return removed


def prune_minority_tips(tree, cell_to_cluster, min_cluster_size=3):
    """Remove tips belonging to very small clusters (outlier cells).

    Args:
        tree: Bio.Phylo tree object (modified in place).
        cell_to_cluster: dict mapping cell name -> cluster id.
        min_cluster_size: clusters with fewer cells are pruned.

    Returns:
        n_pruned: number of tips removed.
    """
    cluster_counts = {}
    for cid in cell_to_cluster.values():
        cluster_counts[cid] = cluster_counts.get(cid, 0) + 1

    small_clusters = {cid for cid, n in cluster_counts.items() if n < min_cluster_size}
    if not small_clusters:
        return 0

    tips_to_prune = [
        tip for tip in tree.get_terminals()
        if tip.name in cell_to_cluster and cell_to_cluster[tip.name] in small_clusters
    ]

    for tip in tips_to_prune:
        tree.prune(tip)

    print(f"  Pruned {len(tips_to_prune)} tips from {len(small_clusters)} "
          f"small clusters (< {min_cluster_size} cells)")
    return len(tips_to_prune)


def prune_tree(nwk_path, labels_path, output_dir,
               remove_outliers=False, branch_threshold=None,
               collapse_clusters=False, min_cluster_size=None):
    """Main pruning pipeline: load tree, apply pruning steps, save result.

    Pruning steps are applied in order:
      0. remove_outliers: remove unassigned and isolated outlier cells.
      1. min_cluster_size: remove tips from clusters smaller than this.
      2. branch_threshold: collapse internal branches shorter than this value.
      3. collapse_clusters: replace each cluster's pure sub-clades with
         single representative leaves.

    Args:
        nwk_path: path to input Newick tree.
        labels_path: path to kmeans_labels.tsv.
        output_dir: directory to save pruned tree.
        remove_outliers: if True, remove unassigned and isolated cells first.
        branch_threshold: branch length threshold for collapsing (None = skip).
        collapse_clusters: if True, collapse each cluster to one representative.
        min_cluster_size: minimum cluster size to keep (None = skip).

    Returns:
        pruned_tree_path: path to the saved pruned tree.
    """
    print(f"\nPruning tree: {nwk_path}")
    tree = load_tree(nwk_path)
    labels_df, cell_to_cluster = load_kmeans_labels(labels_path)

    if remove_outliers:
        removed = remove_outlier_cells(tree, cell_to_cluster)

    if min_cluster_size is not None:
        prune_minority_tips(tree, cell_to_cluster, min_cluster_size)

    if branch_threshold is not None:
        n = prune_short_branches(tree, branch_threshold)
        print(f"  Collapsed {n} internal nodes with branch length < {branch_threshold}")

    if collapse_clusters:
        stats = collapse_cluster_subtrees(tree, cell_to_cluster)
        print(f"  Collapsed pure sub-clades in {len(stats)} clusters:")
        for cid, s in sorted(stats.items()):
            print(f"    cluster_{cid}: {s['n_cells']} cells, "
                  f"{s['n_pure_clades']} pure clades collapsed "
                  f"({s['n_collapsed_tips']} tips -> {s['n_pure_clades']} leaves), "
                  f"{s['n_singletons']} singletons remain")

    suffix_parts = []
    if remove_outliers:
        suffix_parts.append("no_outliers")
    if branch_threshold is not None:
        suffix_parts.append(f"bl{branch_threshold}")
    if collapse_clusters:
        suffix_parts.append("collapsed")
    if min_cluster_size is not None:
        suffix_parts.append(f"min{min_cluster_size}")
    suffix = "_".join(suffix_parts) if suffix_parts else "pruned"

    out_path = os.path.join(output_dir, f"pruned_{suffix}_tree.nwk")
    save_tree(tree, out_path)
    return out_path


def print_summary(input_tsv, method, n_cells, n_bins, best_k, tree_path):
    """Print a summary of the run."""
    print("\n" + "=" * 60)
    print("Run Summary")
    print("=" * 60)
    print(f"  Input file:   {input_tsv}")
    print(f"  DICE method:  {method}")
    print(f"  Cells:        {n_cells}")
    print(f"  Bins:         {n_bins}")
    print(f"  Optimal k:    {best_k}")
    print(f"  Tree file:    {tree_path}")
    print("=" * 60)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run DICE and k-means on MetProject total CN data")
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="Path to total CN input TSV file")
    parser.add_argument("-o", "--output", type=str, default="DICE_SA919",
                        help="Output directory (default: DICE_SA919)")
    parser.add_argument("-m", "--method", type=str, default="balME",
                        choices=["NJ", "uNJ", "balME", "olsME"],
                        help="DICE reconstruction method (default: balME)")
    parser.add_argument("-d", "--dist-type", type=str, default="root",
                        choices=["root", "log", "manhattan", "euclidean"],
                        help="Distance type (default: root)")
    parser.add_argument("-k", "--max-k", type=int, default=15,
                        help="Maximum k for k-means silhouette selection (default: 15)")
    parser.add_argument("-f", "--fastme-path", type=str, default="fastme",
                        help="Path to fastme executable")
    parser.add_argument("--sep", type=str, default="\t",
                        help="Input file separator (default: tab)")
    return parser.parse_args()


def main():
    args = parse_args()

    input_tsv = args.input
    output_dir = args.output
    method = args.method
    run_label = f"{method}_{os.path.splitext(os.path.basename(input_tsv))[0]}"
    dice_output_dir = os.path.join(output_dir, run_label)

    # Step 1: Read total CN profiles into matrix
    cell_names, matrix = read_total_cn_matrix(input_tsv, sep=args.sep)

    # Step 2: Run k-means with silhouette selection
    k_range = range(2, args.max_k + 1)
    best_k, best_labels = run_kmeans_silhouette(matrix, k_range=k_range)

    # Step 3: Save k-means labels
    save_kmeans_labels(cell_names, best_k, best_labels, dice_output_dir)

    # Step 4: Run DICE
    tree_path = run_dice(
        input_tsv, dice_output_dir, method=method,
        dist_type=args.dist_type, fastme_path=args.fastme_path)

    # Step 5: Summary
    print_summary(input_tsv, method, len(cell_names), matrix.shape[1],
                  best_k, tree_path)


if __name__ == "__main__":
    main()
