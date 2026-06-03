#!/usr/bin/env python3
"""
Read Sitka clone labels and tree, and prepare ground-truth file for evaluation.

Usage:
    python read_sitka_labels.py <sitka_csv_gz> [<output_dir>]

Example:
    python read_sitka_labels.py DICE_SA919/SA919_0.5_batch/sitka/SA919_cell_clones_0.5.csv.gz DICE_SA919/SA919_0.5_batch/sitka/

Dependencies:
    pip install biopython pandas
"""

import sys
import os
import pandas as pd
from Bio import Phylo


def read_sitka_labels(sitka_path):
    """Read Sitka cell clone assignments from a .csv.gz file.

    Args:
        sitka_path: Path to .csv.gz with columns cell_id and clone_id.

    Returns:
        DataFrame with columns: cell_id, clone_id
    """
    df = pd.read_csv(sitka_path, compression='gzip')
    print(f"Read {len(df)} cells from {sitka_path}")
    print(f"Columns: {list(df.columns)}")
    print(f"\nNumber of clones: {df['clone_id'].nunique()}")
    print(f"\nClone distribution:")
    clone_counts = df['clone_id'].value_counts().sort_index()
    for clone_id, count in clone_counts.items():
        print(f"  Clone {clone_id}: {count} cells")
    return df


def read_sitka_tree(tree_path):
    """Read and summarize a Sitka Newick tree.

    Args:
        tree_path: Path to tree.newick file.

    Returns:
        Bio.Phylo tree object.
    """
    tree = Phylo.read(tree_path, "newick")
    tips = tree.get_terminals()
    internal = tree.get_nonterminals()
    print(f"\nSitka tree: {tree_path}")
    print(f"  Tips: {len(tips)}")
    print(f"  Internal nodes: {len(internal)}")

    # Check if tree is multifurcating
    max_children = max(len(node.clades) for node in internal)
    if max_children > 2:
        print(f"  Tree is multifurcating (max children per node: {max_children})")
    else:
        print(f"  Tree is binary")

    return tree


def check_label_tree_overlap(sitka_df, tree):
    """Check overlap between Sitka clone labels and tree tip names.

    Args:
        sitka_df: DataFrame with cell_id column.
        tree: Bio.Phylo tree object.
    """
    label_cells = set(sitka_df['cell_id'])
    tree_tips = {tip.name for tip in tree.get_terminals()}

    in_both = label_cells & tree_tips
    in_labels_only = label_cells - tree_tips
    in_tree_only = tree_tips - label_cells

    print(f"\nLabel-tree overlap:")
    print(f"  Cells in both:        {len(in_both)}")
    print(f"  In labels only:       {len(in_labels_only)}")
    print(f"  In tree only:         {len(in_tree_only)}")

    if in_labels_only:
        print(f"  Sample label-only:    {list(in_labels_only)[:3]}")
    if in_tree_only:
        print(f"  Sample tree-only:     {list(in_tree_only)[:3]}")


def trim_tree_to_labels(tree, keep_cells, output_path=None):
    """Trim tree by keeping only tips present in the labels file.

    Args:
        tree: Bio.Phylo tree object.
        keep_cells: set of cell names to keep.
        output_path: Path to save the trimmed tree (optional).

    Returns:
        Trimmed Bio.Phylo tree object.
    """
    all_tips = tree.get_terminals()
    tips_to_remove = [tip for tip in all_tips if tip.name not in keep_cells]

    print(f"\nTrimming tree:")
    print(f"  Original tips:   {len(all_tips)}")
    print(f"  Tips to keep:    {len(all_tips) - len(tips_to_remove)}")
    print(f"  Tips to remove:  {len(tips_to_remove)}")

    # Remove tips not in the labels
    for tip in tips_to_remove:
        tree.prune(tip)

    # Collapse single-child internal nodes
    tree = collapse_singles(tree)

    remaining_tips = tree.get_terminals()
    remaining_internal = tree.get_nonterminals()
    print(f"  Trimmed tips:    {len(remaining_tips)}")
    print(f"  Trimmed internal nodes: {len(remaining_internal)}")

    if output_path:
        Phylo.write(tree, output_path, "newick")
        print(f"  Saved trimmed tree: {output_path}")

    return tree


def collapse_singles(tree):
    """Collapse internal nodes that have only one child."""
    def _collapse(clade):
        # Recurse first
        for child in list(clade.clades):
            _collapse(child)
        # If this node has exactly one child, merge child into parent
        while len(clade.clades) == 1:
            child = clade.clades[0]
            # Merge branch lengths
            if clade.branch_length and child.branch_length:
                child.branch_length += clade.branch_length
            elif clade.branch_length:
                child.branch_length = clade.branch_length
            # Replace this node's attributes with child's
            clade.clades = child.clades
            clade.name = child.name if child.is_terminal() else clade.name

    _collapse(tree.root)
    return tree


def load_kmeans_gt(kmeans_gt_path, sitka_df):
    """Load k-means ground truth and copy it to sitka output for evaluation.

    Uses k-means cluster labels as ground truth (true_clone) instead of
    Sitka's own clone labels.

    Args:
        kmeans_gt_path: Path to kmeans gt.tsv (cell, true_clone, kmeans_cluster).
        sitka_df: DataFrame with cell_id column (used to filter to Sitka cells).

    Returns:
        gt_df: DataFrame with columns cell, true_clone, kmeans_cluster.
    """
    gt_df = pd.read_csv(kmeans_gt_path, sep='\t')
    print(f"\nLoaded k-means ground truth from {kmeans_gt_path}")
    print(f"  Total cells in kmeans gt: {len(gt_df)}")

    # Filter to cells present in Sitka labels
    sitka_cells = set(sitka_df['cell_id'])
    gt_df = gt_df[gt_df['cell'].isin(sitka_cells)].copy()
    print(f"  Cells matching Sitka labels: {len(gt_df)}")

    print(f"\n  K-means clone distribution:")
    clone_counts = gt_df['true_clone'].value_counts().sort_index()
    for clone, count in clone_counts.items():
        print(f"    {clone}: {count} cells")

    return gt_df


def save_gt_for_sitka(gt_df, sitka_df, output_dir):
    """Save ground truth and Sitka labels for evaluation.

    Includes original Sitka clone_id (e.g., A, B, C) alongside k-means labels.

    Args:
        gt_df: DataFrame with cell, true_clone, kmeans_cluster columns.
        sitka_df: DataFrame with cell_id and clone_id columns.
        output_dir: Directory to save output files.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Map cell_id -> clone_id from Sitka labels
    sitka_clone_map = sitka_df.set_index('cell_id')['clone_id']

    # Save sitka_labels.tsv (cell, kmeans_cluster, sitka_clone_id)
    labels_df = gt_df[['cell', 'kmeans_cluster']].copy()
    labels_df['sitka_clone_id'] = labels_df['cell'].map(sitka_clone_map)
    labels_path = os.path.join(output_dir, 'sitka_labels.tsv')
    labels_df.to_csv(labels_path, sep='\t', index=False)
    print(f"\nSaved {labels_path}")

    # Save sitka_gt.tsv (cell, true_clone, kmeans_cluster, sitka_clone_id)
    gt_out = gt_df.copy()
    gt_out['sitka_clone_id'] = gt_out['cell'].map(sitka_clone_map)
    gt_path = os.path.join(output_dir, 'sitka_gt.tsv')
    gt_out.to_csv(gt_path, sep='\t', index=False)
    print(f"Saved {gt_path}")


def main():
    if len(sys.argv) < 2:
        print("Usage: python read_sitka_labels.py <sitka_csv_gz> <kmeans_gt_path> [<output_dir>]")
        sys.exit(1)

    sitka_path = sys.argv[1]
    kmeans_gt_path = sys.argv[2] if len(sys.argv) > 2 else None
    output_dir = sys.argv[3] if len(sys.argv) > 3 else os.path.dirname(sitka_path)

    # Read Sitka clone labels
    sitka_df = read_sitka_labels(sitka_path)

    # Read tree from same directory
    sitka_dir = os.path.dirname(sitka_path)
    tree_path = os.path.join(sitka_dir, 'tree.newick')
    tree = None
    if os.path.exists(tree_path):
        tree = read_sitka_tree(tree_path)
        check_label_tree_overlap(sitka_df, tree)

        # Trim tree to keep only labeled cells
        keep_cells = set(sitka_df['cell_id'])
        trimmed_tree_path = os.path.join(output_dir, 'trimmed_tree.newick')
        trim_tree_to_labels(tree, keep_cells, output_path=trimmed_tree_path)
    else:
        print(f"\nNo tree.newick found in {sitka_dir}")

    # Load k-means ground truth and save for evaluation
    if kmeans_gt_path and os.path.exists(kmeans_gt_path):
        gt_df = load_kmeans_gt(kmeans_gt_path, sitka_df)
        save_gt_for_sitka(gt_df, sitka_df, output_dir)
    else:
        print(f"\nNo k-means gt.tsv provided or not found: {kmeans_gt_path}")


if __name__ == "__main__":
    main()
