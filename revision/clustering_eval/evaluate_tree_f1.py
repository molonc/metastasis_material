#!/usr/bin/env python3
"""
Compute per-clone and overall F1 scores between a DICE phylogenetic tree
and ground-truth clone labels, following the evaluation method from the
DICE paper (Weiner et al.):

    For each ground-truth clone (a subset of cells), identify the clade
    in the reconstructed tree that shows maximum F1 score with respect
    to that clone. Report per-clone and mean F1.

Dependencies:
    pip install biopython pandas
"""

import sys
import os
import pandas as pd
import numpy as np
from Bio import Phylo
from sklearn.metrics import adjusted_rand_score
from scipy.optimize import linear_sum_assignment


def load_tree(nwk_path):
    """Parse a Newick tree file and return a Bio.Phylo tree object."""
    tree = Phylo.read(nwk_path, "newick")
    for tip in tree.get_terminals():
        if tip.name and tip.name.startswith("cell_"):
            tip.name = tip.name[len("cell_"):]
    return tree


def load_ground_truth(gt_path):
    """Load ground-truth labels and k-means clusters from gt.tsv.

    Returns:
        gt_df: DataFrame with columns cell, true_clone, kmeans_cluster
        clone_to_cells: dict mapping clone name -> set of cell names
    """
    df = pd.read_csv(gt_path, sep="\t")
    clone_to_cells = df.groupby("true_clone")["cell"].apply(set).to_dict()
    return df, clone_to_cells


def get_all_clades(tree):
    """Extract the set of leaf names for every internal node (clade) in the tree.

    Returns:
        list of sets, where each set contains the leaf names under one clade.
    """
    clades = []
    for node in tree.find_clades(order="postorder"):
        leaves = {leaf.name for leaf in node.get_terminals()}
        if len(leaves) > 0:
            clades.append(leaves)
    return clades


def compute_f1(pred_set, true_set):
    """Compute precision, recall, and F1 between a predicted clade and a true clone."""
    tp = len(pred_set & true_set)
    if tp == 0:
        return 0.0, 0.0, 0.0
    precision = tp / len(pred_set)
    recall = tp / len(true_set)
    f1 = 2 * precision * recall / (precision + recall)
    return precision, recall, f1


def compute_best_f1_per_clone(clades, clone_to_cells):
    """For each ground-truth clone, find the clade with the maximum F1 score.

    Returns:
        results: list of dicts with clone name, best F1, precision, recall,
                 and the matching clade size.
    """
    results = []
    for clone_name, true_cells in sorted(clone_to_cells.items()):
        best_f1 = 0.0
        best_precision = 0.0
        best_recall = 0.0
        best_clade_size = 0

        for clade_leaves in clades:
            precision, recall, f1 = compute_f1(clade_leaves, true_cells)
            if f1 > best_f1:
                best_f1 = f1
                best_precision = precision
                best_recall = recall
                best_clade_size = len(clade_leaves)

        results.append({
            "clone": clone_name,
            "n_true_cells": len(true_cells),
            "best_clade_size": best_clade_size,
            "precision": best_precision,
            "recall": best_recall,
            "f1": best_f1,
        })

    return results


def get_best_matching_clades(clades, clone_to_cells):
    """Find the best F1-matching clade for each ground-truth clone.

    For each clone, pick the clade that maximizes F1. Then resolve
    overlaps by assigning each cell to the clone whose best clade
    has the highest F1 for that cell.

    Returns:
        list of sets, each containing leaf names for one tree-induced cluster.
        clone_order: list of clone names in the same order.
    """
    all_cells = set()
    for cells in clone_to_cells.values():
        all_cells |= cells

    # Find best clade per clone
    best_clades = {}
    for clone_name, true_cells in clone_to_cells.items():
        best_f1 = 0.0
        best_clade = set()
        for clade_leaves in clades:
            _, _, f1 = compute_f1(clade_leaves, true_cells)
            if f1 > best_f1:
                best_f1 = f1
                best_clade = clade_leaves
        best_clades[clone_name] = (best_clade, best_f1)

    # Assign each cell to the clone whose best-matching clade contains it
    # and has the highest F1 (resolves overlaps)
    clone_order = sorted(clone_to_cells.keys())
    cluster_sets = {c: set() for c in clone_order}

    for cell in all_cells:
        best_clone = None
        best_f1 = -1
        for clone_name in clone_order:
            clade, f1 = best_clades[clone_name]
            if cell in clade and f1 > best_f1:
                best_f1 = f1
                best_clone = clone_name
        if best_clone is None:
            # Cell not in any best clade — assign to clone with nearest clade
            best_clone = clone_order[0]
        cluster_sets[best_clone].add(cell)

    clade_sets = [cluster_sets[c] for c in clone_order]
    return clade_sets, clone_order


def assign_tree_clade_labels(cell_names, clade_sets):
    """Assign each cell a clade label (integer) based on tree-induced clades.

    Returns:
        numpy array of integer labels, one per cell (in cell_names order).
    """
    cell_to_label = {}
    for label, leaves in enumerate(clade_sets):
        for cell in leaves:
            cell_to_label[cell] = label

    labels = np.array([cell_to_label[c] for c in cell_names])
    return labels


def compute_ari_kmeans_vs_tree(gt_df, clades, clone_to_cells):
    """Compute ARI between k-means clusters and tree-induced clades.

    Uses the best F1-matching clades to induce a tree partition, assigns
    each cell a tree-clade label, then computes ARI against k-means labels.

    Returns:
        ari: adjusted rand index
        tree_labels: array of tree-clade assignments
        clade_sets: list of sets of leaf names per tree-induced cluster
        clone_order: list of clone names matching clade_sets order
    """
    cell_names = gt_df["cell"].tolist()
    kmeans_labels = gt_df["kmeans_cluster"].values

    clade_sets, clone_order = get_best_matching_clades(clades, clone_to_cells)
    tree_labels = assign_tree_clade_labels(cell_names, clade_sets)

    ari = adjusted_rand_score(kmeans_labels, tree_labels)
    return ari, tree_labels, clade_sets, clone_order


def print_ari_report(ari, clade_sets, clone_order):
    """Print the ARI report between k-means and tree clades."""
    print(f"Clustering quality: ARI (k-means vs tree clades, k={len(clade_sets)})")
    print("-" * 50)

    for clone_name, clade in zip(clone_order, clade_sets):
        print(f"  Tree clade ({clone_name}): {len(clade)} cells")

    print(f"\n  ARI = {ari:.4f}")

    if ari == 1.0:
        print("  (Perfect agreement)")
    elif ari > 0.8:
        print("  (Strong agreement)")
    elif ari > 0.5:
        print("  (Moderate agreement)")
    else:
        print("  (Weak agreement)")
    print()


def get_bipartitions(tree):
    """Extract the set of bipartitions (splits) from a tree.

    Each bipartition is represented as a frozenset of the smaller side's
    leaf names. Trivial splits (single leaf vs rest) are excluded.

    Returns:
        set of frozensets, each representing one non-trivial bipartition.
    """
    all_leaves = frozenset(leaf.name for leaf in tree.get_terminals())
    bipartitions = set()

    for node in tree.find_clades(order="postorder"):
        if node.is_terminal():
            continue
        clade_leaves = frozenset(leaf.name for leaf in node.get_terminals())
        # Skip the root (whole tree) — it's not an informative split
        if clade_leaves == all_leaves:
            continue
        # Normalize: always store the smaller side
        complement = all_leaves - clade_leaves
        if len(clade_leaves) == 1 or len(complement) == 1:
            continue  # skip trivial splits
        smaller = clade_leaves if len(clade_leaves) <= len(complement) else complement
        bipartitions.add(smaller)

    return bipartitions


def build_reference_tree(clone_to_cells):
    """Build a star-topology reference tree from ground-truth clone labels.

    Each clone forms a clade (subtree). The clades are joined at the root
    in a star topology. Returns a Bio.Phylo tree.
    """
    from Bio.Phylo.BaseTree import Tree, Clade

    clone_clades = []
    for clone_name in sorted(clone_to_cells.keys()):
        cells = sorted(clone_to_cells[clone_name])
        leaf_clades = [Clade(name=c, branch_length=1.0) for c in cells]
        clone_clade = Clade(clades=leaf_clades, branch_length=1.0)
        clone_clades.append(clone_clade)

    root = Clade(clades=clone_clades)
    ref_tree = Tree(root=root)
    return ref_tree


def compute_nrfd(tree, clone_to_cells):
    """Compute Normalized Robinson-Foulds Distance between the inferred tree
    and a reference tree built from ground-truth clone labels.

    NRFD = |S1 symmetric_diff S2| / (|S1| + |S2|)

    where S1, S2 are the sets of non-trivial bipartitions of each tree.
    NRFD = 0 means identical topology, NRFD = 1 means maximally different.

    Returns:
        nrfd: normalized RF distance
        n_shared: number of shared bipartitions
        n_inferred: number of non-trivial bipartitions in inferred tree
        n_reference: number of non-trivial bipartitions in reference tree
    """
    ref_tree = build_reference_tree(clone_to_cells)

    bp_inferred = get_bipartitions(tree)
    bp_reference = get_bipartitions(ref_tree)

    shared = bp_inferred & bp_reference
    sym_diff = len(bp_inferred ^ bp_reference)
    total = len(bp_inferred) + len(bp_reference)

    nrfd = sym_diff / total if total > 0 else 0.0

    return nrfd, len(shared), len(bp_inferred), len(bp_reference)


def print_nrfd_report(nrfd, n_shared, n_inferred, n_reference):
    """Print the Normalized Robinson-Foulds Distance report."""
    print("Tree accuracy: Normalized Robinson-Foulds Distance (NRFD)")
    print("-" * 55)
    print(f"  Bipartitions in inferred tree:  {n_inferred}")
    print(f"  Bipartitions in reference tree: {n_reference}")
    print(f"  Shared bipartitions:            {n_shared}")
    print(f"\n  NRFD = {nrfd:.4f}")

    if nrfd == 0.0:
        print("  (Identical topology)")
    elif nrfd < 0.2:
        print("  (Very similar)")
    elif nrfd < 0.5:
        print("  (Moderately similar)")
    else:
        print("  (Dissimilar)")
    print()


def print_f1_report(results):
    """Print a formatted F1 report."""
    print("\nPer-clone F1 scores (best matching clade):")
    print("-" * 75)
    print(f"{'Clone':<16} {'True cells':>10} {'Clade size':>11} "
          f"{'Precision':>10} {'Recall':>8} {'F1':>8}")
    print("-" * 75)

    f1_scores = []
    for r in results:
        print(f"{r['clone']:<16} {r['n_true_cells']:>10} {r['best_clade_size']:>11} "
              f"{r['precision']:>10.4f} {r['recall']:>8.4f} {r['f1']:>8.4f}")
        f1_scores.append(r["f1"])

    mean_f1 = sum(f1_scores) / len(f1_scores)
    print("-" * 75)
    print(f"{'Mean F1':<16} {'':<10} {'':<11} {'':<10} {'':<8} {mean_f1:>8.4f}")
    print()

    return mean_f1


def find_best_clone_matching(gt_df, clades, clone_to_cells):
    """Find the best one-to-one matching between GT clones and DICE tree clades.

    Uses the Hungarian algorithm on an F1-score cost matrix to find the
    optimal assignment that maximizes total F1 across all clone-clade pairs.

    Steps:
        1. Get tree-induced clade sets (one per GT clone) via best F1 matching.
        2. Build an F1 matrix: rows = GT clones, cols = tree-induced clades.
        3. Run Hungarian algorithm on the negated F1 matrix (minimization).
        4. Produce per-cell label assignments and a confusion matrix.

    Args:
        gt_df: DataFrame with columns cell, true_clone, kmeans_cluster.
        clades: list of sets from get_all_clades().
        clone_to_cells: dict mapping clone name -> set of cell names.

    Returns:
        matching_df: DataFrame with columns gt_clone, matched_dice_clade,
                     n_gt_cells, n_clade_cells, overlap, precision, recall, f1.
        cell_labels_df: DataFrame with columns cell, gt_clone, dice_clade,
                        match (True if gt_clone's matched clade == dice_clade).
        confusion_df: DataFrame confusion matrix (rows=GT, cols=DICE clade).
        accuracy: fraction of cells whose DICE clade matches their GT clone's
                  assigned clade.
    """
    cell_names = gt_df["cell"].tolist()

    # Get tree-induced clade sets and clone ordering
    clade_sets, clone_order = get_best_matching_clades(clades, clone_to_cells)
    n_clones = len(clone_order)

    # Assign each cell a DICE clade label
    tree_labels = assign_tree_clade_labels(cell_names, clade_sets)

    # Build F1 matrix: gt_clones (rows) x dice_clades (cols)
    f1_matrix = np.zeros((n_clones, n_clones))
    for i, gt_clone in enumerate(clone_order):
        gt_cells = clone_to_cells[gt_clone]
        for j in range(n_clones):
            _, _, f1 = compute_f1(clade_sets[j], gt_cells)
            f1_matrix[i, j] = f1

    # Hungarian algorithm: maximize F1 -> minimize negative F1
    row_ind, col_ind = linear_sum_assignment(-f1_matrix)

    # Build the matching result
    gt_to_dice = {}  # gt_clone_name -> dice_clade_index
    matching_rows = []
    for r, c in zip(row_ind, col_ind):
        gt_clone = clone_order[r]
        gt_cells = clone_to_cells[gt_clone]
        dice_clade = clade_sets[c]
        overlap = len(gt_cells & dice_clade)
        prec = overlap / len(dice_clade) if len(dice_clade) > 0 else 0
        rec = overlap / len(gt_cells) if len(gt_cells) > 0 else 0
        f1 = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0

        gt_to_dice[gt_clone] = c
        matching_rows.append({
            "gt_clone": gt_clone,
            "matched_dice_clade": f"dice_clade_{c}",
            "dice_clade_source": clone_order[c],
            "n_gt_cells": len(gt_cells),
            "n_clade_cells": len(dice_clade),
            "overlap": overlap,
            "precision": round(prec, 4),
            "recall": round(rec, 4),
            "f1": round(f1, 4),
        })

    matching_df = pd.DataFrame(matching_rows)

    # Per-cell label assignments
    gt_labels_map = {}
    for clone_name, cells in clone_to_cells.items():
        for c in cells:
            gt_labels_map[c] = clone_name

    cell_label_rows = []
    n_correct = 0
    for i, cell in enumerate(cell_names):
        gt_clone = gt_labels_map[cell]
        dice_clade_idx = tree_labels[i]
        matched_clade_idx = gt_to_dice.get(gt_clone, -1)
        is_match = dice_clade_idx == matched_clade_idx
        if is_match:
            n_correct += 1
        cell_label_rows.append({
            "cell": cell,
            "gt_clone": gt_clone,
            "dice_clade": f"dice_clade_{dice_clade_idx}",
            "dice_clade_source": clone_order[dice_clade_idx],
            "match": is_match,
        })

    cell_labels_df = pd.DataFrame(cell_label_rows)
    accuracy = n_correct / len(cell_names) if len(cell_names) > 0 else 0

    # Confusion matrix: rows = GT clone, cols = DICE clade
    dice_clade_names = [f"dice_clade_{j}" for j in range(n_clones)]
    confusion = np.zeros((n_clones, n_clones), dtype=int)
    for i, cell in enumerate(cell_names):
        gt_clone = gt_labels_map[cell]
        gt_idx = clone_order.index(gt_clone)
        dice_idx = tree_labels[i]
        confusion[gt_idx, dice_idx] += 1

    confusion_df = pd.DataFrame(
        confusion,
        index=clone_order,
        columns=dice_clade_names,
    )
    confusion_df.index.name = "gt_clone"

    return matching_df, cell_labels_df, confusion_df, accuracy


def print_clone_matching_report(matching_df, confusion_df, accuracy):
    """Print the clone matching report."""
    print("Best clone matching (Hungarian algorithm, maximizing F1):")
    print("-" * 95)
    print(f"{'GT Clone':<16} {'Matched Clade':<18} {'Source':<16} "
          f"{'GT cells':>8} {'Clade':>6} {'Overlap':>8} "
          f"{'Prec':>6} {'Rec':>6} {'F1':>6}")
    print("-" * 95)

    for _, r in matching_df.iterrows():
        print(f"{r['gt_clone']:<16} {r['matched_dice_clade']:<18} "
              f"{r['dice_clade_source']:<16} "
              f"{r['n_gt_cells']:>8} {r['n_clade_cells']:>6} "
              f"{r['overlap']:>8} {r['precision']:>6.4f} "
              f"{r['recall']:>6.4f} {r['f1']:>6.4f}")

    mean_f1 = matching_df["f1"].mean()
    print("-" * 95)
    print(f"{'Mean F1':<16} {'':<18} {'':<16} {'':<8} {'':<6} "
          f"{'':<8} {'':<6} {'':<6} {mean_f1:>6.4f}")
    print(f"\n  Matching accuracy: {accuracy:.4f} "
          f"({int(accuracy * len(matching_df) > 0 and sum(matching_df['overlap']))} / "
          f"total cells correctly assigned)")

    if accuracy >= 0.95:
        print("  (Excellent matching)")
    elif accuracy >= 0.8:
        print("  (Good matching)")
    elif accuracy >= 0.5:
        print("  (Moderate matching)")
    else:
        print("  (Poor matching)")

    print(f"\nConfusion matrix (GT clone x DICE clade):")
    print(confusion_df.to_string())
    print()


def save_clone_matching_csv(matching_df, cell_labels_df, confusion_df,
                            accuracy, output_dir):
    """Save clone matching results to CSV files."""
    # Matching summary
    match_path = os.path.join(output_dir, "eval_clone_matching.csv")
    # Append accuracy as last row
    acc_row = pd.DataFrame([{
        "gt_clone": "ACCURACY",
        "matched_dice_clade": "",
        "dice_clade_source": "",
        "n_gt_cells": "",
        "n_clade_cells": "",
        "overlap": "",
        "precision": "",
        "recall": "",
        "f1": round(accuracy, 4),
    }])
    out_df = pd.concat([matching_df, acc_row], ignore_index=True)
    out_df.to_csv(match_path, index=False)
    print(f"  Saved clone matching to {match_path}")

    # Per-cell labels
    cell_path = os.path.join(output_dir, "eval_cell_labels.csv")
    cell_labels_df.to_csv(cell_path, index=False)
    print(f"  Saved per-cell labels to {cell_path}")

    # Confusion matrix
    conf_path = os.path.join(output_dir, "eval_confusion_matrix.csv")
    confusion_df.to_csv(conf_path)
    print(f"  Saved confusion matrix to {conf_path}")


def save_f1_csv(results, mean_f1, output_path):
    """Save per-clone F1 scores and mean F1 to a CSV file."""
    rows = []
    for r in results:
        rows.append({
            "clone": r["clone"],
            "n_true_cells": r["n_true_cells"],
            "best_clade_size": r["best_clade_size"],
            "precision": round(r["precision"], 4),
            "recall": round(r["recall"], 4),
            "f1": round(r["f1"], 4),
        })
    rows.append({
        "clone": "Mean",
        "n_true_cells": "",
        "best_clade_size": "",
        "precision": "",
        "recall": "",
        "f1": round(mean_f1, 4),
    })
    df = pd.DataFrame(rows)
    df.to_csv(output_path, index=False)
    print(f"  Saved F1 scores to {output_path}")


def save_ari_csv(ari, clade_sets, clone_order, output_path):
    """Save ARI result and per-clade sizes to a CSV file."""
    rows = []
    for clone_name, clade in zip(clone_order, clade_sets):
        rows.append({
            "clone": clone_name,
            "tree_clade_size": len(clade),
        })
    rows.append({
        "clone": "ARI",
        "tree_clade_size": round(ari, 4),
    })
    df = pd.DataFrame(rows)
    df.to_csv(output_path, index=False)
    print(f"  Saved ARI results to {output_path}")


def save_nrfd_csv(nrfd, n_shared, n_inferred, n_reference, output_path):
    """Save NRFD result and bipartition counts to a CSV file."""
    rows = [
        {"metric": "bipartitions_inferred", "value": n_inferred},
        {"metric": "bipartitions_reference", "value": n_reference},
        {"metric": "bipartitions_shared", "value": n_shared},
        {"metric": "NRFD", "value": round(nrfd, 4)},
    ]
    df = pd.DataFrame(rows)
    df.to_csv(output_path, index=False)
    print(f"  Saved NRFD results to {output_path}")


def main():
    nwk_path = "output/standard_root_balME_tree.nwk"
    gt_path = "output/gt.tsv"
    output_dir = None

    if len(sys.argv) >= 3:
        nwk_path = sys.argv[1]
        gt_path = sys.argv[2]
    if len(sys.argv) >= 4:
        output_dir = sys.argv[3]

    # Default output_dir: same directory as the tree file
    if output_dir is None:
        output_dir = os.path.dirname(os.path.abspath(nwk_path))
    os.makedirs(output_dir, exist_ok=True)

    print(f"Tree file:         {nwk_path}")
    print(f"Ground truth file: {gt_path}")
    print(f"Output directory:  {output_dir}")

    tree = load_tree(nwk_path)
    gt_df, clone_to_cells = load_ground_truth(gt_path)
    clades = get_all_clades(tree)

    n_leaves = len(tree.get_terminals())
    n_clades = len(clades)
    n_clones = len(clone_to_cells)

    print(f"\nTree leaves:       {n_leaves}")
    print(f"Internal clades:   {n_clades}")
    print(f"Ground-truth clones: {n_clones}")

    # F1 evaluation
    results = compute_best_f1_per_clone(clades, clone_to_cells)
    mean_f1 = print_f1_report(results)
    save_f1_csv(results, mean_f1, os.path.join(output_dir, "eval_f1.csv"))

    # ARI: k-means vs tree clades
    ari, tree_labels, clade_sets, clone_order = compute_ari_kmeans_vs_tree(
        gt_df, clades, clone_to_cells)
    print_ari_report(ari, clade_sets, clone_order)
    save_ari_csv(ari, clade_sets, clone_order,
                 os.path.join(output_dir, "eval_ari.csv"))

    # NRFD: inferred tree vs reference tree from ground-truth clones
    nrfd, n_shared, n_inferred, n_reference = compute_nrfd(tree, clone_to_cells)
    print_nrfd_report(nrfd, n_shared, n_inferred, n_reference)
    save_nrfd_csv(nrfd, n_shared, n_inferred, n_reference,
                  os.path.join(output_dir, "eval_nrfd.csv"))

    # Clone matching: best 1-to-1 assignment between GT clones and DICE clades
    matching_df, cell_labels_df, confusion_df, accuracy = \
        find_best_clone_matching(gt_df, clades, clone_to_cells)
    print_clone_matching_report(matching_df, confusion_df, accuracy)
    save_clone_matching_csv(matching_df, cell_labels_df, confusion_df,
                            accuracy, output_dir)

    # Summary: all metrics in one row
    summary_path = os.path.join(output_dir, "eval_summary.csv")
    df = pd.DataFrame([{
        "tree_file": nwk_path,
        "gt_file": gt_path,
        "n_leaves": n_leaves,
        "n_clades": n_clades,
        "n_clones": n_clones,
        "mean_f1": round(mean_f1, 4),
        "ari": round(ari, 4),
        "nrfd": round(nrfd, 4),
        "matching_accuracy": round(accuracy, 4),
    }])
    df.to_csv(summary_path, index=False)
    print(f"  Saved summary to {summary_path}")

    return mean_f1, ari, nrfd, accuracy


if __name__ == "__main__":
    main()
