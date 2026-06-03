#!/usr/bin/env python3
"""
Generate a synthetic copy number profiles dataset with distinct clonal
populations, then run the full DICE-paper k-means workflow on it.

This creates a profiles.tsv file in DICE input format with 4 clones
(mimicking the T10 dataset structure from the paper: diploid,
hypodiploid, aneuploid A, aneuploid B), runs k-means + silhouette
selection, and computes ARI against known ground-truth labels.

Dependencies:
    pip install numpy pandas scikit-learn
"""

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, adjusted_rand_score
import os

# ================================================================
# Constants
# ================================================================
N_CELLS = 60
N_CHROMS = 3
BINS_PER_CHROM = 50
BIN_SIZE = 1_000_000  # 1 Mbp bins
N_BINS = N_CHROMS * BINS_PER_CHROM

CLONE_SIZES = {"diploid": 15, "hypodiploid": 15, "aneuploid_A": 15, "aneuploid_B": 15}
CLONE_NAMES = list(CLONE_SIZES.keys())


def define_clone_profiles():
    """Define base copy number profiles for each clone."""
    base_profiles = {}

    diploid = np.full(N_BINS, 2, dtype=int)
    base_profiles["diploid"] = diploid.copy()

    hypo = diploid.copy()
    hypo[0:30] = 1       # chr1 bins 0-29: loss
    hypo[50:70] = 1      # chr2 bins 0-19: loss
    hypo[110:150] = 1    # chr3 bins 10-49: loss
    base_profiles["hypodiploid"] = hypo

    aneuA = diploid.copy()
    aneuA[10:40] = 4     # chr1 bins 10-39: gain
    aneuA[50:75] = 3     # chr2 bins 0-24: gain
    base_profiles["aneuploid_A"] = aneuA

    aneuB = diploid.copy()
    aneuB[75:100] = 3    # chr2 bins 25-49: gain
    aneuB[100:130] = 5   # chr3 bins 0-29: gain
    base_profiles["aneuploid_B"] = aneuB

    return base_profiles


def generate_synthetic_dataset(base_profiles, out_tsv="profiles.tsv"):
    """Generate cells with noise around each clone's base profile and write to TSV."""
    print("Generating synthetic CNP dataset ...")

    rows = []
    cell_id = 0
    true_labels = {}

    for clone_name, n_cells in CLONE_SIZES.items():
        base = base_profiles[clone_name]
        for _ in range(n_cells):
            cell_id += 1
            cell_name = f"cell{cell_id}"
            true_labels[cell_name] = clone_name

            profile = base.copy()
            noise_mask = np.random.rand(N_BINS) < 0.05
            noise_vals = np.random.choice([-1, 1], size=N_BINS)
            profile = profile + noise_mask * noise_vals
            profile = np.clip(profile, 0, 8)

            for chrom_idx in range(N_CHROMS):
                chrom_name = f"chr{chrom_idx + 1}"
                for bin_idx in range(BINS_PER_CHROM):
                    global_idx = chrom_idx * BINS_PER_CHROM + bin_idx
                    total_cn = int(profile[global_idx])
                    cn_a = total_cn // 2
                    cn_b = total_cn - cn_a
                    start = bin_idx * BIN_SIZE
                    end = start + BIN_SIZE
                    rows.append([cell_name, chrom_name, start, end, f"{cn_a},{cn_b}"])

    df = pd.DataFrame(rows, columns=["CELL", "chrom", "start", "end", "CN states"])
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"  Wrote {out_tsv}: {N_CELLS} cells × {N_BINS} bins ({len(rows)} rows)")
    print(f"  Clones: {dict(CLONE_SIZES)}\n")

    return true_labels, df


def generate_total_cn_dataset(allele_df, out_tsv="profiles_totalcn.tsv"):
    """Convert allele-specific CN profiles to total CN and write to TSV.

    Changes CN states from 'a,b' format to a single total CN value.
    """
    print("Generating total CN dataset ...")
    df = allele_df.copy()
    split = df["CN states"].str.split(",", expand=True).astype(int)
    df["CN states"] = split[0] + split[1]
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"  Wrote {out_tsv}\n")
    return out_tsv


def read_profiles_matrix(tsv_path):
    """Read allele-specific CN profiles (a,b format) into a cell x bin matrix."""
    print("Reading CNPs into matrix ...")
    df = pd.read_csv(tsv_path, sep="\t")
    df = df.sort_values(["CELL", "chrom", "start"])
    split = df["CN states"].str.split(",", expand=True).astype(int)
    df["cn_total"] = split[0] + split[1]
    pivot = df.pivot_table(index="CELL", columns=["chrom", "start"], values="cn_total", aggfunc="first")
    cell_names = list(pivot.index)
    matrix = pivot.values.astype(float)
    print(f"  Matrix shape: {matrix.shape}\n")
    return cell_names, matrix


def read_total_cn_matrix(tsv_path, sep="\t"):
    """Read total CN profiles (single integer format) into a cell x bin matrix."""
    print(f"Reading total CN profiles from {tsv_path} ...")
    df = pd.read_csv(tsv_path, sep=sep, low_memory=False)
    df = df.sort_values(["CELL", "chrom", "start"])
    df["CN states"] = pd.to_numeric(df["CN states"], errors="coerce")
    pivot = df.pivot_table(index="CELL", columns=["chrom", "start"],
                           values="CN states", aggfunc="first")
    cell_names = list(pivot.index)
    matrix = pivot.values.astype(float)
    nan_count = np.isnan(matrix).sum()
    if nan_count > 0:
        print(f"  Warning: {nan_count} NaN values found, filling with 0")
        matrix = np.nan_to_num(matrix, nan=0.0)
    print(f"  Cells: {len(cell_names)}, Bins: {matrix.shape[1]}")
    print(f"  Matrix shape: {matrix.shape}\n")
    return cell_names, matrix


def run_kmeans_silhouette(matrix, k_range=range(2, 9)):
    """Run k-means for a range of k and select the best by silhouette score."""
    print(f"Running k-means for k = {k_range.start}..{k_range.stop - 1} with Silhouette selection ...\n")
    best_score = -1
    best_k = None
    best_labels = None

    for k in k_range:
        km = KMeans(n_clusters=k, n_init=10, random_state=42)
        labels = km.fit_predict(matrix)
        score = silhouette_score(matrix, labels)
        marker = ""
        if score > best_score:
            best_score = score
            best_k = k
            best_labels = labels
            marker = "  <-- best"
        print(f"  k={k}  silhouette={score:.4f}{marker}")

    print(f"\n  Optimal k = {best_k}\n")

    unique, counts = np.unique(best_labels, return_counts=True)
    for cl, cnt in zip(unique, counts):
        print(f"  Cluster {cl}: {cnt} cells")

    return best_k, best_labels


def compute_ari(cell_names, best_labels, true_labels):
    """Compute Adjusted Rand Index against ground-truth clone labels."""
    print("\nComputing ARI against ground-truth labels ...")
    true_label_vec = [CLONE_NAMES.index(true_labels[c]) for c in cell_names]
    ari = adjusted_rand_score(true_label_vec, best_labels)
    print(f"  ARI = {ari:.4f}")

    if ari == 1.0:
        print("  (Perfect agreement — k-means recovered all ground-truth clones)")
    elif ari > 0.8:
        print("  (Strong agreement)")
    elif ari > 0.5:
        print("  (Moderate agreement)")
    else:
        print("  (Weak agreement)")

    return ari


def print_contingency_table(cell_names, best_labels, true_labels):
    """Print cluster-vs-clone contingency table."""
    print("\n\nCluster-vs-Clone contingency table:")
    print("-" * 55)
    header = f"{'Cluster':>10}" + "".join(f"{c:>14}" for c in CLONE_NAMES)
    print(header)
    print("-" * 55)

    unique = np.unique(best_labels)
    for cl in unique:
        row = f"  {cl:>6}  "
        for clone_name in CLONE_NAMES:
            mask = [
                1
                for i, c in enumerate(cell_names)
                if best_labels[i] == cl and true_labels[c] == clone_name
            ]
            row += f"{sum(mask):>14}"
        print(row)
    print("-" * 55)


def save_ground_truth(cell_names, best_labels, true_labels, out_dir="output"):
    """Save true clone labels and optimal k-means cluster assignments to gt.tsv."""
    os.makedirs(out_dir, exist_ok=True)
    gt_path = os.path.join(out_dir, "gt.tsv")

    df = pd.DataFrame({
        "cell": cell_names,
        "true_clone": [true_labels[c] for c in cell_names],
        "kmeans_cluster": best_labels,
    })
    df.to_csv(gt_path, sep="\t", index=False)
    print(f"\n  Saved ground truth and k-means labels to {gt_path}")


def main():
    np.random.seed(42)
    out_tsv = "profiles.tsv"

    base_profiles = define_clone_profiles()
    true_labels, allele_df = generate_synthetic_dataset(base_profiles, out_tsv)
    total_cn_tsv = generate_total_cn_dataset(allele_df, "profiles_totalcn.tsv")
    cell_names, matrix = read_profiles_matrix(out_tsv)
    best_k, best_labels = run_kmeans_silhouette(matrix)
    compute_ari(cell_names, best_labels, true_labels)
    print_contingency_table(cell_names, best_labels, true_labels)
    save_ground_truth(cell_names, best_labels, true_labels)

    print(f"\n  Output files:")
    print(f"    Allele-specific: {os.path.abspath(out_tsv)}")
    print(f"    Total CN:        {os.path.abspath(total_cn_tsv)}")
    print("  Run DICE:")
    print(f"    dice -i {out_tsv} -o output/ -m balME")
    print(f"    dice -i {total_cn_tsv} -o output/ -m balME -t")


if __name__ == "__main__":
    main()
