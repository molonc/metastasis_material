#!/usr/bin/env python3
"""
Graph clustering with alternative methods (non-Leiden).

Compares Louvain, Label Propagation, Walktrap, Infomap, and FastGreedy
on the same synthetic graph.

Requirements:
    pip install igraph matplotlib numpy
"""

import igraph as ig
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


def build_graph(block_sizes, p_within=0.25, p_between=0.01):
    """Build a synthetic graph with planted community structure (SBM)."""
    n_blocks = len(block_sizes)
    pref_matrix = [[p_within if i == j else p_between
                     for j in range(n_blocks)]
                    for i in range(n_blocks)]
    G = ig.Graph.SBM(sum(block_sizes), pref_matrix, block_sizes)
    print(f"Graph: {G.vcount()} nodes, {G.ecount()} edges\n")
    return G


def run_clustering_methods(G):
    """Run multiple clustering methods and return results."""
    methods = {}

    # Louvain (multilevel modularity optimization)
    louvain = G.community_multilevel()
    methods["Louvain"] = louvain.membership

    # Label Propagation
    lpa = G.community_label_propagation()
    methods["Label Propagation"] = lpa.membership

    # Walktrap (short random walks)
    walktrap = G.community_walktrap().as_clustering()
    methods["Walktrap"] = walktrap.membership

    # Infomap (information-theoretic)
    infomap = G.community_infomap()
    methods["Infomap"] = infomap.membership

    # FastGreedy (greedy modularity optimization)
    fastgreedy = G.community_fastgreedy().as_clustering()
    methods["FastGreedy"] = fastgreedy.membership

    # Edge Betweenness (Girvan-Newman)
    eb = G.community_edge_betweenness().as_clustering()
    methods["Edge Betweenness"] = eb.membership

    for name, membership in methods.items():
        n_clusters = len(set(membership))
        modularity = G.modularity(membership)
        print(f"{name:<20s}  clusters={n_clusters:>3}  modularity={modularity:.4f}")

    return methods


def plot_clustering_methods(G, methods, outfile="other_clustering_methods"):
    """Plot network layouts colored by each clustering method."""
    names = list(methods.keys())
    n_methods = len(names)
    ncols = 3
    nrows = (n_methods + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 5 * nrows))
    axes = axes.flatten()
    layout = G.layout_fruchterman_reingold(seed=42)
    palette = plt.cm.tab20
    n_nodes = G.vcount()

    xs = [layout[i][0] for i in range(n_nodes)]
    ys = [layout[i][1] for i in range(n_nodes)]

    for idx, name in enumerate(names):
        ax = axes[idx]
        membership = methods[name]
        colors = [palette(c % 20) for c in membership]

        for edge in G.get_edgelist():
            x0, y0 = layout[edge[0]]
            x1, y1 = layout[edge[1]]
            ax.plot([x0, x1], [y0, y1], color="lightgrey", lw=0.3, zorder=0)

        ax.scatter(xs, ys, c=colors, s=20, zorder=1, edgecolors="k", linewidths=0.3)
        n_clusters = len(set(membership))
        modularity = G.modularity(membership)
        ax.set_title(f"{name}\n{n_clusters} clusters, Q={modularity:.3f}")
        ax.set_xticks([])
        ax.set_yticks([])

    # Hide unused axes
    for idx in range(n_methods, len(axes)):
        axes[idx].set_visible(False)

    fig.suptitle("Graph Clustering Methods Comparison", fontsize=14, y=1.0)
    plt.tight_layout()
    plt.savefig(f"{outfile}.png", dpi=150, bbox_inches="tight")
    plt.savefig(f"{outfile}.pdf", bbox_inches="tight")
    print(f"\nSaved: {outfile}.png / .pdf")


def main():
    np.random.seed(42)

    # 1. Build graph
    G = build_graph(block_sizes=[40, 35, 30, 25])

    # 2. Run all clustering methods
    methods = run_clustering_methods(G)

    # 3. Plot comparison
    plot_clustering_methods(G, methods)


if __name__ == "__main__":
    main()
