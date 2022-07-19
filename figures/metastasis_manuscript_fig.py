import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


if __name__ == '__main__':
    os.makedirs("figures/output/", exist_ok=True)
    """
    1. Wrangle CSVs: melt to wanted document
    Columns = ['week', 'tumor size', 'id', 'rep', 'site']
    """
    met = pd.read_csv("figures/Fig2_metastatic.csv")
    met = pd.melt(met, id_vars=['weeks'], value_vars=met.columns[1:])
    met = met.assign(id=[i.split("_")[0] for i in met["variable"]])
    met = met.assign(rep=["R"+i.split("_")[1] for i in met["variable"]])
    met = met.assign(site=["metastasis"]*len(met))
    met.columns = ['Weeks', 'ID_rep', 'Tumour Volume', 'ID', 'Replicates', 'Site']
    met = met.dropna()

    nonmet = pd.read_csv("figures/Fig2_nonmetastatic.csv")
    nonmet = pd.melt(nonmet, id_vars=['weeks'], value_vars=nonmet.columns[1:])
    nonmet = nonmet.assign(id=[i.split("_")[0] for i in nonmet["variable"]])
    nonmet = nonmet.assign(rep=["R"+i.split("_")[1] for i in nonmet["variable"]])
    nonmet = nonmet.assign(site=["primary"]*len(nonmet))
    nonmet.columns = ['Weeks', 'ID_rep', 'Tumour Volume', 'ID', 'Replicates', 'Site']
    nonmet = nonmet.dropna()

    # max_tumour_volume = max(met["Tumour Volume"].max(), nonmet["Tumour Volume"].max())
    max_tumour_volume = 1600

    """
    2. Plot Fig 1 Met and nonMet
    """

    fig, ax = plt.subplots(1, 1, figsize=(4.3, 2.5))
    sns.lineplot(data=met, x="Weeks", y="Tumour Volume", hue="ID", style="Replicates",
                 markers=True, dashes=False, linewidth=1, markersize=5, ax=ax)
    plt.axvline(x=14, color='red', linestyle='--', linewidth=1)
    sns.despine()
    plt.legend(bbox_to_anchor=(1, 1), frameon=False, fontsize=7)
    plt.xlabel("Weeks", fontsize=8)
    plt.ylabel("Tumour Volume $mm^3$", fontsize=8)
    plt.xticks(fontsize=7.5)
    plt.yticks(fontsize=7.5)
    plt.xlim(0, 30)
    plt.ylim(0, max_tumour_volume)
    box = ax.get_position()
    ax.set_position([box.x0 + box.height * 0.02, box.y0 + box.height * 0.05, box.width * 0.85, box.height * 0.95])
    plt.title("Tumours with metastasis", weight='bold', fontsize=10)
    plt.savefig("figures/output/fig2_met.png", dpi=500)
    plt.clf()

    fig, ax = plt.subplots(1, 1, figsize=(4.3, 2.5))
    sns.lineplot(data=nonmet, x="Weeks", y="Tumour Volume", hue="ID", style="Replicates",
                 markers=True, dashes=False, linewidth=1, markersize=5, ax=ax)
    plt.axvline(x=14, color='red', linestyle='--', linewidth=1)
    sns.despine()
    plt.legend(bbox_to_anchor=(1.325, 1), frameon=False, fontsize=7)
    plt.xlabel("Weeks", fontsize=8)
    plt.ylabel("Tumour Volume $mm^3$", fontsize=8)
    plt.xticks(fontsize=7.5)
    plt.yticks(fontsize=7.5)
    plt.xlim(0, 30)
    plt.ylim(0, max_tumour_volume)
    box = ax.get_position()
    ax.set_position([box.x0 + box.height * 0.02, box.y0 + box.height * 0.05, box.width * 0.85, box.height * 0.95])
    plt.title("Tumours without metastasis", weight='bold', fontsize=10)
    plt.savefig("figures/output/fig2_nonmet.png", dpi=500)
    plt.clf()

    """
    3. Plot Fig 1, Fig 8 Markers
    """
    fig1_markers = pd.read_csv("figures/Fig2_markers.csv", index_col=0)

    fig, ax = plt.subplots(1, 1, figsize=(2.5, 2.5))
    cmap = ["#5a0064d9", "#fb8811d9"]
    sns.heatmap(fig1_markers, cmap=cmap, cbar=False, linewidth=0.75)
    legend_handles = [Patch(color=cmap[True], label='Positive'),
                      Patch(color=cmap[False], label='Negative')]
    plt.legend(handles=legend_handles, ncol=2, bbox_to_anchor=[0.5, 1.0], loc='lower center',
               frameon=False, fontsize=7.5)
    plt.xticks(rotation=90, fontsize=7.5)
    plt.yticks(rotation=0, fontsize=7.5)
    box = ax.get_position()
    ax.set_position([box.x0 + box.width * 0.2, box.y0 + box.height * 0.2, box.width * 0.85, box.height * 0.85])
    plt.savefig("figures/output/fig2_markers.png", dpi=500)
    plt.clf()

    fig8_markers = pd.read_csv("figures/Fig8_markers.csv", index_col=0)

    fig, ax = plt.subplots(1, 1, figsize=(2.5, 1.5))
    cmap = ["#5a0064d9", "#fb8811d9"]
    sns.heatmap(fig8_markers, cmap=cmap, cbar=False, linewidth=0.75)
    legend_handles = [Patch(color=cmap[True], label='Positive'),
                      Patch(color=cmap[False], label='Negative')]
    plt.legend(handles=legend_handles, ncol=2, bbox_to_anchor=[0.5, 1.0], loc='lower center',
               frameon=False, fontsize=7.5)
    plt.xticks(rotation=90, fontsize=7.5)
    plt.yticks(rotation=0, fontsize=7.5)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.3, box.width, box.height * 0.65])
    plt.savefig("figures/output/fig8_markers.png", dpi=500)
    plt.clf()

    """
    4. Plot EGFR
    """
    fig8_919 = pd.read_csv("figures/Fig8_919_EGFR.csv")
    fig, ax = plt.subplots(1, 1, figsize=(3.8, 2.5))
    sns.boxplot(x="group", y="H-score", data=fig8_919, linewidth=1.5, ax=ax)
    for i, box in enumerate(ax.artists):
        box.set_edgecolor('black')
        box.set_facecolor('white')
    sns.swarmplot(x="group", y="H-score", data=fig8_919, linewidth=1, ax=ax, size=7)
    sns.despine()
    plt.ylim(-10, 200)
    plt.xticks(fontsize=7.5)
    plt.yticks(fontsize=7.5)
    ax.set_xticklabels(["Passage 3 & 4", "Passage 7"])
    plt.xlabel("", fontsize=6)
    plt.ylabel("H-score", fontsize=7.5)
    plt.title("SA919 EGFR H-Score", weight='bold', fontsize=8)
    plt.savefig("figures/output/fig8_919_egfr.png", dpi=500)
    plt.clf()

    fig8_535 = pd.read_csv("figures/Fig8_535_EGFR.csv")
    fig8_535 = fig8_535[fig8_535["group"] == "X4"]
    fig, ax = plt.subplots(1, 1, figsize=(3.8, 2.5))
    sns.boxplot(x="group", y="H-score", data=fig8_535, linewidth=1.5, ax=ax)
    for i, box in enumerate(ax.artists):
        box.set_edgecolor('black')
        box.set_facecolor('white')
    sns.swarmplot(x="group", y="H-score", hue="site", data=fig8_535, linewidth=1, ax=ax, size=7)
    sns.despine()
    plt.ylim(0, 200)
    plt.xticks(fontsize=7.5)
    plt.yticks(fontsize=7.5)
    ax.set_xticklabels(["Passage 4"])
    plt.xlabel("", fontsize=7.5)
    plt.ylabel("H-score", fontsize=7.5)
    plt.legend(bbox_to_anchor=(0.95, 1), frameon=False, fontsize=7.5)
    box = ax.get_position()
    ax.set_position([box.x0 + box.height * 0.02, box.y0 + box.height * 0.05, box.width * 0.85, box.height * 0.95])
    plt.title("SA535 EGFR H-Score", weight='bold', fontsize=8)
    plt.savefig("figures/output/fig8_535_egfr_hscore.png", dpi=500)
    plt.clf()

    fig, ax = plt.subplots(1, 1, figsize=(3.8, 2.5))
    sns.boxplot(x="group", y="Intensity", data=fig8_535, linewidth=1.5, ax=ax)
    for i, box in enumerate(ax.artists):
        box.set_edgecolor('black')
        box.set_facecolor('white')
    sns.swarmplot(x="group", y="Intensity", hue="site", data=fig8_535, linewidth=1, ax=ax, size=7)
    sns.despine()
    plt.ylim(0, 5)
    plt.xticks(fontsize=7.5)
    plt.yticks(fontsize=7.5)
    ax.set_xticklabels(["Passage 4"])
    plt.xlabel("", fontsize=7.5)
    plt.ylabel("Intensity", fontsize=7.5)
    plt.legend(bbox_to_anchor=(0.95, 1), frameon=False, fontsize=7.5)
    box = ax.get_position()
    ax.set_position([box.x0 + box.height * 0.02, box.y0 + box.height * 0.05, box.width * 0.85, box.height * 0.95])
    plt.title("SA535 EGFR Intensity", weight='bold', fontsize=8)
    plt.savefig("figures/output/fig8_535_egfr_intensity.png", dpi=500)
    plt.clf()

    fig, ax = plt.subplots(1, 1, figsize=(3.8, 2.5))
    sns.boxplot(x="group", y="percentage", data=fig8_535, linewidth=1.5, ax=ax)
    for i, box in enumerate(ax.artists):
        box.set_edgecolor('black')
        box.set_facecolor('white')
    sns.swarmplot(x="group", y="percentage", hue="site", data=fig8_535, linewidth=1, ax=ax, size=7)
    sns.despine()
    plt.ylim(0, 100)
    plt.xticks(fontsize=7.5)
    plt.yticks(fontsize=7.5)
    ax.set_xticklabels(["Passage 4"])
    plt.xlabel("", fontsize=7.5)
    plt.ylabel("Percentage", fontsize=7.5)
    plt.legend(bbox_to_anchor=(0.95, 1), frameon=False, fontsize=7.5)
    box = ax.get_position()
    ax.set_position([box.x0 + box.height * 0.02, box.y0 + box.height * 0.05, box.width * 0.85, box.height * 0.95])
    plt.title("SA535 EGFR Percentage", weight='bold', fontsize=8)
    plt.savefig("figures/output/fig8_535_egfr_percentage.png", dpi=500)
    plt.clf()

    """
    5. Plot Ki67
    """
    fig2_ki67 = pd.read_csv("figures/Fig2_ki67.csv")
    fig, ax = plt.subplots(1, 1, figsize=(3.8, 2.5))
    sns.boxplot(x="SA ID", y="Ki67 %", data=fig2_ki67, linewidth=1.5, ax=ax)
    for i, box in enumerate(ax.artists):
        box.set_edgecolor('black')
        box.set_facecolor('white')
    sns.despine()
    plt.ylim(0, 100)
    plt.xticks(rotation=90, fontsize=7.5)
    plt.yticks(fontsize=7.5)
    plt.xlabel("PDX", fontsize=7.5)
    plt.ylabel("Ki67 Percentage", fontsize=7.5)
    box = ax.get_position()
    ax.set_position([box.x0 + box.height * 0.02, box.y0 + box.height * 0.3, box.width, box.height * 0.8])
    plt.savefig("figures/output/fig2_ki67.png", dpi=500)
    plt.clf()

    """
    Old Code
    """
    # fig8_535 = pd.read_csv("figures/Fig8_535_EGFR.csv")
    # fig, ax = plt.subplots(1, 1, figsize=(3.8, 2.5))
    # sns.boxplot(x="group", y="H-score", data=fig8_535, linewidth=1.5, ax=ax)
    # for i, box in enumerate(ax.artists):
    #     box.set_edgecolor('black')
    #     box.set_facecolor('white')
    # sns.swarmplot(x="group", y="H-score", data=fig8_535, linewidth=1, ax=ax, size=5)
    # sns.despine()
    # plt.ylim(0, 200)
    # plt.xticks(fontsize=6)
    # plt.yticks(fontsize=6)
    # ax.set_xticklabels(["Passage 4", "Passage 6 & 7"])
    # plt.xlabel("", fontsize=6)
    # plt.ylabel("H-score", fontsize=6)
    # plt.title("EGFR SA535", weight='bold', fontsize=8)
    # plt.savefig("figures/output/fig8_535_egfr.png", dpi=500)
    # plt.clf()

    # fig8_535 = pd.read_csv("figures/Fig8_535_EGFR.csv")
    # fig, ax = plt.subplots(1, 1, figsize=(3.8, 2.5))
    # for i in range(len(fig8_535["site"])):
    #     if fig8_535["site"][i] != "lung":
    #         fig8_535["site"][i] = "other"
    # fig8_535 = fig8_535.sort_values(by=["site"])
    #
    # sns.boxplot(x="site", y="H-score", data=fig8_535, linewidth=1.5, ax=ax)
    # for i, box in enumerate(ax.artists):
    #     box.set_edgecolor('black')
    #     box.set_facecolor('white')
    # sns.swarmplot(x="site", y="H-score", data=fig8_535, linewidth=1, ax=ax, size=5)
    # sns.despine()
    # plt.ylim(0, 200)
    # plt.xticks(fontsize=6)
    # plt.yticks(fontsize=6)
    # ax.set_xticklabels(["Lung", "Other"])
    # plt.xlabel("", fontsize=6)
    # plt.ylabel("H-score", fontsize=6)
    # plt.title("EGFR SA535", weight='bold', fontsize=8)
    # plt.savefig("figures/output/fig8_535_egfr_lung.png", dpi=500)
    # plt.clf()
    #
    # fig8_535 = pd.read_csv("figures/Fig8_535_EGFR.csv")
    # fig, ax = plt.subplots(1, 1, figsize=(3.8, 2.5))
    # for i in range(len(fig8_535["site"])):
    #     if fig8_535["site"][i] != "lung":
    #         fig8_535["site"][i] = "other"
    # fig8_535 = fig8_535[fig8_535["site"] == "other"]
    # sns.boxplot(x="group", y="H-score", data=fig8_535, linewidth=1.5, ax=ax)
    # for i, box in enumerate(ax.artists):
    #     box.set_edgecolor('black')
    #     box.set_facecolor('white')
    # sns.swarmplot(x="group", y="H-score", data=fig8_535, linewidth=1, ax=ax, size=5)
    # sns.despine()
    # plt.ylim(0, 200)
    # plt.xticks(fontsize=6)
    # plt.yticks(fontsize=6)
    # ax.set_xticklabels(["Passage 4", "Passage 6 & 7"])
    # plt.xlabel("", fontsize=6)
    # plt.ylabel("H-score", fontsize=6)
    # plt.title("EGFR SA535 without lung metastasis", weight='bold', fontsize=8)
    # plt.savefig("figures/output/fig8_535_egfr_nolung.png", dpi=500)
    # plt.clf()
    #