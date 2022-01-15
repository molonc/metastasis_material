#!/bin/sh
# SA535 metastasis
# Tyler corrupt tree encoding version
# echo "Download data"

configPath=/home/htran/storage/datasets/metastasis_results/SA535_new_encoding/SA535_wholedata_v2/config

log_file="${configPath}/corrupt_log"
exec >> $log_file 2>&1 && tail $log_file


# echo "Downloading data"
# snakemake -s $configPath/download_Snakefile --configfile $configPath/snakemake_config.yaml -j 8 2> download.log
# snakemake -s $configPath/download_Snakefile --configfile $configPath/snakemake_config_added_library.yaml -j 8 2> download.log





# echo "Run corrupt tree get features"
# snakemake -s $configPath/corrupt_get_features --configfile $configPath/snakemake_config.yaml -j 8 2> corrupt.log

# echo "Run corrupt tree inference"
# snakemake -s $configPath/corrupt_inference --configfile $configPath/snakemake_config.yaml -j 8 2> corrupt.log


# echo "Run corrupt grow"
# snakemake -s $configPath/corrupt_grow_snakefile --configfile $configPath/snakemake_config.yaml -j 4 2> corrupt.log


# echo "Generate different heatmap plots"
# snakemake -s $configPath/corrupt_plot_heatmap --configfile $configPath/snakemake_config.yaml -j 8 2> corrupt.log


# echo "Generate different heatmap plots"
# snakemake -s $configPath/corrupt_plot_heatmap_leiden --configfile $configPath/snakemake_config.yaml -j 8 2> corrupt.log

# echo "Generate heatmap, trajectory, and prevalences"
# snakemake -s $configPath/generate_trajectory --configfile $configPath/snakemake_config.yaml -j 8 2> corrupt.log

echo "Run CloneAlign download segment data"
snakemake -s $configPath/download_Snakefile_cloneAlign --configfile $configPath/snakemake_config.yaml -j 6 2> clonealign.log

# echo "Run CloneAlign process segment data"
# snakemake -s $configPath/snake_clonealign --configfile $configPath/snakemake_config.yaml -j 8 2> clonealign.log


# echo "Download hmmcopy cnbin data, ascn analysis"
# snakemake -s $configPath/download_Snakefile_schnapp --configfile $configPath/snakemake_config.yaml -j 7 2> ascn.log

echo "Complete"
