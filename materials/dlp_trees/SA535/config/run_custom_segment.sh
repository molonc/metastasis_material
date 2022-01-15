#!/bin/sh
# SA535 metastasis
# Tyler corrupt tree encoding version
# echo "Download data"

configPath=/home/htran/storage/datasets/metastasis_results/SA535_new_encoding/SA535_wholedata_v2/config

log_file="${configPath}/custom_segment_log"
exec >> $log_file 2>&1 && tail $log_file


# echo "Downloading data"
# snakemake -s $configPath/download_Snakefile --configfile $configPath/snakemake_config.yaml -j 8 2> download.log
# snakemake -s $configPath/download_Snakefile --configfile $configPath/snakemake_config_added_library.yaml -j 8 2> download.log

echo "Run custom segment data"
snakemake -s $configPath/download_Snakefile_cloneAlign --configfile $configPath/snakemake_config_clonealign.yaml -j 8 2> clonealign.log

# echo "Run CloneAlign process segment data"
# snakemake -s $configPath/snake_clonealign --configfile $configPath/snakemake_config.yaml -j 8 2> clonealign.log


echo "Complete"
