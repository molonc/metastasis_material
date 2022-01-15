#!/bin/sh
# SA919
# Tyler corrupt tree encoding version
# echo "Download data"

configPath=/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_wholedata/config

log_file="${configPath}/corrupt_log"
exec >> $log_file 2>&1 && tail $log_file

echo "------------------------\n"
# echo "Download data"
# snakemake -s $configPath/download_Snakefile --configfile $configPath/snakemake_config.yaml -j 8 2> download.log




# echo "Run corrupt tree get features"
# snakemake -s $configPath/corrupt_get_features --configfile $configPath/snakemake_config.yaml -j 8 2> corrupt.log

# echo "Run corrupt tree inference"
# snakemake -s $configPath/corrupt_inference --configfile $configPath/snakemake_config.yaml -j 8 2> corrupt.log

# echo "Run download segment data for clonealign"
# snakemake -s $configPath/download_Snakefile_cloneAlign_v2 --configfile $configPath/snakemake_config.yaml -j 8 2> clonealign.log


# echo "Run corrupt grow"
# snakemake -s $configPath/corrupt_Snakefile_nonecells --configfile $configPath/snakemake_config.yaml -j 8 2> corrupt.log
# echo "Run corrupt grow"
# snakemake -s $configPath/corrupt_grow_snakefile --configfile $configPath/snakemake_config.yaml -j 5 2> corrupt.log

# echo "Run corrupt grow"
# snakemake -s $configPath/corrupt_grow_none_cells_snakefile --configfile $configPath/snakemake_config.yaml -j 5 2> corrupt.log



echo "Run clonealign"
snakemake -s $configPath/snake_clonealign --configfile $configPath/snakemake_config.yaml -j 8 2> clonealign.log

# snakemake -s $configPath/corrupt_Snakefile_viz --configfile $configPath/snakemake_config.yaml -j 8 2> corrupt.log


# echo "Generate trajectory"
# snakemake -s $configPath/generate_trajectory --configfile $configPath/snakemake_config.yaml -j 8 2> corrupt.log



# echo "Generate schnapp"
# snakemake -s $configPath/download_Snakefile_schnapp --configfile $configPath/snakemake_config.yaml -j 4 2> corrupt.log


echo "Complete"
