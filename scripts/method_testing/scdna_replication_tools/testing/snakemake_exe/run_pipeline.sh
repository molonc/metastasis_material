#!/bin/sh
# SA535 metastasis project replication timming

configPath=/home/htran/Projects/hakwoo_project/testing_methods/scdna_replication_tools/testing/snakemake_exe


echo "Run preprocessing pipeline"
snakemake -s $configPath/pipeline_snakemake_RT --configfile $configPath/snakemake_config_RT.yaml -j 7 2> rt.log

#echo "Downloading mouse alignment data"
# snakemake -s $configPath/download_Snakefile_mouse --configfile $configPath/snakemake_config_10x.yaml -j 8 2> download_mouse.log

echo "Completed"
