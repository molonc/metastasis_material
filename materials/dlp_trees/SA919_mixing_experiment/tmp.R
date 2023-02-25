
library(dplyr)
input_dir <- '/home/htran/storage/raw_DLP/metastasis_DLP/'
lid <- 'A98166A'
lid <- 'A98166B'
metasamples <- data.table::fread(paste0(input_dir, lid, '/annotation/',lid,'_metrics.csv.gz'))
dim(metasamples)
unique(metasamples$sample_id)
metasamples <- metasamples %>%
  dplyr::filter(sample_id=='AT11367')
metasamples <- metasamples %>%
  dplyr::filter(sample_id=='AT11391')


dim(metasamples)
rm(metasamples)

# Filtering data first
