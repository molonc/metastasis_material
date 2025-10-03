
library(dplyr)
# Uploading data to EGA

# DLP+ in 1 file
# SA919 main experiment 
# SA919 mixing experiment 
# SA535

input_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/dlp_trees/'

sa919_df <- data.table::fread(paste0(input_dir, 'SA919/library_groupings.csv.gz'))
sa919_mixing_df <- data.table::fread(paste0(input_dir, 'SA919_mixing_experiment/Metastasis_Hakwoo_mixing_exp_SA919_results.csv'))
sa535_df <- data.table::fread(paste0(input_dir, 'SA535/library_groupings.csv.gz'))


# library id, sample id, origin, main site, pdxid 
# Bulk RNA-seq
# Small description

head(sa919_df)
cols_use <- c('grouping', 'sample_id', 'mainsite', 'origin', 'pdxid')

sa919_df <- sa919_df %>%
  select(all_of(cols_use)) %>%
  rename(library_id=grouping)
sa919_df$pdxid <- gsub('X0847','X0847-',sa919_df$pdxid)
sa919_df$note <- ''
sa919_df$experiment <- 'main_experiment'
sa919_df$SA_ID <- 'SA919'
dim(sa919_df)


sa535_df <- sa535_df %>%
  select(all_of(cols_use)) %>%
  rename(library_id=grouping)

sa535_df$pdxid <- gsub('X0011_','X0011-',sa535_df$pdxid)
length(unique(sa535_df$library_id))
dup_lib <- sa535_df$library_id[duplicated(sa535_df$library_id)]

sa535_df$note <- ''
sa535_df1 <- sa535_df %>%
  filter(library_id != dup_lib)
dim(sa535_df1)
sa535_df2 <- sa535_df %>%
  filter(library_id == dup_lib)

sa535_df3 <- sa535_df2[2,]
sa535_df4 <- sa535_df2[1,]

sa535_df3$note <- paste0(sa535_df4$sample_id,'_',sa535_df4$origin)

sa535_df <- dplyr::bind_rows(sa535_df1, sa535_df3)
sa535_df$experiment <- 'main_experiment'
sa535_df$SA_ID <- 'SA535'
dim(sa535_df)
View(sa535_df)

View(head(sa919_mixing_df))
selected_cols <- c('library_id','sample_id','mainsite','origin','Tumor_ID')

sa919_mixing_df <- sa919_mixing_df %>%
  select(all_of(selected_cols)) %>%
  rename(pdxid=Tumor_ID) %>%
  mutate(note='', experiment='mixing_experiment')
sa919_mixing_df$SA_ID <- 'SA919'

total_dlp <- dplyr::bind_rows(sa919_df, sa919_mixing_df, sa535_df)
View(total_dlp)
data.table::fwrite(total_dlp, paste0(input_dir, 'metastasis_project_dlp_libraries_metadata.csv'))


sa919_bulk_df <- data.table::fread('/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/metadata/metadata_Hakwoo_bulkRNA_mixing_main_exp.csv')
colnames(sa919_bulk_df)
cols_use <- c('bulk_sid', 'sample_id', 'mainsite', 'origin', 'pdxid','Tumor_ID')
Patient_ID_CID <- 'VBA0847'
sa919_bulk_df <- sa919_bulk_df %>%
  select(all_of(cols_use)) 

sa919_bulk_df1 <- sa919_bulk_df %>%
  filter(pdxid!='SA919')
sa919_bulk_df1$Tumor_ID <- NULL
dim(sa919_bulk_df1)
sa919_bulk_df1$pdxid
sa919_bulk_df1$pdxid <- gsub('X0847','X0847-',sa919_bulk_df1$pdxid)
sa919_bulk_df1$experiment <- 'main_experiment'

sa919_bulk_df2 <- sa919_bulk_df %>%
  filter(pdxid=='SA919')
sa919_bulk_df2$pdxid <- sa919_bulk_df2$Tumor_ID
dim(sa919_bulk_df2)
View(sa919_bulk_df2)
sa919_bulk_df2$Tumor_ID <- NULL
sa919_bulk_df2$experiment <- 'mixing_experiment'

total_bulk <- dplyr::bind_rows(sa919_bulk_df1, sa919_bulk_df2)
View(total_bulk)
total_bulk$SA_ID <- 'SA919'
data.table::fwrite(total_bulk, paste0(input_dir, 'metastasis_project_bulkRNAseq_libraries_metadata.csv'))

