
library(dplyr)
library(data.table)

input_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/revision/jira_tickets_metproj_Hoa/'


## check 535  
df1 <- data.table::fread(paste0(input_dir, 'SA535_jira_tickets_missing_libraries.csv'))
df2 <- data.table::fread(paste0(input_dir, 'metadata_cn_SA535.csv'), header = T)
df3 <- data.table::fread(paste0(input_dir, 'SA535_jira_tickets.csv'))
df5 <- data.table::fread(paste0(input_dir, 'SA535X8_library_groupings.csv'))

head(df1)
head(df2)
head(df3)
head(df5)
df5 <- df5 %>%
  dplyr::rename(library_id=library_labels) %>%
  dplyr::select(library_id, jira_ticket)

df2 <- df2 %>%
  dplyr::select(library_id, jira_ticket)

df <- dplyr::bind_rows(df1, df2, df3, df5)
dim(df)
df$desc <- paste0(df$library_id, df$jira_ticket)
duplicated(df$desc)
df <- df %>% distinct()
data.table::fwrite(df, paste0(input_dir, 'SA535_hmmcopy_libraryIds_jira.csv'))
df <- data.table::fread(paste0(input_dir, 'SA535_hmmcopy_libraryIds_jira.csv'))


atlas <- data.table::fread(paste0(input_dir, 'AX_library_uid_SA535.txt'), header = F)
dim(atlas)
head(atlas)
colnames(atlas) <- c('atlas_AX_id','library_id')
df <- df %>% 
  dplyr::left_join(atlas, by='library_id')
dim(df)
View(df)



# check SA919
df1 <- data.table::fread(paste0(input_dir, 'SA919X7_metadata_Hoa.csv'))
df2 <- data.table::fread(paste0(input_dir, 'ticket_SA919.csv'), header = T)
head(df1)
head(df2)

df1 <- df1 %>%
  dplyr::rename(library_id=!!sym('Library ID_DLP+'), jira_ticket=!!sym('DLP+_jira_ticket')) %>%
  dplyr::select(library_id, jira_ticket)

df2 <- df2 %>%
  dplyr::select(library_id, jira_ticket)

df <- dplyr::bind_rows(df1, df2)
dim(df)
df <- df %>% distinct()


atlas <- data.table::fread(paste0(input_dir, 'AX_library_uid_SA919.txt'), header = F)
dim(atlas)
head(atlas)
colnames(atlas) <- c('atlas_AX_id','library_id')
df <- df %>% 
  dplyr::full_join(atlas, by='library_id')
dim(df)
View(df)

data.table::fwrite(df, paste0(input_dir, 'SA919_hmmcopy_libraryIds_jira.csv'))


df1 <- data.table::fread(paste0(input_dir, 'SA1142_library_groupings.csv.gz'))
head(df1)
df1$nb_filtered_cells
df1 <- df1 %>%
  dplyr::rename(library_id=library_labels) %>%
  dplyr::select(library_id, sample_id)

data.table::fwrite(df1, paste0(input_dir, 'SA1142_hmmcopy_libraryIds.csv'))

df <- tibble::tibble(dataset=c('SA919','SA535'), jira_ticket=c('SC-4094','SC-4095'))

data.table::fwrite(df, paste0(input_dir, 'SA919_SA535_snv_genotyping.csv'))




input_dir <- '/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/revision/jira_tickets_metproj_Hoa/pseudobulk/'


## check 535  
df1 <- data.table::fread(paste0(input_dir, 'ascn_libraryID_SA535.csv'))
df2 <- data.table::fread(paste0(input_dir, 'pseudobulk_lib_ticket_SA535.csv'))
df3 <- data.table::fread(paste0(input_dir, 'pseudobulk_ticket_SA535_v1.csv'))
df4 <- data.table::fread(paste0(input_dir, 'pseudobulk_ticket_SA535.csv'))
df5 <- data.table::fread(paste0(input_dir, 'pseudobulk_tickets_SA535_metastasis.csv'))
df6 <- data.table::fread(paste0(input_dir, 'SA535_pseudobulkticket_SA535.csv'), header = T)


head(df1)
head(df2)
head(df3)
head(df4)
head(df5)
head(df6)
dim(df2)
sum(df2$library_id %in% df3$library_id)
sum(df4$library_id %in% df3$library_id)
sum(df5$library_id %in% df3$library_id)
sum(df6$library_id %in% df3$library_id)
dim(df5)
dim(df3)

df3 <- df3 %>%
  dplyr::select(library_id, jira_ticket)
data.table::fwrite(df3, paste0(input_dir, 'SA919_SA535_pseudobulk.csv'))



