library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)


read_raw_tumour_volume_from_masterfile <- function(){
  ## From master file xlsx, first I manually convert this file into csv file
  input_dir <- '/Users/htran/Documents/storage_tmp/metastasis_metadata/tumour_growth/tumour_growth_processed/'
  master_file <- data.table::fread(paste0(input_dir,'Metastasis_mouse_exp_HL_masterfile.csv'))
  
  # as.numeric(grep('V',colnames(master_file)))
  # "V45" to "V122" 
  # 122-45 + 1 # 78 cols
  cols_used <- c()
  for(i in seq(1, 39, by=1)){
    dt <- paste0('DATE', i)  
    tv <- paste0('TumourVol', i)  
    cols_used <- c(cols_used, dt)
    cols_used <- c(cols_used, tv)
  }
  cols_used
  length(cols_used)
  other_cols <- colnames(master_file)[!colnames(master_file) %in% paste0('V', seq(45, 122, by=1))]   
  other_cols <- stringr::str_replace(other_cols, ' ','_')
  colnames(master_file) <- c(other_cols, cols_used)
  data.table::fwrite(master_file, paste0(input_dir,'Metastasis_mouse_exp_HL_masterfile_processed.csv'))

}

extract_tumour_growth <- function(df, datatag=''){
  
  ## There is no column name for tumour volume in the master file, so I look at master file and annotate them here. 
  cols_used <- c()
  for(i in seq(1, 39, by=1)){
    dt <- paste0('DATE', i)  
    tv <- paste0('TumourVol', i)  
    cols_used <- c(cols_used, dt)
    cols_used <- c(cols_used, tv)
  }
  # cols_used
  length(cols_used)
  
  ## Currently, in master file, data are recorded in wide format, ex: tumour id, date1, vol1, date2, vol2,...
  ## Extracting each record date, and tumour growth volumn for each transplanted id
  ## and put in data frame with long format. 
  datetime_format <- "%m/%d/%y"
  obs_sa <- unique(df$SA_ID)
  obs_sa <- obs_sa[obs_sa!=""]
  
  total <- tibble::tibble()
  meta_tumours <- tibble::tibble()
  
  ## Processing data for each SA id, ex: SA1142, SA919, SA535,...
  for(sa_id in obs_sa){
    tm_df <- df %>%
      dplyr::select(Tumor_ID, SA_ID, Patient_ID, Transplant_Date, all_of(cols_used)) %>%
      dplyr::filter(SA_ID==sa_id)
    print(dim(tm_df))
    
    ## For each line of record is one tumour id, extracting tumour vols for each tumour ids
    for(id in unique(tm_df$Tumor_ID)){
      tmp <- tm_df %>%
        dplyr::filter(Tumor_ID==id)
      # print(dim(tmp))
      transplant_dt <- as.Date(as.character(tmp[1,'Transplant_Date']), datetime_format) # default format in this file is datetime_format <- "%m/%d/%y"
      if(is.na(transplant_dt)){
        transplant_dt <- as.Date(as.character(tmp[1,'Transplant_Date']), "%d/%m/%y") # some exceptions, sometimes date time records are not in the default format above
      }
      
      ## Number of tumour volume records depending on each transplant ids, so remove some empty columns
      empty_columns <- colSums(is.na(tmp) | tmp == "") == nrow(tmp)
      not_empty_columns <- names(empty_columns)[!empty_columns]
      not_empty_columns <- union(not_empty_columns, c('Tumor_ID','SA_ID','Transplant_Date'))
      
      if(sum(grepl('DATE',not_empty_columns))>0){  ## a record of tumour growth exists for this tumour id
        # print('Existing tumour vol record for this tumour id')
        # print(id)
        mt <- tibble::tibble(SA_ID=sa_id, Tumor_ID=id)
        meta_tumours <- dplyr::bind_rows(meta_tumours, mt)
        
        tmp <- tmp %>%
          dplyr::select(all_of(not_empty_columns)) %>%
          as.data.frame()
        
        for(obs_dt in colnames(tmp)[grepl('DATE', colnames(tmp))]){
          i <- as.numeric(gsub('DATE','',obs_dt))
          curr_str <- as.character(tmp[1,paste0('DATE', i)])
          
          ## Sometimes there are exception with date time format as input
          curr_str <- ifelse(grepl('2018',curr_str),gsub('2018','18',curr_str),
                             ifelse(grepl('2019',curr_str),gsub('2019','19',curr_str),curr_str))
          # curr_dt <- ifelse(!is.na(as.Date(curr_str, datetime_format)), as.Date(curr_str, datetime_format),
          #                   as.Date(curr_str, "%d/%m/%y"))
          curr_dt <- as.Date(curr_str, datetime_format)
          # if(is.na(curr_dt2)){
          #   curr_dt <- curr_dt1
          # }  
          # if(is.na(curr_dt)){
          #   # str_dt <- as.character(tmp[1,paste0('DATE', i)])
          #   # str_dt <- gsub('2018','18',str_dt)
          #   # str_dt <- gsub('2019','19',str_dt)
          #   # str_dt <- gsub('2020','20',str_dt)
          #   curr_dt <- as.Date(curr_str, "%d/%m/%y")
          # }
          # print(curr_dt)
          # print(as.numeric(curr_dt - transplant_dt))
          t <- tibble::tibble(Tumor_ID=tmp[1,'Tumor_ID'], SA_ID=tmp[1,'SA_ID'], Transplant_Date=as.character(tmp[1,'Transplant_Date']), 
                              date=as.character(tmp[1,paste0('DATE', i)]), 
                              date_len=as.numeric(curr_dt - transplant_dt),
                              TumourVol=as.character(tmp[1,paste0('TumourVol', i)]))
          total <- dplyr::bind_rows(total, t)
        }
        
        
      }
      
    }
  }
  
  print(dim(total))
  total$week_len <- convert_days2weeks(total$date_len)  ## convert from date length to the weeks time
  total$datetime_issue <- ifelse(total$week_len==-1,'double_check','')  ## There are some cases with negative time, due to the issue with date time format, so we put a flag here. 
  print("Nb issue with date time: ")
  print(sum(total$week_len==-1))
  data.table::fwrite(total, paste0(input_dir,'Metastasis_mouse_exp_tumour_volumes_',datatag, '.csv'))
  data.table::fwrite(meta_tumours, paste0(input_dir,'Metastasis_mouse_exp_tumour_volumes_metadata_',datatag, '.csv'))
  
  
  return(list(meta_tumours=meta_tumours, tumour_growth_vols=total))
  
}

## Main experiment contains different series: SA535, SA919, SA1142, SA605, SA609,...
extract_main_experiment <- function(){
  
  ## From master file xlsx, first I manually convert this file into csv file
  ## See the function: read_raw_tumour_volume_from_masterfile()
  # read_raw_tumour_volume_from_masterfile()
  input_dir <- '/Users/htran/Documents/storage_tmp/metastasis_metadata/tumour_growth/tumour_growth_processed/'
  master_file <- data.table::fread(paste0(input_dir,'Metastasis_mouse_exp_HL_masterfile_processed.csv')) %>% as.data.frame()
  dim(master_file)  
  # master_file <- master_file[!duplicated(master_file$Tumor_ID),]   ## Manually, removing some blank records from this file
  # master_file$Tumor_ID[duplicated(master_file$Tumor_ID)]
  
  ## The list of tumour ids are noted at this file: tumorsize_at_removal_mets_list.csv
  ## Check whether master file contains all different series tumour growth at main experiment
  tumour_main_exp <- data.table::fread(paste0(input_dir,'tumorsize_at_removal_mets_list.csv'))
  # dim(tumour_main_exp)  
  colnames(tumour_main_exp) <- gsub(' ','_',colnames(tumour_main_exp)) # remove space from column names
  length(unique(tumour_main_exp$Tumor_ID))
  sum(unique(tumour_main_exp$Tumor_ID) %in% master_file$Tumor_ID) # good, all tumour ids infos are noted at master file. 
  
  
  ## Check whether master file contains all tumour growth for mixing exp record. 
  main_exp_ls <- unique(tumour_main_exp$Tumor_ID)
  length(unique(main_exp_ls))
  ## Main experiment
  
  # rownames(master_file) <- master_file$Tumor_ID
  master_file1 <- master_file[master_file$Tumor_ID %in% main_exp_ls,]
  dim(master_file1)
  length(unique(master_file1$Tumor_ID))
  # View(master_file1)
  
  ## Extracting tumour volume for each tumour ids and save data as csv file at: Metastasis_mouse_exp_tumour_volumes_main_exp.csv
  res_main <- extract_tumour_growth(master_file1, datatag='main_exp')
  # dim(res_main$meta_tumours)
  
  ## Extracting meta data, we will decide the status 'met', 'non_met' for each tumour id based on these descriptions
  meta_data <- res_main$meta_tumours
  # dim(meta_data)
  cols_use <- c('Tumor_ID','Passage','Source_type','V23','V24','Histology')
  t <- master_file1 %>%
    dplyr::select(all_of(cols_use))
  dim(t)
  meta_data <- meta_data %>%
    dplyr::left_join(t, by='Tumor_ID')
  data.table::fwrite(meta_data, paste0(input_dir, 'Metastasis_mouse_exp_tumour_volumes_metadata_main_exp_annotated.csv'))
  dim(meta_data)
}  


## Mixing experiment for series SA919 (extra experiment, different from main experiment)
extract_mixing_experiment <- function(){
  
  ## From master file xlsx, first I manually convert this file into csv file
  input_dir <- '/Users/htran/Documents/storage_tmp/metastasis_metadata/tumour_growth/tumour_growth_processed/'
  master_file <- data.table::fread(paste0(input_dir,'Metastasis_mouse_exp_HL_masterfile_processed.csv')) %>% as.data.frame()
  dim(master_file)  
  # master_file <- master_file[!duplicated(master_file$Tumor_ID),]
  # master_file$Tumor_ID[duplicated(master_file$Tumor_ID)]
  
  ## The list of tumour ids are noted at this file: tumorsize_at_removal_mets_list.csv
  ## Check whether master file contains all different series tumour growth at main experiment
  tumour_main_exp <- data.table::fread(paste0(input_dir,'tumorsize_at_removal_mets_list.csv'))
  # dim(tumour_main_exp)  
  colnames(tumour_main_exp) <- gsub(' ','_',colnames(tumour_main_exp)) # remove space from column names
  length(unique(tumour_main_exp$Tumor_ID))
  sum(unique(tumour_main_exp$Tumor_ID) %in% master_file$Tumor_ID) # good, all tumour ids infos are noted at master file. 
  
  
  ## Check whether master file contains all tumour growth for mixing exp record. 
  ## Main experiment
  main_exp_ls <- unique(tumour_main_exp$Tumor_ID)
  length(unique(main_exp_ls))
  
  
  ## Mixing experiment SA919
  master_file2 <- master_file[!master_file$Tumor_ID %in% main_exp_ls,]
  dim(master_file2)
  names(master_file)
  master_file2 <- master_file2[grepl('X0847',master_file2$Tumor_ID),] # SA919 series, mouse ids start with X0847
  master_file2 <- master_file2[!grepl('211124',master_file2$Tumor_ID),] # Remove some exceptions, not belong to SA919 mixing experiment
  
  dim(master_file2)
  unique(master_file2$Tumor_ID)
  ## TO DO: define list of mixing exp for master_file2
  # ex: '2164-2112252', '2164-2112253', '21122516'
  
  ## Extracting tumour volume for each tumour ids and save data as csv file at: Metastasis_mouse_exp_tumour_volumes_mixing_exp.csv
  res_mixing <- extract_tumour_growth(master_file2, datatag='mixing_exp')
  dim(res_mixing$meta_tumours)
  meta_data <- res_mixing$meta_tumours
  dim(meta_data)
  cols_use <- c('Tumor_ID','Passage','Source_type','V23','V24','Histology')
  t <- master_file2 %>%
    dplyr::select(all_of(cols_use))
  dim(t)
  meta_data <- meta_data %>%
    dplyr::left_join(t, by='Tumor_ID')
  data.table::fwrite(meta_data, paste0(input_dir, 'Metastasis_mouse_exp_tumour_volumes_metadata_mixing_exp.csv'))
  
  
  ## Extracting meta data, we will decide the status 'met', 'non_met' for each tumour id based on these descriptions
  meta_data <- data.table::fread(paste0(input_dir, 'Metastasis_mouse_exp_tumour_volumes_metadata_mixing_exp.csv'))
  dim(meta_data)
  ## Manual labelling of mixing exp tumour ids
  meta_data <- meta_data %>%
    dplyr::filter(Mixing_exp=='Yes')
  data.table::fwrite(meta_data, paste0(input_dir, 'Metastasis_mouse_exp_tumour_volumes_metadata_mixing_exp_annotated.csv'))
  
  meta_data <- data.table::fread(paste0(input_dir, 'Metastasis_mouse_exp_tumour_volumes_metadata_mixing_exp_annotated.csv'))
  dim(meta_data)
  master_file2 <- master_file2[master_file2$Tumor_ID %in% meta_data$Tumor_ID,]
  dim(master_file2)
  res_mixing <- extract_tumour_growth(master_file2, datatag='mixing_exp')
  
}

