


#quick check for names
names(master_df)

#row_number is used to keep the order of original excel file
master_df$row_number <- c(1:nrow(master_df))+1



## get the names of the variables that only appear once, i.e NOT time related as in tumour growth
anno_varnames <- c("row_number", "experiment_type", "Damian_annotation","DOB",    
                   "Transplant_Date",
                   "V3", "Transplant_site", "Age_at transplant","Strain",
                   "Tumor_ID","Patient_ID",      "SA_ID",  "Passage",
                   "Source_type",     "Treatment_naive/pre","Date_of Dx","Age_at Dx",
                   "Metastasis...17", "Mets_at Dx","Date_of mets Dx", "Mets_site",
"SQL_Barcode", "Cage", "Comments_for surgery", "Wound_dehiscence",
 "V23", "V24", "DLP_ID", "10X_ID",
 "Intrathoracic_Mass", "Histology", "Cell_suspension vol (ul)", "Last_tumor palpation",    
"V31", "Next_tumor palpation","V33",    "Primary_tumor removal",   
 "12_wks from PTR", "16_wks from PTR", "Euthanize", "PTR-Mets", 
 "Time_to mets detection","PET/CT", "Second_surgery",  "8wks_after PTR", 
 "Metastasis...45", "metastatic_site" )

## Check to see if this the tumour growth Date and Volume
ind_var = names(master_df)[!names(master_df)%in% anno_varnames]

# Preprocess the annotation data, 
source(paste0(sourcedir, "master_anno_preProcess_2024.R"))

## get the tumour growth part, here only use row_numbers to uniquely join with annotation
## using the row number is to prevent any accidently un-unique tumour id etc to enter the data
master_df_measure = master_df[, c("row_number", ind_var)]

## Pivot the tumor measure Date from wide to long
master_df_measure_long1 <- master_df_measure[, c("row_number", ind_var[grepl("Date", ind_var, ignore.case = TRUE)])] %>% 
  tidyr::pivot_longer(cols = ind_var[grepl("Date", ind_var, ignore.case = TRUE)], 
                      names_to = c("Seq"),
                      values_to = c("TumorDate"), 
                      names_prefix = "TumorDate_",
                      values_drop_na = TRUE)

master_df_measure_long1$Seq = as.numeric(gsub("Date", "",
                                              master_df_measure_long1$Seq, ignore.case = TRUE))

#change inconsistent data format for Tumor volume to character
#
master_df_measure[, c(ind_var[grepl("TumourVol", ind_var, ignore.case = TRUE)])] = sapply( master_df_measure[, c(ind_var[grepl("TumourVol", ind_var, ignore.case = TRUE)])], as.character)

## Pivot the tumor volume from wide to long
master_df_measure_long2 <- master_df_measure[, c("row_number", 
                                                 ind_var[grepl("TumourVol", ind_var, ignore.case = TRUE)])] %>% 
  tidyr::pivot_longer(cols = ind_var[grepl("TumourVol", ind_var, ignore.case = TRUE)], 
                      names_to = c("Seq"),
                      values_to = c("TumorVol"), 
                      names_prefix = "TumorVol_",
                      values_drop_na = TRUE)

master_df_measure_long2$Seq = as.numeric(gsub("TumourVol", "",
                                              master_df_measure_long2$Seq))

master_df_measure_long <- merge(master_df_measure_long1,
                                 master_df_measure_long2, 
                                 by = c("row_number", "Seq"), 
                                 all = TRUE)


master_df_measure_long$Seq = as.numeric(master_df_measure_long$Seq)


master_df_measure_long <-master_df_measure_long[with(master_df_measure_long,
                                                       order(row_number, Seq)),] 

master_df_measure_long <- merge(master_df_measure_long, master_anno[,c("row_number","Tumor_ID", "Transplant_Date_Date", 
                                                                         "PTR_Date")],
                                 by = "row_number")

# change the tumor date to proper date. some with "%m/%d/%y" format
master_df_measure_long$TumorDate_Date <- as.Date(master_df_measure_long$TumorDate, tryFormats = c("%m/%d/%y"),
                                                 origin = "1970-01-01")

# change the tumor date to proper date. some with "%d/%m/%y" format
for(i in c(1:nrow(master_df_measure_long))){
  if(is.na(master_df_measure_long$TumorDate_Date[i]) &
           !is.na(master_df_measure_long$TumorDate[i])){
             master_df_measure_long$TumorDate_Date[i] <- as.Date(master_df_measure_long$TumorDate[i],
                                                                                tryFormats = c("%d/%m/%Y"),
                                                                 origin = "1970-01-01")
             
           }
}


#check for NA or outlier date and hardcode to correct them if needed
temptemp <- unique(master_df_measure_long[, c("row_number", "Seq", "TumorDate", "TumorDate_Date")])


#change all tumor volume to numerical numbers
master_df_measure_long$TumorVol_num = as.numeric(master_df_measure_long$TumorVol)

#quick check to see which tumor size could not be changed to numerical values
#
temptemp = master_df_measure_long[is.na(master_df_measure_long$TumorVol_num), 
                     c('row_number', "Seq", "TumorVol", "TumorVol_num")]


#missing measurement date or tumor size
temptemp = master_df_measure_long[is.na(master_df_measure_long$TumorDate_Date)|is.na(master_df_measure_long$TumorVol_num),]


#Note, if no tumor measuring date, then it canNOT be decided if the measurement took place before PTR
master_df_measure_long$curve = ifelse(!is.na(master_df_measure_long$TumorDate_Date) &
                                         !is.na(master_df_measure_long$PTR_Date) &
                                         master_df_measure_long$TumorDate_Date <= master_df_measure_long$PTR_Date,
                                       "before PTR", 
                                       ifelse(!is.na(master_df_measure_long$TumorDate_Date) &
                                                              !is.na(master_df_measure_long$PTR_Date) &
                                                              master_df_measure_long$TumorDate_Date > master_df_measure_long$PTR_Date,
                                              "after PTR", NA))
  
# quick check
temp=unique(master_df_measure_long[, c("row_number","Seq", "TumorDate", "TumorDate_Date", "PTR_Date")])

## This is used to 
##quick check to see if the tumor size decreases 
## after discussing with the Sam and other team members, it is possible to have tumor size smaller than 
## the previous measurement due to measurement error. 
## I am not aware of the e-records were checked against the paper or original records. 
master_df_measure_long <- master_df_measure_long %>%
  group_by(row_number) %>%
  mutate(TumorVol_num_lag1 = lag(TumorVol_num, n=1), 
         TumorVol_num_lead1 = lead(TumorVol_num, n=1),
         TumorDate_Date_lag1 = lag(TumorDate_Date, n=1), 
         Seq_lag1 = lag(Seq, n=1),
         same_curve = lag(curve, n= 1) == curve)

master_df_measure_long$Seq_diff = master_df_measure_long$Seq - master_df_measure_long$Seq_lag1
master_df_measure_long$TumorVol_diff = master_df_measure_long$TumorVol_num - master_df_measure_long$TumorVol_num_lag1
master_df_measure_long$TumorDate_diff = master_df_measure_long$TumorDate_Date - master_df_measure_long$TumorDate_Date_lag1

temptemp <- master_df_measure_long[!is.na(master_df_measure_long$Seq_diff) & master_df_measure_long$Seq_diff !=1,
                     c("row_number", "Seq", "Tumor_ID", "Seq_diff", "TumorDate", "TumorVol",
                       "TumorVol_diff", "TumorDate_diff")]

#decreasing size 
temptemp <- master_df_measure_long[(!is.na(master_df_measure_long$TumorVol_diff) & master_df_measure_long$TumorVol_diff<0)|
                                     is.na(master_df_measure_long$TumorDate)|is.na(master_df_measure_long$TumorVol),
                                   
                     c("row_number", "Seq", "Seq_diff", "TumorDate", 
                       "TumorVol_diff")]

temptemp <- master_df_measure_long[!is.na(master_df_measure_long$TumorVol_diff) & master_df_measure_long$TumorVol_diff<0,
                                   c("row_number", "Seq", "Tumor_ID", "Seq_diff", "TumorDate", "TumorVol",
                                     "TumorVol_diff", "TumorDate_diff", "curve", "same_curve")]


temptemp <- merge(temptemp, master_anno[, c("row_number","Cage", "SQL_Barcode", "experiment_type", "Damian_annotation" )], by = "row_number", 
                  all.x = TRUE)

#only focus on main event and first curve event for now
temptemp <- temptemp[temptemp$experiment_type == "main_experiment" & temptemp$curve == "before PTR",]

temptemp <- temptemp[, c( "row_number","Cage",     
              "SQL_Barcode" ,"Tumor_ID",  "TumorDate","TumorVol",
              "TumorVol_diff","curve",     "same_curve", "experiment_type", "Damian_annotation")]
#write.csv(temptemp, file = paste0(resultdir, "datacheck_main_tumoursize.csv"))






temp_long = merge(master_df_measure_long[, !names(master_df_measure_long) %in% 
                                            c("Tumor_ID", "Transplant_Date_Date")], 
                  master_anno, by = c("row_number"))

temp_long$date_len = as.numeric(difftime(temp_long$TumorDate_Date,
                              temp_long$Transplant_Date_Date, units = "days"))


temp_long$growth_rate <- temp_long$TumorVol_diff/as.numeric(temp_long$TumorDate_diff) 


#This is the data used for growth curve.
#missing tumor date will not 
firstcurve_df <- temp_long[!is.na(temp_long$curve) & temp_long$curve == "before PTR",]

#quick check to see if the number of Tumorr_ID is the same as the number of rows originally in the master_df
length(unique(firstcurve_df$Tumor_ID)) == dim(master_df)


