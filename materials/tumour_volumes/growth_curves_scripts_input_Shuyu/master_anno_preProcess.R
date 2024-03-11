
# annotation file
## 
master_anno <- masterfile[, anno_varnames]

names(master_anno) <- gsub("__", "_",gsub(" ", "_", names(master_anno)))
names(master_anno) <- gsub("_$", "", names(master_anno))


# regexpr doesn't work with a specified column in a tibble, so change it to a
# data frame
master_anno <- data.frame(master_anno)



# master_anno <- merge(master_anno, tumour_main_anno[, c("SA_ID", "Tumor_ID")], by = "Tumor_ID", 
#                      all = TRUE)

#temp= unique(master_anno[, c("Tumor_ID", "Patient_ID", "SA_ID.x", "SA_ID.y")])
#check missing SA_ID
unique(master_anno[, c("Patient_ID", "SA_ID")])

#hard coded to get consistent SA_ID
#From ATLAS CID0003463 = SA1019, CID0003465=SA1034 and CID0003466 as you found is SA1020
master_anno$SA_ID <- ifelse(is.na(master_anno$SA_ID) & 
                             master_anno$Patient_ID == "VBA0847", 
                           "SA919",
                           ifelse(is.na(master_anno$SA_ID) & 
                                    master_anno$Patient_ID == "CID0003463", 
                                  "SA1019", 
                                  ifelse(is.na(master_anno$SA_ID) & 
                                           master_anno$Patient_ID == "CID0003465", 
                                         "SA1034", 
                                         ifelse(is.na(master_anno$SA_ID) & 
                                                  master_anno$Patient_ID == "CID0003466",
                                                "SA1020",master_anno$SA_ID))))


master_anno$Patient_ID <- ifelse(is.na(master_anno$Patient_ID) & 
                             master_anno$SA_ID == 'SA919', 
                             "VBA0847", master_anno$Patient_ID)

#check missing SA_ID
unique(master_anno[, c("Patient_ID", "SA_ID")])

# change DOB
master_anno <- date_conversion(master_anno,  varname = "DOB", 
                               replaceVarname = "DOB_Date")

#visual checking
temp=unique(master_anno[, c("row_number", "DOB_Date", "DOB")])

#quick check to make sure there is no UN_NA data is accidentally converted to NA
# should be 0 rows
temp=unique(master_anno[is.na(master_anno$DOB_Date) & !is.na(master_anno$DOB),
                        c("row_number", "DOB_Date", "DOB")])


## Change Transplate_Date
typeof(master_anno$Transplant_Date)

master_anno$Transplant_Date_Date <- as.Date(master_anno$Transplant_Date,
                                            tryFormats = c("%m/%d/%Y","%d/%m/%Y"))



# check if the age at transplante is correct
master_anno[, c("Age_at_transplant_numeric")] <- as.numeric(gsub("wks", "", 
                                                                 master_anno[, c("Age_at_transplant")]))

master_anno$DOB_Date = as.Date(ifelse(is.na(master_anno$DOB_Date) & !is.na(master_anno$Age_at_transplant) &
                                        !is.na(master_anno$Transplant_Date), 
                                      as.Date(master_anno$Transplant_Date - master_anno$Age_at_transplant_numeric*7*24*3600), 
                                      master_anno$DOB_Date), origin = "1970-01-01")


master_anno$mouse_age_transplant <-difftime(master_anno$Transplant_Date_Date, 
                                            master_anno$DOB_Date, units = "weeks")

temp = master_anno[, c("Age_at_transplant", "Age_at_transplant_numeric",
                       "mouse_age_transplant", "DOB", 
                       "DOB_Date", "Transplant_Date")]

temp$diff = temp$Age_at_transplant_numeric - temp$mouse_age_transplant

#CHECK, should be 0 rows
temp[!is.na(temp$diff) & temp$diff > 1,]


master_anno <- date_conversion(master_anno,  varname = "Primary_tumor_removal", 
                               replaceVarname = "PTR_Date")


master_anno <- date_conversion(master_anno,  varname = "Euthanize", 
                               replaceVarname = "Euthanize_Date")

temptemp <- master_anno[, c("row_number", "Euthanize", "Euthanize_Date")]

master_anno$Euthanize_Date <- ifelse(master_anno$Euthanize == "at 1yr old", 
                                     master_anno$DOB_Date+365, 
                                    master_anno$Euthanize_Date)

master_anno$Euthanize_Date <- as.Date(master_anno$Euthanize_Date,
                                      origin = "1970-01-01")
master_anno$SA_ID_Pas = ifelse(!is.na(master_anno$Passage), 
                               apply(cbind(master_anno$SA_ID, "_X", master_anno$Passage), 
                              1, paste0, collapse = ""), master_anno$SA_ID)

