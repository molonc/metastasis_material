
# annotation file
## 
master_anno <- master_df[, anno_varnames]

names(master_anno) <- gsub("__", "_",gsub(" ", "_", names(master_anno)))
names(master_anno) <- gsub("_$", "", names(master_anno))
names(master_anno) <- gsub("/|-", "_", names(master_anno))


# regexpr doesn't work with a specified column in a tibble, so change it to a
# data frame
master_anno <- data.frame(master_anno)



# master_anno <- merge(master_anno, tumour_main_anno[, c("SA_ID", "Tumor_ID")], by = "Tumor_ID", 
#                      all = TRUE)

#temp= unique(master_anno[, c("Tumor_ID", "Patient_ID", "SA_ID.x", "SA_ID.y")])
#check missing SA_ID
unique(master_anno[, c("Patient_ID", "SA_ID")])

#hard coded to get consistent SA_ID
master_anno$SA_ID <- ifelse(master_anno$SA_ID == "SA1142 (need to be checked)",
                            "SA1142", master_anno$SA_ID)


#check missing SA_ID
unique(master_anno[, c("Patient_ID", "SA_ID")])

# change DOB
master_anno <- date_conversion(master_anno,  varname = "DOB", 
                               replaceVarname = "DOB_Date")
master_anno$DOB_Date <- as.Date(master_anno$DOB, tryFormats = c("%d/%m/%Y"),
                                origin = "1970-01-01")


#visual checking
##DOB_date could be wrong as it is impossible to decide if 06/12/2017 is June 12, 2017 or Dec 6th, 2017
## as an example, 
temp=unique(master_anno[, c("row_number", "DOB_Date", "DOB")])

#quick check to make sure there is all non_NA data is converted properly
# should be 0 rows
unique(master_anno[is.na(master_anno$DOB_Date) & !is.na(master_anno$DOB),
                        c("row_number", "DOB_Date", "DOB")])


## Change Transplant_Date
## Same as DOB, could be wrong due to inconsistent date format
typeof(master_anno$Transplant_Date)


master_anno$Transplant_Date_Date <- as.Date(master_anno$Transplant_Date,
                                            tryFormats = c("%m/%d/%Y","%d/%m/%Y"))

#quick visual check 
temp=unique(master_anno[, c("row_number", "Transplant_Date", "Transplant_Date_Date")])


# check if the age at transplant is correct
master_anno[, c("Age_at_transplant_numeric")] <- as.numeric(gsub("wks", "", 
                                                                 master_anno[, c("Age_at_transplant")]))

master_anno$mouse_age_transplant <-difftime(master_anno$Transplant_Date_Date, 
                                            master_anno$DOB_Date, units = "weeks")

temp = master_anno[, c("Age_at_transplant", "Age_at_transplant_numeric",
                       "mouse_age_transplant", "DOB", 
                       "DOB_Date", "Transplant_Date")]

temp$diff = temp$Age_at_transplant_numeric - temp$mouse_age_transplant


#CHECK, should be 0 rows!!
temp[!is.na(temp$diff) & (temp$diff > 1|temp$diff < -1),]


master_anno$PTR_Date <- as.Date(master_anno$Primary_tumor_removal, tryFormats = c("%m/%d/%y"),
                       origin = "1970-01-01")

unique(master_anno[, c("Primary_tumor_removal", "PTR_Date")])

master_anno <- date_conversion(master_anno,  varname = "Euthanize", 
                               replaceVarname = "Euthanize_Date")

master_anno$Euthanize_Date <- as.Date(master_anno$Euthanize,
        tryFormats = c("%m/%d/%y"))

temptemp <- master_anno[, c("row_number", "Euthanize", "Euthanize_Date")]

master_anno$Euthanize_Date <- as.Date(ifelse(master_anno$Euthanize == "at 1yr old" & 
                                       !is.na(master_anno$DOB_Date), 
                                     master_anno$DOB_Date+365, 
                                    master_anno$Euthanize_Date),
                                    origin = "1970-01-01")


master_anno$SA_ID_Pas = ifelse(!is.na(master_anno$Passage), 
                               apply(cbind(master_anno$SA_ID, "_X", master_anno$Passage), 
                              1, paste0, collapse = ""), master_anno$SA_ID)

