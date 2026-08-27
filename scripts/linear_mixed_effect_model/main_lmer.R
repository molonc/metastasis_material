
#import libraries
library(readr)
library(readxl)

# plot libraries
library(ggplot2)

#data manipulate libraries
library(dplyr)

library(tidyr)

library(RColorBrewer)

library(data.table)

# load extra font for publishing
library(extrafont)

#fit lmer model
library(lme4)
library(lmerTest)

maindir <- "C:/Documents/Metastasis/"
sourcedir <- paste0(maindir, "R_files/")
datadir <- paste0(maindir, "Data/")
today <- format(Sys.Date(), format = "%Y%m%d")
resultdir <- paste0(maindir, "Result", today, "/")

# change inconsistent date format to consistent date format. 
# note: need to double check manually to ensure all date is changed properly
date_conversion <- function(data = master_anno, varname = "DOB", 
                            replaceVarname = "DOB_Date"){
  #keep date format consistent
  
  regex = regexpr("\\d{5}", data[, varname])
  regex1 = regexpr("(\\d{2}\\/\\d{2}\\/\\d{4})|(\\d{4}(\\/|-)\\d{2}(\\/|-)\\d{2})", data[, varname])
  
  # temp_df = data.frame(varname = data[, varname], reg1 = regex, 
  #                      reg2 = regex1, index = data[, "index"])
  data$temptemp <- ifelse(is.na(data[, varname]), NA,
                          ifelse(regex>0, as.Date(as.numeric(data[, varname]),
                                                  origin = "1899-12-30"), 
                                 ifelse(regex1 > 0,
                                        as.Date(substr(data[, varname],
                                                       regex1,
                                                       regex1+attr(regex1, "match.length")-1),
                                                tryFormats = c("%m/%d/%Y","%d/%m/%Y", 
                                                               "%Y/%m/%d", "%Y-%m-%d")), 
                                        NA)))
  
  # temp_df = data.frame(varname = data[, varname], reg1 = regex, 
  #                      reg2 = regex1, 
  #                      date_str = substr(data[, varname],
  #                                            regex1,
  #                                            regex1+attr(regex1, "match.length")-1),
  #                      index = data[, "index"], data$temptemp)
  names(data)[names(data)=="temptemp"] = replaceVarname
  
  data[, replaceVarname] = as.Date(data[, replaceVarname],origin = "1970-01-01")
  #data = data[ , -which(names(data) %in% c("temptemp"))]
  return(data)
}

if(!dir.exists(resultdir)) dir.create(resultdir, recursive = TRUE)


## all these files are used only to obtain the main or mixing category. 
# Only main is used for the paper
tumour_main_exp <- read_csv(paste0(datadir, 
                                   "Metastasis_mouse_exp_tumour_volumes_main_exp.csv"))
#to keep the same order as in csv file
tumour_main_exp$index = c(1:nrow(tumour_main_exp))

tumour_main_anno <- read_csv(paste0(datadir, 
                                    "Metastasis_mouse_exp_tumour_volumes_metadata_main_exp_annotated.csv"))

tumour_mixing_exp <- read_csv(paste0(datadir, 
                                     "Metastasis_mouse_exp_tumour_volumes_mixing_exp.csv"))

tumour_mixing_anno <- read_csv(paste0(datadir, 
                                      "Metastasis_mouse_exp_tumour_volumes_metadata_mixing_exp_annotated.csv"))



mix_tumour <- data.frame(Tumor_ID = tumour_mixing_anno$Tumor_ID, category = "Mixing",
                         mainsite = tumour_mixing_anno$mainsite)
main_tumour <- data.frame(Tumor_ID = tumour_main_anno$Tumor_ID, category = "Main", 
                          mainsite = tumour_main_anno$mainsite)

tumour_id_df = rbind(mix_tumour, main_tumour)
## end of ## all these files are used only to obtain the main or mixing category. 


#Damian's category of met and non-met
met_cat<- read_csv(paste0(datadir, "Damian_2024_01_10_Met_Cat.csv"))


# 
# read in here to guarantee a fresh start as I need to keep the order of the data

masterfile <- read_excel(paste0(datadir, 
                                "Metastasis_mouse_exp_HL_masterfile.xlsx"))
#quick check for names
names(masterfile)

#row_number is used to keep the order of original excel file
masterfile$row_number <- c(1:nrow(masterfile))+1

# change a few names to be more meaningful
names(masterfile)[names(masterfile) == "...3"] = "Mouse_ID"

#What is this?
#names(masterfile)[names(masterfile) == "Metastasis...15"] = ""

names(masterfile)[names(masterfile) == "...23"] = "Comments_1"

names(masterfile)[names(masterfile) == "...24"] = "Comments_2"

names(masterfile)[names(masterfile) == "...31"] = "Palpable"

# What date is this one? Column AG, keep the original name as it is for now
# names(masterfile)[names(masterfile) == "...33"] = "Date_??"

anno_varnames <- c("row_number", "DOB","Transplant Date","Mouse_ID","Transplant site",
                   "Age at transplant","Strain",
                   "Tumor ID","Patient ID","SA ID", 
                   "Passage","Source type","Treatment naive/pre",
                   "Date of Dx", "Age at Dx","Metastasis...15", 
                   "Mets at Dx", "Date of mets Dx", "Mets site", 
                   "SQL Barcode", "Cage", "Comments for surgery",
                   "Wound dehiscence", "Comments_1", "Comments_2", 
                   "DLP ID", "10X ID", "Intrathoracic Mass",
                   "Histology", "Cell suspension vol (ul)",
                   "Last tumor palpation",
                   "Palpable", "Next tumor palpation", "...33", 
                   "Primary tumor removal", "12 wks from PTR",
                   "16 wks from PTR","Euthanize",
                   "PTR-Mets","Time to mets detection",
                   "PET/CT", "Second surgery", "8wks after PTR",
                   "Metastasis...43","metastatic site")

ind_var = names(masterfile)[!names(masterfile)%in% anno_varnames]

source(paste0(sourcedir, "master_anno_preProcess.R"))


#check to see if the col_row_number is monotone increasing 
# these variables are used to record tumour size and date in the master data file
regex = regex = regexpr("\\d+", ind_var)
col_row_number = as.numeric(substr(ind_var,
                                   regex,
                                   regex+attr(regex, "match.length")-1))-44
increase_row_number = c(1:(length(col_row_number)))
col_row_number[col_row_number!=increase_row_number]


#rename the column to reflect the tumour measure date and size 
masterfile_measure = masterfile[, c("row_number", ind_var)]
rename_var = ifelse(col_row_number %% 2 == 1, 
                    paste0("Tumour_Date_", (col_row_number+1)/2 %% ((length(col_row_number)/2)+1) ),
                    paste0("Tumour_Size_", col_row_number/2 %% ((length(col_row_number)/2)+1)))
names(masterfile_measure) = c("row_number", rename_var)


masterfile_measure = data.frame(masterfile_measure)
for(i in c(1:length(rename_var))){
  masterfile_measure[, rename_var[i]] = as.character(masterfile_measure[, rename_var[i]])
}


# get the Tumour_date in long format
masterfile_measure_long1 <- masterfile_measure[, c("row_number", rename_var[c(1:(length(rename_var)/2))*2-1])] %>% 
  tidyr::pivot_longer(cols = rename_var[c(1:(length(rename_var)/2))*2-1], 
                      names_to = c("Seq"),
                      values_to = c("Tumour_Date"), 
                      names_prefix = "Tumour_Date_",
                      values_drop_na = TRUE)

## get the Tumour_size in long format
masterfile_measure_long2 <- masterfile_measure[, c("row_number", 
                                                   rename_var[c(1:(length(rename_var)/2))*2])] %>% 
  tidyr::pivot_longer(cols = rename_var[c(1:(length(rename_var)/2))*2], 
                      names_to = c("Seq"),
                      values_to = c("Tumour_Size"), 
                      names_prefix = "Tumour_Size_",
                      values_drop_na = TRUE)

#merge the tumour size and date together
masterfile_measure_long <- merge(masterfile_measure_long1,
                                 masterfile_measure_long2, 
                                 by = c("row_number", "Seq"), 
                                 all = TRUE)
# masterfile_measure_long$impute = "Original"

masterfile_measure_long$Seq = as.numeric(masterfile_measure_long$Seq)


masterfile_measure_long <-masterfile_measure_long[with(masterfile_measure_long,
                                                       order(row_number, Seq)),] 

masterfile_measure_long <- merge(masterfile_measure_long, master_anno[,c("row_number","Tumor_ID",
                                                                         "Transplant_Date",
                                                                         "Transplant_Date_Date", 
                                                                         "PTR_Date")],
                                 by = "row_number")



masterfile_measure_long$Tumour_Size_num = as.numeric(masterfile_measure_long$Tumour_Size)


masterfile_measure_long <- date_conversion( masterfile_measure_long, varname = "Tumour_Date", 
                                            replaceVarname = "Tumour_Date_Date")

#check for NA or outlier date and hardcode to correct them
temptemp <- unique(masterfile_measure_long[, c("row_number", "Seq", "Tumour_Date", "Tumour_Date_Date")])

masterfile_measure_long$Tumour_Date_Date <- as.Date(ifelse(masterfile_measure_long$Tumour_Date_Date == as.Date('2011-11-01'), 
                                                           as.Date('2021-11-01'), masterfile_measure_long$Tumour_Date_Date), 
                                                    origin = "1970-01-01")

masterfile_measure_long$Tumour_Date_Date <- as.Date(ifelse(masterfile_measure_long$Tumour_Date == '29/04/2019', 
                                                           as.Date('29/04/2019', tryFormats = c("%d/%m/%Y")),
                                                           masterfile_measure_long$Tumour_Date_Date), 
                                                    origin = "1970-01-01")


#quick check to see which tumor size could not be changed to numerical values
temptemp = masterfile_measure_long[is.na(masterfile_measure_long$Tumour_Size_num), 
                                   c('row_number', "Seq", "Tumour_Size", "Tumour_Size_num")]


#quick check to see records with missing measurement date or tumor size
temptemp = masterfile_measure_long[is.na(masterfile_measure_long$Tumour_Date_Date)|
                                     is.na(masterfile_measure_long$Tumour_Size_num),]


masterfile_measure_long$curve = ifelse(!is.na(masterfile_measure_long$Tumour_Date_Date) &
                                         !is.na(masterfile_measure_long$PTR_Date) &
                                         masterfile_measure_long$Tumour_Date_Date <= masterfile_measure_long$PTR_Date,
                                       "before PTR", 
                                       ifelse(!is.na(masterfile_measure_long$Tumour_Date_Date) &
                                                !is.na(masterfile_measure_long$PTR_Date) &
                                                masterfile_measure_long$Tumour_Date_Date > masterfile_measure_long$PTR_Date,
                                              "after PTR", NA))

#group by curve (before/after PTR, row_number(each row_number uniquely represents each tumour ))
masterfile_measure_long <- masterfile_measure_long %>%
  group_by(row_number) %>%
  mutate(Tumour_Size_num_lag1 = lag(Tumour_Size_num, n=1), 
         Tumour_Size_num_lead1 = lead(Tumour_Size_num, n=1),
         Tumour_Date_Date_lag1 = lag(Tumour_Date_Date, n=1), 
         Seq_lag1 = lag(Seq, n=1),
         same_curve = lag(curve, n= 1) == curve)

# check the right sequence and tumour size different between two measurements
masterfile_measure_long$Seq_diff = masterfile_measure_long$Seq - masterfile_measure_long$Seq_lag1
masterfile_measure_long$Tumour_Size_diff = masterfile_measure_long$Tumour_Size_num - masterfile_measure_long$Tumour_Size_num_lag1
masterfile_measure_long$Tumour_Date_diff = masterfile_measure_long$Tumour_Date_Date - masterfile_measure_long$Tumour_Date_Date_lag1

masterfile_measure_long = merge(masterfile_measure_long,
                                tumour_id_df[, c("Tumor_ID","category")],
                                by = "Tumor_ID")

#check tumour size difference between measurements. 
temptemp <- masterfile_measure_long[!is.na(masterfile_measure_long$Seq_diff) & masterfile_measure_long$Seq_diff !=1,
                                    c("row_number", "Seq", "Seq_diff", "Tumour_Date", 
                                      "Tumour_Size_diff", "Tumour_Date_diff","category")]

# 
# masterfile_measure_long <- masterfile_measure_long %>%
#   group_by(row_number) %>%
#   mutate(Tumour_Size_num_adj_lag1 = lag(Tumour_Size_num_adj, n=1))
# 
# masterfile_measure_long$Tumour_Size_diff_adj = masterfile_measure_long$Tumour_Size_num_adj - 
#   masterfile_measure_long$Tumour_Size_num_adj_lag1

temp_long = merge(masterfile_measure_long[, !names(masterfile_measure_long) %in% 
                                            c("Tumor_ID", "Transplant_Date_Date")], 
                  master_anno, by = c("row_number"))

temp_long$date_len = as.numeric(difftime(temp_long$Tumour_Date_Date,
                                         temp_long$Transplant_Date_Date, units = "days"))


#temp_long <- merge(temp_long, tumour_id_df, by = "Tumor_ID")

temp_long$growth_rate <- temp_long$Tumour_Size_diff/as.numeric(temp_long$Tumour_Date_diff) 


temp_long <- merge(temp_long, met_cat[, c("Tumor_ID", "Damian","mainsite")],
                   by = "Tumor_ID", all = TRUE)

names(temp_long)[grepl("mainsite", names(temp_long))] = "mainsite_damian"



firstcurve_df <- temp_long[temp_long$curve == "before PTR",]


used_df <- subset(temp_long, category=="Main" & curve == "before PTR" )



# get rid of NA values. With NA values, the predicted values have different length
used_df <- subset(used_df, !is.na(Tumour_Size_num) & !is.na(used_df$mainsite))


used_df$mainsite <- factor(used_df$mainsite, levels = c("non_met", "met"))

used_df$color_SAid <- factor(ifelse(used_df$SA_ID == "SA919", "SA919", 
                                    "Other"), levels = c("SA919", 
                                                         "Other"))


getMatch = function(rexp, str) regmatches(str, regexpr(rexp, str))

getMatch("SA\\d+", used_df$SA_ID)
used_df$SA_ID_new = getMatch("SA\\d+", used_df$SA_ID)
unique(used_df[, c("SA_ID_new", "SA_ID")])

temp = unique(used_df[, c("mainsite_damian", "SA_ID_new")])

id_label <- data.frame(unique(temp[with(temp, order(mainsite_damian, SA_ID_new)),])%>%
                         group_by(mainsite_damian) %>%
                         summarise(color_SAid_new = paste0(SA_ID_new, collapse = ",")) %>%
                         ungroup())

id_label$color_SAid_new = ifelse(id_label$color_SAid_new == "SA1139,SA1146,SA501,SA535,SA604,SA605,SA609,SA919",
                                 "Other \n(SA1139,SA1146,\nSA501,SA535,\nSA604,SA605,\nSA609)",
                                 "Other \n(SA1142,SA1146,\nSA501,SA535,\nSA605,SA609)" )

used_df = merge(used_df, id_label, by = "mainsite_damian")

#used_df <- used_df[ , -which(names(used_df) %in% c( "color_SAid_new"))]

unique(used_df[, c("mainsite_damian", "color_SAid_new", "color_SAid", "SA_ID_new")])

used_df$color_SAid_new = ifelse(used_df$color_SAid == "SA919", "SA919", 
                                used_df$color_SAid_new)

unique(used_df[, c("mainsite_damian", "color_SAid_new", "color_SAid", "SA_ID_new")])

temp_data = used_df[, c("mainsite_damian", "Tumor_ID", "Tumour_Date", "SA_ID", "SA_ID_new",
                        "PTR_Date.x", "date_len", "Tumour_Size_num", "Tumour_Size", "curve",
                        "color_SAid", "color_SAid_new")]


names(temp_data)[names(temp_data) == "PTR_Date.x"] = "PTR_Date"


#Note the temp_data used in that file is actually the output of the above data, 
# it should have the same number of row, and with extra row number for the imported data set

source(paste0(sourcedir, "lmer_subscription.R"))