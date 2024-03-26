
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



maindir <- "C:/Documents/Metastasis_2024/"
sourcedir <- paste0(maindir, "R_files/")
datadir <- paste0(maindir, "Data/")
today <- format(Sys.Date(), format = "%Y%m%d")
resultdir <- paste0(maindir, "Result", today, "/")

source(paste0(sourcedir, "helper_function.R"))
if(!dir.exists(resultdir)) dir.create(resultdir, recursive = TRUE)


master_filename <- "MASTER_FILE_Metastasis_mouse_exp_HL_Damian_updated_14_March_2024.csv"


master_df<- read_csv(paste0(datadir, "March2024Data/",master_filename))

## run the source file to clean up the data 
## note, this is not a file with defined functions. It is to help clarify the stream line of 
## data claning. 
## includesource(paste0(sourcedir, "master_anno_preProcess_2024.R"))
## to handle cleaning for annotation file 
source(paste0(sourcedir, "masterFile_preProcess_2024.R"))

  