# Generate transplant trees diagrams
# 
# Uses transplant and patient data exported from MySQL, and infers parent-daughter relationships from sample names
#
# inputs:  	transplants.csv		transplant details, exported from MySQL table `transplants`
# 	  	patient.csv		patient details, exported from MySQL table `patients`
#
# directories:	inputdir		local directory containing transplants.csv and patient.csv
#		outputdir		local directory for PDFs created
#		 

#install packages if not already done. 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("ggtree")

# Load R packages
library(data.tree)
library(ggtree)
library(rsvg)
library(DiagrammeRsvg)
library(DiagrammeR)

options(stringsAsFactors = FALSE)

# Directories for input and output files
maindir <- path.expand("~/Documents/Aparicio/Tree_Plot/SC-1068/")
inputdir <- paste0(maindir, "Data/")
outputdir <- paste0(maindir, "Results/")
if(!dir.exists(outputdir)) dir.create(outputdir)

#data downloaded from https://atlas.molonc.ca/ Transplants under experiments
transplant_file = "transplants_data_2022_08_03.csv"

# Import transplant data
tsp<-read.csv(paste0(inputdir, transplant_file), header=T, stringsAsFactors=F, na.strings="NULL", colClasses=c("character"))
names(tsp)[names(tsp) == "Transplant.Name"] = 'Transplant_ID'

# series.grown <- unique(tsp$Patient_ID[tsp$Grown>0])
tsp$Drug.Study[is.na(tsp$Drug.Study)] <- ""
tsp$Treated.Control[is.na(tsp$Treated.Control)] <- ""


## Excluding AtiM.Patient.ID empty string.
tsp= subset(tsp, !AtiM.Patient.ID == "")


# Import patient data
#data downloaded from https://atlas.molonc.ca/ Patient Identifiers under Patients & Others
# this data is required to get SAxxx, like SA909 etc.
pat<-read.csv(paste0(inputdir, "patient_identifiers_data_2022_08_04.csv"),header=T,stringsAsFactors=F,na.strings="NULL",colClasses=c("character"))

tsp = merge(tsp, pat[, c("SA.Identifiers", "AtiM.Patient.ID")], by = c("AtiM.Patient.ID"), 
            all.x = TRUE)

# Import filtering data
# Data download from https://airtable.com/shrkb9pHRQI1DwEJb/tbloHKCIrvRJdrWDz
# this dataset is used to filter out all data that is NOT part of Hakwoo 
data_to_select <-read.csv(paste0(inputdir, "Grid view.csv"), header=T, 
                          stringsAsFactors=F,na.strings="NULL",colClasses=c("character"))

temp_data = unique(data_to_select[, names(data_to_select)[names(data_to_select) %in% 
                                    c("ATiM.Patient.ID", "Transplant.Name")]])

names(temp_data) = c("AtiM.Patient.ID", "Transplant_ID")

# quick check to see Transplant_ID is unique.
temp = table(temp_data$Transplant_ID)

# tsp = merge(tsp, pat[, c("SA.Identifiers", "AtiM.Patient.ID")], by = c("AtiM.Patient.ID"), 
#             all.x = TRUE)

# keep the transplant_ID in Hakwoo's data only.
SA_selected = merge(tsp, temp_data, by = c("AtiM.Patient.ID", "Transplant_ID"), 
            all.y = TRUE)

#Sereislist is the vector of unique SA ids. 
# for this SA_selected data set, there is only one. 
serieslist = unique(SA_selected$SA.Identifiers)


# to run tree_plot function, variable "Transplant_ID", "AtiM.Patient.ID", 
# "Passage", "Transplant.Material" are required. The two variables are used to define 
# colors in tree plot. 
for(series in serieslist){
  tree_plot(df = SA_selected, series = "SA919", source = "airtable", resultdir = outputdir)
}
