# install.packages("readxl")
# install.packages("survminer")

library(readxl)
library(survival)
library(ggplot2)
library(survminer)

rootdir = "//Users//shuyufan//Documents//Aparicio//Metastasis"

filenames <- list.files(file.path(rootdir), pattern="Ki67.xlsx",
                        full.names=TRUE)

my_data <- read_excel(filenames, sheet = "Sheet1")

names(my_data) = c("id", "days", "Ki67_high", "Ki67_low")

my_data= my_data[!is.na(my_data$id),]

my_data_hi = my_data[grepl("h", my_data$id), ]
my_data_hi['characteristics'] = "Ki67_high"
my_data_hi = my_data_hi[, c(1,2,3,5)]
names(my_data_hi) = c("id", "days", "status", "characteristics")
                      
my_data_low = my_data[grepl("l", my_data$id), ]
my_data_low['characteristics'] = "Ki67_low"
my_data_low = my_data_low[, c(1,2,4,5)]
names(my_data_low) = c("id", "days", "status", "characteristics")

my_data_long = rbind(my_data_hi, my_data_low)

fit1 = survfit(Surv(days, status) ~characteristics, my_data_long)

survplot_1 = ggsurvplot(
  fit1, 
  data = my_data_long, 
  size = 1,                 # change line size
  palette = 
    c("#E7B800", "#2E9FDF"),# custom color palettes
  conf.int = TRUE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.labs = 
    unique(my_data_long$characteristics),    # Change legend labels
  risk.table.height = 0.3, # Useful to change when you have multiple groups
  ggtheme = theme_bw(),     # Change ggplot2 theme
  # surv.median.line = "hv",  # add the median survival pointer.
  legend.title="",
  pval.method = TRUE
  
)

pdf(paste0(rootdir,"Metastasis_Surv_Ki.pdf"))
print(survplot_1)
dev.off()
# 
# png(paste0(rootdir,"Metastasis_Surv_Ki.png"))
# print(survplot_1)
# dev.off()

coxmodel1 = coxph(Surv(days, status) ~ characteristics, data = my_data_long)
summary(coxmodel1)