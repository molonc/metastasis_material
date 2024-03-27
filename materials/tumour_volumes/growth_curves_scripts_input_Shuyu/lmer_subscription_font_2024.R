library(lmerTest)
library(lme4)
library(rstatix)
library(ggplot2)
library(extrafont)




temp_data <- subset(firstcurve_df, experiment_type == "main_experiment")
patient_id <- read.csv(paste0(datadir, "March2024Data/",'mapping_SAID_PatientID_color.csv'))

#number of unique tumors with main experiment
length(unique(temp_data$Tumor_ID))

#number of unique tumors with main experiment
temp_data <- subset(temp_data,
  Damian_annotation %in% c("non_met", "met"))

temp_data <- merge(temp_data, patient_id, by = "SA_ID")

  
#number of unique tumors with main experiment AND with non_met or met annotation.
## excluding unknown category
  length(unique(temp_data$Tumor_ID))
  
  
temp_data$Damian_annotation <- factor(temp_data$Damian_annotation, 
                                      levels = c("non_met",  "met"))

fit2 <- lmer(TumorVol_num^(1/3) ~ date_len*Damian_annotation+
               (date_len|Tumor_ID), 
             data = temp_data,  REML = FALSE, 
             control = lmerControl(
               optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

summary2_stat = summary(fit2)
round(summary2_stat$coefficients, 2)


temp_data$predm2 <- predict(fit2,re.form=NA)  ## population level
temp_data$predm2fit <- predict(fit2) ## individual level


dim(temp_data[!is.na(temp_data$TumorDate_Date) & !is.na(temp_data$TumorVol_num),])
dim(temp_data)

library(RColorBrewer)
length(unique(temp_data$SA_ID))

# #choose the color
# color_sa = brewer.pal(length(unique(temp_data$SA_ID))+1, "BrBG")
# 
# #look at the color
# plot(x = c(1:length(unique(temp_data$SA_ID))), 
#      y=c(1:length(unique(temp_data$SA_ID))), 
#      col = color_sa, pch= 16)
# #remove the color too faint
# color_sa = color_sa[-5]
# 
# #check again
# plot(x = c(1:length(unique(temp_data$SA_ID))), 
#      y=c(1:length(unique(temp_data$SA_ID))), 
#      col = color_sa, pch= 16)

#get the max time 
max_length = max(temp_data$date_len, na.rm = TRUE)

# get the max Tumor volume
max_size = max(temp_data$TumorVol_num)


color_sa <- patient_id$color_code
names(color_sa) <- patient_id$patient_ID



my_font <- "Helvetica"
lg_pos <- "top"  
thesis_theme <- ggplot2::theme(  
  text = element_text(color="black",size = 8, hjust = 0.5, family=my_font),  
  axis.title.x = element_text(color="black",size=10, hjust = 0.5, family=my_font),  
  axis.title.y = element_text(color="black",size=10, hjust = 0.5, family=my_font),  
  axis.text.x = element_text(color="black",size=7, hjust = 0.5, family=my_font),#, angle = 90
  # axis.text.x = element_blank(),  
  axis.text.y = element_text(color="black",size=7, hjust = 0.5, family=my_font),  
  plot.title = element_text(color="black",size=10, face="bold", hjust=0, family=my_font),  
  legend.title=element_text(color="black",size=7, hjust = 0.5, family=my_font),  
  legend.text=element_text(color="black",size=7, hjust = 0.5, family=my_font),  
  strip.text.x = element_text(color="black",size=10, family=my_font),  
  strip.text.y = element_text(color="black",size=10, family=my_font),  
  legend.spacing.x = unit(0.1, 'mm'),  
  legend.spacing.y = unit(0.1, 'mm'),  
  legend.key.height=unit(0.5,"line"),  
  legend.position = lg_pos,  
  legend.background = element_blank(),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank()  
)

p <- ggplot(data = temp_data, aes(x = date_len, 
                                  y = TumorVol_num^(1/3),
                                  group=Tumor_ID)) +
  geom_point(size = 1.5, alpha=.3) +
  geom_line(linetype = 2, size = 1, alpha=.3) + 
  geom_line(aes(y=predm2fit,group=Tumor_ID, colour = patient_ID), size =1.2,  alpha=.6)+ # note: color by patient id
  scale_x_continuous(limits = c(0, max_length), name = "Days") +
  scale_y_continuous(limits = c(0, max_size^(1/3)), 
                     name = expression(paste(sqrt("Volume", 3)))) +#, trans = "log")#+
  
  labs(colour="Patient ID")+ 
  scale_color_manual(values=color_sa)+
  facet_wrap(vars(Damian_annotation), dir = "h")+ ## note: using horizontal order
  thesis_theme + 
  guides(colour = guide_legend(nrow = 1)) # legends in one row

# theme(
#   axis.text.x = element_text(family = "Helvetica", size = 8),
#   axis.text.y = element_text(family = "Helvetica", size = 8),
#   axis.title.x = element_text(family = "Helvetica", size = 10),
#   axis.title.y = element_text(family = "Helvetica", size = 10),
#   strip.text.x = element_text(family = "Helvetica", size = 12),
#   legend.title=element_text(family = "Helvetica", size = 12),
#   legend.text = element_text(family = "Helvetica", size = 10))

p <- p + geom_line(colour="red",aes(y=predm2, x=date_len), 
                   linewidth = 1.7, alpha=.6,
                   inherit.aes = F)
print(p)


#
ggsave(filename=paste0(resultdir, "tumor_main_cubic_helvetica.png"), 
       height = 3.5,
       width = 4,units = "in")
print(p)
dev.off()

## the device is used to avoid 
# "In grid.Call.graphics(C_points, x$x, x$y, x$pch, x$size) :
#   semi-transparency is not supported on this device: reported only once per page"

ggsave(paste0(resultdir,"tumor_main_cubic_helvetica.svg"),  
       plot = p,  
       height = 3.3,  
       width = 7.5,  
       # useDingbats=F,  
       dpi=200  
)


fit3 <- lmer(TumorVol_num^(1/3) ~ date_len*Damian_annotation+
               (date_len|Tumor_ID), 
             data = subset(temp_data,SA_ID == "SA919"),  REML = FALSE, 
             control = lmerControl(
               optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

summary3_stat = summary(fit3)
round(summary3_stat$coefficients, 3)

SuppTable3 <- temp_data[, c("Tumor_ID","patient_ID", "Damian_annotation",
                                  "PTR_Date", "TumorDate_Date",
                                  "TumorVol_num", "date_len", "color_code")]

write.csv(SuppTable3, paste0(resultdir, "SuppTable3.csv"))
