library(lmerTest)
library(lme4)
library(rstatix)
library(ggplot2)
library(extrafont)


datadir = "C:/Documents/Metastasis/Result20240112/"
temp_data <- read.csv(file =  paste0(datadir, "firstcurve_data.csv"))


temp_data$color_SAid <- factor(ifelse(temp_data$SA_ID == "SA919", "SA919", 
                                    "Other"), levels = c("SA919", 
                                                         "Other"))

fit2 <- lmer(Tumour_Size_num^(1/3) ~ date_len*mainsite+
               (date_len|Tumor_ID), 
             data = temp_data,  REML = FALSE, 
             control = lmerControl(
               optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

summary2_stat = summary(fit2)
round(summary2_stat$coefficients, 2)


temp_data$predm2 <- predict(fit2,re.form=NA)  ## population level
temp_data$predm2fit <- predict(fit2) ## individual level




p <- ggplot(data = temp_data, aes(x = date_len, 
                                y = Tumour_Size_num^(1/3),
                                colour = color_SAid, 
                                group=Tumor_ID)) +
  
  geom_line(colour="grey",aes(y=predm2fit,group=Tumor_ID), size =1.5)+
  geom_point(size = 1.5) +geom_line(linetype = 2, size = 1.5) + #theme(legend.position = "none") +
  scale_x_continuous(limits = c(0, max_length), name = "Days") +
  scale_y_continuous(limits = c(0, max_size^(1/3)), 
                     name = expression(paste(sqrt("Volume", 3)))) +#, trans = "log")#+
  
  labs(colour="Patient ID")+ 
  
  scale_color_manual(values=c("darkorchid1", "lightskyblue4"))+
  
  facet_wrap(vars(mainsite), dir = "v")+
  theme(
    axis.text.x = element_text(family = "Helvetica", size = 12), 
    axis.text.y = element_text(family = "Helvetica", size = 12), 
    axis.title.x = element_text(family = "Helvetica", size = 16),
    axis.title.y = element_text(family = "Helvetica", size = 16),
    strip.text.x = element_text(family = "Helvetica", size = 18), 
    legend.title=element_text(family = "Helvetica", size = 16),
    legend.text = element_text(family = "Helvetica", size = 12)) 

p <- p + geom_line(colour="red",aes(y=predm2, x=date_len), linewidth = 2, inherit.aes = F) 
print(p)

png(filename=paste0(resultdir, "tumor_main_cubic_helvetica.png"),  
    width = 960, height = 960)

print(p)
dev.off()


setEPS()
postscript(paste0(resultdir, "tumor_main_cubic_helvetica.eps"))
print(p)
dev.off()



temp_data_non_met <- subset(temp_data, mainsite == "non_met")
unique(temp_data_non_met$color_SAid_new)
temp_data_non_met$color_SAid_new <- factor(temp_data_non_met$color_SAid_new, 
                                         levels = c("SA919", 
                                                    "Other \n(SA1139,SA1146,\nSA501,SA535,\nSA604,SA605,\nSA609)"))

p <- ggplot(data = temp_data_non_met,
            aes(x = date_len, 
                y = Tumour_Size_num^(1/3),
                colour = color_SAid_new, 
                group=Tumor_ID)) +
  
  geom_line(colour="grey",aes(y=predm2fit,group=Tumor_ID), size =1.5)+
  geom_point(size = 1.1) +geom_line(linetype = 1, size = 0.8) + #theme(legend.position = "none") +
  scale_x_continuous(limits = c(0, max_length), name = "Days") +
  scale_y_continuous(limits = c(0, max_size^(1/3)), 
                     name = expression(paste(sqrt("Volume", 3)))) +#, trans = "log")#+
  
  labs(colour="Patient ID")+ 
  scale_color_manual(values=c("darkorchid1", "lightskyblue4"))+
  ggtitle("Primary tumor growth curve from 8 PDX that didn't develope metastases")+
  theme(
    plot.title = element_text(hjust = 0.5, family = "Helvetica", size = 16),
    axis.text.x = element_text(family = "Helvetica", size = 12), 
    axis.text.y = element_text(family = "Helvetica", size = 12), 
    axis.title.x = element_text(family = "Helvetica", size = 14),
    axis.title.y = element_text(family = "Helvetica", size = 14),
    legend.title=element_text(family = "Helvetica", size = 14),
    legend.text = element_text(family = "Helvetica", size = 12),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "grey", fill=NA, size=1), 
    legend.position="bottom")


p <- p + geom_line(colour="red",aes(y=predm2, x=date_len), linewidth = 2, inherit.aes = F) 
print(p)

ggsave(filename=paste0(resultdir, "tumor_main_cubic_helvetica_non_met.png"), 
       height = 3.5,
       width = 4,units = "in")

ggsave(filename=paste0(resultdir, "tumor_main_cubic_helvetica_non_met.eps"), 
       height = 3.5,
       width = 4,units = "in")




temp_data_met <- subset(temp_data, mainsite == "met")
unique(temp_data_met$color_SAid_new)
temp_data_met$color_SAid_new <- factor(temp_data_met$color_SAid_new, 
                                     levels = c("SA919", 
                                                "Other \n(SA1142,SA1146,\nSA501,SA535,\nSA605,SA609)"))

p <- ggplot(data = temp_data_met,
            aes(x = date_len, 
                y = Tumour_Size_num^(1/3),
                colour = color_SAid_new, 
                group=Tumor_ID)) +
  
  geom_line(colour="grey",aes(y=predm2fit,group=Tumor_ID), size =1.5)+
  geom_point(size = 1.1) +geom_line(linetype = 1, size = 0.8) + #theme(legend.position = "none") +
  scale_x_continuous(limits = c(0, max_length), name = "Days") +
  scale_y_continuous(limits = c(0, max_size^(1/3)), 
                     name = expression(paste(sqrt("Volume", 3)))) +#, trans = "log")#+
  
  labs(colour="Patient ID")+ 
  scale_color_manual(values=c("darkorchid1", "lightskyblue4"))+
  ggtitle("Primary tumor growth curve from 7 PDX that developed metastases")+
  theme(
    plot.title = element_text(hjust = 0.5, family = "Helvetica", size = 16),
    axis.text.x = element_text(family = "Helvetica", size = 12), 
    axis.text.y = element_text(family = "Helvetica", size = 12), 
    axis.title.x = element_text(family = "Helvetica", size = 14),
    axis.title.y = element_text(family = "Helvetica", size = 14),
    legend.title=element_text(family = "Helvetica", size = 14),
    legend.text = element_text(family = "Helvetica", size = 12),
    panel.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(colour = "grey", fill=NA, size=1), 
    legend.position="bottom")


p <- p + geom_line(colour="red",aes(y=predm2, x=date_len), linewidth = 2, inherit.aes = F) 
print(p)

ggsave(filename=paste0(resultdir, "tumor_main_cubic_helvetica_met.png"), 
       height = 3.5,
       width = 4,units = "in")

ggsave(filename=paste0(resultdir, "tumor_main_cubic_helvetica_met.eps"), 
       height = 3.5,
       width = 4,units = "in")