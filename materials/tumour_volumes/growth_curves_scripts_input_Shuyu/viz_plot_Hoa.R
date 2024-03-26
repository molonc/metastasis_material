

input_dir <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/tumour_volumes/growth_curves_scripts_input_Shuyu/'
pts_df <- data.table::fread(paste0(input_dir, 'mapping_SAID_PatientID_color.csv'))

temp_data <- temp_data %>% 
  dplyr::left_join(pts_df, by='SA_ID')


color_sa <- pts_df$color_code
names(color_sa) <- pts_df$patient_ID


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
dev.off()


ggsave(paste0(input_dir,"tumor_main_cubic_helvetica.svg"),  
       plot = p,  
       height = 3.3,  
       width = 7.5,  
       # useDingbats=F,  
       dpi=200  
)

ggsave(paste0(input_dir,"tumor_main_cubic_helvetica.png"),  
       plot = p,  
       height = 3.3,  
       width = 7.5,  
       # useDingbats=F,  
       # type = "cairo-png",
       dpi=200  
)
