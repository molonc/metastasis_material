


# From copy number values, get median values
# Final input like txt files
df <- data.frame(v1=c('a','a','a','b','b','b'),
                 v2=c(1,1,3,2,4,4), stringsAsFactors = F)
df <- df %>% 
  dplyr::select(v1, v2) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(v1) %>% 
  summarize(v2=median(v2), .groups = 'drop')

View(df)

input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_rna/clonealign/SA1035X4XB02879/'

cnv <- read.table(paste0(input_dir,'SA1035X4XB02879.txt'), sep='\t', header=T)
View(head(cnv))
df_cnv <- cnv
df_cnv <- inner_join(df_cnv, annots)

df_cnv <- dplyr::select(df_cnv, ensembl_gene_id, chr, start) %>% 
  group_by(chr) %>% 
  dplyr::mutate(start_order = rank(start)) %>% 
  ungroup() %>% 
  inner_join(df_cnv)

df_cnv$cluster[1:3]
View(head(df_cnv))
segment_df <- readRDS(paste0(input_dir,'SA1035X4XB02879.rds'))
class(segment_df)
View(head(segment_df))
segment_df$chr_desc <- paste0(segment_df$chr,'_',segment_df$start,'_',segment_df$end)
dim(segment_df)
colnames(segment_df)
segment_df$copy_number

segment_df3 <- segment_df %>% 
  dplyr::select(chr, start, end, copy_number) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(chr, start, end) %>% 
  summarize(copy_number=median(copy_number), .groups = 'drop')


segment_df3 <- segment_df %>% 
  dplyr::select(chr, start, end, copy_number) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(chr, start, end) %>% 
  summarize(copy_number=getmode(copy_number), .groups = 'drop')

summary(segment_df3$copy_number)
segment_df3$copy_number[1:3]
# Create the function.
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Create the vector with numbers.
v <- c(2,1,2,3,1,2,3,4,1,5,5,3,2,3)

# Calculate the mode using the user function.
result <- getmode(v)
print(result)


summary(segment_df$copy_number)
is.integer(segment_df$copy_number[1:4])
segment_df$copy_number[1:20]

t <- do.call(is.integer,as.list(segment_df$copy_number))

dim(segment_df2)
segment_df2 <- segment_df3
View(head(segment_df2))
length(unique(segment_df2$chr_desc))
head(df_cnv)

segment_df2 <- dplyr::select(segment_df2, chr, start) %>% 
  group_by(chr) %>% 
  dplyr::mutate(start_order = rank(start)) %>% 
  ungroup() %>% 
  inner_join(segment_df2)

segment_df2$chr <- factor(segment_df2$chr, levels = chr_levels)
summary(segment_df2$copy_number)
segment_df2$copy_number <- as.integer(segment_df2$copy_number)
segment_df2$copy_number <- as.factor(segment_df2$copy_number)

segment_df2 <- drop_na(segment_df2)
levels(segment_df2$copy_number) <- 0:(length(cnv_colors)-1)
segment_df2$cluster <- 'A'

cnv_plot <- dplyr::filter(segment_df2, !chr %in% c("X","Y")) %>% 
  ggplot(aes(x = start_order, y = cluster, fill = copy_number)) +
  geom_raster() +
  facet_wrap(~ chr, scales = "free_x", nrow = 1, switch = "x") +
  #theme(legend.position = "top", axis.text.x = element_blank()) +
  scale_y_discrete(expand = c(0, 0)) +
  #scale_fill_manual(values=cnv_colors, name = "Copy number", guide = 'legend',labels = 0:(length(cnv_colors)-1),drop=FALSE)  +
  #          theme(legend.position = "bottom") +   
  scale_fill_manual(values = cnv_cols, name = "Copy number") +  #, labels = 0:(length(cnv_cols)-1),drop=FALSE) +
  labs(x = "chromosome") +
  theme(strip.background = element_rect(fill = 'white', colour = 'white'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=14), 
        strip.placement = "outside",
        legend.position = "none",
        panel.spacing = unit(c(0.1), 'cm')) +   # MA: was 0.2
  theme(text = element_text(size = 18)) +
  labs(x = "Chromosome", y = "Clone")+
  guides(fill = guide_legend(nrow = 1)) +
  scale_x_continuous(expand = c(0,0))

cnv_plot
main_plot <- cowplot::plot_grid(
  cnv_plot,
  cnv_plot,
  ncol = 1,
  rel_heights = c(2.5, 2.5),
  align = 'v'
)
cowplot::ggdraw(cnv_plot)
library(ggplot2)
main_plot
print(cnv_plot)