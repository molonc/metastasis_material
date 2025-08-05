
df <- tibble::tibble(val=c(3, 4, 3, 3, 3, 4), chr=c('5', '5','7', '7', '10', '10'), 
                     clone=c('B', 'C','B', 'C', 'B', 'C'))       


viz_cnv_profiles_specific_chrs <- function(df){
  chrs <- unique(df$chr)
  chr_levels <- c(as.character(1:23), "X", "Y")
  chr_levels <- chr_levels[chr_levels %in% chrs]
  df$chr <- factor(df$chr, levels = chr_levels)

  cnv_cols <- c('0'='#4880B8', '1'='#A7C9DF','2'='#CCCCCC','3'='#F5CE93','4'='#ED9364',
                '5'='#D2553E','6'='#A42116','7'='#8B1A43','8'='#CB3576','9'='#D06CAD',
                '10'='#C196C4','11'='#D0BAD8')
  vals <- sort(as.numeric(unique(df$cnv)), decreasing=F)
  df$cnv <- as.character(df$cnv)
  df$cnv <- factor(df$cnv, levels = as.character(vals))
  my_font <- 'Helvetica'
  p <- ggplot(df, aes(x=chr, y=clone, fill= cnv)) + 
    geom_tile(colour = "white", size=1.2) + 
    # facet_wrap(Y~., scales = "free") + 
    theme_bw() + 
    scale_fill_manual(values = cnv_cols) + 
    theme(axis.line=element_blank(),
          axis.title = element_text(size=11, family = my_font, colour = 'black'),
          # axis.title.x = element_blank(), 
          # axis.title.y = element_text(size=11), 
          # axis.text.y = element_text(size=11),
          axis.ticks.y = element_blank(),
          axis.text.y = element_text(size=13, family = my_font, colour = 'black'),
          axis.text.x = element_text(size=13, family = my_font, colour = 'black'),
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          panel.grid = element_blank(),
          panel.border = element_blank())+ 
    labs(x='Chromosome',y='Clone',title='') 
  # p
  return(p)
}

p <- viz_cnv_profiles_specific_chrs(stat)
p




stat
df <- stat
