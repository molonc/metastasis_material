library(ggplot2)

p <- ggplot(diamonds, aes(x = carat, y = -0.5)) +
  
  # horizontal boxplots & density plots
  geom_boxplot(aes(fill = cut)) +
  geom_density(aes(x = carat, fill=cut), inherit.aes = F) +
  
  # vertical lines at Q1 / Q2 / Q3
  # stat_boxplot(geom = "vline", aes(xintercept = ..xlower..)) +
  # stat_boxplot(geom = "vline", aes(xintercept = ..xmiddle..)) +
  # stat_boxplot(geom = "vline", aes(xintercept = ..xupper..)) +
  
  # facet_grid(cut ~ .) +
  scale_fill_discrete()
p




# Alternative 1: calculate the box plot's coordinates, & flip them manually before passing the results to ggplot(). Add a density layer in the normal way:


library(dplyr)
library(tidyr)

p.box <- ggplot(diamonds, aes(x = cut, y = carat)) + geom_boxplot()    
p.box.data <- layer_data(p.box) %>%
  select(x, ymin, lower, middle, upper, ymax, outliers) %>%
  mutate(cut = factor(x, labels = levels(diamonds$cut), ordered = TRUE)) %>%
  select(-x)

ggplot(p.box.data) +
  
  # manually plot flipped boxplot
  geom_segment(aes(x = ymin, xend = ymax, y = -0.5, yend = -0.5)) +
  geom_rect(aes(xmin = lower, xmax = upper, ymin = -0.75, ymax = -0.25, fill = cut),
            color = "black") +
  geom_point(data = . %>% unnest(outliers),
             aes(x = outliers, y = -0.5)) +
  
  # vertical lines at Q1 / Q2 / Q3
  geom_vline(data = . %>% select(cut, lower, middle, upper) %>% gather(key, value, -cut),
             aes(xintercept = value)) +
  
  # density plot
  geom_density(data = diamonds, aes(x = carat)) +
  
  facet_grid(cut ~ .) +
  labs(x = "carat") +
  scale_fill_discrete()


# Alternative 2: calculate the density plot's coordinates, 
# & flip them manually before passing the results to ggplot(). 
# Add a box plot layer in the normal way. Flip the whole chart:
p.density <- ggplot(diamonds, aes(x = carat, group = cut)) + geom_density()    
p.density.data <- layer_data(p.density) %>%
  select(x, y, group) %>%
  mutate(cut = factor(group, labels = levels(diamonds$cut), ordered = TRUE)) %>%
  select(-group)
p.density.data <- p.density.data %>%
  rbind(p.density.data %>% 
          group_by(cut) %>% 
          filter(x == min(x)) %>% 
          mutate(y = 0) %>% 
          ungroup())

ggplot(diamonds, aes(x = -0.5, y = carat)) +
  
  # manually flipped density plot
  geom_polygon(data = p.density.data, aes(x = y, y = x, fill = cut), 
               alpha=0.3, color = "black") + #fill = NA, 
  
  # box plot
  # geom_violin(aes(fill = cut), alpha=0.3) + #, group = cut #trim=T
  geom_boxplot(aes(fill = cut), width=0.5, alpha=0.3, outlier.alpha = 0) + 

  # vertical lines at Q1 / Q2 / Q3
  # stat_boxplot(geom = "hline", aes(yintercept = ..lower..)) +
  # stat_boxplot(geom = "hline", aes(yintercept = ..middle..)) +
  # stat_boxplot(geom = "hline", aes(yintercept = ..upper..)) +
  
  # facet_grid(cut ~ .) +
  scale_fill_discrete() +
  coord_flip()








viz_genewise_dispersion <- function(obs_clones=c('B','C')){
  ptittle <- paste0(obs_clones[2],' vs. ', obs_clones[1])
  obs_clones <- paste0('Clone ', obs_clones)
  # color_use <- c('Clone A'='#66C2A5','Clone B'='#FC8D62','Clone C'='#8DA0CB') 
  
  cols_use <- c('cis'='#EE220C','trans'='#0076BA')
  print(cols_use)
  
  df <- diamonds
  unique(df$cut)
  df <- df %>%
    dplyr::filter(cut %in% c('Ideal','Premium')) %>%
    dplyr::mutate(cut=
                    case_when(
                      cut=='Ideal' ~ 'cis',
                      TRUE ~ 'trans'
                    )) %>%
    dplyr::rename(gene_exp=carat)
  df$gene_exp <- df$gene_exp/2
  colnames(df)
  df$cut <- as.factor(df$cut)
  
  p.density <- ggplot(df, aes(x = gene_exp, group = cut)) + geom_density()    
  p.density.data <- layer_data(p.density) %>%
    select(x, y, group) %>%
    mutate(cut = factor(group, labels = levels(df$cut), ordered = TRUE)) %>%
    select(-group)
  p.density.data <- p.density.data %>%
    rbind(p.density.data %>% 
            group_by(cut) %>% 
            filter(x == min(x)) %>% 
            mutate(y = 0) %>% 
            ungroup())
  
  p1 <- ggplot(df, aes(x = -0.5, y = gene_exp)) +
    
    # manually flipped density plot
    geom_polygon(data = p.density.data, aes(x = y, y = x, fill = cut), 
                 alpha=0.3, color = "black") + #fill = NA, 
    
    # box plot
    # geom_violin(aes(fill = cut), alpha=0.3) + #, group = cut #trim=T
    geom_boxplot(aes(fill = cut), width=0.5, alpha=0.3, outlier.alpha = 0) + 
    
    # vertical lines at Q1 / Q2 / Q3
    # stat_boxplot(geom = "hline", aes(yintercept = ..lower..)) +
    # stat_boxplot(geom = "hline", aes(yintercept = ..middle..)) +
    # stat_boxplot(geom = "hline", aes(yintercept = ..upper..)) +
    
    # facet_grid(cut ~ .) +
    scale_fill_manual(values=cols_use) +
    coord_flip() + 
    theme_bw() + 
    theme(#legend.position = 'none',
          panel.grid = element_blank()) + 
    labs(x='', y='', title = ptittle)
  return(p1)
  
}

p1 <- viz_genewise_dispersion(obs_clones=c('A','B'))
p1
p2 <- viz_genewise_dispersion(obs_clones=c('A','C'))
p3 <- viz_genewise_dispersion(obs_clones=c('B','C'))

p_total <- cowplot::plot_grid(p1, p2, p3, nrow = 1)
p_total


viz_violin_plt <- function(obs_clones){
  ptittle <- paste0(obs_clones[2],' vs. ', obs_clones[1])
  # Change colors
  
  
  ToothGrowth$dose <- as.factor(ToothGrowth$dose)
  p <- ggplot(ToothGrowth, aes(x=dose, y=len, color=supp)) +
    geom_violin(position=position_dodge(0.6), alpha=0)+
    geom_jitter(position=position_jitterdodge(0.4)) +
    theme_bw() + 
    theme(#legend.position = 'none',
      panel.grid = element_blank()) + 
    labs(x='', y='', title = ptittle)
    # geom_boxplot(width=1)
  p
  # p+scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  
  return(p)
}

save_dir <- '/Users/hoatran/Documents/BCCRC_projects/hakwoo_project/code/results_bulkRNAseq/SA919_full/figs/'
png(paste0(save_dir,"demo_dispersion.png"), height = 2*150, width=2*600, res = 2*72)
print(p_total)
dev.off()  

