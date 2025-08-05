library(ggplot2)


viz_graph_results <- function(){
  library(ggraph)
  library(igraph)
  library(tidygraph)
  
  input_dir <- "/Users/hoatran/Documents/projects_BCCRC/hakwoo_project/code/metastasis_material/materials/bulkRNAseq/preprocessed_09April2024/"
  save_dir <- paste0(input_dir, 'drivernet_demo/')
  df <- data.table::fread(paste0(save_dir, 'significant_genes_DriverNet_cloneCB_v2.csv'))
  edges_df <- data.table::fread(paste0(save_dir, 'InteractomeFI_2022_wt_Score.csv.gz'))  
  dim(edges_df)  
  head(edges_df)
  dim(df)
  df[1:4,1:5]
  names(df)
  # unique(edges_df$Direction)
  
  # dolphin %>% 
  #   activate(nodes) %>%
  #   mutate(community = as.factor(group_louvain())) %>% 
  #   ggraph() + 
  #   geom_edge_link() + 
  #   geom_node_point(aes(colour = community), size = 5)
  
  
  driver_df <- df %>%
    dplyr::filter(`p-value`<0.05) %>%
    # dplyr::select(gene_symbol, ens_gene_id) %>%
    # dplyr::rename(name=gene_symbol) %>%
    dplyr::mutate(gene_type='driver gene')
  nodes_df <- df %>%
    dplyr::select(gene_symbol, ens_gene_id) %>%
    dplyr::rename(name=gene_symbol) %>%
    dplyr::mutate(gene_type='driver gene', label=name)
  
  cover_genes <- c()
  edges_main <- tibble::tibble()
  for(g in df$gene_symbol){
    gls <- df %>%
      dplyr::filter(gene_symbol==g) %>%
      dplyr::pull(cover_genes)
    gls <- unlist(strsplit(gls,','))
    
    if(length(gls)>0){
      print('______________')
      print(g)
      print(gls)
      cover_genes <- c(cover_genes, gls)
      ed_tmp <- expand.grid(c(g), gls) %>%
        dplyr::rename(x=Var1, y=Var2) %>%
        dplyr::mutate(edge_type='driver_cover_connected')
      edges_main <- dplyr::bind_rows(edges_main, ed_tmp)
    }
  }
  cover_genes <- unique(cover_genes)
  cover_genes <- cover_genes[!cover_genes %in% nodes_df$name]
  nodes_df2 <- tibble::tibble(name=cover_genes, gene_type='cover gene',label='')
  nodes_df <- dplyr::bind_rows(nodes_df, nodes_df2)
  dim(nodes_df)
  dim(edges_main)
  nodes_df <- nodes_df %>%
    dplyr::mutate(size=
                    case_when(gene_type=='driver gene' ~ 5,
                              TRUE ~ 2))
  dim(nodes_df)
  
  
  ## Adding pathways labels here
  ## If nodes are in pathways, coloring nodes by different colors, and edge color between nodes as well
  cls_labels <- sample.int(3, dim(nodes_df)[1], replace = T)
  nodes_df$pathway <- as.character(cls_labels)
  nodes_df <- nodes_df %>%
    dplyr::mutate(
      label=case_when(
        pathway!='1' ~ name,
        TRUE ~ label
      )
    )
  
  
  cols_use <- c('grey', 'green','red')
  names(cols_use) <- unique(nodes_df$pathway)
  nodes_df$pathway <- as.factor(nodes_df$pathway)
  
  # nodes_pw2 <- nodes_df %>%
  #   dplyr::filter(pathway =='2')
  # nodes_pw3 <- nodes_df %>%
  #   dplyr::filter(pathway =='3')
  edges_main <- edges_main %>%
    dplyr::mutate(
      edge_type=case_when(
        x %in% nodes_pw2$name & y %in% nodes_pw2$name ~ 'pathway2',
        x %in% nodes_pw3$name & y %in% nodes_pw3$name ~ 'pathway3',
        TRUE ~ edge_type)
    )
  
  edge_cols_use <- c('grey', 'darkgreen','darkred')
  names(edge_cols_use) <- c('driver_cover_connected','pathway2','pathway3')
  
  edge_cols_use <- c('grey', 'darkgreen','darkred')
  names(edge_cols_use) <- c('driver_cover_connected','pathway2','pathway3')
  
  edges_wd <- c(0.4, 0.8, 0.9)
  names(edges_wd) <- c('driver_cover_connected','pathway2','pathway3')
  
  gt_cols <- c(0.9, 0.1)
  names(gt_cols) <- c("driver gene","cover gene")
  sz <- c(0.9, 0.2)
  names(sz) <- c("driver gene","cover gene")
  fc1 <- c('bold', 'italic')
  names(fc1) <- c("driver gene","cover gene")
  ## get pathway labels here
  # make a tidy graph
  # nodes = NULL, edges = NULL, directed = TRUE, node_key = "name"
  ig <- tidygraph::tbl_graph(nodes = nodes_df, edges = edges_main, directed = FALSE)
  # ig
  # setting theme_graph 
  set_graph_style()
  
  ## To Do: need to check edge directions here 
  # basic plot
  
  ## To do: alpha scale, gene type=signf, not signif, cover genes
  pg <- ggraph(ig, layout='kk') + 
    geom_edge_link(aes(color = edge_type, width=edge_type), alpha=.2) + #, colour='darkgrey'
    geom_node_point(aes(size=size, colour = pathway), alpha=0.3) + #shape = gene_type, 
    geom_node_point(aes(size=size, alpha=gene_type), shape = 1) + #shape = gene_type, , colour='black'
    geom_node_text(aes(label = label), repel=TRUE, max.overlaps = Inf, size=3) + 
    scale_edge_colour_manual(values=edge_cols_use) + 
    scale_edge_width_manual(values=edges_wd) + 
    scale_color_manual(values = cols_use) + 
    scale_alpha_manual(values = gt_cols) +
    theme(legend.position = 'bottom')
  
  # scale_size_continuous() + 
  # scale_size_area(values=sz)
  pg
  # geom_edge_link(aes(color = factor(to), width = log(weight)), alpha = 0.5, 
  #                start_cap = circle(2, 'mm'), end_cap = circle(2, 'mm'))
  # Define colors of locations and characters 
  # V(pg)$color <- "grey20"
  # V(pg)$color[V(pg)$name %in% driver_df$gene_symbol] <- "red"
  
  # BiocManager::install('svglite', ask=F)
  ggsave(  
    filename = paste0(save_dir,"test_BC.svg"),  
    plot = pg,  
    height = 8, width = 12, dpi = 150)
  
  ## to do: significant genes with red color node
  ## size = # connection, with log(size)
  ## pathway genes: including cover genes, with label text
  
  
}



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

# install.packages('BiocManager')
# BiocManager::install("ggraph")
# devtools::install_github('thomasp85/ggraph')
library(ggraph)
library(igraph)
# install.packages("pak")
# pak::pak('thomasp85/ggraph')

test_ggraph <- function(){

  graph <- graph_from_data_frame(highschool)
  ggraph(graph, layout='kk') + 
    geom_edge_link(aes(colour = factor(year))) + 
    geom_node_point()
  
  
  hairball <- graph_from_data_frame(highschool)
  
  # Classify nodes based on popularity gain
  # pop1957 <- degree(delete_edges(hairball, which(E(hairball)$year == 1957)), 
  #                   mode = 'in')
  # pop1958 <- degree(delete_edges(hairball, which(E(hairball)$year == 1958)), 
  #                   mode = 'in')
  # V(hairball)$pop_devel <- ifelse(pop1957 < pop1958, 'increased',
  #                                 ifelse(pop1957 > pop1958, 'decreased', 
  #                                        'unchanged'))
  # V(hairball)$popularity <- pmax(pop1957, pop1958)
  E(hairball)$year <- as.character(E(hairball)$year)
  ggraph(hairball, layout = 'kk') + 
    geom_edge_link(aes(colour = year))
  
  library(dplyr)
  library(tidygraph)
  # BiocManager::install('tidygraph',force = TRUE)
  #add edge type
  nodes <- read.csv("http://users.dimi.uniud.it/~massimo.franceschet/ns/plugandplay/ggraph/dolphin_nodes.csv")
  edges <- read.csv("http://users.dimi.uniud.it/~massimo.franceschet/ns/plugandplay/ggraph/dolphin_edges.csv")
  edges = 
    edges %>% 
    mutate(type = sample(c("love", "friendship"), 
                         nrow(edges), 
                         replace = TRUE) )
  head(nodes)
  ## containing order of node
  head(edges)
  ## x, y coordinates are based on node id - node order
  # make a tidy graph
  # nodes = NULL, edges = NULL, directed = TRUE, node_key = "name"
  dolphin <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)
  dolphin
  
  # setting theme_graph 
  set_graph_style()
  
  # basic plot
  ggraph(dolphin, layout='kk') + 
    geom_edge_link(aes(color = type)) + 
    geom_node_point(aes(shape = sex))+ 
    geom_node_text(aes(label = name), repel=TRUE)
  
}
