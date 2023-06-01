# Combine more channels with the CNV data:
# For now, per chromosome, show:
# 1. CNV data
# 2. GC-corrected read_counts
# 3. Binary-input to corrupt
# 4. Posterior assignment, computed as in how many of the posterior samples, a cell is attached to a bin or its descendent
# Idea: 
# Compile the matrices separately
# Concatenate them
# Design the discrete colour scheme manually
# Use the existing tree-version


# test dataset 
# SA609
# CNV data

library(data.table)
library(phytools)
library(igraph)

sort_chromosome <- function(chr_list) {
  options(scipen=999) # prevent scientific notation
  chr_list[chr_list == 'X'] <- '40' 
  chr_list[chr_list == 'Y'] <- '55' 
  chr_list <- as.numeric(chr_list)
  chr_list <- chr_list[order(chr_list)]
  chr_list[chr_list == '40'] <- 'X' 
  chr_list[chr_list == '55'] <- 'Y' 
  chr_list
}


parse_bin_names <- function(bin_names) {
  # Remove corrupt_tree locus tag if it's there
  bin_names <- gsub('locus_', '', bin_names)
  chr <- gsub('([0-9]+|X|Y)_[0-9]+_[0-9]+', '\\1', bin_names)
  start <- as.numeric(gsub('([0-9]+|X|Y)_([0-9]+)_[0-9]+', '\\2', bin_names))
  end <- as.numeric(gsub('([0-9]+|X|Y)_([0-9]+)_([0-9]+)', '\\3', bin_names))
  data.frame(chr = chr, start = start, end = end, stringsAsFactors = F)
}

fasttable <- function(...) as.data.frame(fread(..., stringsAsFactors = F))

g.get.cells <- function(g=NULL, edge_list_path=NULL) {
  if (is.null(g)) g <- read_ltm_tree(edge_list_path = edge_list_path)
  get_cells(V(g)$name)
}

get_cells <- function(str_arry) str_arry[!grepl('locus_', str_arry)]

g.get.loci <- function(g=NULL, edge_list_path=NULL) {
  if (is.null(g)) g <- read_ltm_tree(edge_list_path = edge_list_path)
  get_loci(V(g)$name)
}

get_loci <- function(str_arry) grep('locus_', str_arry, value = T)

get_decendents <- function(clone_roots, the_graph, min_cell_per_clone=1) {
  decendents <- list()
  for (clone_root in clone_roots) {
    decend_node_ids <- decendents_ids(clone_root, the_graph)
    
    cells_in_clone <- length(decend_node_ids)
    if (cells_in_clone > min_cell_per_clone) {
      decendents[[clone_root]] <- decend_node_ids
    }
  }
  
  return(decendents)
}

decendents_ids <- function(root_id, the_graph) {
  decend_node_ids <- igraph::bfs(
    the_graph,
    root=c(root_id), order=TRUE, neimode='out',
    unreachable=FALSE
  )
  decend_node_ids <- as.numeric(decend_node_ids$order)
  decend_node_ids <- decend_node_ids[!is.na(decend_node_ids)]
  return(decend_node_ids)
}

get_desc_names <- function(the_graph, the_node) {
  desc <- get_decendents(the_node, the_graph=the_graph)
  remove_na_from_list(node_number_to_label(the_graph, desc))
}

setup_graph <- function(g=NULL) {
  V(g)$id <- seq(vcount(g))
  V(g)$name <- seq(vcount(g))
  h <- get_height_dat(g)
  stopifnot(all(h$id == V(g)$name))
  V(g)$height <- as.integer(h$dist)
  g
}

node_number_to_label <- function(the_graph, groups) {
  all_dat <- data.frame(id=V(the_graph)$id, label=V(the_graph)$name, stringsAsFactors = F)
  res <- list()
  for (clade_root in names(groups)) {
    cells_in_clone <- groups[[clade_root]]
    cells <- all_dat$label[match(c(-1, cells_in_clone), all_dat$id)]
    res[[clade_root]] <- cells
  }
  res
}

read_ltm_tree <- function(edge_list) {
  # Find the root
  g <- igraph::graph_from_edgelist(as.matrix(edge_list))
  V(g)$id <- seq(vcount(g))
  return(g)
}

tree_2_edge_list <- function(tree) {
  tree$node.label[1] <- "root"
  node_names <- c(tree$tip.label, tree$node.label)
  
  edges <- tree$edge
  edges <- data.frame(
    source=node_names[edges[, 1]],
    target=node_names[edges[, 2]],
    stringsAsFactors=FALSE
  )
  edges$source <- gsub("cell_", "", edges$source)
  edges$target <- gsub("cell_", "", edges$target)
  
  return(edges)
}

find_root <- function(g) {
  the_root <- NULL
  n_nodes <- vcount(g)
  for (i in seq(n_nodes)) {
    s1 = bfs(g, root=i, order=TRUE, neimode = 'out', unreachable=FALSE)
    if (any(is.na(s1$order)) == FALSE)  {
      the_root <- i
    }
  }
  if (is.null(the_root)) {
    print('Warning! No root was found...')
  }
  the_root
}

remove_na_from_list <- function(some_list) {
  # Why doesn't it keep names?
  tmp <- lapply(names(some_list), function(x) {some_list[[x]][!is.na(some_list[[x]])] })
  names(tmp) <- names(some_list)
  tmp
}

get_height_dat <- function(g) {
  the_root <- find_root(g)
  height_search <- bfs(g, root=c(the_root), order=TRUE, neimode = 'out', unreachable=FALSE, dist = TRUE, rank=TRUE, succ=TRUE, father=TRUE, pred=TRUE)
  res <- data.frame(id=names(height_search$dist), dist=as.numeric(height_search$dist), stringsAsFactors = F)
  res
}

# Returns a binary genotype matrix based on the tree
# Each element i,j in mat, is 1 if cell i has locus j otherwise zero
compute_assignment_vector <- function(tree) {
  g <- read_ltm_tree(tree_2_edge_list(tree))
  cells <- g.get.cells(g)
  loci <- g.get.loci(g)
  mat <- matrix(0, nrow=length(cells), ncol = length(loci))
  rownames(mat) <- cells
  colnames(mat) <- loci
  
  qq <- as.list(1:ncol(mat))
  names(qq) <- colnames(mat)
  vpaths <- igraph::shortest_paths(graph = g, from = 'root',  to = cells, mode = 'out', output = 'vpath')$vpath
  
  for (c_i in 1:length(cells)) {
    pars <- igraph::as_ids(vpaths[[c_i]])
    pars <- pars[-c(1, length(pars))]
    #mat[c_i, match(pars[grepl('locus', pars)], colnames(mat))] <- 1
    mat[c_i, unlist(unname(qq[pars]))] <- 1
  }
  
  stopifnot(sum(mat) < nrow(mat) * ncol(mat))
  mat
}



# The main test
function() {
  mat <- readRDS('~/Desktop/SC-1311/SA609/processed_data/cnv_data.rds')
  binary_mat <- fasttable('/Users/sohrabsalehi/Desktop/SC-1311/SA609/processed_data/bincount/sa609sa666sa777_cnvs_corrupt_no_padding_sans_sphase.csv')
  tree <- read.tree('')
  
  # On shahlab15
  all_mcmc_samples <- '/ssd/sdb1/ssalehi/projects/corrupt21_/corrupt-nextflow/deliverables/upSA609/samples'
  
  
  
  ## Since parsing each tree takes up to 8 seconds, and we 
  # have at least 10,000 of them, do it in parallel
  detectCores() # detect number of cores
  # Go with 10 for now


  # see here:
  # https://www.r-bloggers.com/parallel-r-loops-for-windows-and-linux/
  
  # 1. Divide the MCMC samples into nCore parts
  # 2. Run the whole alg for each subpart and save chunks
  # 3. Combine the chunks
  
  res <-foreach(i=1:10) %do% {
    i
  }
  
}

 
create_dummy_graph <- function() {
  igraph_options(add.vertex.names=T)
  tt <- igraph::make_tree(20)
  V(tt)$name <- seq(vcount(tt))
  
  #tt <- setup_graph(g = tt)
  plot(tt)
  tmp <- igraph::dfs(graph = tt, root = 1, neimode = 'out', unreachable = FALSE, order = TRUE, father = TRUE, order.out = TRUE)
  pn <- igraph::as_ids(tmp$father)
  cn <- igraph::as_ids(tmp$order)
  V(tt)$name
  
  #ff = c(NA, pn[-c(length(pn))])
  ff = pn
  names(ff) <- cn
  ff
  
  oo = igraph::as_ids(seq = tmp$order.out)
  data.frame(cn = cn, 
             pn =  c(NA, pn[-c(1)]))
  
  
  igraph_options(add.vertex.names=F)
  igraph_options()$add.vertex.names
  h<-make_tree(10)
  set.seed(1)
  #h<-set_vertex_attr(h, "name", value=sample(1:10,10))
  gg <- dfs(h, root=1, neimode='out', order=TRUE, father=TRUE,unreachable=FALSE)
  str(gg$father)
  plot(h)
}


# Read the phylo.csv 
parse_phylo_mcmc <- function(phylos.path) {
  # phylos.path <- 'phylo.csv'
  mat <- NULL
  nSamples <- as.numeric(gsub(' phylo.csv', '', system(sprintf('wc -l %s', phylos.path), intern = T))) - 1
  for (i in 1:nSamples) {
    print(i)
    if (i %% 100 == 0) {
      print(sprintf('Processing %d-th tree.', i))
    }
    dat <- fasttable('phylo.csv', skip= (i + 1), nrows=1)$V2
    tree <- read.newick(text=dat)
    tmp <- compute_assignment_vector(tree)
    if (is.null(mat)) {
      mat <- tmp
    } else {
      mat <- mat + tmp
    }
  }
  mat <- mat / nSamples
  saveRDS(object = mat, file = 'posterior_genotype_mat.rds')
  mat
}

parse_phylo_mcmc_old <- function(phylos.path, start, nTrees, temp_path) {
  # phylos.path <- 'phylo.csv'
  mat <- NULL
  dat <- fasttable(file.path(phylos.path ,'phylo.csv'), skip=start, nrows=nTrees)
  nSamples <- nrow(dat)
  for (i in 1:nSamples) {
    if (i %% 10 == 0) {
      print(sprintf('Processing %d-th tree.', i))
    }
    #tree <- read.newick(text=dat$value[i])
    tree <- read.newick(text=dat$V2[i])
    tmp <- compute_assignment_vector(tree)
    if (is.null(mat)) {
      mat <- tmp
    } else {
      mat <- mat + tmp
    }
  }
  mat <- mat / nSamples
  dir.create(temp_path)
  saveRDS(object = mat, file = file.path(temp_path, sprintf('posterior_genotype_mat_%d_%d.rds', start, nTrees)))
  mat
}



#### Minimum working parallel code on Shahlab
library(foreach)
library(doMC)
nCores <- 20
registerDoMC(nCores) 
all_mcmc_samples <- '/ssd/sdb1/ssalehi/projects/corrupt21_/corrupt-nextflow/deliverables/upSA609/samples'
nTrees <- 2
res_list <- foreach(i=1:nCores) %dopar% {
#res_list <- foreach(i=1:10) %do% {
  print(sprintf('Processing batch %d', i))
  parse_phylo_mcmc_old(phylos.path = all_mcmc_samples, start = (i-1)*5 + 2 , nTrees = nTrees, temp_path = '/shahlab/ssalehi/scratch/corrupt/data/temp')
}


mat <- NULL
for (i in 1:length(res_list)) {
  if (is.null(mat))
    mat <- res_list[[i]]
  else 
    mat <- mat + res_list[[i]]
  
  mat <- mat / length(res_list)
}

# TODO: speed things up - we don't want the post processing to take more than the tree-sampling
# Idea1: only compute the paths for one of the siblings and copy for the rest
# Idea2: Ues the decode in nowellpack




# mat_post <- parse_phylo_mcmc(phylos.path = 'phylo.csv')


## Generate some dummy data

# The dummy posterior genotype based on the trees
generate_dummy_post_ass <- function(tree) {
  g <- read_ltm_tree(tree_2_edge_list(tree))
  cells <- g.get.cells(g)
  loci <- g.get.loci(g)
  mat <- matrix(runif(n = length(cells)*length(loci)), nrow=length(cells), ncol = length(loci))
  rownames(mat) <- cells
  colnames(mat) <- loci
  mat
}

# The dummy binary input, based on the filtered.csv
generate_dummy_bin_input <- function(tree, nBins = 200) {
  g <- read_ltm_tree(tree_2_edge_list(tree))
  cells <- g.get.cells(g)
  loci <- g.get.loci(g)
  mat <- matrix(sample(c(0,1), length(cells)*nBins, TRUE), nrow=length(cells), ncol = nBins)
  rownames(mat) <- cells
  colnames(mat) <- sample(loci, nBins)
  mat
}


# The GC-corrected readcounts?
# How does that data look like?


# The CNV
generate_dummy_cn <- function(tree, dt = 'SA922') {
  # between 0 and 11, mostly around 2
  g <- read_ltm_tree(tree_2_edge_list(tree))
  cells <- g.get.cells(g)
  loci <- g.get.loci(g)
  cnv_mat <- readRDS(sprintf('~/Desktop/SC-1311/%s/processed_data/cnv_data.rds', dt))
  probs <- table(as.vector(as.matrix(cnv_mat)))
  mat <- matrix(sample(x = as.numeric(names(probs)), size =  length(cells)*length(loci), prob = unname(probs), replace = TRUE), nrow=length(cells), ncol = length(loci))
  rownames(mat) <- cells
  colnames(mat) <- loci
  mat
}


function() {
  tree <- read.newick('~/Desktop/SC-1311/SA922/processed_data/bincount/10k_10chains/padded/tree.newick')
  mat_post <- generate_dummy_post_ass(tree)
  filtered <- generate_dummy_bin_input(tree, 100)
  mat_cn <- generate_dummy_cn(tree)
  
  # Lets go for 10 cells
  cells <- rownames(mat_post)[1:10] 
  mat_post <- mat_post[rownames(mat_post) %in% cells, ]
  filtered <- filtered[rownames(filtered) %in% cells, ]
  mat_cn <- mat_cn[rownames(mat_cn) %in% cells, ]
  
  
  stopifnot(all(colnames(mat_post) == colnames(filtered)))
  stopifnot(all(colnames(mat_cn) == colnames(filtered)))
  
  
  
  # Interleave by chromosome
  res_mat <- list()
  
  filtered_bins <- parse_bin_names(colnames(filtered))
  filtered_chrs <- sort_chromosome(unique(filtered_bins$chr))
  
  for (chr in filtered_chrs) {
    res_mat[[chr]] <- cbind()
  }
  
  
  # Help me
  
  
}






