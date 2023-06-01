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

remove_na_from_list <- function(some_list) {
  # Why doesn't it keep names?
  tmp <- lapply(names(some_list), function(x) {some_list[[x]][!is.na(some_list[[x]])] })
  names(tmp) <- names(some_list)
  tmp
}

mat <- readRDS('/Users/sohrabsalehi/Desktop/SC-1311/SA609/processed_data/cnv_data.rds')
binary_mat <- fasttable('/Users/sohrabsalehi/Desktop/SC-1311/SA609/processed_data/bincount/sa609sa666sa777_cnvs_corrupt_no_padding_sans_sphase.csv')
tree <- read.tree('')

# On shahlab15
all_mcmc_samples <- '/ssd/sdb1/ssalehi/projects/corrupt21_/corrupt-nextflow/deliverables/upSA609/samples'


dat = fasttable('phylo.csv', skip=20, nrows=1)$V2
tree = read.newick(text=dat)





function() {
  tree = read.newick('/Users/sohrabsalehi/Desktop/SC-1311/SA1000/processed_data/bincount/SuperRun10K_10chains/tree.newick')
  g <- read_ltm_tree(tree_2_edge_list(tree))
  cells <- g.get.cells(g)
  cells <- cells[cells != 'root']
  vpaths = igraph::shortest_paths(graph = g, from = 'root',  to = cells[1], mode = 'out', output = 'vpath')
  
  vpaths <- igraph::shortest_paths(graph = g, from = 'root',  to = cells, mode = 'out', output = 'vpath')
  cos <- vpaths$vpath
  cos <- lapply(1:length(cos), function(i) igraph::as_ids(cos[[i]]))
  igraph::as_ids(cos)
  as.vector()
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
  
  # For each bin, find all its cell descendents and set them to 1
  # locus_i <- 0
  # for (locus in loci) {
  #   # locus = loci[600]
  #   locus_i <- locus_i + 1
  #   desc <- get_desc_names(the_graph = g, the_node = locus)
  #   desc <- unlist(unname(desc))
  #   desc <- desc[!grepl('locus_', desc)]
  #   mat[ match(desc, rownames(mat))  , locus_i] <- 1
  # }
  vpaths <- igraph::shortest_paths(graph = g, from = 'root',  to = cells, mode = 'out', output = 'vpath')$vpath
  for (c_i in 1:length(cells)) {
    pars <- igraph::as_ids(vpaths[[c_i]])
    mat[c_i, match(pars[grepl('locus', pars)], colnames(mat))] <- 1
  }
  
  mat
}


qq = compute_assignment_vector(tree)


# TODO: find a faster one traverse version of this
# 


