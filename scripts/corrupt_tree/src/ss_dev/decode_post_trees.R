library(data.table)
library(phytools)
library(igraph)

fasttable <- function(...) as.data.frame(fread(..., stringsAsFactors = F))

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
  #mat <- mat / nSamples
  
  dir.create(temp_path)
  saveRDS(object = mat, file = file.path(temp_path, sprintf('posterior_genotype_mat_%d_%d.rds', start, nTrees)))
  mat
}

#### Minimum working parallel code on Shahlab
library(foreach)
library(doMC)

# TODO: output Alex's decode code
# TODO: Kevin has access
# TODO: Use all samples... match ntrees


get_chunk_sizes <- function(nWorker, nTask) {
  # nWorker = 38; nTask = 10000
  m <- nTask %/% nWorker
  n <- nTask - m*nWorker
  res <- array(m, nWorker)
  # Redistribute
  if (n > 0) {
    for (i in 1:n) {
      res[i] <- res[i] + 1
    }
  }
  stopifnot(sum(res) == nTask)
  res
}


parallel_corrupt_decode <- function(all_mcmc_samples, nCores = 10, nTotalTrees = 10000) {
  registerDoMC(nCores) 
  nTrees <- get_chunk_sizes(nWorker = nCores, nTask = nTotalTrees)
  
  res_list <- foreach(i=1:nCores) %dopar% {
    #res_list <- foreach(i=1:10) %do% {
    print(sprintf('Processing batch %d', i))
    parse_phylo_mcmc_old(phylos.path = all_mcmc_samples, start = (i-1)*nTrees[i] + 2 , nTrees = nTrees[i], temp_path = '/shahlab/ssalehi/scratch/corrupt/data/temp')
  }
  
  mat <- NULL
  for (i in 1:length(res_list)) {
    if (is.null(mat))
      mat <- res_list[[i]]
    else 
      mat <- mat + res_list[[i]]
  }
  
  mat <- mat / nTotalTrees
  
  output_path <- file.path(dirname(all_mcmc_samples), 'mat.rds')
  saveRDS(output_path)
  print(sprintf('Output at %s', file.path(dirname(all_mcmc_samples), 'mat.rds')))
}

nCores <- 30
all_mcmc_samples <- '/ssd/sdb1/ssalehi/projects/corrupt21_/corrupt-nextflow/deliverables/upSA609/samples'
parallel_corrupt_decode(all_mcmc_samples = all_mcmc_samples, nCores = nCores, nTotalTrees = 10000)

# TODO: speed things up - we don't want the post processing to take more than the tree-sampling
# Idea1: only compute the paths for one of the siblings and copy for the rest
# Idea2: Ues the decode in nowellpack