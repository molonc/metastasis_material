
# v <- NULL
# f <- NULL
# qu <- c()
# p <- NULL
# p_i <- NULL
# cur_indx <- NULL

# global vars
v <<- NULL
f <<- NULL
qu <<- NULL
p <<- NULL
p_i <<- NULL
cur_indx <<- NULL
leaves <<- NULL


make_el <- function(total_nodes = 20, do_plot = FALSE, seed = 1) {
  set.seed(seed)
  g <- igraph::make_tree(total_nodes)
  
  el <- igraph::as_edgelist(g)
  el[, 1] <- paste0('n', el[, 1])
  el[, 2] <- paste0('n', el[, 2])
  if (do_plot) {
    plot(g)
    plot(igraph::graph_from_edgelist(el = el, directed = T))
  } 
  
  el
}

initialise_sdfs <- function(el) {
  leaves <<- setdiff(unique(el[, 2]), unique(el[, 1]))
  nLeavs <- length(leaves)
  nodes <- unique(as.vector(el))
  nNodes <- length(nodes)
  loci <- setdiff(nodes, leaves)
  nLoci <- length(loci)
  
  v <<- array(NA, nNodes) # visited
  f <<- array(NA, nNodes) # father
  qu <<- c(get_root_name(el)) # queue
  p <<- c() # the current path
  p_i <<- matrix(0, nrow = nLeavs, ncol = nLoci)
  
  rownames(p_i) <<- leaves
  colnames(p_i) <<- loci
  loci_dic <<- 1:ncol(p_i)
  names(loci_dic) <<- colnames(p_i)
  
  cell_dic <- 1:nrow(p_i)
  names(cell_dic) <- rownames(p_i)
  
  cur_indx <<- 0
}

get_root_name <- function(edge_list) {
  unique(edge_list[, 1])[!(unique(edge_list[, 1]) %in% unique(edge_list[, 2]))]
}

get_child <- function(el, node) {
  el[, 2][el[,1] == node]
}

get_father <- function(el, node) {
  if (node == 1) return(NA)
  el[, 1][el[ ,2] == node]
}

is_leaf <- function(el, node) {
  return(node %in% leaves)
}


sdfs <- function(qu, el) {
  cur_indx <<- cur_indx + 1
  cN <- qu[length(qu)]; 
  if (length(qu) > 0) qu <- qu[-c(length(qu))]
  fN <- get_father(el = el, node = cN)
  v[cur_indx] <<- cN
  f[cur_indx] <<- fN
  if (is_leaf(el, cN))  {
    p_i[cN - nLoci, p] <<- 1
  } else {
    # Remove all after its father
    if (!is.na(fN)) {
      p <<- p[1:(which(p == fN))]
    }
    p <<- c(p, cN)
    qu <- c(qu, get_child(el, cN))
  }
  if (length(qu) > 0) return(sdfs(qu, el))
  qu
}

# initialise_sdfs(nNodes = nNodes, nLeaves = nLeaves, nLoci = nLoci)
sdfs_fast <- function() {
  cur_indx <<- cur_indx + 1
  if (cur_indx %% 100 == 0) print(cur_indx)
  cN <- qu[length(qu)]; 
  qu <<- qu[-c(length(qu))]
  # Get father
  fN <- el[, 1][el[ ,2] == cN]
  if (length(fN) == 0) fN <- NA
  v[cur_indx] <<- cN
  #f[cur_indx] <<- fN
  if (cN %in% leaves)  {
    p_i[which(rownames(p_i) == cN), loci_dic[p]] <<- 1
  } else {
    # Remove all after its father
    if (!is.na(fN)) {
      p <<- p[1:(which(p == fN))]
    }
    p <<- c(p, cN)
    qu <<- c(qu, el[, 2][el[,1] == cN])
  }
  if (length(qu) > 0) return(sdfs_fast())
  qu
}


sdfs_fast_mem <- function() {
  while (length(qu) > 0) {
    cur_indx <<- cur_indx + 1
    if (cur_indx %% 1000 == 0) print(cur_indx)
    cN <- qu[length(qu)]; 
    qu <<- qu[-c(length(qu))]
    # Get father
    fN <- el[, 1][el[ ,2] == cN]
    if (length(fN) == 0) fN <- NA
    v[cur_indx] <<- cN
    #f[cur_indx] <<- fN
    if (cN %in% leaves)  {
      #p_i[which(rownames(p_i) == cN), unlist(unname(loci_dic[p]))] <<- 1
      p_ii[[cN]] <<- p
    } else {
      # Remove all after its father
      if (!is.na(fN)) {
        p <<- p[1:(which(p == fN))]
      }
      p <<- c(p, cN)
      qu <<- c(qu, el[, 2][el[,1] == cN])
    }
  }

  qu
}



#el <- make_el(total_nodes = 1000, seed = 10)

tree <- read.newick('/Users/sohrabsalehi/Desktop/SC-1311/SA1000/processed_data/bincount/SuperRun10K_10chains/tree.newick')
g <- read_ltm_tree(tree_2_edge_list(tree))
el <- igraph::as_edgelist(g)

# initialise_sdfs(el = el)
# res <- sdfs_fast()



p_ii <- list()
initialise_sdfs(el = el)
res <- sdfs_fast_mem()

for (c_i in 1:length(p_ii)) {
  cell_name <- names(p_ii)[c_i]
  pars <- p_ii[[cell_name]][-c(1)]
  p_i[cell_dic[[cell_name]], unlist(unname(loci_dic[pars]))] <- 1
}



qu
print(v)
print(f)
print(p_i)
p

junk <- function() {
  # Convert a named graph to this version
  # qq = el
  # qq[, 1] <- as.character(el[, 1])
  # qq[, 2] <- as.character(el[, 2])
  # 
  # qq = el
  # qq[, 2] <- paste0('locus', as.character(el[, 2]))
  # qq[, 1] <- paste0('locus', as.character(el[, 1]))
  # 
  # res <- convert_edge_list_to_ape(edge_list_mat = qq, leaf_names = setdiff(unique(qq[, 2]), unique(qq[, 1])))
  # res$modified_edge_list
  # 
  # get_root_name <- function(edge_list) {
  #   unique(edge_list$source)[!(unique(edge_list$source) %in% unique(edge_list$target))]
  # }
  # 
  # convert_edge_list_to_ape <- function(edge_list_mat=NULL, leaf_names) {
  #   options(stringsAsFactors = F)
  #   edge_list <- as.data.frame(edge_list_mat)
  #   colnames(edge_list) <- c('source', 'target')
  #   
  #   root_name <- get_root_name(edge_list)
  #   all_nodes <- unique(c(edge_list$source, edge_list$target))
  #   
  #   internal_nodes <- setdiff(all_nodes, leaf_names)
  #   internal_nodes <- setdiff(internal_nodes, root_name)
  #   
  #   n_leaves <- length(leaf_names)
  #   n_internal <- length(internal_nodes) + 1
  #   
  #   d1 <- data.frame(node_name=leaf_names, node_number=seq(n_leaves))
  #   d2 <- data.frame(node_name=root_name, node_number=n_leaves+1)
  #   # Account for a start where the only loci is the root
  #   if (length(internal_nodes) == 0) {
  #     d <- rbind(d1, d2)
  #   } else {
  #     d3 <- data.frame(node_name=internal_nodes, node_number=seq(n_internal-1))
  #     d3$node_number <- d3$node_number + n_leaves + 1
  #     d <- rbind(d1, d2, d3)
  #   }
  #   
  #   edges <- edge_list
  #   
  #   for (node in all_nodes) {
  #     edges$source[edges$source == node] <- d$node_number[d$node_name == node]
  #     edges$target[edges$target == node] <- d$node_number[d$node_name == node]
  #   }
  #   
  #   edges$source <- as.numeric(edges$source)
  #   edges$target <- as.numeric(edges$target)
  #   
  #   edges <- as.matrix(edges)
  #   colnames(edges) <- NULL
  #   
  #   d <- d[order(d$node_number), ]
  #   tr <- list(edge = edges, tip.label = d$node_name[1:n_leaves], Nnode = n_internal)
  #   class(tr) <- "phylo"
  #   
  #   Nedge <- nrow(tr$edge)
  #   tr$edge.length <- rep(1, Nedge)
  #   list(tree=tr, node_names=d, modified_edge_list=edges)
  # }
  # 
  # g1 <- igraph::graph_from_edgelist(el = qq, directed = T)
  # V(g1)$name
  # as_ids(V(g1))
}