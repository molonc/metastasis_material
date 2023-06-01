# Launch this function in R to convert newick tree to graphml format
# Then load graphml format here in python code 
suppressPackageStartupMessages({
  require("ggtree")
  # require("ape")
  require("optparse")
  require("igraph")
})
option_list <- list(make_option(c("-i", "--inputfile"), type="character", default=NULL, help="input_directory", metavar="character"),
                    make_option(c("-o", "--outputfile"), type="character", default=NULL, help="output_directory", metavar="character"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


tree_2_edge_list <- function(tree) {
  tree$node.label[1] <- "root"
  node_names <- c(tree$tip.label, tree$node.label)
  
  edges <- tree$edge
  edges <- data.frame(source = node_names[edges[, 1]],
                      target = node_names[edges[, 2]],
                      stringsAsFactors = FALSE
  )
  edges$source <- gsub("cell_", "", edges$source)
  edges$target <- gsub("cell_", "", edges$target)
  
  return(edges)
}

read_ltm_tree <- function(edge_list) {
  # Find the root
  g <- igraph::graph_from_edgelist(as.matrix(edge_list))
  V(g)$id <- seq(vcount(g))
  return(g)
}

convert_newick2graphml <- function(input_file, output_file){
  # save_dir <- paste0(results_dir,'graph_cut/')
  # output_file <- paste0(save_dir,'/tree.graphml')
  save_dir <- dirname(output_file)
  if (!file.exists(save_dir)){
    dir.create(save_dir)
  }
  newick <- input_file
  # newick <- paste0(results_dir, 'tree.newick')
  tree <- read.tree(newick)
  edge_list <- tree_2_edge_list(tree)
  g <- read_ltm_tree(edge_list)
  igraph::write.graph(g, output_file,format = 'graphml')

}

print(opt$inputfile)
print(opt$outputfile)
convert_newick2graphml(opt$inputfile, opt$outputfile)

