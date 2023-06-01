suppressPackageStartupMessages({
  require("data.table")
  require("dplyr")
  require("optparse")
  require("ggtree")
  require("gtools")
  require("ape")
  require("ggnet")
  require("network")
  require("tools")
  require("phytools")
  require("igraph")
})



option_list <- list(make_option(c("-i", "--outliers_filtered_fn"), type="character", default=NULL, help="features_mtx_file", metavar="character"),
                    make_option(c("-o", "--output_fn"), type="character", default=NULL, help="output_file", metavar="character"),
                    make_option(c("-f", "--newick_fn"), type="character", default=NULL, help="grown_newick_fn", metavar="character"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


# newick_fn <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_wholedata/corrupt_grow/v2_FN_01/grown_tree.newick'
# output_fn <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_wholedata/corrupt_grow/sub_grown_tree.newick'
# outliers_filtered_fn <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_wholedata/corrupt_grow/filtered_none_cells.csv'
get_subgraph <- function(newick_fn, outliers_filtered_fn, output_fn){
  # Load tree
  # newick <- paste0(results_dir, 'tree.newick')
  tree <- read.tree(newick_fn)
  all_cells <- grep('cell', tree$tip.label, value = T)
  print(length(all_cells))
  
  save_dir <- paste0(dirname(output_fn),'/')
  if (!file.exists(save_dir)){
    dir.create(save_dir)
  }
  # cell_clones <- read.delim(cellclones, stringsAsFactors = FALSE, sep = ",")
  # clones_cellidx <- all_cells[all_cells %in% paste0('cell_',cell_clones$cell_id)]
  # length(clones_cellidx)
  # none_cells <- all_cells[!all_cells %in% clones_cellidx]
  
  filtered_outliers <- data.table::fread(outliers_filtered_fn, stringsAsFactors = F)
  outlier_cells <- unique(filtered_outliers$cells)
  print('DEBUG')
  print(length(outlier_cells))
  
  
  
  tree_cg <- ape::drop.tip(tree, outlier_cells, trim.internal =F, collapse.singles = F)
  # tree_cg
  # root(tree_cg)
  root_name <- ''
  print(root_name %in% tree_cg$node.label)
  for(i in rep(1:length(tree_cg$node.label),1)){
    if(tree_cg$node.label[i]==root_name){
      tree_cg$node.label[i] = 'ROOT'
    }
  }
  print('ROOT' %in% tree_cg$node.label)  # verification
  # ggtree(tree_cg)
  # subtree_fn <- paste0(cg_dir,'subtree.newick')
  # https://www.rdocumentation.org/packages/ape/versions/5.3/topics/write.tree
  write.tree(phy = tree_cg, file = output_fn, tree.names = F)
  
  # loci_ls <- unique(grep("locus", tree_cg$node.label, value = TRUE))
  # print(length(loci_ls))
}

get_subgraph(opt$newick_fn, opt$outliers_filtered_fn, opt$output_fn)





get_subgraph <- function(newick_fn, outliers_filtered_fn, output_fn){
  # Load tree
  newick <- paste0(results_dir, 'corrupt_grow/grown_tree.newick')
  tree <- read.tree(newick)
  all_cells <- grep('cell', tree$tip.label, value = T)
  print(length(all_cells))
  
  cells_clone <- read.csv(paste0(results_dir, 'corrupt_grow/cell_clones.csv'), check.names=F)
  dim(cells_clone)
  save_dir <- paste0(dirname(output_fn),'/')
  if (!file.exists(save_dir)){
    dir.create(save_dir)
  }
  # cell_clones <- read.delim(cellclones, stringsAsFactors = FALSE, sep = ",")
  # clones_cellidx <- all_cells[all_cells %in% paste0('cell_',cell_clones$cell_id)]
  # length(clones_cellidx)
  # none_cells <- all_cells[!all_cells %in% clones_cellidx]
  
  filtered_outliers <- data.table::fread(outliers_filtered_fn, stringsAsFactors = F)
  outlier_cells <- unique(filtered_outliers$cells)
  print('DEBUG')
  print(length(outlier_cells))
  
  ex_clones <- c('I','D','G','K','A')
  ex_clones <- c('L')
  for(c in ex_clones){
    outlier_cells <- cells_clone[cells_clone$clone_id != c,'cell_id']
    print(length(outlier_cells))
    outlier_cells <- paste0('cell_',outlier_cells)
    tree_cg <- ape::drop.tip(tree, outlier_cells, trim.internal =T, collapse.singles = T)
    all_cells <- grep('cell', tree_cg$tip.label, value = T)
    print(length(all_cells))
    # ggtree(tree_cg)
    write.tree(phy = tree_cg, file = paste0(results_dir, 'corrupt_grow/clone_',c,'.newick'), tree.names = F)
    
  }
  
  # tree_cg
  # root(tree_cg)
  root_name <- ''
  View(tree_cg$tip.label)
  print(root_name %in% tree_cg$node.label)
  for(i in rep(1:length(tree_cg$node.label),1)){
    if(tree_cg$node.label[i]==root_name){
      tree_cg$node.label[i] = 'ROOT'
    }
  }
  print('ROOT' %in% tree_cg$node.label)  # verification
  # ggtree(tree_cg)
  # subtree_fn <- paste0(cg_dir,'subtree.newick')
  # https://www.rdocumentation.org/packages/ape/versions/5.3/topics/write.tree
  write.tree(phy = tree_cg, file = paste0(results_dir, 'corrupt_grow/cloneA.newick'), tree.names = F)
  
  # loci_ls <- unique(grep("locus", tree_cg$node.label, value = TRUE))
  # print(length(loci_ls))
}
