source_dir <- '/home/htran/Projects/farhia_project/rscript/dlp/visualize_tree/'
source(paste0(source_dir,'trajectory_utils.R'))

get_clones_subgraph <- function(results_dir, cellclones){
  # Load tree
  newick <- paste0(results_dir, 'tree.newick')
  tree <- read.tree(newick)
  # edge_list <- tree_2_edge_list2(tree)
  # View(head(edge_list))
  
  # g <- read_ltm_tree(edge_list)
  
  # all_cells <- gsub('cell_', '', grep('cell', tree$tip.label, value = T))
  # length(all_cells)
  # cell_clones <- read.delim(cellclones, stringsAsFactors = FALSE, sep = ",")
  # none_cells <- all_cells[!all_cells %in% cell_clones$cell_id]
  # length(none_cells)
  
  
  all_cells <- grep('cell', tree$tip.label, value = T)
  length(all_cells)
  cell_clones <- read.delim(cellclones, stringsAsFactors = FALSE, sep = ",")
  clones_cellidx <- all_cells[all_cells %in% paste0('cell_',cell_clones$cell_id)]
  length(clones_cellidx)
  none_cells <- all_cells[!all_cells %in% clones_cellidx]
  
  
  # none_cells <- none_cells[!none_cells %in% grep('locus_', none_cells, value=T)]
  length(none_cells)
  cg_dir <- paste0(results_dir,'corrupt_grow/')
  if (!file.exists(cg_dir)){
    dir.create(cg_dir)
  }
  
  # edge_list <- edge_list[!edge_list$target %in% none_cells, ]
  # dim(edge_list)
  # edge_list_filtered_path <- paste0(cg_dir,'edges_filtered_SA919.csv')
  # write.csv(edge_list, file = edge_list_filtered_path, quote=F, row.names = F)
  
  # test corrupt grow
  length(none_cells)
  # tree_cg <- ape::drop.tip(tree, none_cells, rooted = T)
  tree_cg <- ape::drop.tip(tree, none_cells, trim.internal =F, collapse.singles = F)
  # tree_cg <- root(tree_cg,1)
  
  # root_name <- get_root_name(edge_list)
  root_name <- ''
  ls_nodes <- unique(c(tree_cg$tip.label, tree_cg$node.label))
  
  root_name %in% tree_cg$node.label
  root_name %in% ls_nodes
  sum('' %in% tree_cg$node.label)
  'ROOT' %in% tree_cg$node.label
  for(i in rep(1:length(tree_cg$node.label),1)){
    if(tree_cg$node.label[i]==root_name){
      tree_cg$node.label[i] = 'ROOT'
    }
  }
  
  # res_cg <- convert_edge_list_to_ape(as.matrix(edge_list))
  
  
  # tree_cg <- subtree
  # tree_cg <- unroot(tree_cg)
  # is.rooted(tree_cg)
  # tree_cg <- as.phylo(tree_cg)
  
  # tree_cg$tip.label <- paste0('cell_',tree_cg$tip.label)
  subtree_fn <- paste0(cg_dir,'subtree1.newick')
  ggtree(tree_cg)
  write.tree(phy = tree_cg, file = paste0(cg_dir,'subtree1.newick'))
  
  write.tree(phy = tree_cg, file = paste0(cg_dir,'subtree1.newick'), tree.names = F)
  write.tree(tree_cg, file = paste0(cg_dir,'subtree2.newick'), append = FALSE, digits = 10, tree.names = T)
  newick_path <- paste0(cg_dir,'subtree1.newick')
  if(file.exists(subtree_fn)){
    return(subtree_fn)
  } else{
    return(NULL)
  }
  
  
  loci_ls1 <- unique(none_cells$loci)
  length(loci_ls1)
  sum(tree_cg$node.label %in% loci_ls1)
  sum(tree_cg$tip.label %in% loci_ls1)
  
  sum(tree$node.label %in% loci_ls1)
  # library(treeio)
  # tree_cg <- tree_subset(tree, "cell_SA919X3XB08939-A98232B-R49-C11", levels_back = 3)
  # sum(tree_cg$tip.label %in% none_cells)
  # newick_clones <- paste0(cg_dir, 'tree.newick.clones')
  # tree_clones <- read.tree(newick_clones)
  # class(tree_clones)
  # 
  # newick <- paste0(results_dir, 'tree.newick')
  # tree <- read.tree(newick)
  # class(tree)
  # the_graph <- read_ltm_tree(tree_2_edge_list(tree_clones))
  # V(the_graph)
  # E(the_graph)
  # p <- ggtree(tree_cg)
  # p
  sum(unique(none_cells$loci) %in% tree$node.label)
  length(tree$node.label)
  tree$node.label[!tree$node.label %in% unique(none_cells$loci)]
}

# get none cells from filtered file 
results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local/'
plot_title <- 'sohrab_version'
cellclones <- paste0(results_dir,'tree_cut_out/cell_clones_if_0.02_af_0.75_p0.75_e0.04.csv')

get_none_cells_cnstate <- function(results_dir, cellclones, tag='sohrab_version'){
  # Get tree with cells in clones only
  subtree_fn <- get_clones_subgraph(results_dir, cellclones)
  cg_dir <- paste0(results_dir,'corrupt_grow/')
  if (!file.exists(cg_dir)){
    dir.create(cg_dir)
  }
  filtered_df <- read.csv(paste0(results_dir, 'filtered.csv'), check.names = F, stringsAsFactors = FALSE)
  cell_clones <- read.delim(cellclones, stringsAsFactors = FALSE, sep = ",")
  filtered_nonecell_df <- filtered_df[!filtered_df$cells %in% paste0('cell_',cell_clones$cell_id),]
  dim(filtered_nonecell_df)
  dim(filtered_df)
  length(unique(filtered_nonecell_df$cells))
  data.table::fwrite(filtered_nonecell_df, paste0(cg_dir, 'filtered_none_cells.csv'), row.names = F, quote = F)
  newick_path <- subtree_fn
  # subtree_fn <- paste0(cg_dir,'tree.newick.clones')
  output_dir <- cg_dir
  
  
}


# run corrupt grow
run_corrupt_grow <- function(project_dir=NULL,input_dir=NULL, output_dir = NULL,newick_path='tree.newick'){
  project_dir <- '/home/htran/Projects/hakwoo_project/corrupt_tree/'
  input_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local_v5/'
  output_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local_v5/corrupt_grow_full/'
  if (!file.exists(output_dir)){
    dir.create(output_dir)
  }
  
  # newick_path <- paste0(input_dir,newick_fn)
  
  corrupt_cnv <- paste0(cg_dir,'filtered_none_cells.csv') 
  none_cells <- read.delim(corrupt_cnv, check.names = FALSE, stringsAsFactors = FALSE, sep = ",")
  # View(none_cells[1:3,1:3])
  corrupt_grow_path <- '/home/htran/storage/install_software/nowellpack/build/install/nowellpack/bin/'
  
  # run with binary file, specify the fpr, fnr
  # cmd_str <- sprintf('%scorrupt-grow --matrix NoisyBinaryCLMatrix --matrix.binaryMatrix %s --matrix.fpRate 0.01 --matrix.fnRate 0.5 --phylo file %s'
  #                    , corrupt_grow_path, corrupt_cnv, newick_path)
  # cmd_str <- sprintf('%scorrupt-grow --matrix NoisyBinaryCLMatrix --matrix.binaryMatrix %s --matrix.fpRate 0.05 --matrix.fnRate 0.5 --phylo file %s'
  #                    , corrupt_grow_path, corrupt_cnv, newick_path)
  
  cmd_str <- sprintf('%scorrupt-grow --matrix NoisyBinaryCLMatrix --matrix.binaryMatrix %s --matrix.fpRate 0.1 --matrix.fnRate 0.2 --phylo file %s'
                     , corrupt_grow_path, corrupt_cnv, newick_path)
  # run with general parameter
  # cmd_str <- sprintf('%scorrupt-grow --matrix ReadOnlyCLMatrix %s --phylo file %s'
  #                    , corrupt_grow_path, corrupt_cnv, newick_path)
  # cmd_str <- gsub('corrupt-grow', corrupt_grow_path, cmd_str)
  print(cmd_str)
  # execute_grow <- system(cmd_str, intern = T)
  # print(execute_grow)
  print("Completed, Voila!!!")
  
  
  corrupt_cnv_df <- read.delim(output_file, check.names = FALSE, stringsAsFactors = FALSE, sep = ",")
  length(tree$node.label)
  locus_ls_3 <- unique(corrupt_cnv_df$loci)
  length(locus_ls_3)
  
  padding1 <- paste0(input_dir, 'bin_cnvs_corrupt_double_padding.csv')  
  padding1_df <- read.delim(padding1, check.names = FALSE, stringsAsFactors = FALSE, sep = ",")
  length(padding1_df$cells)
  
  fcp <- paste0(input_dir, 'total_merged_filtered_states.csv')  
  cp_filtered <- read.csv(fcp, header=T, row.names=1, check.names = F,stringsAsFactors = FALSE)
  
  locus_ls2 <- unique(padding1_df$loci)
  sum(locus_ls2==locus_ls)
  length(locus_ls2)
  setdiff(locus_ls_3, locus_ls2)
  intersect <- intersect(locus_ls_3, locus_ls2)
  length(intersect)
}


