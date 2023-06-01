# Tyler left, right padding version, SA919
# results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local_tyler_left_padding/'
# plot_title <- 'tyler_version'
# cellclones <- paste0(results_dir, 'cell_clones.csv')

# # Sohrab's version, SA919
# results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local/'
# plot_title <- 'sohrab_version'
# cellclones <- paste0(results_dir,'tree_cut_out/cell_clones_if_0.02_af_0.75_p0.75_e0.04.csv')


source_dir <- '/home/htran/Projects/farhia_project/rscript/dlp/visualize_tree/'
source(paste0(source_dir,'trajectory_utils.R'))





# Cut the root short
generate_tree <- function(results_dir,datatag = 'SA919_Tyler') {
  save_dir <- paste0(results_dir,'tree_viz_dream/')
  if (!file.exists(save_dir)){
    dir.create(save_dir)
  }
  
  # Load tree
  newick <- paste0(results_dir, 'tree.newick')
  tree <- read.tree(newick)
  edge_list <- tree_2_edge_list(tree)
  # View(head(edge_list))
  # edge_list_path <- paste0(save_dir,'edges_SA919.csv')
  # write.csv(edge_list,file = edge_list_path,quote=F, row.names = F)
  # ss <- as.matrix(read.csv(edge_list_path, stringsAsFactors = F, header=T))
  g <- read_ltm_tree(edge_list_path)
  
  net <- network(as.matrix(edge_list), directed = T)
  v_names <- network.vertex.names(net)
  res <- convert_edge_list_to_ape(edge_list_path)
  tree <- res$tree
  node_names <- res$node_names
  # head(node_names)
  
  #ggtree(tree) + theme_tree2()
  
  # Get height dat
  h <- get_height_dat(g)
  h <- h[order(h$dist, decreasing = T), ]
  
  # Remove all < rm_threshold to shorten the root
  #rm_threshold <- 70
  # rm_threshold <- 280
  # if (datatag == 'SA532') {
  #   rm_threshold <- 86
  # }
  # rm_threshold <- 20 # nb cells > 3855
  rm_threshold <- 30 
  rem <- h$id[h$dist > rm_threshold]
  # bb <- as.data.frame(ss)
  bb <- edge_list
  bb <- bb[bb$source %in% rem & bb$target %in% rem, ]
  
  # Find a locus to re-attach the root
  first_locus <- grep('locus_', rev(rem), value=T)[1]
  bb <- rbind(bb, data.frame(source = 'root', target = first_locus))
  net <- network(as.matrix(bb), directed = T)
  
  
  #ggnet2(net, directed = T, size=2)
  
  v_names <- network.vertex.names(net)
  length(v_names)
 
  node_names <- res$node_names
  dim(node_names)
  rownames(node_names) <- node_names$node_name
  node_names <- node_names[v_names,]
  meta_sample <- get_clones(node_names, cellclones, results_dir)
  meta_sample <- get_timepoint(meta_sample, results_dir)
  colnames(meta_sample)
  colnames(meta_sample)[which(colnames(meta_sample) == "node_name")] <- "id"
  # sum(meta_sample$id==node_names)==length(v_names)
  meta_sample$size <- rep(0.3,length(meta_sample$id))
  locus_ids <- grep('locus_', meta_sample$id)
  none_cells <- meta_sample$id[meta_sample$clone_id=='None']
  none_cells <- none_cells[!none_cells %in% grep('locus_', none_cells, value=T)]
  meta_sample[none_cells,'size'] <- 0.1
  meta_sample[locus_ids,'size'] <- 0
  # Process root
  rid <- grep('root',meta_sample$id)
  
  meta_sample[rid,'size'] <- 0.4
  meta_sample[rid,'clone_id'] <- 'Root'
  meta_sample[rid,'color'] <- '#FF0000'  
  meta_sample[rid,'color_group'] <- '#4C0000'
  meta_sample[rid,'mainsite'] <- 'Root'
  
  filtered_g <- graph_from_edgelist(as.matrix(bb))
  layoutkk <- fitclone_set_dream_tree_layout(net, filtered_g)
  plot_dream_tree(net, meta_sample, layoutkk, datatag)

}



