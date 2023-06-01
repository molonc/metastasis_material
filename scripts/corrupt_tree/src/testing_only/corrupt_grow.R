corrupt_grow <- function(datatag = NULL, newick_path, padded, post_fix = '', mat = NULL, out_dir = NULL) {
  # The newick file 
  tree <- read.newick(newick_path)
  # Get all the damn loci it has
  all_loci <- c(grep('locus_', tree$tip.label, value = T), grep('locus_', tree$node.label, value = T))
  stopifnot(any(!grep('locus_', all_loci)) == F)
  all_cells <- gsub('cell_', '', grep('cell', tree$tip.label, value = T))
  # Cell that have to be added
  pad_tag <- ifelse(!padded, '_no_padding', '')
  if (is.null(datatag)) {
    delta_mat <- mat
  } else {
    delta_mat <- as.data.frame(fread(sprintf('~/Desktop/SC-1311/%s/processed_data/bincount/%s_cnvs_corrupt%s%s.csv', datatag, tolower(datatag), pad_tag, post_fix), stringsAsFactors = F))
  }
  # Remove cells that are already in the tree
  delta_mat <- delta_mat[!(delta_mat$cells %in% all_cells), ]
  # Keep the shared loci 
  delta_mat$loci <- paste0('locus_', delta_mat$loci)
  # If all tree loci are not present in the file, this is the wrong condition
  stopifnot(all_loci %in% unique(delta_mat$loci))
  delta_mat <- delta_mat[delta_mat$loci %in% all_loci, ]
  stopifnot(all(sort(all_loci) == sort(unique(delta_mat$loci))))
  if (is.null(out_dir)) {
    out_dir <- sprintf('~/Desktop/SC-1311/%s/processed_data/bincount/grow_tree_temp', datatag)
    dir.create(out_dir, showWarnings = F)
    out_path <- file.path(out_dir, paste0(tolower(datatag), '_cnvs_corrupt_', basename(file_path_sans_ext(newick_path)), '.csv'))
  } else {
    out_path <- file.path(out_dir, paste0('cnvs_corrupt_', basename(file_path_sans_ext(newick_path)), '.csv'))
  }
  fwrite(delta_mat, out_path, quote = F, row.names = F)
  # Construct the corrupt-grow command
  cmd_str <- sprintf('corrupt-grow --matrix NoisyBinaryCLMatrix --matrix.binaryMatrix %s --matrix.fpRate 0.1 --matrix.fnRate 0.2 --phylo file %s', out_path, newick_path)
  corrupt_grow_path <- '/Users/sohrabsalehi/projects/nowellpack/build/install/nowellpack/bin/corrupt-grow'
  cmd_str <- gsub('corrupt-grow', corrupt_grow_path, cmd_str)
  print(cmd_str)
  oo <- system(cmd_str, intern = T)
  print(oo)
  result_dir <- yaml.load(oo[2])$outputFolder
  result_dir
}