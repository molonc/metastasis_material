# Hoa Tran
# Put together the subset of cells, lets call them s-cells
# 1. generate corrupt-input matrix (the binary one) for s-cells
# 2. run s-cells through the jitter fix tool
# 3. grab filtered.csv from the jitter-fixed output output.csv if you the .nf file I sent you
# 4. give output.csv to corrupt-grow tool

# 5874 - 4895
# [1] 979 grown cells

# 1. generate corrupt-input matrix (the binary one) for s-cells
get_corrupt_mtx <- function(project_dir=NULL, input_dir=NULL, output_dir = NULL,
                            newick_fn='tree.newick'){
  
  source(paste0(project_dir, "src/testing_only/blob_v2.R"))
  
  # tree <- read.newick(paste0(output_dir,'grown.newick'))
  # The newick file 
  tree <- read.newick(paste0(input_dir,newick_fn))
  length(tree$node.label) # 65
  length(tree$tip.label)
  # Get all the damn loci it has
  all_loci <- c(grep('locus_', tree$tip.label, value = T), grep('locus_', tree$node.label, value = T))
  length(all_loci) # 64
  stopifnot(any(!grep('locus_', all_loci)) == F)
  all_cells <- gsub('cell_', '', grep('cell', tree$tip.label, value = T))
  
  # Read original copy number mtx and get outlier cells
  copynumber <- paste0(input_dir, 'total_merged_filtered_states_original.csv')  
  copy_number <- read.delim(copynumber, check.names = FALSE, stringsAsFactors = FALSE, sep = ",")
  out_fn <- paste0(input_dir, 'corrupt_grow/total_merged_filtered_states_original.csv')
  write.csv(copy_number, out_fn,  row.names = T, quote=F) #get dir
  dim(copy_number)
  cnv_mat <- copy_number[,!(colnames(copy_number) %in% all_cells)]
  dim(cnv_mat)
  write.csv(cnv_mat, paste0(output_dir,'total_merged_filtered_states_outlier_cells.csv'),  row.names = T, quote=F) #get dir
  
  # Convert copy number to binary using new encoding scheme
  output_file <- paste0(output_dir,'bin_cnvs_corrupt_double_padding_outlier_cells.csv')
  pad_mod_cnv_mat(cnv_mat, output_file, pad_left = TRUE)
  if (file.exists(output_file)){
    padding_fn <- output_file
  } else{
    padding_fn <- NULL
  }
  
  return(output_file)
}

run_jitter_fix <- function(project_dir=NULL, input_dir=NULL, output_dir = NULL,newick_fn='tree.newick'){
  
  padding_fn <- get_corrupt_mtx(project_dir, input_dir, output_dir, newick_fn)
  padding_fn <- paste0(output_dir,'bin_cnvs_corrupt_double_padding_outlier_cells.csv')
  
  main_padding_df <- read.delim(paste0(input_dir,'bin_cnvs_corrupt_double_padding.csv'), 
                                check.names = FALSE, stringsAsFactors = FALSE, sep = ",")
  cells_clone <- read.delim(paste0(input_dir, 'cell_clones.csv'), 
                            stringsAsFactors = FALSE, sep=",")
  none_padding_df <- main_padding_df[main_padding_df$cells %ni% cells_clone$cell_id,]
  dim(none_padding_df)
  length(unique(none_padding_df$cells))
  unique(cells_clone$clone_id)
  data.table::fwrite(none_padding_df, 
                     paste0(output_dir,'bin_cnvs_corrupt_double_padding_none_cells.csv'), 
                     row.names = F, quote = F)  
  
  outlier_df <- read.delim(paste0(input_dir,'corrupt_grow/bin_cnvs_corrupt_double_padding_outlier_cells.csv'), 
                                check.names = FALSE, stringsAsFactors = FALSE, sep = ",")
  length(unique(outlier_df$cells))
  total_unassign <- rbind(none_padding_df,outlier_df)
  total_unassign_fn <- paste0(output_dir,'bin_cnvs_corrupt_double_padding_total_unassign_cells.csv')
  data.table::fwrite(total_unassign, total_unassign_fn, 
                     row.names = F, quote = F)  
  
  
  nf_dir <- '/home/htran/storage/install_software/corrupt-nextflow/'
  cmd_str <- sprintf('%snextflow run %sbinary-infer-pipeline-short-jitter-no-filter.nf -resume --tipInclusionProbabilities %s', 
                     nf_dir, nf_dir, total_unassign_fn)
  print(cmd_str)
  # excute_oo <- system2(cmd_str)
  # print(excute_oo)
  # current_dir <- getwd()
  out_jitter_dir <- sprintf('deliverables/binary-infer-pipeline-short-jitter-no-filter/')
  
  cmd_str2 <- sprintf('mv %soutput.csv %sbin_cnvs_corrupt_double_padding_jitter.csv',
                      out_jitter_dir, output_dir)
  print(cmd_str2)
  # excute_oo2 <- system(cmd_str2, intern = T)
  # print(excute_oo2)
  padding_jitter <- paste0(output_dir,'bin_cnvs_corrupt_double_padding_jitter.csv')
  padding_jitter_filtered <- paste0(output_dir,'bin_cnvs_corrupt_double_padding_jitter_filtered.csv')
  # jitter_df <- read.delim(padding_jitter_filtered, check.names = FALSE, stringsAsFactors = FALSE, sep = ",")
  # sum(unique(jitter_df$loci) %in% all_loci)
  # length(unique(jitter_df$loci))
  
  if(file.exists(padding_jitter)){
    padding_jitter_df <- read.delim(padding_jitter, check.names = FALSE, stringsAsFactors = FALSE, sep = ",")
    length(unique(padding_jitter_df$loci))
    # pos <- grep('locus_', padding_jitter_df$loci, value = F)
    # padding_jitter_df[pos,'loci'] <- gsub('locus_', '',padding_jitter_df[pos,'loci'])
    sum(unique(padding_jitter_df$loci) %in% all_loci)
    if(sum(unique(padding_jitter_df$loci) %in% all_loci)==length(all_loci)){
      padding_jitter_df <- padding_jitter_df[padding_jitter_df$loci %in% all_loci,]
      padding_jitter_df$loci[1:3]
      all_loci[1:3]
      dim(padding_jitter_df)
      if(!is.null(dim(padding_jitter_df))){
        data.table::fwrite(padding_jitter_df, padding_jitter_filtered, row.names = F, quote = F)  
      }
      # ls_cells <- gsub('cell_', '',padding_jitter_df$cells)
      # ls <- ls_cells[ls_cells %in% all_cells]
      # length(unique(ls))
      # test_cell <- 'SA919X3XB08939-A98232B-R48-C25'
      # test_cell %in% padding_jitter_df$cells
      # test_cell %in% all_cells
      # all_cells[1:3]
    }
    return(padding_jitter_filtered)
  } else{
    return(NULL)
  }  
}

results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_filtered/'
corrupt_cnv <- paste0(results_dir,'corrupt_grow/filtered_outliers.csv') 
newick_path <- paste0(results_dir,'corrupt_grow/grown_nonecells.newick')


run_corrupt_grow <- function(project_dir=NULL,input_dir=NULL, output_dir = NULL,newick_fn='tree.newick'){
  project_dir <- '/home/htran/Projects/hakwoo_project/corrupt_tree/'
  input_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local_v5/'
  output_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local_v5/corrupt_grow_full/'
  if (!file.exists(output_dir)){
    dir.create(output_dir)
  }
  corrupt_cnv <- run_jitter_fix(project_dir, input_dir, output_dir, newick_fn)
  newick_path <- paste0(input_dir,newick_fn)
  
  corrupt_cnv <- paste0(output_dir,'bin_cnvs_corrupt_double_padding_jitter_filtered.csv') 
  corrupt_grow_path <- '/home/htran/storage/install_software/nowellpack/build/install/nowellpack/bin/'
  
  # run with binary file, specify the fpr, fnr
  cmd_str <- sprintf('%scorrupt-grow --matrix NoisyBinaryCLMatrix --matrix.binaryMatrix %s --matrix.fpRate 0.01 --matrix.fnRate 0.5 --phylo file %s'
                     , corrupt_grow_path, corrupt_cnv, newick_path)
  cmd_str <- sprintf('%scorrupt-grow --matrix NoisyBinaryCLMatrix --matrix.binaryMatrix %s --matrix.fpRate 0.05 --matrix.fnRate 0.5 --phylo file %s'
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

# fstr <- paste0(input_dir, 'straighten_output.csv') 
# fstr_df <- read.delim(fstr, check.names = FALSE, stringsAsFactors = FALSE, sep = ",")
# View(head(fstr_df))
# fs <- unique(fstr_df$loci)
# 
# gstr <- paste0(output_dir, 'bin_cnvs_corrupt_double_padding_jitter.csv') 
# gstr_df <- read.delim(fstr, check.names = FALSE, stringsAsFactors = FALSE, sep = ",")
# View(head(gstr_df))
# gs <- unique(gstr_df$loci)
# df <- setdiff(fs, gs)
# length(df)
# length(intersect(fs,gs))
# 
# total_straighten <- rbind(fstr_df, gstr_df)
# fn <- paste0(output_dir,'total_straighten_output.csv')
# data.table::fwrite(total_straighten, fn, row.names = F, quote = F)





# filtered_df <- read.delim(paste0(results_dir,'filtered.csv'), check.names = FALSE, stringsAsFactors = FALSE, sep = ",")
# View(head(filtered_df[,1:3]))
# t <- grep('locus_16_', filtered_df$loci, value = T)
# 
# input_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_2252_local_v2/'
# padding1 <- paste0(input_dir, 'bin_cnvs_corrupt_double_padding.csv')  
# padding1_df <- read.delim(padding1, check.names = FALSE, stringsAsFactors = FALSE, sep = ",")
# dim(padding1_df)
# 
# padding2 <- paste0(input_dir,'corrupt_grow/bin_cnvs_corrupt_double_padding_jitter.csv')  
# padding2_df <- read.delim(padding2, check.names = FALSE, stringsAsFactors = FALSE, sep = ",")
# dim(padding2_df)
# View(head(padding2_df))
# pos <- grep('locus_', padding2_df$loci, value = F)
# padding2_df[pos,'loci'] <- gsub('locus_', '',padding2_df[pos,'loci'])
# loci1 <- unique(padding1_df$loci)
# loci2 <- unique(padding2_df$loci)
# setdiff(paste0('locus_',loci1), loci2) 
# length(intersect(loci1, loci2))
# 
# straighten <- paste0(input_dir, 'straighten_output.csv')  
# straighten_df <- read.delim(straighten, check.names = FALSE, stringsAsFactors = FALSE, sep = ",")
# dim(straighten_df)
# setdiff(straighten_df$loci,padding2_df$loci)
# 
# sum(all_loci %in% straighten_df$loci)
# 
# 
