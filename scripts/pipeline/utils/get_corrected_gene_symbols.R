library(dplyr)
# Issue: gene symbols that biologist provide can contain alias, not real symbol that used in computation
# ref_df : dataframe with a colname 'gene_symb'
# obs_genes: the list of uncorrected genes
# reference_set_name: name that you want to save output as

get_corrected_gene_symbol <- function(ref_df, reference_set_name='', obs_genes=NULL){
  ref_dir <- 'yourdir where you keep the mapping file/'
  database_genes <- data.table::fread(paste0(ref_dir, 'custom_ens_genes_symbols.txt'), sep='\t') %>% as.data.frame()
  # dim(database_genes)
  database_genes$corrected_gene_symb <- database_genes$`Approved symbol`
  database_genes$previous_symbols <- database_genes$`Previous symbols`
  database_genes$alias <- database_genes$`Alias symbols`
  
  meta_genes <- data.table::fread('yourdir/Symbol_ensembl.csv') %>% as.data.frame() # your set of common scRNAseq genes
  ## Or using the database
  # meta_genes <- annotables::grch38 %>% 
  #   dplyr::select(gene_symb = symbol)
  # meta_genes <- annotables::grch37 %>% 
  #   dplyr::select(gene_symb = symbol)
  
  
  if(is.null(obs_genes)){
    obs_genes <- ref_df %>%
      dplyr::filter(!gene_symb %in% meta_genes$Symbol) %>%
      dplyr::pull(gene_symb)
    obs_genes <- unique(obs_genes)
  }
  
  anno <- tibble::tibble()
  for(t in obs_genes){
    # print(t)
    re <- grepl(toupper(t),database_genes$alias)
    if(sum(re==TRUE)>0){
      tmp <- database_genes[re,c('corrected_gene_symb','previous_symbols','alias')]
      tmp$gene_symb <- t
      anno <- dplyr::bind_rows(anno, tmp)
    }
    re1 <- grepl(toupper(t),database_genes$previous_symbols)
    if(sum(re1==TRUE)>0){
      tmp1 <- database_genes[re1,c('corrected_gene_symb','previous_symbols','alias')]
      tmp1$gene_symb <- t
      anno <- dplyr::bind_rows(anno, tmp1)
    }
  }
  out <- tibble::tibble()
  for(i in seq(1:dim(anno)[1])){
    # gs <- unlist(strsplit(anno$previous_symbols[i],', '))
    gs <- unlist(strsplit(anno$alias[i],', '))
    if(length(gs)>0){
      tmp <- tibble::tibble(corrected_gene_symb=rep(anno$corrected_gene_symb[i],length(gs)),gene_symb=gs)
      out <- dplyr::bind_rows(out, tmp)  
    }
    
  }
  for(i in seq(1:dim(anno)[1])){
    # gs <- unlist(strsplit(anno$previous_symbols[i],', '))
    gs <- unlist(strsplit(anno$previous_symbols[i],', '))
    if(length(gs)>0){
      tmp <- tibble::tibble(corrected_gene_symb=rep(anno$corrected_gene_symb[i],length(gs)),gene_symb=gs)
      out <- dplyr::bind_rows(out, tmp)  
    }
    
  }
  print(dim(out))
  out <- out %>%
    dplyr::filter(gene_symb %in% toupper(obs_genes))
  # correct_mapping <- summary(as.factor(out$gene_symb))
  # correct_mapping <- correct_mapping[correct_mapping==1] # only one mapping between alias and correct gene name
  # out1 <- out %>%
  #   dplyr::filter(gene_symb %in% names(correct_mapping))
  # out2 <- out %>%
  #   dplyr::filter(!gene_symb %in% names(correct_mapping))
  # out1
  ref_df <- ref_df %>%
    dplyr::filter(!gene_symb %in% out$gene_symb)
  # dim(ref_df)
  out <- out %>%
    dplyr::select(corrected_gene_symb) %>%
    dplyr::rename(gene_symb=corrected_gene_symb)
  ref_df <- dplyr::bind_rows(ref_df, out)
  dim(ref_df)
  if(reference_set_name!=''){
    data.table::fwrite(ref_df, paste0(ref_dir, reference_set_name,'.csv'))  
  }
  
  return(ref_df)
  
}