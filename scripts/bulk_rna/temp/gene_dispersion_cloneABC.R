## Full version: dddd


main(){
  
  # From gene_dispersion_cloneA.R
  load_input_data()
  
  # ## Loading cnv file
  obs_clones <- c('A','B')
  cnv <- get_cis_trans_gene_type(obs_clones)
  dim(cnv)
}