






get_stat_test <- function(){
  
  # Get list of bulk RNA-seq genes: 19K genes: genome genes
  meta_genes <- data.table::fread(paste0(script_dir, 'SA919_DE_analysis_DESeq2_Hoa_09April2024/meta_genes_SA919.csv.gz'))  
  
  # Get # cis, trans genes for each DE comparison
  
  # Get bootstrap test for number of occurrence
  
}
get_bootstrap_stat <- function(our_obs, ref_set, genome_genes, nsamples=1000){
  resamples <- lapply(1:nsamples, function(i) sample(genome_genes, size=length(our_obs), replace = T))
  
  our <- sum(our_obs %in% ref_set)
  # occurs <- sapply(resamples, count_occurance)
  occurs <- sapply(resamples, function(s) {return(sum(s %in% ref_set))})
  names(occurs) <- paste0('R',seq(1:length(occurs)))
  occurs['our'] <- our
  # length(occurs)
  r <- rank(occurs) #,ties.method = "max"
  pval <- (nsamples + 1 - r['our'])/nsamples
  return(list(CI=r['our'],pval=pval))
}


get_bootstrap_stat <- function(our_obs, ref_set, genome_genes, nsamples=1000){
  resamples <- lapply(1:nsamples, function(i) sample(genome_genes, size=length(our_obs), replace = T))
  
  our <- sum(our_obs %in% ref_set)
  # occurs <- sapply(resamples, count_occurance)
  occurs <- sapply(resamples, function(s) {return(sum(s %in% ref_set))})
  names(occurs) <- paste0('R',seq(1:length(occurs)))
  occurs['our'] <- our
  # length(occurs)
  r <- rank(occurs) #,ties.method = "max"
  pval <- (nsamples + 1 - r['our'])/nsamples
  return(list(CI=r['our'],pval=pval))
}

