

# Calculate FDR values 
# Get comparison files 
# Load significant_genes csv files 
# Divide into in-cis, in-trans genes 
# Bootstrap function 
# Record confident interval and pvalue 
# Do it for 4 different mapping
# Save output to 1 csv file 

library(dplyr)
ref_dif <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/Symbol_ensembl.csv'

# yaxis: logFC values 
# xaxis : log mean expres sctransform A - B 


pair_groups_v2 <- pair_groups
pair_groups_v2$datatag <- ifelse(grepl('X',pair_groups_v2$desc),'SA535_CX5461',
                                 ifelse(!grepl('X',pair_groups_v2$desc),paste0(pair_groups_v2$datatag,'_CISPLATIN'),pair_groups_v2$datatag))

ref_genes_desc <- 'CF' # 553 genes 
pair_groups_v2 <- get_confident_interval(ref_genes_desc, pair_groups_v2)
cf_incis <- pair_groups_v2 %>%
  dplyr::group_by(datatag) %>%
  dplyr::summarise(avg_incis = round(mean(CF_incis_pval),4),
                   sd_incis=round(sd(CF_incis_pval),4),
                   nb_sigf = sum(CF_incis_pval<0.05)) %>%
  ungroup()

# datatag          avg_incis sd_incis nb_sigf
# * <chr>                <dbl>    <dbl>   <int>
#   1 SA1035_CISPLATIN    0.245    0.291        2  need attention
# 2 SA535_CISPLATIN     0.0022   0.0059       7
# 3 SA535_CX5461        0.0004   0.0007       6
# 4 SA609_CISPLATIN     0        0            6

unique(pair_groups_v2$datatag)
pair_groups_v2[pair_groups_v2$datatag=='SA1035_CISPLATIN' & pair_groups_v2$CF_incis_pval>=0.05,'desc']

cf_intrans <- pair_groups_v2 %>%
  dplyr::group_by(datatag) %>%
  dplyr::summarise(avg_intrans = round(mean(CF_intrans_pval),4),
                   sd_intrans=round(sd(CF_intrans_pval),4),
                   nb_sigf = sum(CF_intrans_pval<0.05)) %>%
  ungroup()

# datatag          avg_intrans sd_intrans nb_sigf
# * <chr>                  <dbl>      <dbl>   <int>
#   1 SA1035_CISPLATIN       0          0           4
# 2 SA535_CISPLATIN        0.038      0.100       6, one case is not significant
# 3 SA535_CX5461           0          0           6
# 4 SA609_CISPLATIN        0          0           6
pair_groups_v2[pair_groups_v2$datatag=='SA535_CISPLATIN' & pair_groups_v2$CF_intrans_pval>=0.05,c('CF_intrans_pval','desc')]

ref_genes_desc <- 'BS' # 983 genes 
pair_groups_v2 <- get_confident_interval(ref_genes_desc, pair_groups_v2)
cf_incis_bs <- pair_groups_v2 %>%
  dplyr::group_by(datatag) %>%
  dplyr::summarise(avg_incis = round(mean(BS_incis_pval),4),
                   sd_incis=round(sd(BS_incis_pval),4),
                   nb_sigf = sum(BS_incis_pval<0.05)) %>%
  ungroup()
# datatag          avg_incis sd_incis nb_sigf
# * <chr>                <dbl>    <dbl>   <int>
#   1 SA1035_CISPLATIN    0.0011   0.0022       4
# 2 SA535_CISPLATIN     0        0            7
# 3 SA535_CX5461        0        0            6
# 4 SA609_CISPLATIN     0        0            6

cf_intrans_bs <- pair_groups_v2 %>%
  dplyr::group_by(datatag) %>%
  dplyr::summarise(avg_intrans = round(mean(BS_intrans_pval),4),
                   sd_intrans=round(sd(BS_intrans_pval),4),
                   nb_sigf = sum(BS_intrans_pval<0.05)) %>%
  ungroup()
# datatag          avg_intrans sd_intrans nb_sigf
# * <chr>                  <dbl>      <dbl>   <int>
#   1 SA1035_CISPLATIN           0          0       4
# 2 SA535_CISPLATIN            0          0       7
# 3 SA535_CX5461               0          0       6
# 4 SA609_CISPLATIN            0          0       6

ref_genes_desc <- 'COSMIC'
#126 cosmic genes
pair_groups_v2 <- get_confident_interval(ref_genes_desc, pair_groups_v2)
cf_incis_cosmic <- pair_groups_v2 %>%
  dplyr::group_by(datatag) %>%
  dplyr::summarise(avg_incis = round(mean(COSMIC_incis_pval),4),
                   sd_incis=round(sd(COSMIC_incis_pval),4),
                   nb_sigf = sum(COSMIC_incis_pval<0.05)) %>%
  ungroup()

# datatag          avg_incis sd_incis nb_sigf
# * <chr>                <dbl>    <dbl>   <int>
#   1 SA1035_CISPLATIN     0.832   0.178        0
# 2 SA535_CISPLATIN      0.718   0.194        0
# 3 SA535_CX5461         0.904   0.101        0
# 4 SA609_CISPLATIN      0.975   0.0317       0


cf_intrans_cosmic <- pair_groups_v2 %>%
  dplyr::group_by(datatag) %>%
  dplyr::summarise(avg_intrans = round(mean(COSMIC_intrans_pval),4),
                   sd_intrans=round(sd(COSMIC_intrans_pval),4),
                   nb_sigf = sum(COSMIC_intrans_pval<0.05)) %>%
  ungroup()

# datatag          avg_intrans sd_intrans nb_sigf
# * <chr>                  <dbl>      <dbl>   <int>
#   1 SA1035_CISPLATIN       1.00      0.0005       0
# 2 SA535_CISPLATIN        0.953     0.0957       0
# 3 SA535_CX5461           0.958     0.0758       0
# 4 SA609_CISPLATIN        0.982     0.0362       0

ref_genes_desc <- 'cis_res'  # 91 genes
pair_groups_v2 <- get_confident_interval(ref_genes_desc, pair_groups_v2)
cf_incis_cisr <- pair_groups_v2 %>%
  dplyr::group_by(datatag) %>%
  dplyr::summarise(avg_incis = round(mean(cis_res_incis_pval),4),
                   sd_incis=round(sd(cis_res_incis_pval),4),
                   nb_sigf = sum(cis_res_incis_pval<0.05)) %>%
  ungroup()

cf_intrans_cisr <- pair_groups_v2 %>%
  dplyr::group_by(datatag) %>%
  dplyr::summarise(avg_intrans = round(mean(cis_res_intrans_pval),4),
                   sd_intrans=round(sd(cis_res_intrans_pval),4),
                   nb_sigf = sum(cis_res_intrans_pval<0.05)) %>%
  ungroup()


# CF 553 genes z

# Objective: Fig 5.9 Gene set membership cis and trans differential expression between resistant and sensitive clones, 
# whether obtained genes set membership follow random organization - null hypothesis or significantly correlated with reference genes set? 
# We use null hypothesis test to compare our observed DE genes list with 1000 random populations of genes with same sample size from total 33514 human symbol genes. 
# We reject null hypothesis in case our observed data are significant different from random organization (the confidence level 
# of our observed genes is greater than 95% confidence interval - reflect a significant level of p-value is less than 0.05)

# Genes membership - mapped to 553 PanCancer Core Fitness reference genes 
# almost tests are significant, except the case of SA1035: 2 tests are significant:  Res(X7:H) vs Sen(X7:E), Res(X8:H) vs Sen(X8:E) with p-value: 0.01, 0.001, 
# and 2 cases we can not reject null hypothesis: Res(X7:G) vs Sen(X7:E), Res(X7:G) vs Sen(X7:E) with p-value are
# 0.37, 0.59. The percentage of in cis genes in these 2 cases are lower than 10% of total DE genes, so it is easy to explain output of null hypothesis tests. 


# Genes membership - mapped to 983 Broad Sanger reference genes: 
# In cis genes: all DE comparisons are significant different from null hypothesis with 95% confidence interval. 
# In trans genes: all DE comparisons are significant different from null hypothesis with 95% confidence interval. 


# Genes membership - mapped to 126 cosmic reference genes: 
# In cis genes: can not reject null hypothesis with 95% confidence interval for all null hypothesis tests in 3 series of data. 
# In trans genes: can not reject null hypothesis with 95% confidence interval for null hypothesis tests in 3 series of data. 


# Genes membership - mapped to 91 cisplatin resistance reference genes: 
# In cis genes: can not reject null hypothesis with 95% confidence interval for null hypothesis tests in 3 series of data. 
# In trans genes: can not reject null hypothesis with 95% confidence interval for null hypothesis tests in 3 series of data. 



# Genes membership - mapped to 553 PanCancer Core Fitness reference genes 
# In cis genes: almost tests are significant with p-value < 0.05 (SA609: mean of p-value: 0, 
#                                                                 SA535: Cisplatin: mean of p-value:0.002 ,SA535:CX5461: mean p-value: 0) 
# .Except the case of SA1035: 2 tests are significant:  Res(X7:H) vs Sen(X7:E), Res(X8:H) vs Sen(X8:E) with p-value: 0.01, 0.001, 
# and 2 cases we can not reject null hypothesis: Res(X7:G) vs Sen(X7:E), Res(X7:G) vs Sen(X7:E) with p-value are
# 0.37, 0.59. The percentage of in cis genes in these 2 cases are lower than 10% of total DE genes, so it is easy to explain output of null hypothesis tests. 

# In trans genes: almost tests are significant with p-value < 0.05 (SA609: mean of p-value: 0, 
                                                                  # SA1035: mean of p-value: 0, 
#                                                                 SA535: Cisplatin: mean of p-value:0.038 ,
                                                                  # SA535:CX5461: mean p-value: 0) 
# .Except the case of SA535:Cisplatin: SA535_UUTTT_R_UUUUU_J with p-value: 0.266



# Genes membership - mapped to 983 Broad Sanger reference genes: 
# In cis genes: all DE comparisons are significant different from null hypothesis with 95% confidence interval. 
# (SA609: mean of p-value: 0, 
  # SA1035: mean of p-value: 0.0011, 
  # SA535: Cisplatin: mean of p-value:0 ,
  # SA535:CX5461: mean p-value: 0) 

# In trans genes: all DE comparisons are significant different from null hypothesis with 95% confidence interval. 
# (SA609: mean of p-value: 0, 
# SA1035: mean of p-value: 0, 
# SA535: Cisplatin: mean of p-value:0 ,
# SA535:CX5461: mean p-value: 0) 


# Genes membership - mapped to 126 cosmic reference genes: 
# In cis genes: can not reject null hypothesis with 95% confidence interval for all permutation tests in 3 series of data. 
# Gene membership are in 90% of confidence interval for series SA609, SA535:CX5461. And 70% of confidence interval for series: SA1035, SA535: Cisplatin
# (SA609: mean of p-value: 0.97, 
# SA1035: mean of p-value: 0.83, 
# SA535: Cisplatin: mean of p-value:0.71 ,
# SA535:CX5461: mean p-value: 0.90) 

# In trans genes: can not reject null hypothesis with 95% confidence interval for all permutation tests in 3 series of data. 
# Gene membership are in 90% of confidence interval for series SA609, SA1035, SA535:Cisplatin, SA535:CX5461.
# (SA609: mean of p-value: 0.98, 
# SA1035: mean of p-value: 1.0, 
# SA535: Cisplatin: mean of p-value:0.95 ,
# SA535:CX5461: mean p-value: 0.95) 


# Genes membership - mapped to 91 cisplatin resistance reference genes: 
# In cis genes: can not reject null hypothesis with 95% confidence interval for all permutation tests in 3 series of data. 
# Gene membership are in 90% of confidence interval for series SA609, SA1035, SA535:Cisplatin, SA535:CX5461.
# (SA609: mean of p-value: 0.98, 
# SA1035: mean of p-value: 0.94, 
# SA535: Cisplatin: mean of p-value:0.99 ,
# SA535:CX5461: mean p-value: 1) 

# In trans genes: can not reject null hypothesis with 95% confidence interval for all permutation tests in 3 series of data. 
# (SA609: mean of p-value: 1.0, 
# SA1035: mean of p-value: 1.0, 
# SA535: Cisplatin: mean of p-value:0.96 ,
# SA535:CX5461: mean p-value: 1.0) 


get_confident_interval <- function(){
  ## get DE analysis, and list of cis, trans genes
  # get number of cis genes overlapping with genome regions
  # get random genes, and check how many genes fall into cis category 
  # get bootstrap values 
  
}


get_confident_interval <- function(ref_genes_desc, pair_groups=NULL){
  ref_dif <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/'
  base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
  input_dir <- paste0(base_dir,'SA535_total_rna_v2/SA535-v6/')
  if(is.null(pair_groups)){
    pair_groups_fn <- paste0(input_dir,'comparisons_drug_res.csv')
    id609 <- paste0(base_dir,'SA609_rna/deg_analysis/SA609-v6/')
    id1035 <- paste0(base_dir,'SA1035_rna/deg_analysis/SA1035-v6/')
    id535 <- paste0(base_dir,'SA535_total_rna_v2/SA535-v6/')
    pair_groups <- read.csv(pair_groups_fn, header=T, check.names=F, stringsAsFactors=F)
    # View(pair_groups)
    print(dim(pair_groups))
    # length(unique(pair_groups$desc))
    # unique(pair_groups$datatag)
    pair_groups$result_dir <- ifelse(pair_groups$datatag=='SA609',paste0(id609,pair_groups$desc,'/'),
                                     ifelse(pair_groups$datatag=="SA1035",paste0(id1035,pair_groups$desc,'/'),
                                            paste0(id535,pair_groups$desc,'/')))
    
  }
  
  if(ref_genes_desc=='BS'){ # Broad Sanger
    bs_genes <- read.csv(paste0(ref_dif, '41467_2019_13805_MOESM6_ESM_integrated_broad_sanger_essential_genes.csv'))
    ess_genes <- intersect(bs_genes$Broad, bs_genes$Sanger)
    ref_set <- unique(ess_genes)
  }else if(ref_genes_desc=='cis_res'){ # cisplatin resistance genes
    ref_genes <- read.table(paste0(ref_dif, 'cisplatin_resistance_genes.txt'),sep='\t',header = T, check.names = F, stringsAsFactors = F)
    ref_set <- unique(ref_genes$gene_symbol)
  }else if(ref_genes_desc=='COSMIC'){
    ref_genes <- read.table(paste0(ref_dif, 'oncogene_cosmic.txt'), sep = '\t', check.names = F, stringsAsFactors = F)
    ref_set <- unique(ref_genes$V1)
  }else if(ref_genes_desc=='CF'){ # core fitness genes 
    cancer_ref_genes_df <- read.csv(paste0(ref_dif,'Behan_CFgenes.csv'), stringsAsFactors=F, check.names = F)
    colnames(cancer_ref_genes_df)[which(colnames(cancer_ref_genes_df) == "ADAM PanCancer Core-Fitness genes")] <- "PanCancer_Fitness_genes"
    ref_set <- unique(cancer_ref_genes_df$PanCancer_Fitness_genes)
  }else{
    stop('Input reference mapping list!!!')
  }
  print(ref_genes_desc)
  print(length(ref_set))
  
  
  genome_genes_df <- read.csv(paste0(ref_dif, 'Symbol_ensembl.csv'), check.names = F, stringsAsFactors = F)
  # dim(genome_genes_df)
  genome_genes <- genome_genes_df$Symbol
  genome_genes <- unique(genome_genes)
  print(length(genome_genes))
  
  
  minLogFC <- 0.25
  FDR_cutoff <- 0.01
  pValueThrs <- 0.05
  
  # Mirela just added new group
  pair_groups <- pair_groups %>%
    dplyr::filter(desc != 'SA535_UUTTTTT_S_T_UXXXX_U')
  rownames(pair_groups) <- pair_groups$desc
  # print(dim(pair_groups))
  pair_groups[,paste0(ref_genes_desc,'_incis_pval')] <- NA
  pair_groups[,paste0(ref_genes_desc,'_intrans_pval')] <- NA
  for(de in pair_groups$desc){
    signif_genes <- read.csv(paste0(pair_groups[de,'result_dir'],'signif_genes.csv'), check.names = F, stringsAsFactors=F)
    signif_genes <- signif_genes[abs(signif_genes$logFC)>minLogFC & signif_genes$FDR<FDR_cutoff  & signif_genes$PValue<pValueThrs,]
    # signif_genes$is_broad_sanger_gene <- signif_genes$gene_symbol %in% ess_genes
    print(dim(signif_genes))
    
    incis <- signif_genes %>% 
      dplyr::filter(!is.na(classified_gene_dlp)) %>% 
      dplyr::pull(gene_symbol)
    
    intrans <- signif_genes %>% 
      dplyr::filter(is.na(classified_gene_dlp)) %>% 
      dplyr::pull(gene_symbol)
    
    
    print(length(incis))
    print(length(intrans))
    
    res_cis <- get_bootstrap_stat(incis, ref_set, genome_genes) #nsamples=1000
    res_trans <- get_bootstrap_stat(intrans, ref_set, genome_genes) #nsamples=1000
    pair_groups[de,paste0(ref_genes_desc,'_incis_pval')] <- res_cis$pval
    pair_groups[de,paste0(ref_genes_desc,'_intrans_pval')] <- res_trans$pval
    
  }
  print(ref_genes_desc)
  print('incis genes')
  print(summary(pair_groups[,paste0(ref_genes_desc,'_incis_pval')]))
  
  print('intrans genes')
  print(summary(pair_groups[,paste0(ref_genes_desc,'_intrans_pval')]))
  
  return(pair_groups)
}


count_occurance <- function(s){
  return(sum(s %in% ref_set))
}

get_bootstrap_stat <- function(our_obs, ref_set, genome_genes, nsamples=1000){
  resamples <- lapply(1:nsamples, function(i) sample(genome_genes, size=length(our_obs), replace = T))
  
  our <- sum(our_obs %in% ref_set)
  occurs <- sapply(resamples, count_occurance)
  names(occurs) <- paste0('R',seq(1:length(occurs)))
  occurs['our'] <- our
  # length(occurs)
  r <- rank(occurs) #,ties.method = "max"
  pval <- (nsamples + 1 - r['our'])/nsamples
  return(list(CI=r['our'],pval=pval))
}




