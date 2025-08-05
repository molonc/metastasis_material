suppressPackageStartupMessages({
  library(tidyverse)
  library(SingleCellExperiment)
  # library(scater)
  library(data.table)
  library(methods)
  # library(scran)
  library(parallel)
  # library(feather)
  # library(annotables)
  library(dplyr)
  library(ggplot2)
  # library(snpEnrichment)
  # library(scMerge)
  options(dplyr.summarise.inform = FALSE)
  options(tidyverse.quiet = TRUE)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(IRanges)
})
# initial.options <- commandArgs(trailingOnly = FALSE)
# file.arg.name <- "--file="
# script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
# script.basename <- dirname(script.name)
# # print(script.basename)
# source(paste0(script.basename, "/utils.R"))
# 
# 
# parser <- ArgumentParser(description = "Get reads data for each clone")
# parser$add_argument('--output_file', type = 'character', metavar = 'character',
#                     help="Output path for gene CN table (feather file).")
# parser$add_argument('--cellclones', type = 'character', default=NULL, metavar = 'FILE',
#                     help="cellclones file")
# parser$add_argument('--obs_clones', type = 'character', default=NULL, metavar = 'character',
#                     help="obs_clones")
# parser$add_argument('--input_dir', type = 'character', default=NULL, metavar = 'character',
#                     help="input_dir")
# parser$add_argument('--datatag', type = 'character', default=NULL, metavar = 'character',
#                     help="datatag")
# parser$add_argument('--library_grouping', type = 'character', default=NULL, default='library_groupings.csv', metavar = 'FILE',
#                     help="library_groupings.csv file")
# 
# 
# args <- parser$parse_args()
# # 
# output_file <- args$output_file
# print(output_file)
# 
# print(args$library_grouping)
# print(args$cellclones)
# # print(args$input_dir)
# if(!file.exists(args$segment_file) | !file.exists(args$library_grouping) 
#    | !file.exists(args$cellclones)){
#   stop('Input files do not exist, pls double check')
# }
# cellclones_fn <- args$cellclones
# input_dir <- args$input_dir
# obs_clones <- args$obs_clones

calc_mode <- function(x) {
  keys <- unique(x)
  keys[which.max(tabulate(match(x, keys)))]
}

get_library_id <- function(cell_ids, cores_use=2) {
  
  labels <- mclapply(strsplit(cell_ids, "-"), function(x) {
    return(x[2])
  }, mc.cores = cores_use)
  return(as.character(labels))
  
}
get_gene_coordinates <- function(segments){ #, min_pc_overlap=0.1
  # segments$chr <- paste0('chr',chrs_df$chr)
  # library(org.Hs.eg.db)
  # library(IRanges)
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  g <- genes(txdb, single.strand.genes.only=FALSE)
  
  segments_gr <- makeGRangesFromDataFrame(segments, keep.extra.columns = TRUE)
  overlaps <- findOverlaps(g, segments_gr, ignore.strand = TRUE)
  overlapping_ranges <- pintersect(g[queryHits(overlaps)], segments_gr[subjectHits(overlaps)], ignore.strand = TRUE, drop.nohit.ranges = TRUE)
  
  percentage_overlap <- width(overlapping_ranges) / width(segments_gr[subjectHits(overlaps)])
  # wd_overlap <- width(overlapping_ranges)
  # percentage_overlap <- width(overlapping_ranges) / width(g[queryHits(overlaps)])
  
  percentage_overlap_sum <- sapply(percentage_overlap, sum)
  # wd_overlap_sum <- sapply(wd_overlap, sum)
  gene_cn <- tibble(entrezgene = names(g)[queryHits(overlaps)], percent_overlap=percentage_overlap_sum) #wd_overlap=wd_overlap_sum
  # print(dim(gene_cn))
  # length(unique(gene_cn$entrezgene))
  # t <- gene_cn[duplicated(gene_cn$entrezgene),]
  # dim(t)
  # head(t)
  # t1 <- gene_cn[gene_cn$entrezgene=='10000',]
  gene_cn <- gene_cn %>%
    cbind(mcols(segments_gr[subjectHits(overlaps)]))
  # head(gene_cn)
  # dim(gene_cn)
  # length(unique(gene_cn$entrezgene))
  ext_rows <- gene_cn %>% # Remove the cases where a given gene that map to several regions.
    dplyr::group_by(entrezgene, cluster) %>%
    dplyr::summarise(percent_overlap=max(percent_overlap))%>% 
    dplyr::mutate(desc=paste0(entrezgene, cluster, percent_overlap))%>%
    dplyr::ungroup()%>%
    dplyr::pull(desc)
  # length(ext_rows)
  gene_cn <- gene_cn %>%
    dplyr::mutate(desc=paste0(entrezgene, cluster, percent_overlap))%>%
    dplyr::filter(desc %in% ext_rows)%>%
    dplyr::select(-desc)%>%
    dplyr::mutate(percent_overlap=round(percent_overlap, 4))
  
  
  # dim(gene_cn)
  gene_cn$ensembl_gene_id <- mapIds(org.Hs.eg.db,
                                    keys = gene_cn$entrezgene,
                                    column="ENSEMBL",
                                    keytype="ENTREZID",
                                    multiVals="first") #"CharacterList"
  # ens_mapped <- mapIds(org.Hs.eg.db,
  #             keys = gene_cn$entrezgene,
  #             column="ENSEMBL",
  #             keytype="ENTREZID",
  #             multiVals="list")
  
  ## Other way to do it
  # chrs <- c(as.character(1:22), "X")
  # t <- annotables::grch38 %>%
  #   dplyr::select(ensembl_gene_id = ensgene, entrezgene=entrez)  %>%
  #   dplyr::filter(entrezgene %in% gene_cn$entrezgene) #  & chr %in% chrs
  # #   inner_join(deg_df) %>%
  
  
  # ens_mapped <- unlist(ens_mapped)
  # ens_mapped_out <- data.frame(entrezgene=names(ens_mapped), ensembl_gene_id=ens_mapped)
  # ens_mapped_out <- ens_mapped_out %>%
  #   dplyr::filter(!is.na(ensembl_gene_id))
  # gene_cn <-  gene_cn %>% inner_join(ens_mapped_out, by='entrezgene')
  # print(dim(gene_cn))
  
  # Filter for complete entries
  gene_cn <- gene_cn %>%
    as.data.frame %>%
    drop_na()
  
  # gene_cn <- gene_cn %>%
  #   dplyr::filter(percent_overlap>min_pc_overlap)
  # print(dim(gene_cn))
  # gene_cn <- annotables::grch37 %>% 
  #   dplyr::select(ensembl_gene_id = ensgene, symbol) %>% 
  #   inner_join(gene_cn)
  # print(dim(gene_cn))
  # gene_cn$percent_overlap <- round(gene_cn$percent_overlap, 4)
  return(gene_cn)
}
check_input_data_files <- function(cells_id, input_dir){
  obs_libs <- get_library_id(cells_id)
  obs_libs <- unique(obs_libs)
  missing_libs <- c()
  for(l in obs_libs){
    # lib_fn <- paste0(input_dir,l,'/',l,'/hmmcopy/',l,'_reads.csv.gz') # csv.gz or csv
    lib_fn1 <- paste0(input_dir,l,'/hmmcopy/','reads.csv.gz') # csv.gz or csv
    # lib_fn <- paste0(download_dir,l,'/',l,'/hmmcopy/','reads.csv.gz') # csv.gz or csv
    lib_fn2 <- paste0(input_dir,l,'/hmmcopy/',l,'_reads.csv.gz') # csv.gz or csv
    lib_fn <- ifelse(file.exists(lib_fn1),lib_fn1,
                     ifelse(file.exists(lib_fn2),lib_fn2,''))
    if(lib_fn==''){
      print(l)
      missing_libs <- c(missing_libs, l)
      # print('Do not exist hmmcopy reads for library, double check input data')
    } 
  }  
  if(!is.null(missing_libs)){
    print('Do not exist hmmcopy reads for library, double check input data: ')  
    print(paste(missing_libs, collapse = '; '))
  }
  
  return(missing_libs)
}
get_mapping_genes <- function(segments){
  if(!grepl('chr',segments$chr[1])){
    segments <- segments %>%
      dplyr::mutate(chr=paste0("chr", chr))
  }
  # unassigned_clone_lbs <- paste0('clone_',c('None','Unassigned', 'Un','un'))
  segments <- segments %>%
    # dplyr::filter(clone_id =='clone_D') %>%
    # dplyr::filter(!clone_id %in% unassigned_clone_lbs) %>%
    dplyr::rename(cluster=clone_id)#%>%
  # dplyr::select(chr,start,end,cluster,copy_number,pct_pure)
  gene_cn <- get_gene_coordinates(segments)
  return(gene_cn)
}
# ex: obs_clones <- c('I,'B-D') 
get_reads_per_clone <- function(input_dir, cellclones_fn, output_file, 
                                obs_clones=NULL, tag=NULL, cores_use=NULL){
  if(is.null(cores_use)){
    cores_use <- 5
  }
  nbcores <- detectCores()
  if(cores_use > nbcores){
    cores_use <- nbcores
  }
  cores_use
  output_dir <- paste0(dirname(output_file),'/')
  if(!file.exists(output_dir)){
    dir.create(output_dir)
  }
  print("Read clone labels file")
  # input_dir <- paste0(input_dir,'/')
  print(input_dir)
  
  # cell_clones <- read.delim(cellclones_fn, stringsAsFactors = FALSE, sep = ",")
  cell_clones <- data.table::fread(cellclones_fn)
  print("DEBUG 1")
  print(dim(cell_clones))
  print("Remove unassigned cells")
  
  cell_clones <- cell_clones %>% 
    dplyr::filter(!clone_id %in% c('Un','unassigned','Unassigned','None','un'))
  print(dim(cell_clones))
  print(summary(as.factor(cell_clones$clone_id)))
  
  missing_libs <- check_input_data_files(cell_clones$cell_id, input_dir)
  if(!is.null(missing_libs)){
    df <- tibble::tibble(library_id=missing_libs)
    data.table::fwrite(df, paste0(input_dir, 'missing_libraries.csv'))
    stop("Missing data, check input data first")
  }
  
  
  # cell_clones <- cell_clones %>% 
  #   dplyr::mutate(clone_id = case_when(clone_id=='A' ~ 'R',  # the most updated csv file R --> A
  #                                      TRUE ~ clone_id))
  # ref_dif <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/'
  # reference_genes_df <- data.table::fread(paste0(ref_dif,'Symbol_ensembl.csv')) %>% as.data.frame()
  # reference_genes_df <- reference_genes_df %>%
  #   dplyr::rename(ensembl_gene_id=Ensembl)
  # obs_clones <- 'T'
  chrs <- c(as.character(1:22), "X")  # Remove Y from analysis
  
  ## Reference genes from annotables package
  r38 <- annotables::grch38 %>% 
    dplyr::select(ensembl_gene_id = ensgene, gene_symbol = symbol, chr) %>%
    dplyr::filter(chr %in% chrs) %>% 
    dplyr::select(-chr)
  # dim(r38)
  
  cls <- unique(cell_clones$clone_id)
  if(is.null(obs_clones)){
    print("Get reads data for each clone")
    obs_clones <- cls
  }else{
    obs_clones <- union(obs_clones, cls) # take into account all clones in datasets. 
    # obs_clones <- obs_clones[obs_clones %in% cls]   
  }
  if(is.null(tag)){
    tag <- paste(obs_clones, sep='_')
  }  
  print(obs_clones)
  # cols_use <- c('chr', 'start', 'end', 'width', 'gc', 'map', 'reads','cell_id')
  total_reads <- tibble::tibble()
  mapped_clones <- tibble::tibble()
  for(c in obs_clones){
    cs <- as.character(unlist(strsplit(c, '-')))
    cs <- cs[cs %in% cls]
    if(length(cs)!=sum(cs %in% cls)){
      print('Double check input data!!!')
    }
    cells_id <- cell_clones %>%
      dplyr::filter(clone_id %in% cs) %>%
      dplyr::pull(cell_id)
    if(length(cells_id)>2000){ #sampling data, in case of a big clone, just extracting 2000 random cells from this clone
      cells_id <- sample(cells_id, 2000)
    }
    print('__________________________________')
    print(paste0('Observed clone: ',c,' and nb cells: ', length(cells_id)))
    obs_libs <- get_library_id(cells_id)
    obs_libs <- unique(obs_libs)
    print(paste0('Clone: ',c,' and list of observed libraries: '))
    print(paste(obs_libs, collapse = ' '))
    obs_reads <- list()
    for(l in obs_libs){
      # lib_fn <- paste0(input_dir,l,'/',l,'/hmmcopy/',l,'_reads.csv.gz') # csv.gz or csv
      lib_fn1 <- paste0(input_dir,l,'/hmmcopy/','reads.csv.gz') # csv.gz or csv
      # lib_fn <- paste0(download_dir,l,'/',l,'/hmmcopy/','reads.csv.gz') # csv.gz or csv
      lib_fn2 <- paste0(input_dir,l,'/hmmcopy/',l,'_reads.csv.gz') # csv.gz or csv
      lib_fn <- ifelse(file.exists(lib_fn1),lib_fn1,
                       ifelse(file.exists(lib_fn2),lib_fn2,''))
      if(lib_fn==''){
        stop('Do not exist hmmcopy reads for library, double check input data')
      } 
      
      if(file.exists(lib_fn)){
        tmp <- data.table::fread(lib_fn)
        # View(head(tmp))
        # dim(tmp)
        tmp <- tmp %>%
          dplyr::filter(cell_id %in% cells_id)
        tmp <- tmp %>%
          dplyr::select(chr, start, end, width, gc, map, reads, cell_id, state) # 
        print(l)
        print(dim(tmp))
        if(dim(tmp)[1]>0){
          obs_reads[[l]] <- tmp  
        }else{
          print(paste0('***Warning:  Do not have any output for clone: ',c))
        }
        
      }else{
        print(paste0('***ERROR:  Clone: ',c, ' and library: ', l))
        # print('Do not exist hmmcopy reads for library, double check input data')
        stop('Do not exist hmmcopy reads for library, double check input data')
      }
      
    }
    reads_df <- as.data.frame(dplyr::bind_rows(obs_reads))
    print(dim(reads_df))
    # data.table::fwrite(reads_df,file=paste0(output_dir,'clone_',c,'.csv'), row.names=F, quote=F, sep=',')
    ## Version 1
    # stat_reads_df <- reads_df %>% 
    #   dplyr::group_by(chr, start, end) %>% 
    #   dplyr::summarise(reads=sum(reads),
    #                    gc=mean(gc),
    #                    map=mean(map))%>% 
    #   dplyr::ungroup()%>%
    #   dplyr::select(chr, start, end, gc, map, reads)
    
    ## Version 2
    stat_reads_df <- reads_df %>% 
      dplyr::group_by(chr, start, end) %>% 
      dplyr::summarise(reads=sum(reads, na.rm=T),
                       gc=round(mean(gc, na.rm=T), 4),
                       avg_mappability=round(mean(map, na.rm=T), 4),
                       cn_mean=round(mean(state, na.rm=T), 4),
                       cn_min=min(state, na.rm=T),
                       cn_max=max(state, na.rm=T),
                       cn_median=median(state, na.rm=T),
                       cn_mode=calc_mode(state),
                       pct_pure=round(sum(cn_mode == state)/n(), 4),
                       n_cells=n()
      )%>% 
      dplyr::ungroup()#%>%
    # dplyr::select(chr, start, end, cn_median, cn_mode, avg_mappability, cn_mean, cn_min, cn_max, reads, gc)
    
    print(dim(stat_reads_df))
    # print(summary(stat_reads_df$pct_pure))
    # print(sum(stat_reads_df$pct_pure>=0.6))
    # print(sum(stat_reads_df$avg_mappability<0.9))
    # colnames(stat_reads_df)
    # stat_reads_df <- stat_reads_df %>%
    #   dplyr::select(chr, start, end, width, gc, map, reads)
    
    # datatag <- args$datatag #'SA535'
    stat_reads_df$clone_id <- c
    # data.table::fwrite(stat_reads_df,file=paste0(output_dir,'reads_clone_',c,'.csv'))
    total_reads <- dplyr::bind_rows(total_reads, stat_reads_df)
    gene_cn <- get_mapping_genes(stat_reads_df)
    print(dim(gene_cn))
    gene_cn$clone_id <- c
    # data.table::fwrite(gene_cn,file=paste0(output_dir,'mapped_gene_clone_',c,'.csv.gz'))
    gene_cn <- gene_cn %>%
      dplyr::group_by(ensembl_gene_id, clone_id) %>%
      dplyr::summarise(avg_mappability=mean(avg_mappability),
                       cn_median=mean(cn_median),pct_pure=mean(pct_pure)) 
    
    mapped_clones <- dplyr::bind_rows(mapped_clones, gene_cn)
    gene_cn2 <- gene_cn %>%
      dplyr::select(-clone_id) %>%
      dplyr::rename(!!sym(paste0('avg_mappability_',c)):=avg_mappability,
                    !!sym(paste0('cn_median_',c)):=cn_median,
                    !!sym(paste0('pct_pure_',c)):=pct_pure)
    # dim(gene_cn2)
    # data.table::fwrite(gene_cn2,file=paste0(output_dir,'filtered_mapped_gene_clone_',c,'.csv.gz'))
    r38 <- r38 %>% left_join(gene_cn2, by='ensembl_gene_id')
    # dim(r38)
  }  
  
  data.table::fwrite(total_reads,file=paste0(output_dir,'total_reads_',tag,'.csv.gz'))
  data.table::fwrite(mapped_clones,file=paste0(output_dir,'mapped_clones_',tag,'.csv.gz'))
  r38 <- r38 %>% tidyr::drop_na()
  data.table::fwrite(r38,file=paste0(output_dir,'mapped_wholedata_',tag,'_grch38.csv.gz'))
  
  mapped_clones <- mapped_clones %>%
    dplyr::select(-avg_mappability, -pct_pure) %>%
    tidyr::pivot_wider(names_from='clone_id', values_from='cn_median')
  
  data.table::fwrite(mapped_clones,file=paste0(output_dir,'mapped_wholedata_',tag,'_v2.csv.gz'))
  
}



# # 6206 bins
# cellclones_fn <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cell_clones/SA609_cell_clones_metadata.csv'
# # cellclones_fn <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cell_clones/SA609_cell_clones.csv'
# input_dir <- '/home/htran/storage/raw_DLP/drug_resistance_DLP/SA609/'
# output_dir <- '/home/htran/storage/raw_DLP/drug_resistance_DLP/SA609/reads_clones_v3/'
# output_file <- paste0(output_dir,'out.csv')
# obs_clones <- c('B','R','H','C')
# tag <- 'SA609'
# get_reads_per_clone(input_dir, cellclones_fn, output_file, obs_clones, tag, cores_use=NULL)
# 
# 
# df <- data.table::fread(paste0(output_dir,'mapped_wholedata_SA609_grch38.csv.gz'))
# dim(df)
# colnames(df)



cellclones_fn <- '/home/htran/Projects/hakwoo_project/metastasis_material/materials/dlp_trees/SA535/cell_clones.csv.gz'
input_dir <- '/home/htran/storage/raw_DLP/metastasis_DLP/SA535/'
output_dir <- '/home/htran/storage/raw_DLP/metastasis_DLP/SA535/reads_clones_v3/'
output_file <- paste0(output_dir,'out.csv')
obs_clones <- NULL
tag <- 'SA535'
get_reads_per_clone(input_dir, cellclones_fn, output_file, obs_clones, tag, cores_use=NULL)


# paste0(output_dir,'mapped_wholedata_',tag,'_v2.csv.gz')

# /home/htran/anaconda3/envs/sisyphus38/bin/python jira_ticket_from_file.py --input_library_fn /home/htran/storage/raw_DLP/metastasis_DLP/SA535/missing_libraries.csv --output_dir /home/htran/storage/raw_DLP/metastasis_DLP/SA535/
