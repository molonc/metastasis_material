suppressPackageStartupMessages({
  library(tidyverse)
  library(SingleCellExperiment)
  # library(scater)
  library(data.table)
  library(methods)
  # library(scran)
  library(parallel)
  # library(feather)
  library(annotables)
  library(dplyr)
  library(ggplot2)
  # library(snpEnrichment)
  # library(fgsea)
  # library(scMerge)
})

# library(extrafont)
# font_import(prompt=F) # import all your fonts
# fonts()
# Based on plottinglist function from schnapps package, Marc William
# xaxis_order = "genome_position"
# CNbins is a dataframe with columns: start, end, chr, state

# tag: site Metastasis, Primary
viz_CN_profile <- function(CNbins, cn_change, save_dir, obs_clone, 
                           de_desc='', annotate_genes=T, tag='', obs_chr=NULL, add_legend=T){
  options(ggrepel.max.overlaps = Inf)
  legend_pos <- "none"
  if(add_legend){
    legend_pos <- "top"
  }
  CNbins <- CNbins %>%
    dplyr::select(c(obs_clone,'chr_desc'))
  
  # cnv <- median_cnv %>%
  #   pivot_longer(!chr_desc, names_to = "state", values_to = "A")
  
  CNbins <- CNbins %>%
    dplyr::rename(state=obs_clone)
  
  print(dim(CNbins))
  CNbins <- get_chr_infos(CNbins)
  
  obs_chrs = as.character(c(paste0(1:22), "X"))
  CNbins <- CNbins %>%
    dplyr::filter(chr %in% obs_chrs)
  
  plottitle <- paste0("SA919: Clone ",as.character(strsplit(obs_clone, "_")[[1]][1]), ' - ',tag) #"Median CN profile of ",
  pointsize = 0.2
  alphaval = 1
  maxCN = 10
  statecol = "state"
  y_axis_trans = "identity"
  xaxis_order = "genome_position"
  pl <- prepare_data(CNbins, maxCN = 20)
  statecolpal <- schnapps::scCNstate_cols()
  CNbin2 <- pl$CNbins 
  CNbin2$copy <- round(CNbin2$state)
  CNbin2$state <- round(CNbin2$state)
  CNbin2 <- CNbin2 %>%
    dplyr::mutate(state = ifelse(state >= 11, "11+", paste0(state))) %>%
    dplyr::mutate(state = factor(paste0(state), levels = c(paste0(seq(0, 10, 1)), "11+")))
  
  
  
  gCN <-  ggplot2::ggplot(CNbin2, ggplot2::aes(x = idx, y = copy)) +
    ggplot2::geom_vline(xintercept = pl$chrbreaks, col = "grey90", alpha = 0.3) +
    # ggrastr::geom_point_rast(ggplot2::aes_string(col = statecol), size = pointsize, alpha = alphaval) +
    ggplot2::geom_jitter(ggplot2::aes_string(col = statecol), size = pointsize, alpha = alphaval, width = 0.01, height = 0.1) +
    # ggplot2::geom_point(ggplot2::aes_string(col = statecol), size = pointsize, alpha = alphaval) +  #
    ggplot2::scale_color_manual(name = "CNV",
                                breaks = names(statecolpal),
                                labels = names(statecolpal),
                                values = statecolpal,
                                drop = FALSE) +
    ggplot2::theme(text=element_text(family="Arial"),
                   axis.title.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_text(color='black', hjust = 0.5, size=9, family="Arial"),
                   plot.title = ggplot2::element_text(color='black', hjust = 0.5, size=10, family="Arial"),
                   axis.ticks.y = ggplot2::element_blank(),
                   legend.position = "none") +
    ggplot2::scale_x_continuous(breaks = pl$chrticks, labels = pl$chrlabels, expand = c(0, 0), limits = c(pl$minidx, pl$maxidx)) + #,guide = ggplot2::guide_axis(check.overlap = TRUE)) +
    ggplot2::scale_y_continuous(breaks = seq(0, maxCN, 2), limits = c(0, maxCN),trans = y_axis_trans) +
    # ggplot2::xlab() +
    # ggplot2::ylab() +
    ggplot2::labs(x="Chromosome",y="Median CNV", title=plottitle) + 
    cowplot::theme_cowplot() +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow=1, #, byrow = TRUE
                                                   override.aes = list(alpha=1, size = 0.9, shape = 15))) +
    ggplot2::theme(legend.title = ggplot2::element_text(color='black', hjust = 0.5, size=8, family="Arial"), 
                   legend.position = legend_pos)
  
  
  # Just testing, add label
  if(annotate_genes){
    cn_change$chr_desc <- as.character(cn_change$chr_desc)
    cn_change <- cn_change[!duplicated(cn_change$symbol),]
    df_annot <- CNbin2 %>% inner_join(cn_change, by="chr_desc")
    print(dim(df_annot))
    if(!is.null(obs_chr)){
      df_annot <- df_annot %>%
        dplyr::filter(chr %in% as.character(obs_chr))
    }
    gCN <- gCN + ggrepel::geom_text_repel(data = df_annot, aes(label = symbol), size = 3,
                                          nudge_x = .15,
                                          box.padding = 0.5,
                                          nudge_y = .1, max.overlaps = Inf,
                                          min.segment.length = 0,  # draw segment lines, not matter how short they are
                                          color='black', segment.alpha = 0.6)  #segment.size = 0.05, 
    
  }
  
  # gCN
  # png(paste0(save_dir, obs_clone, '_de_',de_desc,'_mediancn.png'), height = 2*250, width=2*1000, res = 2*72)
  # print(gCN)
  # dev.off()
  res <- list(pl=pl, cnplt=gCN, obs_clone=obs_clone, de_desc=de_desc, tag=tag, cn_change=cn_change)
  # saveRDS(res, paste0(save_dir, obs_clone, '_de_',de_desc,'_mediancn.rds'))
  return(res)
  
}
get_related_pathways <- function(cn_change, save_dir, datatag, pathway=NULL){
  print(dim(cn_change))
  cn_change <- cn_change[!duplicated(cn_change$symbol),]
  rownames(cn_change) <- cn_change$symbol
  if(is.null(pathway)){
    gmt_fn="/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/pathway_set/h.all.v7.0.symbols.gmt"
    hm <- fgsea::gmtPathways(gmt_fn) 
    pathway <- dplyr::bind_rows(lapply(names(hm), function(x) {
      pt <- data.frame(gene_id=as.character(hm[[x]]), pathway=rep(x, length(hm[[x]])), stringsAsFactors = F) 
      return(pt)
    }))
    # print(dim(pathway))
  }
  for(g in cn_change$symbol){
    if(g %in% pathway$gene_id){
      cn_change[g,'pathway'] <- paste(gsub('HALLMARK_','',pathway[pathway$gene_id==g,'pathway']), collapse='--')
    }else{
      cn_change[g,'pathway'] <- NA
    }
    
  }
  cn_change <- cn_change %>%
    dplyr::filter(!is.na(pathway))
  # print(cn_change$pathway)
  write.csv(cn_change, paste0(save_dir,datatag,'_pathway.csv'), quote = F, row.names = F)
  
  return(cn_change)
} 
## Get enriched pathways using gprofiler
## First, install package gprofiler2
## obs_genes_symb: vector of genes symbols
## datatag: for naming, ex: 'SA919'
## save_data: save data or return a data frame

get_gprofiler_pathways_obsgenes <- function(obs_genes_symb, save_dir, datatag, 
                                            custom_id=NULL, pathway_fn=NULL, save_data=F){
  library(gprofiler2)
  if(is.null(pathway_fn)){
    # pathway_fn = 'yourdir/h.all.v7.0.symbols.gmt'  
    pathway_fn="/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/pathway_set/h.all.v7.0.symbols.gmt"
  }
  
  ref_set <- fgsea::gmtPathways(pathway_fn)
  # for(s in names(ref_set)){ ## just quick check the size of each set
  #   print(s)
  #   print(length(ref_set[[s]]))
  # }
  if(is.null(custom_id)){
    custom_id <- gprofiler2::upload_GMT_file(pathway_fn)  ## need to generate one id for each usage
  }
  stat <- NULL
  ## correction_method: one of 'fdr', 'gSCS', 'bonferroni'
  ## using the common use method correction_method='gSCS'
  gostres <- gprofiler2::gost(list(obs_genes_symb), organism = custom_id, correction_method='gSCS')
  ths <- 0.1
  if(!is.null(gostres$result)){
    stat <- gostres$result
    cols_use <- c('p_value','intersection_size','precision','recall','term_id')
    stat <- stat %>%
      dplyr::select(all_of(cols_use)) %>%
      dplyr::rename(reference_set=term_id, nb_signif_genes=intersection_size) %>%
      dplyr::filter(p_value<ths) # just to be sure
    
    
    # Get pathway genes 
    for(i in seq(nrow(stat))){
      pw_set <- stat$reference_set[i]
      ref_genes <- ref_set[[pw_set]]
      # obs_genes <- deg_df$gene_symbol
      intersect_genes <- intersect(obs_genes_symb, ref_genes)
      stat$signif_genes[i] <- paste(intersect_genes, collapse=',')
    }
    if(save_data){
      added_time <- gsub(':','',format(Sys.time(), "%Y%b%d_%X"))
      data.table::fwrite(stat, paste0(save_dir, 'gprofiler_pathways_',added_time,'.csv.gz'))  
    }
  }  
  return(stat)
  
} 
viz_cn_change <- function(median_cnv, obs_clone1, obs_clone2, save_dir){
  if(!dir.exists(save_dir)){
    dir.create(save_dir)  
  }
  
  
  de_desc <- paste0(obs_clone1[1],'_vs_',obs_clone2[1])
  cnv_mat <- median_cnv[,c(obs_clone1[1], obs_clone2[1], "chr_desc")]
  # var_genes <- apply(cnv_mat[,!colnames(cnv_mat) %in% c("chr_desc")], 1, var)
  # cnv_mat <- cnv_mat[var_genes > 0,]
  cnv_mat <- cnv_mat[abs(cnv_mat[,obs_clone1[1]]-cnv_mat[,obs_clone2[1]])>=1,]
  print(dim(cnv_mat))
  # View(cnv_mat)
  print('debug 1')
  cnv_mat <- cnv_mat[!duplicated(cnv_mat$chr_desc),]
  # head(cnv_mat)
  cn_label <- data.frame(chr_desc=unique(cnv_mat$chr_desc))
  cn_label <- get_chr_infos(cn_label)
  cn_change <- get_gene_coordinates_v2(cn_label, 0.6) 
  print(dim(cn_change))
  print('debug 2')
  cn_change <- cn_change %>% left_join(cnv_mat, by=c('chr_desc'))
  if(nrow(cn_change)>0){
    write.csv(cn_change, paste0(save_dir, de_desc,'_CN_change.csv'), quote=F, row.names=F)
  }else{
    if(nrow(cnv_mat)>0){
      write.csv(cnv_mat, paste0(save_dir, de_desc,'_cnv_mat.csv'), quote=F, row.names=F)
    }
  }
    # head(cn_change)
  print('debug 3')
  res_cl1 <- viz_CN_profile(median_cnv, cn_change, save_dir, obs_clone1[1], de_desc, annotate_genes=F, obs_clone1[2], NULL, T)
  # res_cl1$cnplt
  res_cl2 <- viz_CN_profile(median_cnv, cn_change, save_dir, obs_clone2[1], de_desc, annotate_genes=T, obs_clone2[2], NULL, F)
  
  res_cl22 <- viz_CN_profile(median_cnv, cn_change, save_dir, obs_clone2[1], de_desc, annotate_genes=F, obs_clone2[2], NULL, F)
  # res_cl2$cnplt
  
  p_total <- cowplot::plot_grid(res_cl1$cnplt, res_cl2$cnplt, ncol=1, labels =c('a','b'), rel_heights = c(1.1,1)) #, align = 'v'
  png(paste0(save_dir, de_desc,'_medianCN.png'), height = 2*600, width=2*1100, res = 2*72)
  print(p_total)
  dev.off()
  
  
  p_total2 <- cowplot::plot_grid(res_cl1$cnplt, res_cl22$cnplt, ncol=1, labels =c('a','b'), rel_heights = c(1.17,1)) #, align = 'v'
  png(paste0(save_dir, de_desc,'_medianCN_without_anno.png'), height = 2*400, width=2*1100, res = 2*72)
  print(p_total2)
  dev.off()
  
  # ggsave(
  #   filename = paste0(save_dir, de_desc,'_medianCN_without_anno_arial.pdf'),
  #   plot = p_total2,
  #   height = 4,
  #   width = 11,
  #   useDingbats=F)
  print('debug 4')
  df <- as.data.frame(res_cl1$cn_change)
  
  # pathway_fn="/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/pathway_set/h.all.v7.0.symbols.gmt"
  # custom_id <- gprofiler2::upload_GMT_file(pathway_fn)
  custom_id <- 'gp__aFc2_raJn_V2I'
  if(nrow(df)>0){
    get_gprofiler_pathways_obsgenes(unique(res_cl1$cn_change$symbol), save_dir, datatag, 
                                                custom_id, pathway_fn=NULL, save_data=T)
    get_related_pathways(df, save_dir, de_desc, NULL)
  }
    
}
prepare_data <- function(CNbins, maxCN = 20){
  binsize <- CNbins$end[1] - CNbins$start[1] + 1
  obs_chrs = as.character(c(paste0(1:22), "X"))
  bins <- schnapps::getBins(chrom.lengths = schnapps::hg19_chrlength, binsize = binsize) %>%
    dplyr::filter(chr %in% intersect(unique(CNbins$chr), obs_chrs)) %>%
    dplyr::mutate(idx = 1:dplyr::n())
  
  
  CNbins <- dplyr::full_join(bins, CNbins) %>%
    # dplyr::filter(!is.na(copy)) %>%
    dplyr::filter(!is.na(state)) %>%
    dplyr::filter(chr %in% obs_chrs) %>%
    # dplyr::mutate(copy = ifelse(copy > maxCN, maxCN, copy)) %>%
    dplyr::mutate(state = ifelse(state > maxCN, maxCN, state)) %>%
    dplyr::mutate(idxs = forcats::fct_reorder(factor(idx), idx)) %>%
    dplyr::mutate(CNs = forcats::fct_reorder(ifelse(is.na(state), NA,
                                                    paste0("CN", state)), state))
  
  #get breaks - first index of each chromosome
  chrbreaks <- bins %>%
    dplyr::filter(chr %in% unique(CNbins$chr)) %>%
    dplyr::group_by(chr) %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::pull(idx)
  
  #get ticks - median bin of each chromosome
  chrticks <- bins %>%
    dplyr::filter(chr %in% unique(CNbins$chr)) %>%
    dplyr::group_by(chr) %>%
    dplyr::summarise(idx = round(median(idx))) %>%
    dplyr::pull(idx)
  
  chrlabels <- gtools::mixedsort(unique(CNbins$chr))
  minidx <- min(bins$idx)
  maxidx <- max(bins$idx)
  pl <- list(CNbins = CNbins, chrbreaks = chrbreaks, chrticks = chrticks, chrlabels = chrlabels, minidx = minidx, maxidx = maxidx)
  return(pl)
}


get_chr_infos <- function(median_cnv) {
  chr_desc <- as.character(median_cnv$chr_desc)
  obs_chrs = as.character(c(paste0(1:22), "X"))
  
  chrs <- sapply(strsplit(chr_desc, "_"), function(x) {
    return(x[1])
  })
  starts <- sapply(strsplit(chr_desc, "_"), function(x) {
    return(x[2])
  })
  ends <- sapply(strsplit(chr_desc, "_"), function(x) {
    return(x[3])
  })
  median_cnv$chr <- as.character(chrs)
  median_cnv$start <- as.numeric(starts)
  median_cnv$end <- as.numeric(ends)
  median_cnv <- median_cnv %>% 
    dplyr::filter(chr %in% obs_chrs)
  # return(list(chr=as.character(chrs),start=as.character(starts),end=as.character(ends)))
  return(median_cnv)
}
calc_mode <- function(x) {
  keys <- unique(x)
  keys[which.max(tabulate(match(x, keys)))]
}
get_gene_coordinates_v2 <- function(segments, min_pc_overlap=0.1){
  segments$chr <- paste0('chr',segments$chr)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(org.Hs.eg.db)
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  g <- genes(txdb, single.strand.genes.only=FALSE)
  
  segments_gr <- makeGRangesFromDataFrame(segments, keep.extra.columns = TRUE)
  overlaps <- findOverlaps(g, segments_gr, ignore.strand = TRUE)
  overlapping_ranges <- pintersect(g[queryHits(overlaps)], segments_gr[subjectHits(overlaps)], ignore.strand = TRUE, drop.nohit.ranges = TRUE)
  
  percentage_overlap <- width(overlapping_ranges) / width(segments_gr[subjectHits(overlaps)])
  percentage_overlap_sum <- sapply(percentage_overlap, sum)
  
  gene_cn <- tibble(entrezgene = names(g)[queryHits(overlaps)], percent_overlap=percentage_overlap_sum)
  gene_cn$ensembl_gene_id <- mapIds(org.Hs.eg.db, 
                                    keys = gene_cn$entrezgene, 
                                    column="ENSEMBL", 
                                    keytype="ENTREZID", 
                                    multiVals="first")
  
  gene_cn <- gene_cn %>%
    cbind(mcols(segments_gr[subjectHits(overlaps)]))
  
  
  # Filter for complete entries
  gene_cn <- gene_cn %>%
    as.data.frame %>%
    drop_na()
  
  gene_cn <- gene_cn %>%
    dplyr::filter(percent_overlap>min_pc_overlap) %>%
    dplyr::select(-entrezgene)
  
  print(dim(gene_cn))
  gene_cn <- annotables::grch37 %>%
    dplyr::select(ensembl_gene_id = ensgene, symbol) %>%
    inner_join(gene_cn)
  print(dim(gene_cn))
  
  return(gene_cn)
}
get_gene_coordinates <- function(segments){ #, min_pc_overlap=0.1
  # segments$chr <- paste0('chr',chrs_df$chr)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(org.Hs.eg.db)
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  g <- genes(txdb, single.strand.genes.only=FALSE)
  
  segments_gr <- makeGRangesFromDataFrame(segments, keep.extra.columns = TRUE)
  overlaps <- findOverlaps(g, segments_gr, ignore.strand = TRUE)
  overlapping_ranges <- pintersect(g[queryHits(overlaps)], segments_gr[subjectHits(overlaps)], ignore.strand = TRUE, drop.nohit.ranges = TRUE)
  
  percentage_overlap <- width(overlapping_ranges) / width(segments_gr[subjectHits(overlaps)])
  percentage_overlap_sum <- sapply(percentage_overlap, sum)
  
  gene_cn <- tibble(entrezgene = names(g)[queryHits(overlaps)], percent_overlap=percentage_overlap_sum)
  gene_cn$ensembl_gene_id <- mapIds(org.Hs.eg.db, 
                                    keys = gene_cn$entrezgene, 
                                    column="ENSEMBL", 
                                    keytype="ENTREZID", 
                                    multiVals="first")
  
  gene_cn <- gene_cn %>%
    cbind(mcols(segments_gr[subjectHits(overlaps)]))
  
  
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
  print(dim(gene_cn))
  
  return(gene_cn)
}
get_library_id <- function(cell_ids) {
  
  labels <- lapply(strsplit(cell_ids, "-"), function(x) {
    return(x[2])
  })
  return(as.character(labels))
}

split_raw_segments_by_sample <- function(save_dir, raw_segments=NULL){
  # save_dir_2 <- paste0(save_dir,'gene_cn_samples/')
  # if (!file.exists(save_dir_2)){
  #   dir.create(save_dir_2)
  # }
  
  if(is.null(raw_segments)){
    raw_segments <- readRDS(paste0(save_dir,'raw_segments.rds'))
  }
  dim(raw_segments)
  sample_ids <- get_sample_id(raw_segments$cell_names, cores_use=5)
  samples <- unique(sample_ids)
  raw_segments$sample_id <- sample_ids
  for(s in samples){
    print(s)
    raw_segments_tmp <- raw_segments[raw_segments$sample_id==s,]
    # print(dim(raw_segments_tmp))
    # raw_segments_tmp <- raw_segments_tmp[ , colnames(raw_segments_tmp) != "sample_id"]
    # raw_segments_tmp <- raw_segments_tmp[ , !(colnames(raw_segments_tmp) %in% c("sample_id"))]
    
    raw_segments_tmp <- subset(raw_segments_tmp, select = c(chr,start,end,copy_number,median,multiplier,cell_names))
    print(dim(raw_segments_tmp))
    save_dir_sample <- paste0(save_dir,s,'/')
    if (!file.exists(save_dir_sample)){
      dir.create(save_dir_sample)
    }
    saveRDS(raw_segments_tmp, file = paste0(save_dir_sample, s, '.rds'))
  }
}

get_sample_id <- function(cell_ids) {
  labels <- lapply(strsplit(cell_ids, "-"), function(x) {
    return(x[1])
  })
  return(as.character(labels))
}

split_data_by_sample <- function(save_dir, gene_cn_summarized2=NULL){
  # save_dir_2 <- paste0(save_dir,'gene_cn_samples/')
  # if (!file.exists(save_dir_2)){
  #   dir.create(save_dir_2)
  # }
  
  if(is.null(gene_cn_summarized2)){
    gene_cn_summarized2 <- feather::read_feather(paste0(save_dir,'gene_cn_summarized2.feather'))
  }
  sample_ids <- get_sample_id(gene_cn_summarized2$cell_names, cores_use=8)
  samples <- unique(sample_ids)
  gene_cn_summarized2$sample_id <- sample_ids
  for(s in samples){
    print(s)
    gene_cn_tmp <- gene_cn_summarized2[gene_cn_summarized2$sample_id==s,]
    print(dim(gene_cn_tmp))
    gene_cn_tmp <- gene_cn_tmp[ , !(names(gene_cn_tmp) == "sample_id")]
    feather::write_feather(gene_cn_tmp, paste0(save_dir_2,s,'.feather'))
  }
}



get_library_metadata <- function(input_dir){
  grouping_df <- read.csv(paste0(input_dir,'library_groupings.csv'),header=T, check.names=F, stringsAsFactors = FALSE)
  # colnames(grouping_df)[which(names(grouping_df) == "passage")] <- "time"
  # colnames(grouping_df)[which(names(grouping_df) == "treatment_st")] <- "treatmentSt"
  # 
  # grouping_df <- grouping_df[,c("library_labels","treatmentSt","time"), drop=F]
  # grouping_df$library_labels <- as.factor(grouping_df$library_labels)
  return(grouping_df)
}

# entrezgene, ensembl_gene_id
# https://stackoverflow.com/questions/62140483/how-to-interpret-dplyr-message-summarise-regrouping-output-by-x-override
get_gene_cn_summary <- function(unique_genes, base_name, save_dir, gene_cn, output_file, cores_use=8){
  
  # if(is.null(gene_cn)){
  #   gene_cn <- readRDS(paste0(save_dir,'gene_cn.rds'))
  # }
  print("Get gene_cn_summarized")
  # gene_cn_summarized <- mclapply(unique_genes, function(g) {
  #   dplyr::filter(gene_cn, ensembl_gene_id == g)  %>%
  #     dplyr::group_by(ensembl_gene_id, cell_names, time, cluster, copy_number) %>%
  #     dplyr::summarise(percent_overlap=sum(percent_overlap),
  #                      n_parts=n()) %>% 
  #     ungroup()
  # }, mc.cores = cores_use) %>% bind_rows()
  
  # TO DO: add time column to list
  gene_cn_summarized <- mclapply2(unique_genes, function(g) {
    dplyr::filter(gene_cn, ensembl_gene_id == g)  %>%
      dplyr::group_by(ensembl_gene_id, cell_names, time, cluster, copy_number) %>%
      dplyr::summarise(percent_overlap=sum(percent_overlap),
                       n_parts=n(), .groups = 'drop') %>% 
      ungroup()
  }, mc.cores = cores_use) %>% bind_rows()
  # feather::write_feather(gene_cn_summarized, paste0(save_dir,base_name,'_gene_cn_summarized.feather'))
  
  
  print("Get gene_cn_summarized2")
  cat("Summarizing to gene specific CN...\n")
  # gene_cn_summarized2 <- mclapply(unique_genes, function(g) {
  #   dplyr::filter(gene_cn_summarized, ensembl_gene_id == g)  %>%
  #     dplyr::group_by(ensembl_gene_id, cell_names, time, cluster) %>%
  #     dplyr::summarise(copy_number_mean=sum(copy_number*percent_overlap)/sum(percent_overlap),
  #                      copy_number_min=min(copy_number),
  #                      copy_number_max=max(copy_number),
  #                      copy_number_mode=copy_number[which.max(percent_overlap)][1],
  #                      n_parts=sum(n_parts)) %>% 
  #     ungroup()
  # }, mc.cores = cores_use) %>% bind_rows()
  
  
  gene_cn_summarized2 <- mclapply2(unique_genes, function(g) {
    dplyr::filter(gene_cn_summarized, ensembl_gene_id == g)  %>%
      dplyr::group_by(ensembl_gene_id, cell_names, time, cluster) %>%
      dplyr::summarise(copy_number_mean=sum(copy_number*percent_overlap)/sum(percent_overlap),
                       copy_number_min=min(copy_number),
                       copy_number_max=max(copy_number),
                       copy_number_mode=copy_number[which.max(percent_overlap)][1],
                       n_parts=sum(n_parts), .groups = 'drop') %>% 
      ungroup()
  }, mc.cores = cores_use) %>% bind_rows()
  
  print("Save output to feather files")
  # feather::write_feather(gene_cn_summarized2, paste0(save_dir,'gene_cn_summarized2.feather'))
  feather::write_feather(gene_cn_summarized2, output_file)
  
  
  cat("Completed.\n")
}



preprocess_data <- function(sce, raw_cnvs, base_name, save_dir, saveOutput=TRUE,
                            min_counts_per_gene=50, min_counts_per_cell=100, max_copy_number=7){
  
  print(paste0("Debug: raw_cnvs: ",dim(raw_cnvs)[1],'  ',dim(raw_cnvs)[2]))
  ## Important -- select autosomes only
  
  autosomal_genes <- annotables::grch37 %>% 
    dplyr::select(ensembl_gene_id = ensgene, chr) %>% 
    dplyr::filter(chr %in% as.character(1:22)) %>% 
    .$ensembl_gene_id
  
  
  # print(paste0("Debug: autosomal_genes",length(autosomal_genes)))
  
  ## Filter CNV data
  cnv <- dplyr::filter(raw_cnvs, use_gene) %>%
    dplyr::rename(clone = cluster,
                  copy_number=median_cnmode) %>% 
    dplyr::select(ensembl_gene_id, clone, copy_number) %>% 
    # dplyr::filter(clone %in% present_clones) %>%
    spread(clone, copy_number)
  
  cnv_mat <- cnv %>%
    as.data.frame %>%
    column_to_rownames("ensembl_gene_id") %>%
    as.matrix
  
  print("Debug: cnv_mat ")
  print(colnames(cnv_mat))
  
  common_genes <- intersect(intersect(rownames(cnv_mat), rownames(sce)), autosomal_genes)
  
  print(paste0("Debug: nb common_genes: ",length(common_genes)))
  sce <- sce[common_genes,]
  cnv_mat <- cnv_mat[common_genes,]
  
  print(paste0("Debug: cnv_mat: ",dim(cnv_mat)[1]," ",dim(cnv_mat)[2]))
  print(paste0("Debug: sce: ",dim(sce)[1]," ",dim(sce)[2]))
  if(saveOutput){
    saveRDS(sce, file = paste0(save_dir,base_name,'_sce_filtered.rds'))
    saveRDS(cnv_mat, file = paste0(save_dir,base_name,'_cnv_mat.rds'))
  }
  
  # One more layer of preprocessing cnv data
  # Possible error: use_gene vectors are different for each clone
  # rv_genes <- c()
  # for(c in rownames(cnv_mat)){
  #   if(sum(is.na(cnv_mat[c,]))>0){
  #     rv_genes <- c(rv_genes,c)
  #   }
  # }
  # print(paste0("Debug: rv_genes in cnv_mat: ",length(rv_genes)))
  # cnv_mat <- cnv_mat[!rownames(cnv_mat) %in% rv_genes,]
  # print(dim(cnv_mat))
  # sce <- sce[rownames(cnv_mat),]
  # print(dim(sce))
  
  processed_data <- clonealign:::preprocess_for_clonealign(sce, cnv_mat, 
                                              min_counts_per_gene, 
                                              min_counts_per_cell, TRUE, 10, max_copy_number, TRUE)
  
  # processed_data <- filter_data(processed_data$gene_expression_data, processed_data$copy_number_data)
  
  if(saveOutput){
    saveRDS(processed_data, file = paste0(save_dir,base_name,'_processed_data.rds'))
  }
  # print(paste0("Debug: processed_data ",length(processed_data)))
  return(processed_data)
}

# save_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_rna/clonealign/SA1035X4XB02879/'
# base_name <- 'SA1035X4XB02879'
# processed_data <- readRDS(paste0(save_dir,base_name,'_processed_data.rds'))
# dim(processed_data$copy_number_data)
# sce <- readRDS(paste0(save_dir,base_name,'_sce_filtered.rds'))
# cnv_mat <- readRDS(paste0(save_dir,base_name,'_cnv_mat.rds'))
# dim(sce)
# dim(cnv_mat)

# One more layer of filtering in case clonealign dont work well
filter_data <- function(gene_expression_data, cnv_mat){
  var_log_counts <- colVars(as.matrix(log2(gene_expression_data+1)))
  gex_var_quantile <- 0.5
  var_quantile <- quantile(var_log_counts, probs = as.numeric(gex_var_quantile),na.rm=T)
  # # rowData(sce)$total_counts > 0 & 
  print(var_quantile)
  genes <- colnames(gene_expression_data)
  genes <- genes[var_log_counts > var_quantile]
  print(paste0("Nb common genes after twice filtering: ",length(genes)))
  gene_expression_data <- gene_expression_data[,colnames(gene_expression_data) %in% genes]
  cnv_mat <- cnv_mat[rownames(cnv_mat) %in% colnames(gene_expression_data),]
  list(gene_expression_data = gene_expression_data, copy_number_data = cnv_mat, 
       retained_cells = rownames(gene_expression_data), 
       retained_genes = colnames(gene_expression_data))
}
# max_copy_number=11

filtering_10x <- function(sce, gex_var_quantile=0.5){
  rowData(sce)$var_log_counts <- rowVars(as.matrix(logcounts(sce)))
  var_quantile <- quantile(rowData(sce)$var_log_counts, probs = as.numeric(gex_var_quantile))
  # rowData(sce)$total_counts > 0 & 
  sce <- sce[rowData(sce)$var_log_counts > var_quantile,]
  usual_suspects <- grepl("^HSP|^FOS|^JUN|^MALAT1|^UBC", rowData(sce)$Symbol)
  sce <- sce[!usual_suspects,]
  return(sce)
}

run_CloneAlign <- function(sce, processed_data, output_file, base_name, save_dir='', 
                           saveOutput=T, max_copy_number=6, max_iteration=50, clone_probs_thrs=0.5){
  print(paste0("Dim sce expression data: ", dim(processed_data$gene_expression_data)))
  cat(paste0("\n Dim expression data: ", dim(processed_data$gene_expression_data)[1],' ',
                                         dim(processed_data$gene_expression_data)[2]), 
                                         file = paste0(save_dir,"clonealign_log.txt"), append = TRUE)
  
  print(dim(processed_data$copy_number_data))
  print("Run CloneAlign")
  ca <- run_clonealign(processed_data$gene_expression_data, 
                       processed_data$copy_number_data,
                       initial_shrinks = c(2,5,10), #args$initial_shrink,
                       n_repeats = 3,
                       mc_samples = 1,
                       learning_rate = 0.07,  #0.07
                       max_iter = max_iteration,  #5e2
                       data_init_mu = TRUE,
                       saturation_threshold = max_copy_number,  # in fitness paper: saturation threshold = 6
                       clone_call_probability = 0.9)
  
  
  sce_filtered <- sce[, rownames(processed_data$gene_expression_data)]
  print(paste0("Debug: sce_filtered",dim(sce_filtered)))
  cat(paste0("\n Debug: sce_filtered",dim(sce_filtered)), file = paste0(save_dir,"clonealign_log.txt"), append = TRUE)
  
  clone_fit <- data.frame(
    cell_id=rownames(processed_data$gene_expression_data),
    id=sce_filtered$library_label,
    clone=ca$clone,
    stringsAsFactors = FALSE
  )
  
  ca$clone_fit <- clone_fit
  ca$cnv_mat <- processed_data$copy_number_data
  ca$sce_filtered <- sce_filtered
  # ca$cnv_mat <- cnv_mat_remapped
  
  # Save results
  # saveRDS(ca, args$outfname)
  print("Save output of CloneAlign")
  # saveRDS(ca, paste0(save_dir,'cloneAlign_result_counts_v3.rds'))
  saveRDS(ca, output_file)
  
  print("Generate stat plots")
  #paste0(save_dir,base_name,"_clonealign_convergence.png")
  png(paste0(save_dir,base_name,"_convergence.png"), height = 2*300, width=2*400,res = 2*72)
  pc <- qplot(seq_along(ca$convergence_info$elbo), ca$convergence_info$elbo, 
              geom = c("point", "line")) + labs(x = "Iteration", y = "ELBO")
  print(pc)
  dev.off()
  
  png(paste0(save_dir,base_name,"_clonealign_heatmap_clone_probs.png"), height = 2*600, width=2*700,res = 2*72)
  p <- pheatmap::pheatmap(ca$ml_params$clone_probs)
  print(p)
  dev.off()
  cat("Completed.\n")
  cat(paste0("\n Completed: ",base_name), file = paste0(save_dir,"clonealign_log.txt"), append = TRUE)
  
  get_clones_assignment(ca, base_name, save_dir, clone_probs_thrs)
  
}  


get_clones_assignment <- function(ca, base_name, save_dir, clone_probs=0.5){
  
  # save_dir <- paste0(dirname(output_file),'/')
  # if (!file.exists(save_dir)){
  #   dir.create(save_dir)
  # }
  # print(save_dir)
  # base_name <- basename(output_file)
  # base_name <- gsub('.rds','',base_name)
  if(is.null(ca)){
    print('Load clone align output from folder...')
    ca <- readRDS(paste0(save_dir, base_name, '_ca.rds'))
    # ca <- readRDS(paste0(save_dir, base_name,'/', base_name, '_ca.rds'))
  }
  cat(paste0("\n Get clone assignment: ",base_name), file = paste0(save_dir,"clonealign_log.txt"), append = TRUE)
  print(paste0('With prob = 0.4: ',sum(ca$ml_params$clone_probs>0.4)))  
  cat(paste0('\n ',base_name,'  with prob >= 0.4: ',sum(ca$ml_params$clone_probs>=0.4)), file = paste0(save_dir,"clonealign_log.txt"), append = TRUE)
  
  print(paste0('With prob = 0.45: ',sum(ca$ml_params$clone_probs>0.45)))
  cat(paste0('\n ',base_name,'  with prob >= 0.45: ',sum(ca$ml_params$clone_probs>=0.45)), file = paste0(save_dir,"clonealign_log.txt"), append = TRUE)
  
  print(paste0('With prob = 0.5: ',sum(ca$ml_params$clone_probs>0.5)))
  cat(paste0('\n ',base_name,'  with prob >= 0.5: ',sum(ca$ml_params$clone_probs>=0.5)), file = paste0(save_dir,"clonealign_log.txt"), append = TRUE)
  
  print(length(ca$ml_params$clone_probs))
  
  
  cat(paste0('\n ',base_name,'  Nb 10x cells: ',dim(ca$sce_filtered)), file = paste0(save_dir,"clonealign_log.txt"), append = TRUE)
  cat(paste0('\n ',base_name,' Clone Probability Threshold: ',clone_probs), file = paste0(save_dir,"clonealign_log.txt"), append = TRUE)
  
  # ca <- clonealign:::recompute_clone_assignment(ca, clone_probs)
  clones <- recompute_clone_assignment_v2(ca, clone_probs)
  ca$clone <- clones
  cat(paste0('\n ',base_name,' CloneAlign output: ',summary(as.factor(clones))), file = paste0(save_dir,"clonealign_log.txt"), append = TRUE)
 
  # ca1 <- clonealign:::recompute_clone_assignment(ca, 0.45)
  
  df <- ca$clone_fit
  df$clone <- ca$clone
  # rownames(df) <- df$cell_id
  # df_clones <- df[df$clone!='unassigned',]
  # dim(df_clones)
  # unique(df_clones$clone)
  # rownames(df_clones) <- df_clones$cell_id
  
  # write.csv(df_clones, output_file, row.names=F, quote=F) #paste0(save_dir,'cells_cloneAlign.csv')
  write.csv(df, paste0(save_dir, base_name,'_clones.csv'), row.names=F, quote=F)
  print("Save output of clone assignment into folder")
}

# Take into account the case where cells can belong to 2 clones
recompute_clone_assignment_v2 <- function(ca, clone_assignment_probability = 0.95) 
{
  clone_names <- colnames(ca$ml_params$clone_probs)
  clones <- apply(ca$ml_params$clone_probs, 1, function(r) {
    if (max(r) < clone_assignment_probability) {
      return("unassigned")
    }
    # pos_max <- which.max(r)
    idx <- which(r==max(r), arr.ind=TRUE)
    # if(length(idx)>1){
    #   print(idx)
    # }
    cl <- unique(clone_names[idx])
    cl1 <- paste(cl, collapse="_")
    return(cl1)
    # return(clone_names[which.max(r)])
  })
  return(clones)
  # ca$clone <- clones
  # ca
}

# base_name is sample_id
plot_prevalences <- function(df_clones, cell_clones, save_dir, base_name){
  # cell_clones <- read.delim(cellclones, stringsAsFactors = FALSE, check.names=F, sep = ",")
  cell_clones$sample <- get_sample_id(cell_clones$cell_id, 3)
  cell_clones <- cell_clones[cell_clones$sample==base_name,]
  
  counts <- table(cell_clones$clone_id)
  dlp_prop <- data.frame(clone=names(counts), prop=as.numeric(counts))
  pd <- ggplot(data=dlp_prop, aes(x=clone, y=prop)) +
    geom_bar(stat="identity", fill="steelblue", width=0.4)+
    theme_minimal()
  pd <- pd + labs(title=paste0("DLP clones ",base_name), 
           x="Clones", y = "# cells")
  
  df_clones <- df_clones[df_clones$clone!='unassigned',]
  counts1 <- table(df_clones$clone)
  drna_prop <- data.frame(clone=names(counts1), prop=as.numeric(counts1))
  pr <- ggplot(data=drna_prop, aes(x=clone, y=prop)) +
    geom_bar(stat="identity", fill="steelblue", width=0.4)+
    theme_minimal()
  pr <- pr + labs(title=paste0("CloneAlign ",base_name), 
                  x="Clones", y = "# cells")
  
  p <- cowplot::plot_grid(pd, pr, align = "hv", ncol = 2)
  png(paste0(save_dir,base_name,"_clones_proportion.png"), height = 2*250, width=2*600,res = 2*72)
  print(p)
  dev.off()
}

plot_clones_rd <- function(ca, input_file, save_dir){
  df <- ca$clone_fit
  df$clone <- ca$clone
  rownames(df) <- df$cell_id
  processed_data <- readRDS(paste0(save_dir,'processed_data.rds'))
  sce_filtered <- ca$sce_filtered
  sce_filtered$cell_type <- 'unassigned'
  dim(sce_filtered)
  cells_use <- intersect(colnames(sce_filtered),df$cell_id)
  length(cells_use)
  
  sce_filtered <- sce_filtered[,cells_use]
  sce_filtered$cell_type <- df[colnames(sce_filtered),'clone']
  sce_filtered <- runPCA(sce_filtered, ntop = 1000, ncomponents = 50, exprs_values = "counts")
  sce_filtered <- runTSNE(sce_filtered, dimred = "PCA", n_dimred = 50, ncomponents = 2)
  dim(sce_filtered)
  sce_filtered1 <- sce_filtered[,sce_filtered$cell_type!='unassigned']
  
  sce_filtered2 <- sce_filtered[,sce_filtered$cell_type=='unassigned']
  sce_filtered2$treatment_status <- as.factor(sce_filtered2$treatment_status)
  summary(sce_filtered2$treatment_status)
  
  dim(sce_filtered)
  dim(sce_filtered1)
  sce_filtered1 <- runPCA(sce_filtered1, ntop = 1000, ncomponents = 50, exprs_values = "counts")
  sce_filtered1 <- runTSNE(sce_filtered1, dimred = "PCA", n_dimred = 50, ncomponents = 2)
  
  saveRDS(sce_filtered1, paste0(save_dir,'sce_filtered_assigned_v2.rds'))
  saveRDS(sce_filtered, paste0(save_dir,'sce_filtered_celltype_v2.rds'))
  p <- plotReducedDim(sce_filtered, dimred="TSNE",colour_by = 'cell_type')
  # p
  
  p1 <- plotReducedDim(sce_filtered, dimred="TSNE",colour_by = 'treatment_status')
  # p1
  
  p2 <- plotReducedDim(sce_filtered1, dimred="TSNE",colour_by = 'cell_type')
  meta_data_ls <- readRDS(paste0(input_dir,'tree_viz_dream/SA1035_Tyler_meta_data_ls.rds'))
  meta_data_ls$clone_meta
  colorcode <- c('#FFD92F','#8DA0CB','#E78AC3')
  names(colorcode) <- c('C','D','E')
  p2 <- p2 + scale_fill_manual(values = colorcode, guide = FALSE)
  # p2
  p3 <- plotReducedDim(sce_filtered1, dimred="PCA",colour_by = 'treatment_status')
  # p3
  library(cowplot)
  ptotal <- plot_grid(plots = list(p2, p3), nrow=1)
  png(paste0(save_dir,"clonealign_treatmentSt_PCA.png"), height = 2*400, width=2*1120,res = 2*72)
  print(ptotal)
  dev.off()
  
  ptsne <- plot_grid(plots = list(p, p1), nrow=1)
  png(paste0(save_dir,"clonealign_treatmentSt_TSNE_all.png"), height = 2*400, width=2*550,res = 2*72)
  print(p2)
  dev.off()
  
  
}


scMerge_remove_batch_effects <- function(sce, save_dir, base_name){
  library(BiocParallel)
  # exprs_mat = SummarizedExperiment::assay(sce, 'counts')
  # print(dim(exprs_mat))
  # set.seed(1)
  
  # Manual computation of scSEG list
  # param = BiocParallel::MulticoreParam(workers = 4, progressbar = FALSE)
  # segIndx_df = scMerge::scSEGIndex(exprs_mat = exprs_mat, BPPARAM = param)
  # ## Closing the parallelisation
  # BiocParallel::register(BPPARAM = BiocParallel::SerialParam())
  # segIndx_df$gene_ens <- rownames(segIndx_df)
  # seg = segIndx_df$gene_ens
  data("segList", package = "scMerge") 
  seg = segList$human$human_scSEG
  length(seg)
  kmeansK = c(3,3)
  sce$batch <- sce$batch_info
  sce1 <- SingleCellExperiment(assays = list(counts = as.matrix(counts(sce)),
                                             logcounts=as.matrix(logcounts(sce))), 
                                      colData = colData(sce))
  if(length(unique(sce$batch))==1){
    # saveRDS(sce, paste0(save_dir, base_name,'_corrected.rds'))
    return(sce)
  } else{
    marker = NULL
    replicate_prop <- 0.8
    scMerge_res <- scMerge::scMerge(
      sce_combine = sce1, 
      ctl = seg,
      kmeansK = kmeansK,
      marker = marker,
      assay_name = "scMerge_res",
      replicate_prop = replicate_prop,
      cell_type = NULL, # unsupervised
      verbose=T,
      BPPARAM = MulticoreParam(workers = 4))
    saveRDS(scMerge_res, paste0(save_dir, base_name,'_corrected.rds'))
    
  }
    return(scMerge_res)
}

plot_prevalences_v2 <- function(df_clones, cellclones, save_dir, base_name){
  cell_clones <- read.delim(cellclones, stringsAsFactors = FALSE, check.names=F, sep = ",")
  cell_clones$sample <- get_sample_id(cell_clones$cell_id, 1)
  cell_clones <- cell_clones[cell_clones$sample==base_name,]
  df_clones <- df_clones[df_clones$clone!='unassigned',]
  counts <- table(cell_clones$clone_id)
  dlp_prop <- data.frame(clone=names(counts), prop=as.numeric(counts))
  pd <- ggplot(data=dlp_prop, aes(x=clone, y=prop)) +
    geom_bar(stat="identity", fill="steelblue", width=0.4)+
    theme_minimal()
  pd <- pd + labs(title=paste0("DLP clones ",base_name), 
                  x="Clones", y = "# cells")
  
  counts1 <- table(df_clones$clone)
  drna_prop <- data.frame(clone=names(counts1), prop=as.numeric(counts1))
  pr <- ggplot(data=drna_prop, aes(x=clone, y=prop)) +
    geom_bar(stat="identity", fill="steelblue", width=0.4)+
    theme_minimal()
  pr <- pr + labs(title=paste0("CloneAlign ",base_name), 
                  x="Clones", y = "# cells")
  
  p <- cowplot::plot_grid(pd, pr, align = "hv", ncol = 2)
  png(paste0(save_dir,base_name,"_clones_proportion.png"), height = 2*250, width=2*480,res = 2*72)
  print(p)
  dev.off()
}


# raw_cnvs: cluster_cn_summarized
# var_thres=0.05: remove genes without any variance in all clones
get_median_genotype <- function(raw_cnvs, save_dir, base_name, var_thres=0.05){
  print(dim(raw_cnvs))
  
  chrs <- c(as.character(1:22), "X")  # Remove Y from analysis
  autosomal_genes <- annotables::grch37 %>% 
    dplyr::select(ensembl_gene_id = ensgene, chr) %>% 
    dplyr::filter(chr %in% chrs) %>% 
    .$ensembl_gene_id
  # print(paste0("Debug: autosomal_genes",length(autosomal_genes)))
  
  ## Filter CNV data
  cnv <- raw_cnvs %>%
    dplyr::rename(clone = cluster,
                  copy_number=median_cn) %>% 
    dplyr::select(ensembl_gene_id, clone, copy_number) %>% 
    # dplyr::filter(clone %in% present_clones) %>%
    spread(clone, copy_number)
  
  print(dim(cnv))
  
  cnv_mat <- cnv %>%
    as.data.frame %>%
    column_to_rownames("ensembl_gene_id") %>%
    as.matrix
  
  print("Debug: cnv_mat ")
  print(colnames(cnv_mat))
  common_genes <- intersect(rownames(cnv_mat), autosomal_genes)
  cnv_mat <- cnv_mat[common_genes,]
  print(dim(cnv_mat))
  
  # Get variant genes
  # var_genes <- apply(cnv_mat, 1, var)
  # cnv_mat <- cnv_mat[var_genes >= var_thres,]
  cnv_mat <- as.data.frame(cnv_mat)
  saveRDS(cnv_mat, paste0(save_dir,base_name,'_cnv_mat.rds'))
  
  
  cnv_mat$ensembl_gene_id <- rownames(cnv_mat)
  write.csv(cnv_mat, paste0(save_dir,base_name,'_cnv_mat.csv'), row.names = F, quote=F)
  
  cnv_mat <- cnv_mat %>%
    pivot_longer(!ensembl_gene_id, names_to = "clone", values_to = "cnv")
  
  print(dim(cnv_mat))
  write.csv(cnv_mat, paste0(save_dir,base_name,'_cnv_mat_pivot.csv'), row.names = F, quote=F)
  # saveRDS(cnv_mat, paste0(save_dir,'cnv_mat_pivot.rds'))
  
}

get_median_genotype_v2 <- function(datatag, results_dir, copynumber_fn=NULL, cellclone_fn=NULL,
                                   calcul_distance=F, distance_type='Manhattan'){
  
  if(is.null(copynumber_fn)){
    copynumber_fn <- paste0(results_dir,'total_merged_filtered_states.csv')
  }
  if(is.null(cellclone_fn)){
    cellclone_fn <- paste0(results_dir,'cell_clones.csv')  
  }
  save_dir <- paste0(results_dir,'CN_profile/')
  if(!dir.exists(save_dir)){
    dir.create(save_dir)
  }
  
  copynumber <- read.csv(copynumber_fn, header=T, 
                         row.names = 1, check.names = F,stringsAsFactors = FALSE)
  
  # cell_clones contain 2 columns of cell_id, and clone_id
  # ex:           cell_id                      clone_id
  # 1   SA535X4XB02498-A98163A-R09-C11          C
  cell_clones <- read.csv(cellclone_fn, check.names = F,stringsAsFactors = FALSE)
  
  print(dim(copynumber))
  cells_use <- colnames(copynumber)[colnames(copynumber) %in% cell_clones$cell_id]
  copynumber <- copynumber[,cells_use]
  copynumber$chr_desc <- rownames(copynumber)
  cnv <- copynumber %>%
    pivot_longer(!chr_desc, names_to = "cell_id", values_to = "copy_number")
  
  print(dim(cnv))
  cnv <- cnv %>% inner_join(cell_clones, by=c("cell_id"))
  
  # get median genotype 
  print("Get median genotype")
  stat_cnv <- cnv %>%
    dplyr::group_by(clone_id, chr_desc) %>%
    dplyr::summarise(median_cn=median(copy_number)) #,mode_cn=calc_mode(copy_number)
  
  
  median_cnv <- stat_cnv %>%
    dplyr::select(chr_desc, median_cn, clone_id) %>%
    group_by(clone_id) 
  
  
  median_cnv <- median_cnv %>%
    pivot_wider(names_from = clone_id, values_from = median_cn)
  
  median_cnv <- as.data.frame(median_cnv)
  write.csv(median_cnv, paste0(save_dir,'median_cnv.csv'), quote = F, row.names = F)
  median_cnv <- median_cnv %>%
    tibble::column_to_rownames('chr_desc')
  
  # median_cnv <- median_cnv %>% select(-c(chr_desc))  
  # chr_infos <- get_chr_infos(median_cnv$chr_desc)
  # median_cnv$chr <- chr_infos$chr
  # median_cnv$start <- chr_infos$start
  # median_cnv$end <- chr_infos$end
  # head(median_cnv)
  median_cnv <- median_cnv[unique(chrs$chr_desc),]
  # dim(median_cnv) # SA609
  # rownames(median_cnv)[1]
  # ls_features <- c()
  # for(i in seq(dim(median_cnv)[1])){
  #   if(sum(is.na(median_cnv[i,]))==0){
  #     ls_features <- c(ls_features, i)
  #   }
  # }
  # length(ls_features)
  # median_cnv <- median_cnv[ls_features,]
  # dim(median_cnv)
  res <- compute_dist_mat(median_cnv, save_dir, use_hamming = F)
  head(res$dist_to_median)
  p <- plot_heatmap_genotype(res$dist_to_median, distance_type, datatag, save_dir)
  res$hm <- p
  dim(res$out_mtx)
  # if(calcul_distance){
  #   if(distance_type=='Hamming'){
  #     dis_mtx <- compute_dist_mat(median_cnv, save_dir, use_hamming = T)
  #     p <- plot_heatmap_genotype(dis_mtx, distance_type, datatag, save_dir)
  #   }else if(distance_type=='Manhattan'){
  #     dis_mtx <- compute_dist_mat(median_cnv, save_dir, use_hamming = F)
  #     head(dis_mtx)
  #     dim(dis_mtx)
  #     p <- plot_heatmap_genotype(dis_mtx, distance_type, datatag, save_dir)
  #   }else{ # get results for both case
  #     dis_mtx1 <- compute_dist_mat(median_cnv, save_dir, use_hamming = T)
  #     p <- plot_heatmap_genotype(dis_mtx1, 'Hamming', datatag, save_dir)
  #     dis_mtx2 <- compute_dist_mat(median_cnv, save_dir, use_hamming = F)
  #     p <- plot_heatmap_genotype(dis_mtx2, 'Manhattan', datatag, save_dir)
  #   }
  #   
  # }
  return(res)
}
detect_CN_change <- function(de, median_cnv, meta_cells, save_dir){
  meta_cells <- meta_cells %>%
    dplyr::mutate(
      ms = case_when(
        mainsite =="Metastasis" ~ "Met",
        TRUE ~ "Pri"
      ),
      desc=paste0(ms,': clone ',clone_id,' (',nb_cells,')')
    )
 
  for(i in seq(1:nrow(de))){
    cl1 <- de$cl1[i]
    cl2 <- de$cl2[i]
    t1 <-  de$t1[i]
    t2 <-  de$t2[i]
    # desc1 <- paste0(tolower(t1)," cells (n=",meta_cells[meta_cells$clone_id==cl1 & meta_cells$mainsite==t1,'nb_cells'],")")
    # desc2 <- paste0(tolower(t2)," cells (n=",meta_cells[meta_cells$clone_id==cl2 & meta_cells$mainsite==t2,'nb_cells'],")")
    # View(meta_cells)
    desc1 <- meta_cells %>%
      dplyr::filter(clone_id==cl1) %>% # & mainsite==t1
      dplyr::pull(desc)
    desc2 <- meta_cells %>%
      dplyr::filter(clone_id==cl2) %>% # & mainsite==t2
      dplyr::pull(desc)
    print(desc1)
    print(desc2)
    desc1 <- paste(desc1, collapse = '; ')
    desc2 <- paste(desc2, collapse = '; ')
    # viz_cn_change(median_cnv, c(paste0(cl1,'_',t1),desc1), c(paste0(cl2,'_',t2),desc2), 
    #               paste0(save_dir,cl1,'_',t1,'_vs_',cl2,'_',t2,'/'))
    viz_cn_change(median_cnv, c(paste0(cl1,'_',t1),desc1), c(paste0(cl2,'_',t2),desc2), 
                  paste0(save_dir,cl1,'_',t1,'_vs_',cl2,'_',t2,'/'))
    # obs_clone1 <- c(paste0(cl1,'_',t1),desc1)
    # obs_clone2 <- c(paste0(cl2,'_',t2),desc2)
    # "K_Primary" %in% colnames(median_cnv)
  }
}

detect_CN_change_v1 <- function(de, median_cnv, meta_cells, save_dir){
  for(i in seq(1:nrow(de))){
    cl1 <- de$cl1[i]
    cl2 <- de$cl2[i]
    t1 <-  de$t1[i]
    t2 <-  de$t2[i]
    
    desc1 <- paste0(tolower(t1)," (n=",meta_cells[meta_cells$clone_id==cl1 & meta_cells$mainsite==t1,'nb_cells'],")")
    desc2 <- paste0(tolower(t2)," cells (n=",meta_cells[meta_cells$clone_id==cl2 & meta_cells$mainsite==t2,'nb_cells'],")")
    print(desc1)
    print(desc2)
    viz_cn_change(median_cnv, c(paste0(cl1,'_',t1),desc1), c(paste0(cl2,'_',t2),desc2), 
                  paste0(save_dir,cl1,'_',t1,'_vs_',cl2,'_',t2,'/'))
    
  }
}
get_median_genotype_v3 <- function(copynumber_fn, 
                                   datatag, save_dir,
                                   cellclone_fn=NULL, library_grouping_fn=NULL){
  if(is.null(cellclone_fn)){
    cellclone_fn <- paste0(results_dir,'cell_clones.csv')  
  }
  if(is.null(library_grouping_fn)){
    library_grouping_fn <- paste0(results_dir,'library_groupings.csv')  
  }
  if(!dir.exists(save_dir)){
    dir.create(save_dir)
  }
  
  # copynumber <- read.csv(copynumber_fn, header=T, row.names = 1, check.names = F,stringsAsFactors = FALSE)
  copynumber <- as.data.frame(data.table::fread(copynumber_fn))
  rownames(copynumber) <- copynumber$V1
  copynumber$V1 <- NULL
  # cell_clones contain 2 columns of cell_id, and clone_id
  # ex:           cell_id                      clone_id
  # 1   SA535X4XB02498-A98163A-R09-C11          C
  cell_clones <- data.table::fread(cellclone_fn) %>% as.data.frame()
  dim(cell_clones)
  metasample <- data.table::fread(library_grouping_fn) %>% as.data.frame()
  dim(metasample)
  head(metasample)
  cell_clones <- cell_clones %>%
    dplyr::filter(!clone_id %in% c('None','unassigned'))
  
  
  metasample <- metasample %>%
    dplyr::rename(library_id=grouping) %>%
    dplyr::select(library_id, mainsite, sample_id)
  print(dim(copynumber))
  
  
  cell_clones$library_id <- get_library_id(cell_clones$cell_id)
  cell_clones$sample_id <- get_sample_id(cell_clones$cell_id)
  cell_clones <- cell_clones %>% left_join(metasample, by=c("library_id","sample_id"))
  meta_cells <- cell_clones
  dim(meta_cells)
  dim(meta_cells)
  meta_cells <- meta_cells %>%
    dplyr::group_by(clone_id, mainsite) %>%
    dplyr::summarise(nb_cells=n()) %>%
    dplyr::ungroup()
  # View(meta_cells)
  
  meta_cells <- meta_cells %>%
    dplyr::filter(nb_cells>=20 & clone_id !='None')
  dim(meta_cells)
  colnames(meta_cells)
  write.csv(meta_cells, paste0(save_dir,'meta_cells.csv'), quote = F, row.names = F)
  
  copynumber <- copynumber[,colnames(copynumber)[colnames(copynumber) %in% cell_clones$cell_id]]
  copynumber$chr_desc <- rownames(copynumber)
  cnv <- copynumber %>%
    pivot_longer(!chr_desc, names_to = "cell_id", values_to = "copy_number")
  
  print(dim(cnv))
  cnv <- cnv %>% inner_join(cell_clones, by=c("cell_id"))
  length(unique(cnv$cell_id))
  dim(cnv)
  colnames(cnv)
  cnv <- cnv %>% inner_join(meta_cells, by=c("clone_id","mainsite"))
  dim(cnv)
  
  cnv$clone_label <- paste0(cnv$clone_id,'_',cnv$mainsite)
  length(unique(cnv$clone_label))
  
  # get median genotype 
  print("Get median genotype")
  # stat_cnv <- cnv %>%
  #   dplyr::group_by(clone_id, chr_desc) %>%
  #   dplyr::summarise(median_cn=median(copy_number),
  #                    mode_cn=calc_mode(copy_number))
  # 
  stat_cnv <- cnv %>%
    dplyr::group_by(clone_label, chr_desc) %>%
    dplyr::summarise(median_cn=median(copy_number)) %>%
    dplyr::ungroup()
  dim(stat_cnv)
  # median_cnv <- stat_cnv %>%
  #   dplyr::select(chr_desc, median_cn, clone_id, clone_label) %>%
  #   group_by(clone_label) 
  head(stat_cnv)
  
  median_cnv_pivot <- stat_cnv %>%
    pivot_wider(names_from = clone_label, values_from = median_cn)
  
  median_cnv_pivot <- as.data.frame(median_cnv_pivot)
  write.csv(median_cnv_pivot, paste0(save_dir,'median_cnv.csv'), quote = F, row.names = F)
  data.table::fwrite(stat_cnv, paste0(save_dir,'stat_median_cnv.csv'), quote = F, row.names = F)
  
}  


