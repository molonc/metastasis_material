require("yaml")
require("tools")
require("phytools")
require("igraph")
require("RColorBrewer")
require("outliers")
require("data.table")
require("dplyr")
require("optparse")


# source("/datadrive/tree_viz/ss_viz/general_utils.R")
# source("/datadrive/tree_viz/ss_viz/tree_and_heatmap_utils.R")

newick_2_edge_list <- function(newick_path) {
  funny_str <- generate_random_funny_string()
  res_path = paste0(file_path_sans_ext(newick_path), funny_str,  '.csv')
  the_tree = read.newick(newick_path)
  edges = the_tree$edge
  the_tree$node.label[1] <- 'root'
  the_names <- c(the_tree$tip.label, the_tree$node.label)
  named_edges <- data.frame(source=the_names[edges[, 1]], target=the_names[edges[, 2]], stringsAsFactors = F)
  named_edges[1, ]$source
  named_edges$source <- gsub('cell_', '', named_edges$source)
  named_edges$target <- gsub('cell_', '', named_edges$target)
  write.csv(named_edges, res_path, row.names = F)
  res_path
}

fitclone_plot_tree_heatmap <- function(g, datatag, edge_list_path, mat, tag='', out_dir=NULL, just_tree_loci=NULL) {
    fitclone_plot_tree_heatmap_more(g = g, datatag = datatag, edge_list_path = edge_list_path, mat=mat, tag = tag, out_dir = out_dir, just_tree_loci=just_tree_loci)
}

fitclone_plot_tree_heatmap_more <- function(g, datatag, edge_list_path, mat, tag='', out_dir=NULL, aug_cut=NULL, chr_filter=NULL, timepoint_filter=NULL, just_tree_loci=NULL, save_file=TRUE, chr_font_size = 3.8) {
    if (is.null(g)) {
        g <- read_ltm_tree(edge_list_path)
    }
    
    if (is.null(out_dir)) {
        out_dir <- sprintf('~/Desktop/SC-1311/%s/plots', datatag)
        dir.create(out_dir, showWarnings = F)
    }
    
    tcc <- get_tcc(g)
    n_total_cells <- sum(tcc$Freq)
    
    print('Generating heatmap object...')
    pp <- tree_and_heatmap(g = g, datatag = datatag, edge_list_path = edge_list_path, 
                           corrupt_input_mat=mat,
                           the_cut = NULL, ts=NULL, outpath = NULL, tag = '', aug_cut = aug_cut, 
                           chr_filter=chr_filter, time_point_filter = timepoint_filter, 
                           just_tree_loci = just_tree_loci, chr_font_size = chr_font_size)
    p_res <- pp
    n_total_cells <- 0
    if (save_file) {
        pp <- pp + 
            ggtitle(get_title_str(datatag, n_total_cells), subtitle = get_ID_from_edge_list_path(edge_list_path)) + 
            theme(legend.text=element_text(size=15), 
                        plot.title = element_text(size = 28), 
                        plot.subtitle=element_text(size=20, face="italic", color="black"))
        
        coef = 1
        height = 8*coef
        width = 14 * coef
        sp <- exp_path; if (is.null(sp)) sp <- out_dir; 
        out_path <- file.path(sp, get_file_name_postfix(datatag, edge_list_path, tag, '.png'))
        print('Saving the heatmap object...')
        ggsave(plot = pp, filename = out_path,  width = width, height = height, units = "in", limitsize = FALSE)
    } else {
        out_path = NA
    }
    list(out_path=out_path, p=p_res)
}

tree_and_heatmap <- function(g, datatag, edge_list_path, the_cut=NULL, aug_cut=NULL, ts=NULL, outpath=NULL, tag='', clone_dic=NULL, time_point_filter=NULL, not_in_cut=FALSE, corrupt_input_mat=NULL, just_tree_loci = NULL, add_tp_indicators = FALSE, chr_filter=NULL, cells_to_keep=NULL, chr_font_size = 3.8) {
    if (is.null(g)) {
        g <- read_ltm_tree(edge_list_path)
    }
    
    if (is.null(outpath)) {
        outpath <- file.path('~/Desktop/SC-1311/', datatag, '/phylogeny/temp/', paste0('tree_', tag, '_', generate_random_str(), '.png'))
    }
    
    dir.create(dirname(outpath), recursive = T, showWarnings = F)
    
    # Load CN data
    if (!is.null(corrupt_input_mat)) {
        mat = corrupt_input_mat
        rownames(mat) <- gsub('locus_', '', rownames(mat))
        mat <- as.data.frame(mat)
    } else {
        mat <- load_new_cn_data(datatag)
    }
    
    if (!is.null(chr_filter)) {
        chr <- gsub('([0-9]+|X|Y)_[0-9]+_[0-9]+', '\\1', rownames(mat))
        mat <- mat[chr %in% chr_filter, ]
    }
    
    mat <- sort_mat_by_bins(the_mat = mat)  
    mat <- remove_pad_bins(mat)
    
    # @Sohrab: update for the aug_cut
    if (!is.null(the_cut)) {
        desc <- get_decendents(the_cut, the_graph = g, min_cell_per_clone = 1)
        res <- retrieve_time_point_per_clade(g, desc, datatag, FALSE)
        
        if (not_in_cut) {
            clusterCut <- ltm_clust_to_cuttree(res, seq_along(names(res)))
            mat <- mat[, !(colnames(mat) %in% unlist(res))]
        } else {
            mat <- mat[, colnames(mat) %in% names(clusterCut)]
            mat <- mat[, match(names(clusterCut), colnames(mat))]  
        }
    }
    
    if (!is.null(aug_cut)) {
        #TODO
    }
    
    if (!is.null(time_point_filter)) {
        mat <- mat[, parse_cell_names(colnames(mat)) %in% time_point_filter]
    }
     
    # The tree
    res <- convert_edge_list_to_ape(edge_list_path)
    tree <- res$tree
    node_names <- res$node_names
    the_height <- find_tail_height(g = g, edge_list_path = edge_list_path)$h_cut

    if (!is.null(the_height)) {
        tree_res <- trim_tree_before_height(tree, node_names, g, the_height, 20)
        tree <- tree_res$tree
    }
 
    # Add the heatmap
    sub_mat <- heatmap_mat_from_bin_cellID(mat, update_colnames = FALSE)
    
    if (!is.null(just_tree_loci)) {
        if (just_tree_loci == TRUE) {
            the_loci <- grep('locus', node_names$node_name, value = T)
            the_loci <- gsub('locus_', '', the_loci)
            sub_mat <- sub_mat[, colnames(sub_mat) %in% the_loci]
        }
    }
    
    colnames_filter <- strip_chr_names(colnames(sub_mat))
    ss <- which(colnames_filter!='')
    colnames_filter <- data.frame(from = colnames(sub_mat)[ss], new_label=colnames_filter[ss])
    
    # Add white borders between chromosomes
    genotype = as.tibble(sub_mat)
    white_border_size <- 0
    if (!is.null(just_tree_loci)) white_border_size <- 0
        
    if (white_border_size > 0) {
        for (cname in colnames_filter$from[-c(1)]) {
            for (i in 1:white_border_size) {
                new_name = gsub('([0-9]+|X|Y)_(.*)', paste0('\\1_', i-1, '_\\2'), cname)
                genotype <- add_column(genotype, !!new_name := NA, .before = cname)
            }
        }
    }
    
    # Convert to character format
    d = dim(genotype)
    cn = colnames(genotype)
    rn = rownames(sub_mat)
    genotype = as.matrix(genotype)
    genotype <- as.character(genotype)
    genotype <- matrix(data=genotype, nrow = d[1], ncol = d[2])
    colnames(genotype) <- cn
    rownames(genotype) <- rn
    
    colnames_level = unique(colnames(genotype))
    
    # Add colours for CNV value
    cn_colours <- legacy_colours_for_CN_heatmap()
    names(cn_colours) <- paste0(seq_along(cn_colours)-1)
    
    tip.loci <- grep('locus', tree$tip.label, value= T)
    tree <- drop.tip(tree, tip.loci)
    
    # Plot
    if (!is.null(aug_cut)) {
        xtree <- groupOTU(tree, .node=aug_cut)
        p <- ggtree(xtree, aes(color=group))
    } else {
        p <- ggtree(tree)
    }
    
    offset <- 10
    
    breaks <- factor(x = names(cn_colours), levels = names(cn_colours), ordered = T)
    
    p <- sos_heat(p = p, data = genotype, offset=offset, width=8, colnames=TRUE, color=NULL, colnames_level=colnames_level, colnames_filter=colnames_filter, font.size = chr_font_size, colnames_offset_y = -40) + 
        scale_fill_manual(breaks=breaks, values=cn_colours, guide = guide_legend(nrow=1, direction='horizontal', label.position = 'bottom'))
    
    # Add the cluster annotations
    if (!is.null(the_cut) | !is.null(aug_cut)) {
        tcc <- get_tcc(g)
        ts <- get_edge_data(g, tcc, datatag, edge_list_path)
        the_cut <- get_the_cut(datatag, edge_list_path, g=g)
        aug_cut <- add_siblings_cut(datatag = datatag, edge_list_path = edge_list_path, the_cut = the_cut, the_g = g)
        new_clades <- setdiff(names(aug_cut), unlist(the_cut))
        candi_edge <- get_candi_edge(the_cut = the_cut, aug_cut = aug_cut, ts = ts, tcc = tcc, new_clades = new_clades)
        
        if (datatag == 'SA666') {
            # Handles the missing clone
            clone_dic <- cosmic_clone_dic(datatag = 'SA1000', edge_list_path = get_edge_list_path_for_datatag('SA1000'))
        } else {
            clone_dic <- universal_clone_dic(ts = ts, aug_cut = aug_cut, new_clades = new_clades, tcc = tcc)
        }
        
        p <- annoate_tree_aug(p = p, candiate_edges = candi_edge, g = g, aug_cut = aug_cut, the_height = the_height, clone_dic = clone_dic, draw_bar = FALSE)
    }
    
    ppp <- p + theme(legend.position="top")
    ppp
}


# ## Launch
# option_list <- list(make_option(c("-n", "--newick"), type="character", default=NULL, help="input_newick", metavar="character"),
#                     make_option(c("-t", "--datatag"), type="character", default=NULL, help="library_id", metavar="character")          
# )
# 
# opt_parser <- OptionParser(option_list=option_list)
# opt <- parse_args(opt_parser)
# edge_list_path <- newick_2_edge_list(opt$newick)
# g <- read_ltm_tree(edge_list_path)
# fitclone_plot_tree_heatmap(g, opt$datatag, edge_list)


