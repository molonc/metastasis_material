#!/usr/bin/env Rscript

suppressMessages(library(ape))
suppressMessages(library(argparse))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggtree))
suppressMessages(library(gtools))
suppressMessages(library(RColorBrewer))

cn_colours <- structure(
    c(
        "#3182BD", "#9ECAE1", "#CCCCCC", "#FDCC8A", "#FC8D59", "#E34A33",
        "#B30000", "#980043", "#DD1C77", "#DF65B0", "#C994C7", "#D4B9DA"
    ),
    names=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11+")
)

clone_palette_20 <- c(
    "#be5f72", "#d74058", "#dc4229", "#a6552c", "#df956f", "#e47a33",
    "#d49f34", "#836e2c", "#b2ad5a", "#92b539", "#4c7d38", "#4dc041",
    "#5dba7f", "#47b8c3", "#6280ca", "#7b57db", "#ce8bd1", "#934f94",
    "#cb48cb", "#d74391"
)
clone_none_black <- "#1B1B1B"

get_args <- function() {
    p <- ArgumentParser(description="Plot cell copynumber heatmap with tree")
    p$add_argument("--outdir", "-od", help="output dir")
    p$add_argument("--tree", "-t", help="cell newick tree file")
    p$add_argument("--copynumber", "-cn", help="cell copynumber tsv file")
    p$add_argument("--pdf", "-o", help="output plot pdf")

    p$add_argument("-c", "--clones", help="cell clone tsv file")
    p$add_argument(
        "-n", "--normalize-ploidy", action="store_true",
        help="normalize all cell ploidy to 2"
    )
    p$add_argument(
        "-b", "--branch-lengths", type="integer",
        help="set all tree branch lengths to this value"
    )
    p$add_argument("--grouping_file", default=NULL, help="file with sample grouping")
    p$add_argument("--datatag", default='SA', type="character", help="datatag")

    return(p$parse_args())
}

read_tsv <- function(fn, ...) {
    df <- read.delim(fn, check.names=FALSE, stringsAsFactors=FALSE, sep=",", ...)
    return(df)
}

calc_state_mode <- function(states) {
    state_levels <- unique(states)
    state_mode <- state_levels[
        which.max(tabulate(match(states, state_levels)))
    ]
    return(state_mode)
}

normalize_cell_ploidy <- function(copynumber) {
    cell_ids <- colnames(copynumber)
    cell_ids <- cell_ids[!(cell_ids %in% c("chr", "start", "end", "width"))]

    for(cell_id in cell_ids) {
        state_mode <- calc_state_mode(copynumber[[cell_id]])
        copynumber[[cell_id]] <- as.integer(ceiling(
            copynumber[[cell_id]] / (state_mode / 2)
        ))
    }
    return(copynumber)
}

format_tree <- function(tree, brlen) {
    locus_tips <- grep('locus', tree$tip.label, value=TRUE)
    tree <- drop.tip(tree, locus_tips)

    if (!is.null(brlen)) {
        tree <- compute.brlen(tree, brlen)
    }

    tree$tip.label <- gsub('cell_', '', tree$tip.label)

    return(tree)
}

get_clone_members <- function(clones) {
    clone_members <- list()
    for(c in unique(clones$clone_id)) {
        if(c != "None") {
            clone_members[[c]] <- clones[clones$clone_id == c, "cell_id"]
        }
    }
    return(clone_members)
}

make_clone_palette <- function(levels) {
    if (length(levels) <= 12) {
        pal <- brewer.pal(max(length(levels), 3), "Set3")
    } else if (length(levels) <= 20) {
        pal <- clone_palette_20
    } else {
        pal <- colorRampPalette(brewer.pal(8, "Set2"))(length(levels))
        print("WARNING: more clones than palette can accomodate!")
    }
    names(pal) <- levels
    pal <- pal[levels]
    return(pal)
}

make_tree_ggplot <- function(tree, clones) {
    if(!is.null(clones)) {
        clone_members <- get_clone_members(clones)
        tree <- groupOTU(tree, clone_members)

        clone_levels <- mixedsort(unique(clones$clone_id))
        clone_pal <- make_clone_palette(clone_levels)
        clone_pal[["0"]] <- clone_none_black

        tree_aes <- aes(x, y, colour=group)
    } else {
        tree_aes <- aes(x, y)
    }

    p <- ggplot(tree, tree_aes) +
        geom_tree(size=0.25) +
        coord_cartesian(expand=FALSE) +
        ylim(0.5, length(tree$tip.label) + 0.5) +
        theme_void()

    if(!is.null(clones)) {
        p <- p + scale_colour_manual(values=clone_pal)
    }

    return(p)
}

make_discrete_palette <- function(pal_name, levels) {
    if (length(levels) <= 8) {
        pal <- brewer.pal(max(length(levels), 3), pal_name)
        } else if (length(levels) <= 12) {
        pal <- brewer.pal(max(length(levels), 3), "Set3")
    } else if (length(levels) <= 20) {
        pal <- clone_palette_20
    } else {
        pal <- clone_palette_20
        print("WARNING: more clones than palette can accomodate!")
    }
    names(pal) <- levels
    pal <- pal[levels]
    return(pal)
}

format_copynumber_values <- function(copynumber) {
    copynumber[copynumber > 11] <- 11
    for(col in colnames(copynumber)) {
        values <- as.character(copynumber[, col])
        values[values == "11"] <- "11+"
        copynumber[, col] <- values
    }
    return(copynumber)
}

space_copynumber_columns <- function(copynumber, spacer_cols) {
    chroms <- sapply(strsplit(colnames(copynumber), ":"), function(x) x[[1]])
    spacer <- as.data.frame(matrix(
        data=NA, nrow=nrow(copynumber), ncol=spacer_cols
    ))
    chrom_copynumber_dfs <- list()
    for(chrom in mixedsort(unique(chroms))) {
        chrom_copynumber <- copynumber[, chroms == chrom, drop=FALSE]
        chrom_copynumber_dfs <- c(chrom_copynumber_dfs, list(chrom_copynumber))
        chrom_copynumber_dfs <- c(chrom_copynumber_dfs, list(spacer))
    }
    chrom_copynumber_dfs[length(chrom_copynumber_dfs)] <- NULL
    copynumber <- do.call(cbind, chrom_copynumber_dfs)

    return(copynumber)
}

get_ordered_cell_ids <- function(tree_plot_dat) {
    return(rev(arrange(tree_plot_dat[tree_plot_dat$isTip, ], y)$label))
}

format_copynumber <- function(copynumber, tree_plot_dat, spacer_cols=20) {
    if (!("chr" %in% colnames(copynumber))) {
        loci <- sapply(rownames(copynumber), strsplit, "_")
        copynumber$chr <- unname(sapply(loci, '[[', 1))
        copynumber$start <- as.numeric(unname(sapply(loci, '[[', 2)))
        copynumber$end <- as.numeric(unname(sapply(loci, '[[', 3)))
        copynumber$width <- (copynumber$end - copynumber$start + 1)
    }
    copynumber$chr <- gsub("chr", "", copynumber$chr)
    copynumber <- arrange(copynumber, as.numeric(chr), chr, start)

    rownames(copynumber) <- paste0(
        copynumber$chr, ":", copynumber$start, ":", copynumber$end
    )
    copynumber <- subset(copynumber, select=-c(chr, start, end, width))
    copynumber <- as.data.frame(t(copynumber))

    copynumber <- copynumber[get_ordered_cell_ids(tree_plot_dat), ]

    copynumber <- format_copynumber_values(copynumber)
    copynumber <- space_copynumber_columns(copynumber, spacer_cols)

    return(copynumber)
}


format_copynumber_v2 <- function(copynumber, clones, spacer_cols=20) {
    if (!("chr" %in% colnames(copynumber))) {
        loci <- sapply(rownames(copynumber), strsplit, "_")
        copynumber$chr <- unname(sapply(loci, '[[', 1))
        copynumber$start <- as.numeric(unname(sapply(loci, '[[', 2)))
        copynumber$end <- as.numeric(unname(sapply(loci, '[[', 3)))
        copynumber$width <- (copynumber$end - copynumber$start + 1)
    }
    copynumber$chr <- gsub("chr", "", copynumber$chr)
    copynumber$chr_lb <- copynumber$chr
    copynumber$chr_lb <- ifelse(copynumber$chr_lb=="X",23,
                             ifelse(copynumber$chr_lb=="Y",24,copynumber$chr_lb))
    
    # copynumber <- arrange(copynumber, as.numeric(chr), chr, start)
    copynumber <- arrange(copynumber, as.numeric(chr_lb), chr, start)
    
    rownames(copynumber) <- paste0(
        copynumber$chr, ":", copynumber$start, ":", copynumber$end
    )
    copynumber <- subset(copynumber, select=-c(chr_lb, chr, start, end, width))
    copynumber <- as.data.frame(t(copynumber))
    
    # copynumber <- copynumber[get_ordered_cell_ids(tree_plot_dat), ]
    copynumber <- copynumber[clones$cell_id, ]
    
    copynumber <- format_copynumber_values(copynumber)
    copynumber <- space_copynumber_columns(copynumber, spacer_cols)
    
    return(copynumber)
}


format_clones <- function(clones, tree_plot_dat) {
    tree_cells <- get_ordered_cell_ids(tree_plot_dat)
    clones <- merge(clones, data.frame(cell_id=tree_cells), all=TRUE)
    clones[is.na(clones$clone_id), "clone_id"] <- "None"

    clone_counts <- clones %>% group_by(clone_id) %>% summarise(count=n())

    for(i in 1:nrow(clone_counts)) {
        clone <- unlist(clone_counts[i, "clone_id"], use.names=FALSE)
        clone_count <- unlist(clone_counts[i, "count"], use.names=FALSE)
        clone_label <- paste0(clone, " (", clone_count, ")")
        clones[clones$clone_id == clone, "clone_label"] <- clone_label
    }

    rownames(clones) <- clones$cell_id
    clones <- clones[tree_cells, ]

    return(clones)
}

# Hoa format clone, without tree
format_clones_v2 <- function(clones, copynumber, grouping_file, order_type='mainsite') {
    cells_use <- colnames(copynumber)
    # dim(clones)
    # dim(copynumber)
    # summary(as.factor(clones$clone_id))
    clones <- merge(clones, data.frame(cell_id=cells_use), all=TRUE)
    clones[is.na(clones$clone_id), "clone_id"] <- "None"
    
    clone_counts <- clones %>% group_by(clone_id) %>% summarise(count=n())
    
    for(i in 1:nrow(clone_counts)) {
        clone <- unlist(clone_counts[i, "clone_id"], use.names=FALSE)
        clone_count <- unlist(clone_counts[i, "count"], use.names=FALSE)
        clone_label <- paste0(clone, " (", clone_count, ")")
        clones[clones$clone_id == clone, "clone_label"] <- clone_label
    }
    
    rownames(clones) <- clones$cell_id
    
    clones <- clones[cells_use, ]
    clones$library_id <- get_library_id(cells_use)
    clones$sample_id <- get_sample_id(cells_use)
    groupings <- read_tsv(grouping_file)
    colnames(groupings)[which(names(groupings) == "grouping")] <- "library_id"
    clones <- left_join(clones, groupings, by = c("library_id","sample_id"))
    
    if(order_type=='pdxid'){
        print('Ordering data by pdxid')
        clones <- clones[order(clones$pdxid, clones$clone_id, decreasing = F), ]
        if(length(unique(clones$pdxid))==1){
            clones <- clones[order(clones$clone_id, decreasing = F), ]
        }
    } else if(order_type=='mainsite'){
        print('Ordering data by main site')
        clones <- clones[order(clones$mainsite, clones$clone_id, decreasing = T), ]
    } else if(order_type=='origin'){
        print('Ordering data by organ')
        clones <- clones[order(clones$origin, clones$clone_id, decreasing = F), ]
    } else{
        clones <- clones[order(clones$clone_id, decreasing = F), ]
    }
    
    return(clones)
}


make_corrupt_tree_heatmap <- function(tree_ggplot) {
    tree_annot_func = AnnotationFunction(
        fun=function(index) {
            pushViewport(viewport(height=1))
            grid.draw(ggplotGrob(tree_ggplot)$grobs[[5]])
            popViewport()
        },
        var_import=list(tree_ggplot=tree_ggplot),
        width=unit(4, "cm"),
        which="row"
    )
    tree_annot <- HeatmapAnnotation(
        tree=tree_annot_func, which="row", show_annotation_name=FALSE
    )

    n_cells <- sum(tree_ggplot$data$isTip)
    tree_hm <- Heatmap(matrix(nc=0, nr=n_cells), left_annotation=tree_annot)

    return(tree_hm)
}

find_largest_contiguous_group <- function(x) {
    starts <- c(1, which(diff(x) != 1 & diff(x) != 0) + 1)
    ends <- c(starts[-1] - 1, length(x))
    largest <- which.max(ends - starts + 1)
    return(x[starts[largest]:ends[largest]])
}

get_clone_label_pos <- function(clones) {
    clone_label_pos <- list()
    for(clone in unique(clones$clone_id)) {
        if(!grepl("None", clone)) {
            clone_idx <- which(clones$clone_id == clone)
            clone_idx <- find_largest_contiguous_group(clone_idx)
            clone_label_pos[[as.character(clone)]] <-
                as.integer(round(mean(clone_idx)))
        }
    }
    return(clone_label_pos)
}

get_library_id <- function(cell_ids) {
    labels <- sapply(strsplit(cell_ids, "-"), function(x) {
        return(x[2])
    })
    return(as.character(labels))
}

get_sample_id <- function(cell_ids) {
    labels <- sapply(strsplit(cell_ids, "-"), function(x) {
        return(x[1])
    })
    return(as.character(labels))
}

# get_groupings <- function(cell_ids, grouping_file) {
#     library_labels <- get_library_id(cell_ids)
#     sample_ids <- get_sample_id(cell_ids)
#     groupings <- read_tsv(grouping_file)
# 
#     tmp <- as.data.frame(library_id=library_labels, sample_id=sample_ids, stringsAsFactors=FALSE)
#     tmp <- left_join(tmp, groupings, by = c("library_id","sample_id"))
#     return(tmp$grouping)
# }

get_library_grouping <- function(cell_ids, grouping_file) {
    library_ids <- get_library_id(cell_ids)
    groupings <- read_tsv(grouping_file)
    colnames(groupings)[which(names(groupings) == "grouping")] <- "library_id"
    sample_ids <- get_sample_id(cell_ids)
    tmp <- data.frame(cell_id=cell_ids, library_id=library_ids, sample_id=sample_ids, stringsAsFactors=FALSE)
    tmp <- left_join(tmp, groupings, by = c("library_id","sample_id"))
    return(tmp)
}

make_left_annot <- function(copynumber, clones, grouping_file) {
    annot_colours <- list()
    annot_cols <- 3

    library_labels <- get_library_id(rownames(copynumber))
    library_levels <- mixedsort(unique(library_labels))
    annot_colours$Sample <- make_discrete_palette("Set2", library_levels)
    library_legend_rows <- 10

    grouping_labels <- NULL
    annot_colours$Groupings <- NULL
    grouping_legend <- NULL
    if (!is.null(grouping_file)) {
        # grouping_labels <- get_groupings(rownames(copynumber), grouping_file)
        grouping_labels <- get_library_id(rownames(copynumber))
        grouping_levels <- mixedsort(unique(grouping_labels))
        annot_colours$Groupings <- make_discrete_palette("Set2", grouping_levels)
        grouping_legend <- list(nrow=10)
        annot_cols <- 4
    }

    if(!is.null(clones)) {
        clone_levels <- unique(clones$clone_label)
        clone_level_none <- clone_levels[grepl("None", clone_levels)]
        clone_levels <- mixedsort(clone_levels[!grepl("None", clone_levels)])

        clone_pal <- make_clone_palette(clone_levels)
        if(length(clone_level_none > 0)) {
            clone_pal[[clone_level_none]] <- clone_none_black
        }
        annot_colours$Clone <- clone_pal

        clone_label_generator <- function(index) {
            clone_label_pos <- get_clone_label_pos(clones)
            y_pos <- 1 - unlist(clone_label_pos) / nrow(clones)
            grid.text(
                names(clone_label_pos), 0.5, y_pos,
                just=c("centre", "centre")
            )
        }

        clone_legend_rows <- 10
        if(length(clone_levels) > 10) {
            clone_legend_rows <- round(sqrt(length(clone_levels) * 4))
        }

        left_annot <- HeatmapAnnotation(
            Clone=clones$clone_label, clone_label=clone_label_generator,
            Sample=library_labels, Groupings=grouping_labels,
            col=annot_colours, show_annotation_name=c(TRUE, FALSE, TRUE),
            which="row", annotation_width=unit(rep(0.4, annot_cols), "cm"),
            annotation_legend_param=list(
                Clone=list(nrow=clone_legend_rows),
                Sample=list(nrow=library_legend_rows),
                Groupings=grouping_legend
            )
        )
    } else {
        left_annot <- HeatmapAnnotation(
            Sample=library_labels, col=annot_colours,
            which="row", simple_anno_size=unit(0.4, "cm"),
            annotation_legend_param=list(
                Sample=list(nrow=library_legend_rows)
            )
        )
    }

    return(left_annot)
}


make_left_annot_v2 <- function(copynumber, clones, grouping_file) {
    annot_colours <- list()
    annot_cols <- 3
    
    library_labels <- get_library_id(rownames(copynumber))
    library_levels <- mixedsort(unique(library_labels))
    # annot_colours$Sample <- make_discrete_palette("Set2", library_levels)
    # library_legend_rows <- 10
    
    grouping_labels <- NULL
    annot_colours$Groupings <- NULL
    grouping_legend <- NULL
    if (!is.null(grouping_file)) {
        # grouping_labels <- get_groupings(rownames(copynumber), grouping_file)
        grouping_labels <- get_library_id(rownames(copynumber))
        grouping_levels <- mixedsort(unique(grouping_labels))
        # annot_colours$Groupings <- make_discrete_palette("Set2", grouping_levels)
        grouping_legend <- list(nrow=10)
        annot_cols <- 5
        
        lib_group_meta <- get_library_grouping(rownames(copynumber), grouping_file)
        mainsite_labels <- lib_group_meta$mainsite
        mainsite_levels <- mixedsort(unique(mainsite_labels))
        annot_colours$MainSite <- make_discrete_palette("Set1", mainsite_levels)
        library_legend_rows <- 10
        pdx_labels <- lib_group_meta$pdxid
        pdx_levels <- mixedsort(unique(pdx_labels))
        annot_colours$Pdx <- make_discrete_palette("Accent", pdx_levels)
        # library_legend_rows <- 10
    }
    
    if(!is.null(clones)) {
        clone_levels <- unique(clones$clone_label)
        clone_level_none <- clone_levels[grepl("None", clone_levels)]
        clone_levels <- mixedsort(clone_levels[!grepl("None", clone_levels)])
        
        clone_pal <- make_clone_palette(clone_levels)
        if(length(clone_level_none > 0)) {
            clone_pal[[clone_level_none]] <- clone_none_black
        }
        annot_colours$Clone <- clone_pal
        
        clone_label_generator <- function(index) {
            clone_label_pos <- get_clone_label_pos(clones)
            y_pos <- 1 - unlist(clone_label_pos) / nrow(clones)
            grid.text(
                names(clone_label_pos), 0.5, y_pos,
                just=c("centre", "centre")
            )
        }
        
        clone_legend_rows <- 10
        if(length(clone_levels) > 10) {
            clone_legend_rows <- round(sqrt(length(clone_levels) * 4))
        }
        # Hoa, modification
        annot_cols <- 4  # remove lib idx
        left_annot <- HeatmapAnnotation(
            Clone=clones$clone_label, clone_label=clone_label_generator,
            MainSite=mainsite_labels, Pdx=pdx_labels,
            col=annot_colours, show_annotation_name=c(TRUE, FALSE, TRUE),
            which="row", annotation_width=unit(rep(0.4, annot_cols), "cm"),
            annotation_legend_param=list(
                Clone=list(nrow=clone_legend_rows),
                MainSite=list(nrow=library_legend_rows),
                Pdx=list(nrow=library_legend_rows)
            )
        )
        # left_annot <- HeatmapAnnotation(
        #     Clone=clones$clone_label, clone_label=clone_label_generator,
        #     MainSite=mainsite_labels, Pdx=pdx_labels, Groupings=grouping_labels,
        #     col=annot_colours, show_annotation_name=c(TRUE, FALSE, TRUE),
        #     which="row", annotation_width=unit(rep(0.4, annot_cols), "cm"),
        #     annotation_legend_param=list(
        #         Clone=list(nrow=clone_legend_rows),
        #         MainSite=list(nrow=library_legend_rows),
        #         Pdx=list(nrow=library_legend_rows),
        #         Groupings=grouping_legend
        #     )
        # )
    } else {
        left_annot <- HeatmapAnnotation(
            MainSite=mainsite_labels, col=annot_colours,
            which="row", simple_anno_size=unit(0.4, "cm"),
            annotation_legend_param=list(
                MainSite=list(nrow=library_legend_rows)
            )
        )
    }
    
    return(left_annot)
}


# Add origin columns
make_left_annot_v3 <- function(copynumber, clones, grouping_file) {
    annot_colours <- list()
    annot_cols <- 3
    
    library_labels <- get_library_id(rownames(copynumber))
    library_levels <- mixedsort(unique(library_labels))
    # annot_colours$Sample <- make_discrete_palette("Set2", library_levels)
    # library_legend_rows <- 10
    
    grouping_labels <- NULL
    annot_colours$Groupings <- NULL
    grouping_legend <- NULL
    if (!is.null(grouping_file)) {
        # grouping_labels <- get_groupings(rownames(copynumber), grouping_file)
        grouping_labels <- get_library_id(rownames(copynumber))
        grouping_levels <- mixedsort(unique(grouping_labels))
        # annot_colours$Groupings <- make_discrete_palette("Set2", grouping_levels)
        grouping_legend <- list(nrow=10)
        annot_cols <- 5
        
        lib_group_meta <- get_library_grouping(rownames(copynumber), grouping_file)
        mainsite_labels <- lib_group_meta$mainsite
        mainsite_levels <- mixedsort(unique(mainsite_labels))
        annot_colours$MainSite <- make_discrete_palette("Set1", mainsite_levels)
        
        origin_labels <- lib_group_meta$origin
        origin_levels <- mixedsort(unique(origin_labels))
        annot_colours$Origin <- make_discrete_palette("Pastel1", origin_levels)
        
        library_legend_rows <- 10
        pdx_labels <- lib_group_meta$pdxid
        pdx_levels <- mixedsort(unique(pdx_labels))
        annot_colours$Pdx <- make_discrete_palette("Accent", pdx_levels)
        # library_legend_rows <- 10
    }
    
    if(!is.null(clones)) {
        clone_levels <- unique(clones$clone_label)
        clone_level_none <- clone_levels[grepl("None", clone_levels)]
        clone_levels <- mixedsort(clone_levels[!grepl("None", clone_levels)])
        
        clone_pal <- make_clone_palette(clone_levels)
        if(length(clone_level_none > 0)) {
            clone_pal[[clone_level_none]] <- clone_none_black
        }
        annot_colours$Clone <- clone_pal
        
        clone_label_generator <- function(index) {
            clone_label_pos <- get_clone_label_pos(clones)
            y_pos <- 1 - unlist(clone_label_pos) / nrow(clones)
            grid.text(
                names(clone_label_pos), 0.5, y_pos,
                just=c("centre", "centre")
            )
        }
        
        clone_legend_rows <- 10
        if(length(clone_levels) > 10) {
            clone_legend_rows <- round(sqrt(length(clone_levels) * 4))
        }
        # Hoa, modification
        # annot_cols <- 4  # remove lib idx
        annot_cols <- 5
        left_annot <- HeatmapAnnotation(
            Clone=clones$clone_label, clone_label=clone_label_generator,
            MainSite=mainsite_labels, Origin=origin_labels, Pdx=pdx_labels,
            col=annot_colours, show_annotation_name=c(TRUE, FALSE, TRUE, TRUE, TRUE),
            which="row", annotation_width=unit(rep(0.5, annot_cols), "cm"),
            annotation_legend_param=list(
                Clone=list(nrow=clone_legend_rows),
                MainSite=list(nrow=library_legend_rows),
                Origin=list(nrow=library_legend_rows),
                Pdx=list(nrow=library_legend_rows)
            )
        )
        # left_annot <- HeatmapAnnotation(
        #     Clone=clones$clone_label, clone_label=clone_label_generator,
        #     MainSite=mainsite_labels, Pdx=pdx_labels, Groupings=grouping_labels,
        #     col=annot_colours, show_annotation_name=c(TRUE, FALSE, TRUE),
        #     which="row", annotation_width=unit(rep(0.4, annot_cols), "cm"),
        #     annotation_legend_param=list(
        #         Clone=list(nrow=clone_legend_rows),
        #         MainSite=list(nrow=library_legend_rows),
        #         Pdx=list(nrow=library_legend_rows),
        #         Groupings=grouping_legend
        #     )
        # )
    } else {
        left_annot <- HeatmapAnnotation(
            MainSite=mainsite_labels, col=annot_colours,
            which="row", simple_anno_size=unit(0.4, "cm"),
            annotation_legend_param=list(
                MainSite=list(nrow=library_legend_rows)
            )
        )
    }
    
    return(left_annot)
}


get_chrom_label_pos <- function(copynumber) {
    chrom_label_pos <- c()
    chroms <- sapply(strsplit(colnames(copynumber), ":"), function(x) x[[1]])
    uniq_chroms <- c(as.character(1:22), "X", "Y")
    for(chrom in uniq_chroms) {
        chrom_idx <- which(chroms == chrom)
        chrom_label_pos[[chrom]] <- as.integer(round(mean(chrom_idx)))
    }
    return(chrom_label_pos)
}

# From ComplexHeatmap, needed for modified anno_mark
recycle_gp = function(gp, n = 1) {
    for(i in seq_along(gp)) {
        x = gp[[i]]
        gp[[i]] = c(rep(x, floor(n/length(x))), x[seq_len(n %% length(x))])
    }
    return(gp)
}

# From ComplexHeatmap, modified
anno_mark = function(at, labels, which = c("column", "row"),
    side = ifelse(which == "column", "top", "right"),
    lines_gp = gpar(), labels_gp = gpar(), padding = 0.5,
    link_width = unit(5, "mm"), link_height = link_width,
    link_gp = lines_gp,
    extend = unit(0, "mm")) {

    which = match.arg(which)[1]

    if(!is.numeric(at)) {
        stop_wrap(paste0("`at` should be numeric ", which, " index corresponding to the matrix."))
    }

    n = length(at)
    link_gp = recycle_gp(link_gp, n)
    labels_gp = recycle_gp(labels_gp, n)
    labels2index = structure(seq_along(at), names = labels)
    at2labels = structure(labels, names = at)

    if(length(extend) == 1) extend = rep(extend, 2)
    if(length(extend) > 2) extend = extend[1:2]
    if(!inherits(extend, "unit")) extend = unit(extend, "npc")

    height = link_width + max_text_width(labels, gp = labels_gp)
    width = unit(1, "npc")

    .pos = NULL
    .scale = NULL

    column_fun = function(index) {
        n = length(index)

        # adjust at and labels
        at = intersect(index, at)
        if(length(at) == 0) {
            return(NULL)
        }
        labels = at2labels[as.character(at)]

        # labels_gp = subset_gp(labels_gp, labels2index[labels])
        # link_gp = subset_gp(link_gp, labels2index[labels])

        if(is.null(.scale)) {
            .scale = c(0.5, n+0.5)
        }
        pushViewport(viewport(yscale = c(0, 1), xscale = .scale))
        if(inherits(extend, "unit")) extend = convertWidth(extend, "native", valueOnly = TRUE)
        # text_height = convertWidth(grobHeight(textGrob(labels, gp = labels_gp))*(1+padding), "native", valueOnly = TRUE)
        text_height = convertWidth(grobWidth(textGrob(labels, gp = labels_gp))*(2+padding), "native", valueOnly = TRUE)
        if(is.null(.pos)) {
            i2 = which(index %in% at)
            pos = i2 # position of rows
        } else {
            pos = .pos[which(index %in% at)]
        }
        h1 = pos - text_height*0.5
        h2 = pos + text_height*0.5
        pos_adjusted = smartAlign(h1, h2, c(.scale[1] - extend[1], .scale[2] + extend[2]))
        h = (pos_adjusted[, 1] + pos_adjusted[, 2])/2

        n2 = length(labels)
        # grid.text(labels, h, rep(max_text_width(labels, gp = labels_gp), n2), default.units = "native", gp = labels_gp, rot = 0, just = "center")
        grid.text(labels, h, rep(grobHeight(textGrob(labels, gp = labels_gp)), n2), default.units = "native", gp = labels_gp, rot = 0, just = "center")
        link_height = link_height - unit(1, "mm")
        grid.segments(pos, unit(rep(1, n2), "npc"), pos, unit(1, "npc")-rep(link_height*(1/3), n2), default.units = "native", gp = link_gp)
        grid.segments(pos, unit(1, "npc")-rep(link_height*(1/3), n2), h, unit(1, "npc")-rep(link_height*(2/3), n2), default.units = "native", gp = link_gp)
        grid.segments(h, unit(1, "npc")-rep(link_height*(2/3), n2), h, unit(1, "npc")-rep(link_height, n2), default.units = "native", gp = link_gp)
        upViewport()
    }

    fun = column_fun

    anno = AnnotationFunction(
        fun = fun,
        fun_name = "anno_mark",
        which = which,
        width = width,
        height = height,
        n = -1,
        var_import = list(at, labels2index, at2labels, link_gp, labels_gp, padding, .pos, .scale,
            side, link_width, link_height, extend),
        show_name = FALSE
    )

    # anno@subset_rule$at = subset_by_intersect

    anno@subsetable = TRUE
    return(anno)
}

make_bottom_annot <- function(copynumber) {
    chrom_label_pos <- get_chrom_label_pos(copynumber)
    bottom_annot <- HeatmapAnnotation(chrom_labels=anno_mark(
        at=chrom_label_pos,
        labels=names(chrom_label_pos),
        side="bottom",
        padding=0.5, extend=0.01
    ), show_annotation_name=FALSE)
    return(bottom_annot)
}

make_copynumber_heatmap <- function(copynumber, clones, grouping_file) {
    copynumber_hm <- Heatmap(
        name="Copy Number",
        as.matrix(copynumber),
        col=cn_colours,
        na_col="white",
        show_row_names=FALSE,
        cluster_rows=FALSE,
        cluster_columns=FALSE,
        show_column_names=FALSE,
        bottom_annotation=make_bottom_annot(copynumber),
        left_annotation=make_left_annot_v2(copynumber, clones, grouping_file),
        heatmap_legend_param=list(nrow=4),
        use_raster=TRUE,
        raster_quality=5
    )
    return(copynumber_hm)
}


make_copynumber_heatmap_v2 <- function(copynumber, clones, grouping_file) {
    copynumber_hm <- Heatmap(
        name="Copy Number",
        as.matrix(copynumber),
        col=cn_colours,
        na_col="white",
        show_row_names=FALSE,
        cluster_rows=FALSE,
        cluster_columns=FALSE,
        show_column_names=FALSE,
        bottom_annotation=make_bottom_annot(copynumber),
        left_annotation=make_left_annot_v3(copynumber, clones, grouping_file),
        heatmap_legend_param=list(nrow=4),
        use_raster=TRUE,
        raster_quality=5
    )
    return(copynumber_hm)
}
make_cell_copynumber_tree_heatmap <- function(tree, copynumber, clones,
                                              brlen, grouping_file) {
    tree <- format_tree(tree, brlen)

    tree_ggplot <- make_tree_ggplot(tree, clones)
    copynumber <- format_copynumber(copynumber, tree_ggplot$data)
    if(!is.null(clones)) {
        clones <- format_clones(clones, tree_ggplot$data)
    }

    tree_hm <- make_corrupt_tree_heatmap(tree_ggplot)
    # copynumber_hm <- make_copynumber_heatmap(copynumber, clones, grouping_file)
    copynumber_hm <- make_copynumber_heatmap_v2(copynumber, clones, grouping_file)
    h <- tree_hm + copynumber_hm

    draw(
         h,
         padding=unit(c(10, 2, 2, 2), "mm"),
         annotation_legend_side="right",
         heatmap_legend_side="bottom"
    )
}

# ords <- c('pdx','mainsite','organ')
# result_dir <- "/home/htran/storage/datasets/metastasis_results/SA535_new_encoding/Tyler_2361/"
# copynumber <- read.csv(paste0(result_dir,'total_merged_filtered_states.csv'), header=T, row.names=1, check.names = F,stringsAsFactors = FALSE)
# clones <- read_tsv(paste0(result_dir,'cell_clones.csv'))
# grouping_file <- paste0(result_dir,'library_groupings.csv')
make_cell_copynumber_heatmap_ordering <- function(copynumber, clones, 
                                                  grouping_file, output_dir, 
                                                  order_type='pdxid', datatag='',
                                                  ht=550,wd=1400) {
    
    save_dir <- paste0(output_dir,"/heatmap/")
    if (!file.exists(save_dir)){
        dir.create(save_dir)
    }
    if(!is.null(clones)) {
        clones <- format_clones_v2(clones, copynumber, grouping_file, order_type)
    }
    print("DEBUG 3")
    copynumber <- format_copynumber_v2(copynumber, clones)
    print("DEBUG 4")
    copynumber_hm <- make_copynumber_heatmap_v2(copynumber, clones, grouping_file)
    print("DEBUG 5")
    # padding of the plot. Elements correspond to bottom, left, top, right paddings.
    png(paste0(save_dir,datatag,'_',order_type,'_hm.png'), height = 2*(ht+250), width=2*wd, res = 2*72)
    # png(paste0(result_dir,'tree_viz_mainsite.png'), height = 2*800, width=2*1400, res = 2*72)
    draw(
        copynumber_hm,
        padding=unit(c(6, 4, 2, 2), "mm"),
        annotation_legend_side="right",
        heatmap_legend_side="bottom"
    )
    dev.off()
}

generate_heatmap <- function(){
    argv <- get_args()
    print("DEBUG 1")
    print(argv$copynumber)
    copynumber <- read.csv(argv$copynumber, header=T, row.names=1, check.names = F,stringsAsFactors = FALSE)
    # copynumber <- read_tsv(argv$copynumber)
    print(dim(copynumber))
    
    if(!is.null(argv$clones)) {
        clones <- read_tsv(argv$clones)
    } else {
        clones <- NULL
    }
    copynumber <- copynumber[,colnames(copynumber) %in% clones$cell_id]
    ords <- c("pdxid","mainsite","origin","clone_id")
    for(ot in ords){
        print("DEBUG 2")
        print(argv$grouping_file)
        print(argv$outdir)
        print(dim(clones))
        make_cell_copynumber_heatmap_ordering(copynumber, clones,
                                              argv$grouping_file, argv$outdir, 
                                              order_type=ot, argv$datatag)
    }
    
}

generate_heatmap_v2 <- function(){
    argv <- get_args()
    print("DEBUG 1")
    print(argv$copynumber)
    copynumber <- read.csv(argv$copynumber, header=T, row.names=1, check.names = F,stringsAsFactors = FALSE)
    # copynumber <- read_tsv(argv$copynumber)
    print(dim(copynumber))
    
    if(!is.null(argv$clones)) {
        clones <- read_tsv(argv$clones)
    } else {
        clones <- NULL
    }
    
    copynumber <- copynumber[,colnames(copynumber) %in% clones$cell_id]
    cell_ids <- colnames(copynumber)
    meta_grouping <- get_library_grouping(cell_ids, argv$grouping_file)
    # for(mt in unique(meta_grouping$mainsite)){
    #     print(mt)
    #     cells_use <- meta_grouping[meta_grouping$mainsite==mt,'cell_id']
    #     copynumber_tmp <- copynumber[,cells_use]
    #     print(dim(copynumber_tmp))
    #     make_cell_copynumber_heatmap_ordering(copynumber_tmp, clones,
    #                                           argv$grouping_file, argv$outdir, 
    #                                           "clone_id", paste0(mt,'_',argv$datatag),
    #                                           ht=round(ncol(copynumber_tmp)/6),wd=1400)
    # }
    
    for(pd in unique(meta_grouping$pdxid)){
        print(pd)
        cells_use <- meta_grouping[meta_grouping$pdxid==pd,'cell_id']
        copynumber_tmp <- copynumber[,cells_use]
        print(dim(copynumber_tmp))
        make_cell_copynumber_heatmap_ordering(copynumber_tmp, clones,
                                              argv$grouping_file, argv$outdir, 
                                              "clone_id", paste0(pd,'_',argv$datatag),
                                              ht=round(ncol(copynumber_tmp)/6),wd=1400)
    }
    for(ori in unique(meta_grouping$origin)){
        print(ori)
        cells_use <- meta_grouping[meta_grouping$origin==ori,'cell_id']
        copynumber_tmp <- copynumber[,cells_use]
        print(dim(copynumber_tmp))
        make_cell_copynumber_heatmap_ordering(copynumber_tmp, clones,
                                              argv$grouping_file, argv$outdir, 
                                              "clone_id", paste0(ori,'_',argv$datatag),
                                              ht=round(ncol(copynumber_tmp)/6),wd=1400)
    }
   
    
}

main <- function() {
    argv <- get_args()
    
    tree <- read.tree(argv$tree)
    print(argv$tree)
    print(argv$copynumber)
    copynumber <- read.csv(argv$copynumber, header=T, row.names=1, check.names = F,stringsAsFactors = FALSE)
    # copynumber <- read_tsv(argv$copynumber)
    
    if(argv$normalize_ploidy) {
        print(paste(
            "WARNING: ceiling() will be applied to non-integer copynumber",
            "during normalization"
        ))
        copynumber <- normalize_cell_ploidy(copynumber)
    }

    if(!is.null(argv$clones)) {
        clones <- read_tsv(argv$clones)
    } else {
        clones <- NULL
    }

    # pdf(argv$pdf, width=10)
    png(argv$pdf, height = 2*800, width=2*1400, res = 2*72)
    make_cell_copynumber_tree_heatmap(
        tree, copynumber, clones, argv$branch_lengths, argv$grouping_file
    )
    dev.off()
}

generate_heatmap_v2()
# generate_heatmap()
# main()  # comment back