#!/usr/bin/env Rscript
library(argparse)
library(dplyr)
library(readr)
library(tidyr)


get_args <- function() {
    p <- ArgumentParser()

    p$add_argument("--states", "-s", help = "CN states file")
    p$add_argument("--binary_features", "-b",help = "jitter fixed binary features matrix csv file")

    p$add_argument("--out", "-o", help = "output binary features csv file")

    return(p$parse_args())
}

create_state_matrix <- function(states) {
    bin_names <- paste("locus", states$chr, states$start, states$end, sep="_")

    states <- select(states, -chr, -start, -end, -width)
    # bin_names <- states$loci
    states <- as.matrix(states)

    rownames(states) <- bin_names
    colnames(states) <- paste("cell", colnames(states), sep = "_")

    return(states)
}

create_feature_matrix <- function(features) {
    features <- tidyr::spread(features, cells, tipInclusionProbabilities)

    loci <- features$loci
    features <- select(features, -loci)
    features <- as.matrix(features)
    rownames(features) <- loci

    return(features)
}

# For Tyler version, there are some left padding bin with negative values
parse_bin_names <- function(bin_names) {
    # Remove corrupt_tree locus tag if it's there
    bin_names <- gsub('locus_', '', bin_names)
    chr <- gsub('([0-9]+|X|Y)_(-1|[0-9]+)_[0-9]+', '\\1', bin_names)
    start <- as.numeric(gsub('([0-9]+|X|Y)_(-1|[0-9]+)_[0-9]+', '\\2', bin_names))
    end <- as.numeric(gsub('([0-9]+|X|Y)_(-1|[0-9]+)_([0-9]+)', '\\3', bin_names))
    bins <- data.frame(
        chr = chr, start = start, end = end, stringsAsFactors = FALSE
    )
    return(bins)
}

# parse_bin_names <- function(bin_names) {
#     # Remove corrupt_tree locus tag if it's there
#     bin_names <- gsub('locus_', '', bin_names)
#     chr <- gsub('([0-9]+|X|Y)_[0-9]+_[0-9]+', '\\1', bin_names)
#     start <- as.numeric(gsub('([0-9]+|X|Y)_([0-9]+)_[0-9]+', '\\2', bin_names))
#     end <- as.numeric(gsub('([0-9]+|X|Y)_([0-9]+)_([0-9]+)', '\\3', bin_names))
#     bins <- data.frame(
#         chr = chr, start = start, end = end, stringsAsFactors = FALSE
#     )
#     return(bins)
# }

sort_mat_by_bins <- function(the_mat) {
    bin_info <- parse_bin_names(rownames(the_mat))

    # Sort the matrix by their chromosome (just inside the chromosome)
    bin_info$chr[bin_info$chr == 'X'] <- '40' 
    bin_info$chr[bin_info$chr == 'Y'] <- '55' 
    bin_info$chr <- as.numeric(bin_info$chr)

    the_mat <- the_mat[order(bin_info$chr, bin_info$start), ]
    
    return(the_mat)
}

value_mode <- function(x) {
    uniq_x <- unique(x)
    val <- uniq_x[which.max(tabulate(match(x, uniq_x)))]
    return(val)
}

compute_state_ploidy_ratios <- function(states) {
    ploidies <- apply(states, 2, value_mode)
    ratios <- round(sweep(states, MARGIN=2, ploidies, "/"), 5)

    return(ratios)
}

compute_state_ploidy_ratio_freq <- function(ratios, locus_idx, cell_idx) {

    ratio_freq <- as.data.frame(
        table(as.vector(as.matrix(ratios[locus_idx, cell_idx]))),
        stringsAsFactors = FALSE
    )

    ratio_freq$freq <- ratio_freq$Freq/sum(ratio_freq$Freq)
    ratio_freq$ratio <- as.numeric(ratio_freq$Var1)
    ratio_freq <- ratio_freq[, c("ratio", "freq")]

    ratio_freq <- ratio_freq[order(ratio_freq$freq, decreasing = TRUE), ]

    return(ratio_freq)
}

create_sync_cnv_feature <- function(state_ploidy_ratios, ratio, locus_idx,
                                    bin_info, i) {

    # if (locus_idx == 3609) {
        # print(state_ploidy_ratios[locus_idx, ])
        # print(which(state_ploidy_ratios[locus_idx, ] == ratio))
        # print(paste(ratio))
    # }

    new_feature <- state_ploidy_ratios[1, , drop = FALSE]
    new_feature[, ] <- 0

    new_start <- bin_info$start[locus_idx] + i * 100
    rownames(new_feature) <- paste(
        "locus", bin_info$chr[locus_idx], new_start, new_start + 1, sep = "_"
    )

    new_feature[, which(state_ploidy_ratios[locus_idx, ] == ratio)] <- 1

    return(new_feature)
}

format_new_features <- function(new_features) {
    # sort new_features by rows again to put the new loci in their place...
    new_features <- sort_mat_by_bins(new_features)

    new_features <- as.data.frame(new_features)
    new_features$loci <- rownames(new_features)
    rownames(new_features) <- NULL

    new_features <- tidyr::gather(
        new_features, key = 'cells', value = 'tipInclusionProbabilities', -loci
    )
    feature_cols <- c("cells", "loci", "tipInclusionProbabilities")
    new_features <- new_features[, feature_cols]

    return(new_features)
}

insert_sync_cnv_features <- function(states, features,
                                  feature_freq_threshold = 0.6) {

    states <- sort_mat_by_bins(states)
    features <- sort_mat_by_bins(features)
    states <- states[, colnames(features)]

    bin_info <- parse_bin_names(rownames(states))
    state_ploidy_ratios <- compute_state_ploidy_ratios(states)

    loci_feature_freq <- rowMeans(features)
    # print(loci_feature_freq["locus_20_3000001_3500000"])
    loci_of_interest  <- which(unname(
        loci_feature_freq >= feature_freq_threshold
    ))

    my_loci <- which("locus_20_3000001_3500000" == names(loci_feature_freq))
    # print(my_loci %in% loci_of_interest)
    # print(my_loci)

    new_features <- features

    for (l in loci_of_interest) {
        state_locus_idx <- which(
            rownames(state_ploidy_ratios) == rownames(features)[l]
        )
        cell_idx <- which(features[l, ] == 1)

        state_ploidy_ratio_freq <- compute_state_ploidy_ratio_freq(
            state_ploidy_ratios, state_locus_idx, cell_idx
        )

        # if(l == my_loci) {
            # print(state_locus_idx)
            # print(state_ploidy_ratio_freq)
        # }

        if (nrow(state_ploidy_ratio_freq) == 1) {next}

        for (i in 2:nrow(state_ploidy_ratio_freq)) {

            if (state_ploidy_ratio_freq[i, "freq"] < 0.05) {break}
        
            new_feature <- create_sync_cnv_feature(
                state_ploidy_ratios, state_ploidy_ratio_freq[i, "ratio"],
                state_locus_idx, bin_info, i
            )

            stopifnot(sum(new_feature) > 0)
            # if (l == my_loci) {
                # print(paste("new feature", i, sum(new_feature)))
            # }

            new_features <- rbind(new_features, new_feature)
        }
    }

    new_features <- format_new_features(new_features)

    return(new_features)
}


argv <- get_args()

states <- read.csv(argv$states, check.names = F, stringsAsFactors = FALSE)
features <- read.csv(argv$binary_features, check.names = F, stringsAsFactors = FALSE)

# states <- read_tsv(argv$states, col_types = cols(
#     .default = col_integer(),
#     chr = col_character()
# ))
# features <- read_csv(argv$binary_features, col_types = cols(
#     .default = col_character(),
#     tipInclusionProbabilities = col_double()
# ))

# st <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local_tyler/total_merged_filtered_states.csv'
# st <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local_tyler/bin_cnvs_corrupt_double_padding.csv'
# 
# ft <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local_tyler/corrupt_tree_straightened_features.csv'
# states <- read.csv(st, check.names = F, stringsAsFactors = FALSE)
# features <- read.csv(ft, check.names = F, stringsAsFactors = FALSE)
states <- create_state_matrix(states)
features <- create_feature_matrix(features)

new_features <- insert_sync_cnv_features(states, features)
# output_fn <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local_tyler/corrupt_tree_sync-cnv_features.csv'
# write_csv(new_features, output_fn)
write_csv(new_features, argv$out)

if (!is.null(warnings())) {
    print(warnings())
}
