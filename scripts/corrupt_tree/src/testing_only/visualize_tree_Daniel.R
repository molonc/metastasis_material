# read in the files
reads <- fread("*_reads.csv.gz")
metrics <- read.csv("*_metrics.csv.gz")
# for proper chr sorting
reads$chr <- factor(reads$chr, levels = c(1:22, "X", "Y"))
# Building the hierarchical tree, skip if you don't care about the sorting
slice <- reads[, c("chr", "start", "cell_id", "state")]
slice$pos <- paste0(slice$chr, ":", slice$start)
wide <- spread(slice[, -c(1:2)], pos, state)
rownames(wide) <- wide$cell_id
wide$cell_id <- NULL
cluster <- hclust(dist(wide), method = "ward.D")
slice$cell_id <- factor(slice$cell_id, level = rownames(wide)[cluster$ord])
cols <- c(rev(brewer.pal(n = 3, "Blues"))[1:2], "#CCCCCC", tail(brewer.pal(n = 8, "OrRd"), 6))
cols <- c(cols, cols[c(9, 9, 9)])
names(cols) <- 1:12
png(file = "bulk_hmmcopy_heatmap.png", 2400, 1600)
g <- ggplot(slice, aes(start, cell_id, fill = as.factor(state))) + geom_tile() + scale_fill_manual(values = cols) + facet_grid(~chr, scales = "free", space = "free", switch = "x") + scale_x_continuous(expand = c(0, 0), breaks = NULL) + theme(panel.spacing = unit(0.1, "lines"))
print(g)
dev.off()


reads <- fread("*_reads.csv.gz")
lim <- max(quantile(tmp$copy, na.rm = TRUE, 0.99), 4)
tmp <- subset(reads, cell_id == "SOME ID OF YOUR CHOICE")
ggplot(tmp, aes(start, copy * 2, col = as.factor(state))) + geom_point(size = 0.5) + 
  facet_grid(cell_id ~ chr, s, space = "free_x", switch = "x") + 
  scale_x_continuous(expand = c(0, 0), breaks = NULL) + theme(limit = c(0, lim)) + 
  scale_col_manual(values = cols)



