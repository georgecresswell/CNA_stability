options(rlib_downstream_check = FALSE)
library(ggplot2)
library(reshape2)
library(ggtree)
library(circlize)
library(ComplexHeatmap)
library(magrittr)
library(cowplot)
theme_set(theme_cowplot())
library(GenomicRanges)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)
genome = BSgenome.Hsapiens.UCSC.hg38

# Only take cells with the same ploidy
ploidy_increase_lim = 1

# Minimum number of bins to be used in the analysis
min_bin_num = 2

# Expected ploidy
PLOIDY = 2

# Read in the reads
scCNA = read.csv("POTs_EP1_reads.csv", stringsAsFactors = F)

# Make as a factor
scCNA$chr = factor(scCNA$chr, levels = c(1:22,"X","Y"))

# Read in the metrics
scCNA_metrics = read.csv("POTs_EP1_metrics.csv", stringsAsFactors = F)

# Strict filtering of good quality non-S-phase cells
good_cells = scCNA_metrics$cell_id[scCNA_metrics$quality > 0.95 & scCNA_metrics$total_reads > 250000 & 
                                     scCNA_metrics$is_s_phase_prob < 0.1]

# Get those cells
scCNA_good_cells = scCNA[scCNA$cell_id %in% good_cells,]

# Sort matrix by the cells then location
scCNA_good_cells = scCNA_good_cells[order(scCNA_good_cells$cell_id, scCNA_good_cells$chr, 
                                          scCNA_good_cells$start),]

# Only take bins which are ideal
scCNA_good_cells = scCNA_good_cells[scCNA_good_cells$ideal == "True",]

# Make a bin id
scCNA_good_cells$bin_id = factor(paste0(scCNA_good_cells$chr,"_",
                                        scCNA_good_cells$start,"_",scCNA_good_cells$end),
                                 levels = unique(paste0(scCNA_good_cells$chr,"_"
                                                        ,scCNA_good_cells$start,"_",scCNA_good_cells$end)))

# Make that a matrix
scCNA_mat = reshape2::dcast(scCNA_good_cells, cell_id ~ bin_id, value.var = "state")

# Take only those with no NAs (i.e. ideal in some but not others)
scCNA_mat = scCNA_mat[,apply(scCNA_mat, 2, function(i) length(which(is.na(i)))==0)]

# Make cell id the rownames
rownames(scCNA_mat) = scCNA_mat$cell_id
scCNA_mat$cell_id = NULL

# Transpose
scCNA_mat = t(scCNA_mat)

# Do the ploidy filter
scCNA_mat = scCNA_mat[,colMeans(scCNA_mat) < (mean(colMeans(scCNA_mat)) + ploidy_increase_lim)]

# Split out row names by chr, start, end
locations_split = strsplit(rownames(scCNA_mat), "_")

# Add location columns
chr   = unlist(lapply(locations_split, function(i) i[1]))
start = as.numeric(unlist(lapply(locations_split, function(i) i[2])))
end   = as.numeric(unlist(lapply(locations_split, function(i) i[3])))

# Make a dataframe
chr_pos = data.frame(chr = chr,
                     start = start,
                     end = end)

# Make an identifier
ids = sapply(1:nrow(scCNA_mat), function(i) paste0(c(chr[i], scCNA_mat[i,]), collapse = "_"))

# Collapse
run_len = rle(ids)

# Cum summed
cumulative_len = cumsum(run_len$lengths)

# Create a matrix
res = data.frame(chrom = chr[cumulative_len],
                 start = start[c(1, cumulative_len[-length(cumulative_len)]+1)],
                 end = end[cumulative_len],
                 size = run_len$lengths,
                 scCNA_mat[cumulative_len,])

# Filter for size
res = res[res$size >= min_bin_num,]
res$size = NULL

# Rearrange
res = cbind(res[,1:3], reshape2::melt(res[,4:ncol(res)]))
res = res[,c(4,1:3,5)]
colnames(res) = c("sample_id", "chrom", "start", "end", "cn")

# Remove redundant part of the cell name
res$sample_id = gsub("[.]POTs", "", unlist(lapply(strsplit(as.character(res$sample_id), split = "_"), 
                                                  function(i) paste0(i[6], collapse = "_"))))

options(scipen = 999)

# Write it out
write.table(res, file = paste0("SW620_DLP_regions.tsv"),
            sep = "\t", row.names = F, quote = F)

######################################################################################
##### A few moments later... #########################################################
######################################################################################

# Read in the MEDICC2 tree
tree = try(read.tree(paste0("SW620_DLP_regions_final_tree.new")))

# Collect the tree info in a df
breakdown = data.frame(A = tree$edge[,1],
                       B = tree$edge[,2],
                       Edge_Length = tree$edge.length,
                       A_Label = c(tree$tip.label, tree$node.label)[tree$edge[,1]],
                       B_Label = c(tree$tip.label, tree$node.label)[tree$edge[,2]])

# Branch lengths
med_bl = median(breakdown[breakdown$A_Label!="","Edge_Length"])
men_bl = signif(mean(breakdown[breakdown$A_Label!="","Edge_Length"]), digits = 2)

# Basic tree structure
p = ggtree(tree, size = 1) + geom_tiplab(hjust = -0.5, size = 4) +
           geom_tippoint(col = "black", size = 3) +
           ggtitle(paste0("SW620 (branch length=",med_bl,"/",men_bl,")")) + 
           theme(title=element_text(size=20, colour = "gray20"),
                 plot.title = element_text(hjust = 0.5))

# Label only branches of a certain size
p = p + geom_label(data = p$data[p$data$branch.length > 10,], aes(x = branch, label = branch.length),
           label.size = 0, nudge_y = 0.2)

# Alt version
p = ggtree(tree, size = 1) + geom_tiplab(hjust = -0.5, size = 0) +
  geom_tippoint(col = "black", fill = c(rep("#75bbfd", times = length(tree$tip.label)-1), "white"), 
                size = 3, shape = 21) +
  ggtitle(paste0("SW620")) + 
  theme_tree2() +
  xlab("") +
  theme(title=element_text(size=20, colour = "gray20"),
        axis.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5)) 

ggsave(plot = p, filename = "SW620_single_cell_tree.png", 
       width = 5, height = 5)

######################################################################################
##### Now we will plot a heatmap #####################################################
######################################################################################

# Remove redundant part of the cell names
colnames(scCNA_mat) = gsub("-POTs", "", unlist(lapply(strsplit(as.character(colnames(scCNA_mat)), split = "_"), 
                                                      function(i) paste0(i[6], collapse = "_"))))

# Make them factors
chr.ends = cumsum(table(factor(chr, levels = c(1:22,"X","Y"))))

# Bottom annotation
bot.column.anno = columnAnnotation(link = anno_mark(at = round(c(0, chr.ends[-24]) +
                                   rle(chr)$lengths / 2),
                                   labels = c(1:22, "X", "Y"),
                                   side = "bottom",
                                   link_width = unit(0.25, "cm")))

pdf("SW620_single_cell_heatmap.pdf", 
    width = 7, height = 4)

# Make a heatmap ordered by phylogeny
Heatmap(t(scCNA_mat[,get_taxa_name(p)[-length(get_taxa_name(p))]]), 
        cluster_rows = F, cluster_columns = F, 
        show_column_names = F,
        show_row_names = F,
        border = T,
        use_raster = F,
        column_title = "SW620",
        name = "CN",
        bottom_annotation = bot.column.anno,
        col = colorRamp2(c(-1, 0, PLOIDY, PLOIDY*2, PLOIDY*4), 
                         c("grey", "#0504aa", "white", "#ff474c", "#e50000")),
        heatmap_legend_param = list(at = c(0, 1, PLOIDY, PLOIDY*2, PLOIDY*4), 
                                    labels = paste0(c(0, 1, PLOIDY, PLOIDY*2, paste0(PLOIDY*4,"+")))))

#Add lines
for(boundary in chr.ends[-24] / nrow(scCNA_mat)) {
  
  #Add the lines
  decorate_heatmap_body("CN", {
    grid.lines(c(boundary, boundary), c(0, 1), gp = gpar(lty = "dotted", lwd = 1))
  })
  
}

dev.off()

# Special function for writing out the supplementary data
outputHeatmapData = function(mat, filename) {
  
  res = cbind(do.call(rbind, strsplit(rownames(mat), split = "_")), mat)
  
  colnames(res)[1:3] = c("chr", "start", "end")
  
  rownames(res) = NULL

  write.csv(res, file = filename, quote = F, row.names = F)
  
}

# Perform it
outputHeatmapData(scCNA_mat, 
                  filename = "SW620_single_cell_heatmap_data.csv")

######################################################################################
##### Now we will calculate timed tree events ########################################
######################################################################################

# Collect the tree info in a df
breakdown = data.frame(A = tree$edge[,1],
                       B = tree$edge[,2],
                       Edge_Length = tree$edge.length,
                       A_Label = c(tree$tip.label, tree$node.label)[tree$edge[,1]],
                       B_Label = c(tree$tip.label, tree$node.label)[tree$edge[,2]])

# CN profiles as calculated by MEDICC2
cn_profiles = read.table(paste0("SW620_DLP_regions_final_cn_profiles.tsv"),
                         stringsAsFactors = F, header = T)

# Get those edges that include the tip
tip_edges = breakdown[!(grepl("internal_", breakdown$B_Label) | grepl("diploid", breakdown$B_Label)),]

# Calculate the cna difference of each
tip_cnas = apply(tip_edges, 1, function(i) cn_profiles[cn_profiles$sample_id==i["B_Label"],"cn"] - 
                   cn_profiles[cn_profiles$sample_id==i["A_Label"],"cn"])

# Count gains and losses, each gain and loss is only counted as one
tip_gain = apply(tip_cnas > 0, 1, sum)
tip_loss = apply(tip_cnas < 0, 1, sum)

# Make a dataframe
tip_events = data.frame(cn_profiles[1:length(tip_loss),2:4], loss = tip_loss, gain = tip_gain)

# Now make the windows across the genome
windows = data.frame(do.call(rbind, strsplit(rownames(scCNA_mat), split = "_")))

# Various format measures
colnames(windows) = c("chr", "start", "end")
windows$chr = factor(windows$chr, levels = c(1:22,"X","Y"))
windows$name  = paste0("win",1:nrow(windows))
windows$start = as.numeric(as.character(windows$start))
windows$end   = as.numeric(as.character(windows$end))
windows_gr    = GRanges(seqnames = paste0("chr", windows$chr), 
                        ranges = IRanges(start = windows$start, end = windows$end))

# Make the tip events into a GRanges object
tip_events_gr = GRanges(seqnames = tip_events$chrom, ranges = IRanges(start = tip_events$start, 
                                                                      end = tip_events$end))

# Start with NAs
res = data.frame(matrix(NA, nrow = length(windows_gr), ncol = 2))

# Where do the windows and the regions overlap?
overlaps = findOverlaps(windows_gr, tip_events_gr)

# Add them in 
for(i in 1:length(overlaps)) {
  r = queryHits(overlaps)[i]
  res[r,] = tip_events[subjectHits(overlaps)[i],4:ncol(tip_events)]
}

# Take the column names
colnames(res) = colnames(tip_events)[4:ncol(tip_events)]

# Combine with the windows
data_per_window_summary = cbind(windows, res)

# Recurrence plot in the style of Timon Heide
chr_centers =
  filter(D3GB::GRCh38.bands) %>%
  split(., .$chr) %>%
  sapply(function(x) max(x$end)/2)

centromer_position =
  filter(D3GB::GRCh38.bands, score=="acen") %>%
  split(., .$chr) %>%
  sapply(function(x) mean(range(c(x$end,x$start))))
  
chr_centers        = chr_centers[1:22]
centromer_position = centromer_position[1:22]

ylm = 8

plot_df=
  data_per_window_summary %>%
  mutate(
    chr = gsub("chr", "", chr), 
    chr = factor(chr, levels = as.character(1:22))
  ) %>%
  filter(!is.na(chr)) %>%
  mutate(pos=(start+end)/2-chr_centers[as.character(chr)])

plot_cna_recurence = ggplot(data = plot_df, aes(x=pos, y=gain)) +
  geom_area(aes(alpha="tip", fill="gains")) +
  geom_area(aes(y=-loss, alpha="tip", fill="losses")) +
  geom_hline(yintercept=0) +
  geom_vline(data=data.frame(chr=names(centromer_position), pos=centromer_position-chr_centers),
  aes(xintercept=pos), linetype=3, size=0.5, alpha=0.2) +
  facet_grid(~ factor(chr, levels = as.character(1:22)), scales="free_x", space="free_x") +
  theme(strip.text.x=element_text(vjust=-102), strip.background.x=element_blank()) +
  scale_alpha_manual(values=c("clonal"=1, "intermediate"=0.5, "tip"=0.8)) +
  scale_fill_manual(values=c("losses"="#3c73a8","gains"="#ec2d01")) +
  scale_x_continuous(breaks=c(0), labels="") +
  scale_y_continuous(breaks=seq(-ylm, ylm, by=2),
  labels=abs(seq(-ylm, ylm, by=2)),
  limits=c(-ylm,ylm)) +
  xlab("Chromosomes") + ylab("Number of cells") +
  guides(fill="none", alpha="none")

# switch of clipping for the strip texts:
plot_cna_recurence = ggplotGrob(plot_cna_recurence)
for (i in which(grepl("strip-t", plot_cna_recurence$layout$name))){
  plot_cna_recurence$grobs[[i]]$layout$clip = "off"
}
ggsave(plot = plot_cna_recurence, 
       file = "SW620_tip_frequencies_single_cell.png",
       width = 12.7, height = 5)
