#!/applications/R/R-3.3.2/bin/Rscript

# Generate read counts for each peak and ranLoc (with equivalent width
# distribution to peaks) and plot read counts against locus width

# /applications/R/R-3.3.2/bin/Rscript ./REC8_peaks_ranLoc_read_counts_RPKM.R REC8_HA_Rep2_ChIP REC8_MYC_Rep1_input PE

#ChIPname <- "REC8_HA_Rep2_ChIP"
#inputname <- "REC8_MYC_Rep1_input"
#libType <- "PE"

args <- commandArgs(trailingOnly = TRUE)
ChIPname <- args[1]
inputname <- args[2]
libType <- args[3]

library(GenomicAlignments)
library(ShortRead)
library(rtracklayer)
library(regioneR)

chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chrStart <- c(1, 1, 1, 1, 1)
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))

bamDir <- "/home/ajt200/analysis/REC8_pooled/"
peakDir <- "/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep1_input_p0.001_q0.01/"

## ChIP
## Load BAM and create RangedData object
#if(libType == "PE") {
#  ChIP_readGAlignments <- readGAlignmentPairs(paste0(bamDir, ChIPname,
#                                                     "_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam"))
#} else if(libType == "SE") {
#  ChIP_readGAlignments <- readGAlignments(paste0(bamDir, ChIPname,
#                                                 "_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam"))
#} else {
#  stop("Library type neither SE nor PE")
#}
#ChIP_ranges <- ranges(ChIP_readGAlignments)
#ChIP_chrs <- as.character(as.data.frame(seqnames(ChIP_readGAlignments))[,1])
#ChIP_strand <- as.character(strand(ChIP_readGAlignments))
#ChIP_ranged <- RangedData(space = ChIP_chrs,
#                          ranges = ChIP_ranges,
#                          strand = ChIP_strand) 
#ChIP_ranged <- ChIP_ranged[space(ChIP_ranged) != "chloroplast" &
#                           space(ChIP_ranged) != "mitochondria",]
#names(ChIP_ranged) <- sub("", "Chr", names(ChIP_ranged))
#save(ChIP_ranged, file = paste0(ChIPname, "_RangedData.RData"))
load(paste0(ChIPname, "_RangedData.RData"))

# Calculate library size
ChIP_size <- length(space(ChIP_ranged))
# Calculate "per million" scaling factor
ChIP_RPM_scaling_factor <- ChIP_size/1e+06
# Convert ranged data into GRanges object
ChIP_rangedGR <- as(ChIP_ranged, "GRanges")


## input
## Load BAM and create RangedData object
#if(libType == "PE") {
#  input_readGAlignments <- readGAlignmentPairs(paste0(bamDir, inputname,
#                                                      "_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam"))
#} else if(libType == "SE") {
#  input_readGAlignments <- readGAlignments(paste0(bamDir, inputname,
#                                                  "_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam"))
#} else {
#  stop("Library type neither SE nor PE")
#}
#input_ranges <- ranges(input_readGAlignments)
#input_chrs <- as.character(as.data.frame(seqnames(input_readGAlignments))[,1])
#input_strand <- as.character(strand(input_readGAlignments))
#input_ranged <- RangedData(space = input_chrs,
#                           ranges = input_ranges,
#                           strand = input_strand) 
#input_ranged <- input_ranged[space(input_ranged) != "chloroplast" &
#                             space(input_ranged) != "mitochondria",]
#names(input_ranged) <- sub("", "Chr", names(input_ranged))
#save(input_ranged, file = paste0(inputname, "_RangedData.RData"))
load(paste0(inputname, "_RangedData.RData"))

# Calculate library size
input_size <- length(space(input_ranged))
# Calculate "per million" scaling factor
input_RPM_scaling_factor <- input_size/1e+06
# Convert ranged data into GRanges object
input_rangedGR <- as(input_ranged, "GRanges")

# Load features
load(paste0(peakDir, ChIPname, "_rangerPeaksGR_arm_mergedOverlaps_noMinWidth.RData"))
load(paste0(peakDir, ChIPname, "_rangerPeaksGR_peri_mergedOverlaps_noMinWidth.RData"))
rangerPeaksGR_arm_mergedOverlaps <- sortSeqlevels(rangerPeaksGR_arm_mergedOverlaps)
rangerPeaksGR_peri_mergedOverlaps <- sortSeqlevels(rangerPeaksGR_peri_mergedOverlaps)
peaksGR <- sort(c(rangerPeaksGR_arm_mergedOverlaps, rangerPeaksGR_peri_mergedOverlaps))

## Calculate log2(ChIP_RPM/input_RPM) at features and random loci
# Create random loci of the same number and size distribution as for features
set.seed(8349235)
ranLocGR <- randomizeRegions(peaksGR,
                             genome = genome,
                             per.chromosome = TRUE,
                             allow.overlaps = TRUE)

## ChIP
# Calculate RPM and RPKM for each feature and random locus
peaks_ChIP_RPM <- countOverlaps(peaksGR, ChIP_rangedGR)/
                  ChIP_RPM_scaling_factor
peaks_ChIP_RPKM <- peaks_ChIP_RPM/
                   (width(peaksGR)/1e+03)
save(peaks_ChIP_RPM,
     file = paste0(ChIPname, "_RPM_at_", ChIPname, "_peaks.RData"))
save(peaks_ChIP_RPKM,
     file = paste0(ChIPname, "_RPKM_at_", ChIPname, "_peaks.RData"))

ranLoc_ChIP_RPM <- countOverlaps(ranLocGR, ChIP_rangedGR)/
                   ChIP_RPM_scaling_factor
ranLoc_ChIP_RPKM <- ranLoc_ChIP_RPM/
                    (width(ranLocGR)/1e+03)
save(ranLoc_ChIP_RPM,
     file = paste0(ChIPname, "_RPM_at_", ChIPname, "_peaks_ranLoc.RData"))
save(ranLoc_ChIP_RPKM,
     file = paste0(ChIPname, "_RPKM_at_", ChIPname, "_peaks_ranLoc.RData"))

## input
# Calculate RPM and RPKM for each feature and random locus
peaks_input_RPM <- countOverlaps(peaksGR, input_rangedGR)/
                   input_RPM_scaling_factor
peaks_input_RPKM <- peaks_input_RPM/
                    (width(peaksGR)/1e+03)
save(peaks_input_RPM,
     file = paste0(inputname, "_RPM_at_", ChIPname, "_peaks.RData"))
save(peaks_input_RPKM,
     file = paste0(inputname, "_RPKM_at_", ChIPname, "_peaks.RData"))

ranLoc_input_RPM <- countOverlaps(ranLocGR, input_rangedGR)/
                    input_RPM_scaling_factor
ranLoc_input_RPKM <- ranLoc_input_RPM/
                     (width(ranLocGR)/1e+03)
save(ranLoc_input_RPM,
     file = paste0(inputname, "_RPM_at_", ChIPname, "_peaks_ranLoc.RData"))
save(ranLoc_input_RPKM,
     file = paste0(inputname, "_RPKM_at_", ChIPname, "_peaks_ranLoc.RData"))

## Calculate log2(ChIP/input) RPM and RPKM ratios
peaks_log2_ChIP_input_RPM <- log2((peaks_ChIP_RPM+1)/(peaks_input_RPM+1))
ranLoc_log2_ChIP_input_RPM <- log2((ranLoc_ChIP_RPM+1)/(ranLoc_input_RPM+1))

peaks_log2_ChIP_input_RPKM <- log2((peaks_ChIP_RPKM+1)/(peaks_input_RPKM+1))
ranLoc_log2_ChIP_input_RPKM <- log2((ranLoc_ChIP_RPKM+1)/(ranLoc_input_RPKM+1))

## Plot feature width histogram, and log2(ChIP/input) RPKM and RPM vs
# width and cumulative fraction of loci (features and random loci)
pdf(paste0(ChIPname, "_peaks_width_hist_and_log2ChIPinput_RPKM_RPM_vs_locus_width_and_ecdf.pdf"),
    height = 8, width = 12)
par(mfrow = c(2, 3),
    mar = c(6, 6, 2, 2),
    mgp = c(4, 1.5, 0))

## log2(ChIP RPKM/input RPKM)
# feature width histogram
hist(width(peaksGR), 
     breaks = 250,
     col = "grey50",
     border = NA,
     lwd = 2,
     xlab = "Peak width (bp)",
     ylab = "Peaks",
     main = "",
     cex.lab = 2, cex.axis = 2)
abline(v = mean(width(peaksGR)),
       col = "red", lty = 2, lwd = 2)

# log2(ChIP RPKM/input RPKM) vs random locus and feature width
plot(x = ranLoc_log2_ChIP_input_RPKM,
     y = width(ranLocGR),
     type = "p", pch = 16, cex = 0.25, col = "black",
     log = "y",
     xlim = c(min(c(ranLoc_log2_ChIP_input_RPKM,
                    peaks_log2_ChIP_input_RPKM)),
              max(c(ranLoc_log2_ChIP_input_RPKM,
                    peaks_log2_ChIP_input_RPKM))),
     ylim = c(10, 10000),
     xaxt = "n", yaxt = "n",
     xlab = "", ylab = "")
axis(side = 1,
     at = seq(round(min(c(ranLoc_log2_ChIP_input_RPKM,
                          peaks_log2_ChIP_input_RPKM))),
              round(max(c(ranLoc_log2_ChIP_input_RPKM,
                          peaks_log2_ChIP_input_RPKM))),
              by = 2),
     lwd.tick = 2,
     cex.axis = 2)
axis(side = 2,
     at = c(seq(10, 90, by = 10),
            seq(100, 900, by = 100),
            seq(1000, 10000, by = 1000)),
     labels = c("10", rep("", 8),
                expression("10"^"2"), rep("", 8),
                expression("10"^"3"), rep("", 8),
                expression("10"^"4")),
     lwd.tick = 2,
     cex.axis = 2)
par(new = T)
plot(x = peaks_log2_ChIP_input_RPKM,
     y = width(peaksGR),
     type = "p", pch = 16, cex = 0.25, col = "red",
     log = "y",
     xlim = c(min(c(ranLoc_log2_ChIP_input_RPKM,
                    peaks_log2_ChIP_input_RPKM)),
              max(c(ranLoc_log2_ChIP_input_RPKM,
                    peaks_log2_ChIP_input_RPKM))),
     ylim = c(10, 10000),
     xaxt = "n", yaxt = "n",
     xlab = expression("Log"[2]*"(ChIP RPKM/input RPKM)"),
     ylab = "Locus width (bp)",
     main = bquote(.(ChIPname)~"peaks"),
     cex.lab = 2, cex.main = 2)
abline(h = mean(width(peaksGR)),
       col = "red", lty = 2, lwd = 2)
legend("topleft",
       legend = c(as.expression(bquote("Peaks"~
                                       italic("r"[s])~
                                       "="~
                                       .(round(cor(peaks_log2_ChIP_input_RPKM,
                                                   width(peaksGR),
                                                   method = "spearman"),
                                               digits = 2)))),
                  as.expression(bquote("Random loci"~
                                       italic("r"[s])~
                                       "="~
                                       .(round(cor(ranLoc_log2_ChIP_input_RPKM,
                                                   width(ranLocGR),
                                                   method = "spearman"),
                                               digits = 2))))),
       col = "white",
       text.col = c("red", "black"),
       ncol = 1, cex = 1, lwd = 1, bty = "n") 
box(lwd = 2)

# log2(ChIP RPKM/input RPKM) vs cumulative fraction of loci
plot(ecdf(ranLoc_log2_ChIP_input_RPKM),
     pch = 19, cex = 0.5, col = "black",
     xlim = c(min(c(ranLoc_log2_ChIP_input_RPKM,
                    peaks_log2_ChIP_input_RPKM)),
              max(c(ranLoc_log2_ChIP_input_RPKM,
                    peaks_log2_ChIP_input_RPKM))),
     xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", main = "")
axis(side = 1,
     at = seq(round(min(c(ranLoc_log2_ChIP_input_RPKM,
                          peaks_log2_ChIP_input_RPKM))),
              round(max(c(ranLoc_log2_ChIP_input_RPKM,
                          peaks_log2_ChIP_input_RPKM))),
              by = 2),
     lwd.tick = 2,
     cex.axis = 2)
axis(side = 2,
     at = seq(0, 1, by = 0.25),
     labels = c("0", "", "0.5", "", "1"),
     lwd.tick = 2,
     cex.axis = 2)
par(new = T)
plot(ecdf(peaks_log2_ChIP_input_RPKM),
     pch = 19, cex = 0.5, col = "red",
     xlim = c(min(c(ranLoc_log2_ChIP_input_RPKM,
                    peaks_log2_ChIP_input_RPKM)),
              max(c(ranLoc_log2_ChIP_input_RPKM,
                    peaks_log2_ChIP_input_RPKM))),
     xaxt = "n", yaxt = "n",
     xlab = expression("Log"[2]*"(ChIP RPKM/input RPKM)"),
     ylab = "Cumulative fraction of loci",
     main = "",
     cex.lab = 2)
legend("left",
       legend = c("Peaks", "Random loci"),
       col = "white",
       text.col = c("red", "black"),
       ncol = 1, cex = 1, lwd = 1, bty = "n") 
box(lwd = 2)


## log2(ChIP RPM/input RPM)
# feature width histogram
hist(width(peaksGR), 
     breaks = 250,
     col = "grey50",
     border = NA,
     lwd = 2,
     xlab = "Peak width (bp)",
     ylab = "Peaks",
     main = "",
     cex.lab = 2, cex.axis = 2)
abline(v = mean(width(peaksGR)),
       col = "red", lty = 2, lwd = 2)

# log2(ChIP RPM/input RPM) vs random locus and feature width
plot(x = ranLoc_log2_ChIP_input_RPM,
     y = width(ranLocGR),
     type = "p", pch = 16, cex = 0.25, col = "black",
     log = "y",
     xlim = c(min(c(ranLoc_log2_ChIP_input_RPM,
                    peaks_log2_ChIP_input_RPM)),
              max(c(ranLoc_log2_ChIP_input_RPM,
                    peaks_log2_ChIP_input_RPM))),
     ylim = c(10, 10000),
     xaxt = "n", yaxt = "n",
     xlab = "", ylab = "")
axis(side = 1,
     at = seq(round(min(c(ranLoc_log2_ChIP_input_RPM,
                          peaks_log2_ChIP_input_RPM))),
              round(max(c(ranLoc_log2_ChIP_input_RPM,
                          peaks_log2_ChIP_input_RPM))),
              by = 1),
     lwd.tick = 2,
     cex.axis = 2)
axis(side = 2,
     at = c(seq(10, 90, by = 10),
            seq(100, 900, by = 100),
            seq(1000, 10000, by = 1000)),
     labels = c("10", rep("", 8),
                expression("10"^"2"), rep("", 8),
                expression("10"^"3"), rep("", 8),
                expression("10"^"4")),
     lwd.tick = 2,
     cex.axis = 2)
par(new = T)
plot(x = peaks_log2_ChIP_input_RPM,
     y = width(peaksGR),
     type = "p", pch = 16, cex = 0.25, col = "red",
     log = "y",
     xlim = c(min(c(ranLoc_log2_ChIP_input_RPM,
                    peaks_log2_ChIP_input_RPM)),
              max(c(ranLoc_log2_ChIP_input_RPM,
                    peaks_log2_ChIP_input_RPM))),
     ylim = c(10, 10000),
     xaxt = "n", yaxt = "n",
     xlab = expression("Log"[2]*"(ChIP RPM/input RPM)"),
     ylab = "Locus width (bp)",
     main = bquote(.(ChIPname)~"peaks"),
     cex.lab = 2, cex.main = 2)
abline(h = mean(width(peaksGR)),
       col = "red", lty = 2, lwd = 2)
legend("topleft",
       legend = c(as.expression(bquote("Peaks"~
                                       italic("r"[s])~
                                       "="~
                                       .(round(cor(peaks_log2_ChIP_input_RPM,
                                                   width(peaksGR),
                                                   method = "spearman"),
                                               digits = 2)))),
                  as.expression(bquote("Random loci"~
                                       italic("r"[s])~
                                       "="~
                                       .(round(cor(ranLoc_log2_ChIP_input_RPM,
                                                   width(ranLocGR),
                                                   method = "spearman"),
                                               digits = 2))))),
       col = "white",
       text.col = c("red", "black"),
       ncol = 1, cex = 1, lwd = 1, bty = "n")
box(lwd = 2)

# log2(ChIP RPM/input RPM) vs cumulative fraction of loci
plot(ecdf(ranLoc_log2_ChIP_input_RPM),
     pch = 19, cex = 4, col = "black",
     xlim = c(min(c(ranLoc_log2_ChIP_input_RPM,
                    peaks_log2_ChIP_input_RPM)),
              max(c(ranLoc_log2_ChIP_input_RPM,
                    peaks_log2_ChIP_input_RPM))),
     xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", main = "")
axis(side = 1,
     at = seq(round(min(c(ranLoc_log2_ChIP_input_RPM,
                          peaks_log2_ChIP_input_RPM))),
              round(max(c(ranLoc_log2_ChIP_input_RPM,
                          peaks_log2_ChIP_input_RPM))),
              by = 1),
     lwd.tick = 2,
     cex.axis = 2)
axis(side = 2,
     at = seq(0, 1, by = 0.25),
     labels = c("0", "", "0.5", "", "1"),
     lwd.tick = 2,
     cex.axis = 2)
par(new = T)
plot(ecdf(peaks_log2_ChIP_input_RPM),
     pch = 19, cex = 4, col = "red",
     xlim = c(min(c(ranLoc_log2_ChIP_input_RPM,
                    peaks_log2_ChIP_input_RPM)),
              max(c(ranLoc_log2_ChIP_input_RPM,
                    peaks_log2_ChIP_input_RPM))),
     xaxt = "n", yaxt = "n",
     xlab = expression("Log"[2]*"(ChIP RPM/input RPM)"),
     ylab = "Cumulative fraction of loci",
     main = "",
     cex.lab = 2)
legend("left",
       legend = c("Peaks", "Random loci"),
       col = "white",
       text.col = c("red", "black"),
       ncol = 1, cex = 1, lwd = 1, bty = "n") 
box(lwd = 2)
dev.off()

### ChIP RPKM
## feature width histogram
#hist(width(peaksGR), 
#     breaks = 250,
#     col = "grey50",
#     border = NA,
#     lwd = 2,
#     xlab = "Peak width (bp)",
#     ylab = "Peaks",
#     main = "",
#     cex.lab = 2, cex.axis = 2)
#abline(v = mean(width(peaksGR)),
#       col = "red", lty = 2, lwd = 2)
#
## ChIP RPKM vs random locus and feature width
#plot(x = ranLoc_ChIP_RPKM+1,
#     y = width(ranLocGR),
#     type = "p", pch = 16, cex = 0.25, col = "black",
#     log = "y",
#     xlim = c(1, 1000),
#     ylim = c(10, 10000),
#     xaxt = "n", yaxt = "n",
#     xlab = "", ylab = "")
#axis(side = 1,
#     at = c(1:10,
#            seq(20, 90, by = 10),
#            seq(100, 1000, by = 100)),
#     labels = c("1", rep("", 8),
#                "10", rep("", 8),
#                expression("10"^"2"), rep("", 8),
#                expression("10"^"3")),
#     lwd.tick = 2,
#     cex.axis = 2)
#axis(side = 2,
#     at = c(seq(10, 90, by = 10),
#            seq(100, 900, by = 100),
#            seq(1000, 10000, by = 1000)),
#     labels = c("10", rep("", 8),
#                expression("10"^"2"), rep("", 8),
#                expression("10"^"3"), rep("", 8),
#                expression("10"^"4")),
#     lwd.tick = 2,
#     cex.axis = 2)
#par(new = T)
#plot(x = peaks_ChIP_RPKM+1,
#     y = width(peaksGR),
#     type = "p", pch = 16, cex = 0.25, col = "red",
#     log = "y",
#     xlim = c(1, 1000),
#     ylim = c(10, 10000),
#     xaxt = "n", yaxt = "n",
#     xlab = "ChIP-seq coverage (RPKM+1)",
#     ylab = "Locus width (bp)",
#     main = bquote(.(ChIPname)~"peaks"),
#     cex.lab = 2, cex.main = 2)
#abline(h = mean(width(peaksGR)),
#       col = "red", lty = 2, lwd = 2)
#legend("topleft",
#       legend = c("Peaks", "Random loci"),
#       col = "white",
#       text.col = c("red", "black"),
#       ncol = 1, cex = 1, lwd = 1, bty = "n") 
#box(lwd = 2)
#
## ChIP RPKM vs cumulative fraction of loci
#plot(ecdf(ranLoc_ChIP_RPKM+1),
#     pch = 19, cex = 0.5, col = "black",
#     xlim = c(1, 1000),
#     xaxt = "n", yaxt = "n",
#     xlab = "", ylab = "", main = "")
#axis(side = 1,
#     at = c(1:10,
#            seq(20, 90, by = 10),
#            seq(100, 1000, by = 100)),
#     labels = c("1", rep("", 8),
#                "10", rep("", 8),
#                expression("10"^"2"), rep("", 8),
#                expression("10"^"3")),
#     lwd.tick = 2,
#     cex.axis = 2)
#axis(side = 2,
#     at = seq(0, 1, by = 0.25),
#     labels = c("0", "", "0.5", "", "1"),
#     lwd.tick = 2,
#     cex.axis = 2)
#par(new = T)
#plot(ecdf(peaks_ChIP_RPKM+1),
#     pch = 19, cex = 0.5, col = "red",
#     xlim = c(0, 1000),
#     xaxt = "n", yaxt = "n",
#     xlab = "ChIP-seq coverage (RPKM+1)",
#     ylab = "Cumulative fraction of loci",
#     main = "",
#     cex.lab = 2)
#legend("topleft",
#       legend = c("Peaks", "Random loci"),
#       col = "white",
#       text.col = c("red", "black"),
#       ncol = 1, cex = 1, lwd = 1, bty = "n") 
#box(lwd = 2)
#
### ChIP RPM
## feature width histogram
#hist(width(peaksGR), 
#     breaks = 250,
#     col = "grey50",
#     border = NA,
#     lwd = 2,
#     xlab = "Peak width (bp)",
#     ylab = "Peaks",
#     main = "",
#     cex.lab = 2, cex.axis = 2)
#abline(v = mean(width(peaksGR)),
#       col = "red", lty = 2, lwd = 2)
#
## ChIP RPM vs random locus and feature width
#plot(x = ranLoc_ChIP_RPM+1,
#     y = width(ranLocGR),
#     type = "p", pch = 16, cex = 0.25, col = "black",
#     log = "y",
#     xlim = c(1, 100),
#     ylim = c(10, 10000),
#     xaxt = "n", yaxt = "n",
#     xlab = "", ylab = "")
#axis(side = 1,
#     at = c(1:10,
#            seq(20, 100, by = 10)),
#     labels = c("1", rep("", 8),
#                "10", rep("", 8),
#                expression("10"^"2")),
#     lwd.tick = 2,
#     cex.axis = 2)
#axis(side = 2,
#     at = c(seq(10, 90, by = 10),
#            seq(100, 900, by = 100),
#            seq(1000, 10000, by = 1000)),
#     labels = c("10", rep("", 8),
#                expression("10"^"2"), rep("", 8),
#                expression("10"^"3"), rep("", 8),
#                expression("10"^"4")),
#     lwd.tick = 2,
#     cex.axis = 2)
#par(new = T)
#plot(x = peaks_ChIP_RPM+1,
#     y = width(peaksGR),
#     type = "p", pch = 16, cex = 0.25, col = "red",
#     log = "y",
#     xlim = c(1, 100),
#     ylim = c(10, 10000),
#     xaxt = "n", yaxt = "n",
#     xlab = "ChIP-seq coverage (RPM+1)",
#     ylab = "Locus width (bp)",
#     main = bquote(.(ChIPname)~"peaks"),
#     cex.lab = 2, cex.main = 2)
#abline(h = mean(width(peaksGR)),
#       col = "red", lty = 2, lwd = 2)
#legend("topleft",
#       legend = c("Peaks", "Random loci"),
#       col = "white",
#       text.col = c("red", "black"),
#       ncol = 1, cex = 1, lwd = 1, bty = "n") 
#box(lwd = 2)
#
## ChIP RPM vs cumulative fraction of loci
#plot(ecdf(ranLoc_ChIP_RPM+1),
#     pch = 19, cex = 0.5, col = "black",
#     xlim = c(1, 100),
#     xaxt = "n", yaxt = "n",
#     xlab = "", ylab = "", main = "")
#axis(side = 1,
#     at = c(1:10,
#            seq(20, 100, by = 10)),
#     labels = c("1", rep("", 8),
#                "10", rep("", 8),
#                expression("10"^"2")),
#     lwd.tick = 2,
#     cex.axis = 2)
#axis(side = 2,
#     at = seq(0, 1, by = 0.25),
#     labels = c("0", "", "0.5", "", "1"),
#     lwd.tick = 2,
#     cex.axis = 2)
#par(new = T)
#plot(ecdf(peaks_ChIP_RPM+1),
#     pch = 19, cex = 0.5, col = "red",
#     xlim = c(0, 100),
#     xaxt = "n", yaxt = "n",
#     xlab = "ChIP-seq coverage (RPM+1)",
#     ylab = "Cumulative fraction of loci",
#     main = "",
#     cex.lab = 2)
#legend("topleft",
#       legend = c("Peaks", "Random loci"),
#       col = "white",
#       text.col = c("red", "black"),
#       ncol = 1, cex = 1, lwd = 1, bty = "n") 
#box(lwd = 2)
#dev.off()
