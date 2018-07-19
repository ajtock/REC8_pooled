#!/applications/R/R-3.3.2/bin/Rscript

# Usage:
# ./peak_width_histogram.R REC8_HA_Rep1

library(GenomicRanges)

args <- commandArgs(trailingOnly = T)
libName <- args[1]

inDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/"

breaks <- 250

# Arm peaks
load(paste0(inDir, libName,
            "_armrangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData"))
load(paste0(inDir, libName,
            "_armrangerPeaksGR_minuslog10_p0.001_q0.01_qval_sorted_noMinWidth.RData"))

pdf(paste0(inDir, libName, "_armrangerPeaksGR_hist.pdf"),
    height = 4, width = 8)
par(mfrow = c(1, 2),
    mar = c(4, 4, 2, 2),
    mgp = c(3, 1, 0))
hist(width(armrangerPeaksGRmergedOverlaps),
     breaks = breaks,
     col = "grey60",
     border = NA,
     lwd = 2,
     xlab = paste0(libName, " arm peak width (bp)"),
     ylab = "Peaks",
     main = "Overlaps merged",
     cex.main = 1,
     cex.lab = 1,
     cex.axis = 1)
abline(v = mean(width(armrangerPeaksGRmergedOverlaps)),
       col = "red", lty = 2, lwd = 1)
hist(width(armrangerPeaksGR),
     breaks = breaks,
     col = "grey60",
     border = NA,
     lwd = 2,
     xlab = paste0(libName, " arm peak width (bp)"),
     ylab = "Peaks",
     main = "Overlaps unmerged",
     cex.main = 1,
     cex.lab = 1,
     cex.axis = 1)
abline(v = mean(width(armrangerPeaksGR)),
       col = "red", lty = 2, lwd = 1)
dev.off()

# Pericentromeric peaks
load(paste0(inDir, libName,
            "_perirangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData"))
load(paste0(inDir, libName,
            "_perirangerPeaksGR_minuslog10_p0.001_q0.01_qval_sorted_noMinWidth.RData"))

pdf(paste0(inDir, libName, "_perirangerPeaksGR_hist.pdf"),
    height = 4, width = 8)
par(mfrow = c(1, 2),
    mar = c(4, 4, 2, 2),
    mgp = c(3, 1, 0))
hist(width(perirangerPeaksGRmergedOverlaps),
     breaks = breaks,
     col = "grey60",
     border = NA,
     lwd = 2,
     xlab = paste0(libName, " peri peak width (bp)"),
     ylab = "Peaks",
     main = "Overlaps merged",
     cex.main = 1,
     cex.lab = 1,
     cex.axis = 1)
abline(v = mean(width(perirangerPeaksGRmergedOverlaps)),
       col = "red", lty = 2, lwd = 1)
hist(width(perirangerPeaksGR),
     breaks = breaks,
     col = "grey60",
     border = NA,
     lwd = 2,
     xlab = paste0(libName, " peri peak width (bp)"),
     ylab = "Peaks",
     main = "Overlaps unmerged",
     cex.main = 1,
     cex.lab = 1,
     cex.axis = 1)
abline(v = mean(width(perirangerPeaksGR)),
       col = "red", lty = 2, lwd = 1)
dev.off()

