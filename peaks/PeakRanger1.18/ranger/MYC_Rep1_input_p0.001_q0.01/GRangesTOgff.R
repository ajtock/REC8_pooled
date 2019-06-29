#!/applications/R/R-3.3.2/bin/Rscript

# Convert GRanges objects to gff files

# Usage:
# ./GRangesTOgff.R REC8_HA_Rep1

library(GenomicRanges)

args <- commandArgs(trailingOnly = T)
libName <- args[1]

inDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/"

load(paste0(inDir, libName, "_armrangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData"))
load(paste0(inDir, libName, "_perirangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData"))

armPeaks <- data.frame(armrangerPeaksGRmergedOverlaps)
periPeaks <- data.frame(perirangerPeaksGRmergedOverlaps)

armrangerPeaksGRmergedOverlaps <- NULL
perirangerPeaksGRmergedOverlaps <- NULL

armPeaksgff <- cbind(armPeaks[,1], rep(".", length(armPeaks[,1])), rep(paste0(libName, "_peak")), armPeaks[,2], armPeaks[,3], rep(".", length(armPeaks[,1])), rep(".", length(armPeaks[,1])), rep(".", length(armPeaks[,1])), rep(".", length(armPeaks[,1])))
colnames(armPeaksgff) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
write.table(armPeaksgff, file = paste0(inDir, libName, "_ranger_armPeaks_mergedOverlaps.gff"), row.names = F, col.names = F, quote = F, sep = "\t")

periPeaksgff <- cbind(periPeaks[,1], rep(".", length(periPeaks[,1])), rep(paste0(libName, "_peak")), periPeaks[,2], periPeaks[,3], rep(".", length(periPeaks[,1])), rep(".", length(periPeaks[,1])), rep(".", length(periPeaks[,1])), rep(".", length(periPeaks[,1])))
colnames(periPeaksgff) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
write.table(periPeaksgff, file = paste0(inDir, libName, "_ranger_periPeaks_mergedOverlaps.gff"), row.names = F, col.names = F, quote = F, sep = "\t")


