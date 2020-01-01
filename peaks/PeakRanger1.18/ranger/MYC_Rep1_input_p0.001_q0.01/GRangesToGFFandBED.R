#!/applications/R/R-3.3.2/bin/Rscript

# Convert GRanges objects to GFF files

# Usage:
# ./GRangesToGFFandBED.R REC8_HA_Rep2_ChIP

library(GenomicRanges)

args <- commandArgs(trailingOnly = T)
libName <- args[1]

load(paste0(libName, "_rangerPeaksGR_arm_mergedOverlaps_noMinWidth.RData"))
load(paste0(libName, "_rangerPeaksGR_peri_mergedOverlaps_noMinWidth.RData"))

armPeaks <- data.frame(rangerPeaksGR_arm_mergedOverlaps)
periPeaks <- data.frame(rangerPeaksGR_peri_mergedOverlaps)

rangerPeaksGR_arm_mergedOverlaps <- NULL
rangerPeaksGR_peri_mergedOverlaps <- NULL

#armPeaksGFF <- cbind(armPeaks[,1], rep(".", length(armPeaks[,1])), rep(paste0(libName, "_peak")), armPeaks[,2], armPeaks[,3], rep(".", length(armPeaks[,1])), rep(".", length(armPeaks[,1])), rep(".", length(armPeaks[,1])), rep(".", length(armPeaks[,1])))
#colnames(armPeaksGFF) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
#write.table(armPeaksGFF, file = paste0(libName, "_rangerPeaksGR_arm_mergedOverlaps.gff"), row.names = F, col.names = F, quote = F, sep = "\t")
#
#periPeaksGFF <- cbind(periPeaks[,1], rep(".", length(periPeaks[,1])), rep(paste0(libName, "_peak")), periPeaks[,2], periPeaks[,3], rep(".", length(periPeaks[,1])), rep(".", length(periPeaks[,1])), rep(".", length(periPeaks[,1])), rep(".", length(periPeaks[,1])))
#colnames(periPeaksGFF) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
#write.table(periPeaksGFF, file = paste0(libName, "_rangerPeaksGR_peri_mergedOverlaps.gff"), row.names = F, col.names = F, quote = F, sep = "\t")

armPeaksBED <- data.frame(chr = as.character(armPeaks$seqnames),
                          start = as.integer(armPeaks$start-1),
                          end = as.integer(armPeaks$end),
                          stringsAsFactors = F)
write.table(armPeaksBED,
            file = paste0(libName, "_ranger_peaks_in_chromosome_arms.bed"),
            sep = "\t", quote = F,
            row.names = F, col.names = F)

periPeaksBED <- data.frame(chr = as.character(periPeaks$seqnames),
                           start = as.integer(periPeaks$start-1),
                           end = as.integer(periPeaks$end),
                           stringsAsFactors = F)
write.table(periPeaksBED,
            file = paste0(libName, "_ranger_peaks_in_centromeres_and_pericentromeres.bed"),
            sep = "\t", quote = F,
            row.names = F, col.names = F)
