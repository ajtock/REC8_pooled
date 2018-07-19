#!/applications/R/R-3.3.2/bin/Rscript

# Separate arm and pericentromeric peaks
# Convert to GRanges objects
# (with overlapping peaks unmerged or merged)
# sort by decreasing -log10(qval) for use with weeder2

# Usage:
# peaks_peakranger_ranger_minuslog10_p0.001_q0.01_TreadsNormCreads_sort_by_minuslog10_qval_arm_peri_fmt_peakSummits.R REC8_HA_Rep1

args <- commandArgs(trailingOnly = TRUE)
libName <- args[1]

library(GenomicRanges)
library(dplyr)

# Definitions
peakDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/"
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)

rangerPeaks <- read.table(paste0(peakDir,
                                 libName,
                                 "_peaks_peakranger_ranger_p0.001_q0.01_TreadsNormCreads.narrowPeak"))
rangerPeaks <- cbind(rangerPeaks[,1:3],
                     rangerPeaks[,7:10])
colnames(rangerPeaks) <- c("chr", "start0based", "end",
                           "sigval", "pval", "qval", "summit0based")
rangerPeaks <- as.data.frame(cbind(rangerPeaks$chr,
                                   rangerPeaks$start0based+1,
                                   rangerPeaks$end,
                                   rangerPeaks$sigval,
                                   rangerPeaks$pval,
                                   rangerPeaks$qval,
                                   rangerPeaks$summit0based))
colnames(rangerPeaks) <- c("chr", "start", "end",
                           "sigval", "pval", "qval", "summit0based")

rangerPeaks.arms <- NULL
rangerPeaks.peri <- NULL
for(i in 1:5) {
  chr.rangerPeaks <- rangerPeaks[rangerPeaks$chr == i,]
  chr.rangerPeaks.arms <- chr.rangerPeaks %>%
    filter(end < pericenStart[i] | start > pericenEnd[i])
  chr.rangerPeaks.peri <- chr.rangerPeaks %>%
    filter(end > pericenStart[i] & start < pericenEnd[i])
  rangerPeaks.arms <- rbind(rangerPeaks.arms,
                            chr.rangerPeaks.arms)
  rangerPeaks.peri <- rbind(rangerPeaks.peri,
                            chr.rangerPeaks.peri)
}
print("Unfiltered arm peaks:")
print(dim(rangerPeaks.arms)[[1]])
print("Unfiltered pericentromeric and centromeric peaks:")
print(dim(rangerPeaks.peri)[[1]])
print("Unfiltered peaks:")
print(dim(rangerPeaks.arms)[[1]] + dim(rangerPeaks.peri)[[1]])

write.table(rangerPeaks.arms,
            file = paste0(peakDir,
                          libName,
                          "_arm_peaks_peakranger_ranger_minuslog10_p0.001_q0.01.txt"))
write.table(rangerPeaks.peri,
            file = paste0(peakDir,
                          libName,
                          "_peri_peaks_peakranger_ranger_minuslog10_p0.001_q0.01.txt"))

### Create GRanges objects sorted by decreasing -log10(qval)
# arm
armrangerPeaksGR <- sort(GRanges(seqnames = paste0("Chr", rangerPeaks.arms$chr),
                                 ranges = IRanges(start = rangerPeaks.arms$start,
                                                  end = rangerPeaks.arms$end),
                                 strand = "+",
                                 sigval = rangerPeaks.arms$sigval,
                                 pval = rangerPeaks.arms$pval,
                                 qval = rangerPeaks.arms$qval,
                                 summit0based = rangerPeaks.arms$summit0based),
                         by = ~ qval, decreasing = T)
save(armrangerPeaksGR,
     file = paste0(peakDir,
                   libName,
                   "_armrangerPeaksGR_minuslog10_p0.001_q0.01_qval_sorted_noMinWidth.RData"))
# Merge overlapping peaks
armrangerPeaksGRmergedOverlaps <- reduce(armrangerPeaksGR)
save(armrangerPeaksGRmergedOverlaps,
     file = paste0(peakDir,
                   libName,
                   "_armrangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData"))
# Remove peaks < 50 bp
armrangerPeaksGR_50bpMinWidth <- armrangerPeaksGR[width(armrangerPeaksGR) >= 50]
save(armrangerPeaksGR_50bpMinWidth,
     file = paste0(peakDir,
                   libName,
                   "_armrangerPeaksGR_minuslog10_p0.001_q0.01_qval_sorted_50bpMinWidth.RData"))
armrangerPeaksGRmergedOverlaps_50bpMinWidth <- armrangerPeaksGRmergedOverlaps[width(armrangerPeaksGRmergedOverlaps) >= 50]
save(armrangerPeaksGRmergedOverlaps_50bpMinWidth,
     file = paste0(peakDir,
                   libName,
                   "_armrangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_50bpMinWidth.RData"))

# peri
perirangerPeaksGR <- sort(GRanges(seqnames = paste0("Chr", rangerPeaks.peri$chr),
                                  ranges = IRanges(start = rangerPeaks.peri$start,
                                                   end = rangerPeaks.peri$end),
                                  strand = "+",
                                  sigval = rangerPeaks.peri$sigval,
                                  pval = rangerPeaks.peri$pval,
                                  qval = rangerPeaks.peri$qval,
                                  summit0based = rangerPeaks.peri$summit0based),
                          by = ~ qval, decreasing = T)
save(perirangerPeaksGR,
     file = paste0(peakDir,
                   libName,
                   "_perirangerPeaksGR_minuslog10_p0.001_q0.01_qval_sorted_noMinWidth.RData"))
# Merge overlapping peaks
perirangerPeaksGRmergedOverlaps <- reduce(perirangerPeaksGR)
save(perirangerPeaksGRmergedOverlaps,
     file = paste0(peakDir,
                   libName,
                   "_perirangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData"))
# Remove peaks < 50 bp
perirangerPeaksGR_50bpMinWidth <- perirangerPeaksGR[width(perirangerPeaksGR) >= 50]
save(perirangerPeaksGR_50bpMinWidth,
     file = paste0(peakDir,
                   libName,
                   "_perirangerPeaksGR_minuslog10_p0.001_q0.01_qval_sorted_50bpMinWidth.RData"))
perirangerPeaksGRmergedOverlaps_50bpMinWidth <- perirangerPeaksGRmergedOverlaps[width(perirangerPeaksGRmergedOverlaps) >= 50]
save(perirangerPeaksGRmergedOverlaps_50bpMinWidth,
     file = paste0(peakDir,
                   libName,
                   "_perirangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_50bpMinWidth.RData"))

