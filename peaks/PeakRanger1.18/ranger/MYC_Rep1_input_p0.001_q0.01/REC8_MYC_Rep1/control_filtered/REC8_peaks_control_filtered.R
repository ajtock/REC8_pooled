#!/applications/R/R-3.3.2/bin/Rscript

# Remove REC8 peaks that overlap negative-control-library peaks by >= 1 bp
# Convert to GRanges objects
# (with overlapping peaks unmerged or merged)
# sort by decreasing -log10(qval) for use with weeder2

# Usage:
# ./REC8_peaks_control_filtered.R REC8_MYC_Rep1_ChIP

#ChIP_prefix <- "REC8_MYC_Rep1_ChIP"

args <- commandArgs(trailingOnly = TRUE)
ChIP_prefix <- args[1]

library(GenomicRanges)
library(dplyr)

REC8_inDir <- "/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep1_input_p0.001_q0.01/"
control_inDir <- "/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep1_input_p0.05_q0.1/"

rangerPeaks_arm <- read.table(paste0(REC8_inDir, ChIP_prefix,
                                     "_rangerPeaks_arm.txt"))
rangerPeaks_peri <- read.table(paste0(REC8_inDir, ChIP_prefix,
                                      "_rangerPeaks_peri.txt"))
controlMYCPeaks_arm <- read.table(paste0(control_inDir, "control_MYC_Rep1_ChIP",
                                         "_rangerPeaks_arm.txt"))
controlMYCPeaks_peri <- read.table(paste0(control_inDir, "control_MYC_Rep1_ChIP",
                                          "_rangerPeaks_peri.txt"))
controlHAPeaks_arm <- read.table(paste0(control_inDir, "control_HA_Rep1_ChIP",
                                        "_rangerPeaks_arm.txt"))
controlHAPeaks_peri <- read.table(paste0(control_inDir, "control_HA_Rep1_ChIP",
                                         "_rangerPeaks_peri.txt"))

# arm
# Combine controlMYCPeaks and controlHAPeaks
controlPeaksGR_arm <- sort(GRanges(seqnames = paste0("Chr",
                                                     c(controlMYCPeaks_arm$chr,
                                                       controlHAPeaks_arm$chr)),
                                   ranges = IRanges(start = c(controlMYCPeaks_arm$start,
                                                              controlHAPeaks_arm$start),
                                                    end = c(controlMYCPeaks_arm$end,
                                                            controlHAPeaks_arm$end)),
                                   strand = "*"))

# Create GRanges objects sorted by decreasing -log10(qval)
rangerPeaksGR_arm <- sort(GRanges(seqnames = paste0("Chr", rangerPeaks_arm$chr),
                                  ranges = IRanges(start = rangerPeaks_arm$start,
                                                   end = rangerPeaks_arm$end),
                                  strand = "*",
                                  sigval = rangerPeaks_arm$sigval,
                                  pval = rangerPeaks_arm$pval,
                                  qval = rangerPeaks_arm$qval,
                                  summit0based = rangerPeaks_arm$summit0based),
                          by = ~ qval, decreasing = T)

# Remove REC8 peaks overlapping combined control peaks 
rangerPeaksGR_arm_control_overlaps <- findOverlaps(query = controlPeaksGR_arm,
                                                   subject = rangerPeaksGR_arm,
                                                   ignore.strand = TRUE,
                                                   select = "all")
rangerPeaksGR_arm_control_filtered <- rangerPeaksGR_arm[-subjectHits(rangerPeaksGR_arm_control_overlaps)]
save(rangerPeaksGR_arm_control_filtered,
     file = paste0("./", ChIP_prefix,
                   "_rangerPeaksGR_arm_control_filtered_minuslog10Qsorted_noMinWidth.RData"))
# Merge overlapping peaks
rangerPeaksGR_arm_control_filtered_mergedOverlaps <- reduce(rangerPeaksGR_arm_control_filtered)
print(rangerPeaksGR_arm_control_filtered_mergedOverlaps)
save(rangerPeaksGR_arm_control_filtered_mergedOverlaps,
     file = paste0("./", ChIP_prefix,
                   "_rangerPeaksGR_arm_control_filtered_mergedOverlaps_noMinWidth.RData"))
# Remove peaks < 50 bp
rangerPeaksGR_arm_control_filtered_50bpMinWidth <- rangerPeaksGR_arm_control_filtered[width(rangerPeaksGR_arm_control_filtered) >= 50]
save(rangerPeaksGR_arm_control_filtered_50bpMinWidth,
     file = paste0("./", ChIP_prefix,
                   "_rangerPeaksGR_arm_control_filtered_minuslog10Qsorted_50bpMinWidth.RData"))
rangerPeaksGR_arm_control_filtered_mergedOverlaps_50bpMinWidth <- rangerPeaksGR_arm_control_filtered_mergedOverlaps[width(rangerPeaksGR_arm_control_filtered_mergedOverlaps) >= 50]
save(rangerPeaksGR_arm_control_filtered_mergedOverlaps_50bpMinWidth,
     file = paste0("./", ChIP_prefix,
                   "_rangerPeaksGR_arm_control_filtered_mergedOverlaps_50bpMinWidth.RData"))

# peri
# Combine controlMYCPeaks and controlHAPeaks
controlPeaksGR_peri <- sort(GRanges(seqnames = paste0("Chr",
                                                      c(controlMYCPeaks_peri$chr,
                                                        controlHAPeaks_peri$chr)),
                                    ranges = IRanges(start = c(controlMYCPeaks_peri$start,
                                                               controlHAPeaks_peri$start),
                                                     end = c(controlMYCPeaks_peri$end,
                                                             controlHAPeaks_peri$end)),
                                    strand = "*"))

# Create GRanges objects sorted by decreasing -log10(qval)
rangerPeaksGR_peri <- sort(GRanges(seqnames = paste0("Chr", rangerPeaks_peri$chr),
                                   ranges = IRanges(start = rangerPeaks_peri$start,
                                                    end = rangerPeaks_peri$end),
                                   strand = "*",
                                   sigval = rangerPeaks_peri$sigval,
                                   pval = rangerPeaks_peri$pval,
                                   qval = rangerPeaks_peri$qval,
                                   summit0based = rangerPeaks_peri$summit0based),
                           by = ~ qval, decreasing = T)

# Remove REC8 peaks overlapping combined control peaks 
rangerPeaksGR_peri_control_overlaps <- findOverlaps(query = controlPeaksGR_peri,
                                                    subject = rangerPeaksGR_peri,
                                                    ignore.strand = TRUE,
                                                    select = "all")
rangerPeaksGR_peri_control_filtered <- rangerPeaksGR_peri[-subjectHits(rangerPeaksGR_peri_control_overlaps)]
save(rangerPeaksGR_peri_control_filtered,
     file = paste0("./", ChIP_prefix,
                   "_rangerPeaksGR_peri_control_filtered_minuslog10Qsorted_noMinWidth.RData"))
# Merge overlapping peaks
rangerPeaksGR_peri_control_filtered_mergedOverlaps <- reduce(rangerPeaksGR_peri_control_filtered)
print(rangerPeaksGR_peri_control_filtered_mergedOverlaps)
save(rangerPeaksGR_peri_control_filtered_mergedOverlaps,
     file = paste0("./", ChIP_prefix,
                   "_rangerPeaksGR_peri_control_filtered_mergedOverlaps_noMinWidth.RData"))
# Remove peaks < 50 bp
rangerPeaksGR_peri_control_filtered_50bpMinWidth <- rangerPeaksGR_peri_control_filtered[width(rangerPeaksGR_peri_control_filtered) >= 50]
save(rangerPeaksGR_peri_control_filtered_50bpMinWidth,
     file = paste0("./", ChIP_prefix,
                   "_rangerPeaksGR_peri_control_filtered_minuslog10Qsorted_50bpMinWidth.RData"))
rangerPeaksGR_peri_control_filtered_mergedOverlaps_50bpMinWidth <- rangerPeaksGR_peri_control_filtered_mergedOverlaps[width(rangerPeaksGR_peri_control_filtered_mergedOverlaps) >= 50]
save(rangerPeaksGR_peri_control_filtered_mergedOverlaps_50bpMinWidth,
     file = paste0("./", ChIP_prefix,
                   "_rangerPeaksGR_peri_control_filtered_mergedOverlaps_50bpMinWidth.RData"))
