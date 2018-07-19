####################################################################################
# Generate peak density genome profiles, with cumulative chromosomal coordinates   #
####################################################################################

# Usage:
# Rscript peak_density_genomeProfiles_all_REC8_HA_Rep1_peaks.R 10000 10kb

library(segmentSeq)

args <- commandArgs(trailingOnly = TRUE)
winSize <- as.numeric(args[1])
winName <- as.character(args[2])
peakSet <- "REC8_HA_Rep1_peaks"

inDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/"
outDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/genome_wide/peak_density_genomeProfiles/"

# Import peaks as GRanges object
load(paste0(inDir,
            "REC8_HA_Rep1_armrangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData"))
load(paste0(inDir,
            "REC8_HA_Rep1_perirangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData"))
peaksGR <- sort(c(armrangerPeaksGRmergedOverlaps, perirangerPeaksGRmergedOverlaps))
print(peaksGR)

chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
# pericentromeric regions are as defined in Supplemental Table S26 of Ziolkowski et al. (2017) Genes Dev. 31
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)

################################
# make cumulative genomes      #
################################

sumchr <- cumsum(c(0, chrLens))
print(sumchr)
sumchr_tot <- sumchr[length(sumchr)]
print(sumchr_tot)

centromeres <- sapply(seq_along(centromeres), function(x) {
  centromeres[x] + sumchr[x]
})
print(centromeres)
pericenStart <- sapply(seq_along(pericenStart), function(x) {
  pericenStart[x] + sumchr[x]
})
print(pericenStart)
pericenEnd <- sapply(seq_along(pericenEnd), function(x) {
  pericenEnd[x] + sumchr[x]
})
print(pericenEnd)

# Calculate peak density
cumPeakDensity <- NULL
for(i in 1:5) {
  # define windows as GRanges object
  seqWindows <- seq(1, chrLens[i], by = winSize)
  seqWindows <- c(seqWindows, chrLens[i])
  cumWindows <- seqWindows + sumchr[i]
  windowsGR <- GRanges(seqnames = chrs[i],
                       ranges = IRanges(start = seqWindows, width = winSize),
                       strand = "*")
  # count peaks within windows
  chrPeaks <- peaksGR[seqnames(peaksGR) == chrs[i]]
  winPeaks <- countOverlaps(windowsGR, chrPeaks)
  chrDat <- cbind(cumWindows, winPeaks)
  write.table(chrDat,
              file = paste0(outDir, peakSet, "_density_chr", i, "_", winName, ".txt"))
  cumPeakDensity <- rbind(cumPeakDensity, chrDat)
}
write.table(cumPeakDensity,
            file = paste0(outDir, peakSet, "_density_genome_", winName, ".txt"))

chrPeakDensity <- mclapply(1:5, function(i) {
  read.table(paste0(outDir, peakSet, "_density_chr", i, "_", winName, ".txt"))
}, mc.cores = 5)

# smooth density values with MA filter
test <- seq(1, 1000, by = 1)
j = 100
ma <- rep(1, test[j])/test[j]

filt_peakDensity <- NULL
filt_peakDensity_noNA <- NULL
for(i in 1:5) {
  filt_chrPeakDensity <- stats::filter(chrPeakDensity[[i]][,2], ma)
  which_na <- which(is.na(filt_chrPeakDensity) == TRUE)
  left_na <- which_na[which(which_na < 100)]
  left_val <- filt_chrPeakDensity[left_na[length(left_na)]+1]
  filt_chrPeakDensity[left_na] <- left_val
  right_na <- which_na[which(which_na > 100)]
  right_val <- filt_chrPeakDensity[right_na[1]-1]
  filt_chrPeakDensity[right_na] <- right_val
  filt_chrPeakDensity_noNA <- filt_chrPeakDensity[!is.na(filt_chrPeakDensity)]
  filt_chrPeakDensity <- cbind(chrPeakDensity[[i]][,1], filt_chrPeakDensity)
  filt_chrPeakDensity_noNA <- cbind(chrPeakDensity[[i]][,1], filt_chrPeakDensity_noNA)
  write.table(filt_chrPeakDensity,
              file = paste0(outDir, "filt_", peakSet, "_density_chr", i, "_", winName, ".txt"))
  write.table(filt_chrPeakDensity_noNA,
              file = paste0(outDir, "filt_noNA_", peakSet, "_density_chr", i, "_", winName, ".txt"))
  filt_peakDensity <- rbind(filt_peakDensity, filt_chrPeakDensity)
  filt_peakDensity_noNA <- rbind(filt_peakDensity_noNA, filt_chrPeakDensity_noNA)
}
write.table(filt_peakDensity,
            file = paste0(outDir, "filt_", peakSet, "_density_genome_", winName, ".txt"))
write.table(filt_peakDensity_noNA,
            file = paste0(outDir, "filt_noNA_", peakSet, "_density_genome_", winName, ".txt"))

