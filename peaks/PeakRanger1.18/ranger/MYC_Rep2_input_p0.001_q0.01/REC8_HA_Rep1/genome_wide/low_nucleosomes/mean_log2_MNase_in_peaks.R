#!/applications/R/R-3.3.2/bin/Rscript

# Calculate mean Z-score-standardised log2-transformed library-size-normalised MNase
# coverage between the start and end of each peak
# Select peaks with MNase coverage <= 0

library(segmentSeq)

inDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/"
outDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/genome_wide/low_nucleosomes/"

# Chromosome definitions
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")

# Import peaks as GRanges object
load(paste0(inDir,
            "REC8_HA_Rep1_armrangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData"))
load(paste0(inDir,
            "REC8_HA_Rep1_perirangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData"))
peaksGR <- sort(c(armrangerPeaksGRmergedOverlaps, perirangerPeaksGRmergedOverlaps))
strand(peaksGR) <- "*"
print("***********peaks***********")
print(peaksGR)

MNase_norm <- read.table("/projects/ajt200/BAM_masters/nucleosomes/WT/coverage/nakedDNA_untrimmed_input/log2ChIPinput/log2wtNucNakedDNAuntrimmed_norm_allchrs_coverage_coord_tab.bed")
allCovPeaks <- NULL
for(i in 1:5) {
  print(i)
  chrCov <- MNase_norm[MNase_norm[,1] == chrs[i],]
  chrCov <- chrCov[,4]
  print(i)
  covCoords <- seq(1, length(chrCov), by = 1)
  covGR <- GRanges(seqnames = chrs[i],
                   ranges = IRanges(start = covCoords,
                                    width = 1),
                   strand = "*")
  print(i)
  chrPeaksGR <- peaksGR[seqnames(peaksGR) == chrs[i]]
  peaksOverlaps <- getOverlaps(chrPeaksGR,
                               covGR,
                               whichOverlaps = T)
  MNase_norm_peakCov <- sapply(peaksOverlaps,
                               function(x) mean(chrCov[x]))
  print(i)
  chrPeaks <- cbind(rep(chrs[i], length(chrPeaksGR)),
                    start(chrPeaksGR),
                    end(chrPeaksGR),
                    MNase_norm_peakCov)
  allCovPeaks <- rbind(allCovPeaks, chrPeaks)
}
colnames(allCovPeaks) <- c("Chr", "Start", "End", "MNase_norm_peakCov")
write.table(allCovPeaks,
            file = paste0(outDir,
                          "REC8_HA_Rep1_peaks_mean_log2_MNase_cov.txt"),
             sep = "\t", row.names = F, quote = F)

allCovPeaks <- read.table(paste0(outDir,
                                 "REC8_HA_Rep1_peaks_mean_log2_MNase_cov.txt"),
                          header = T)

# Select peaks with MNase coverage <= 0
lowNucPeaks <- allCovPeaks[allCovPeaks$MNase_norm_peakCov <= 0,]                      
write.table(lowNucPeaks,
            file = paste0(outDir,
                          "REC8_HA_Rep1_peaks_mean_log2_MNase_cov_lessThanOrEqualToZero.txt"),
            sep = "\t", row.names = F, quote = F)


