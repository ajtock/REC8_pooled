#!/applications/R/R-3.3.2/bin/Rscript

# Calculate mean library-size-normalised wt and kss
# coverage between the start and end of each peak
# Select peaks with greatest change in kss (top 10%)

library(segmentSeq)

inDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/"
outDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/genome_wide/REC8_HA_Rep1_wt_kss_diff/"

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

wt_norm <- read.table("/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab.bed",
                      colClasses = c(NA, NA, "NULL", NA))
kss_norm <- read.table("/home/meiosis/ajt200/analysis/180622_Chris_lambing_ChIP_REC8_HA_Col_kss/kss/coverage/common_input_MYC_Rep2/log2ChIPinput/log2_kss_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab.bed",
                       colClasses = c(NA, NA, "NULL", NA))

allCovPeaks <- NULL
for(i in 1:5) {
  print(i)
  wt_chrCov <- wt_norm[wt_norm[,1] == chrs[i],]
  wt_chrCov <- wt_chrCov[,3]
  kss_chrCov <- kss_norm[kss_norm[,1] == chrs[i],]
  kss_chrCov <- kss_chrCov[,3]
  print(i)
  covCoords <- seq(1, length(wt_chrCov), by = 1)
  covGR <- GRanges(seqnames = chrs[i],
                   ranges = IRanges(start = covCoords,
                                    width = 1),
                   strand = "*")
  print(i)
  chrPeaksGR <- peaksGR[seqnames(peaksGR) == chrs[i]]
  peaksOverlaps <- getOverlaps(chrPeaksGR,
                               covGR,
                               whichOverlaps = T)
  wt_norm_peakCov <- sapply(peaksOverlaps,
                            function(x) mean(wt_chrCov[x]))
  kss_norm_peakCov <- sapply(peaksOverlaps,
                             function(x) mean(kss_chrCov[x]))
  print(i)
  chrPeaks <- cbind(rep(chrs[i], length(chrPeaksGR)),
                    start(chrPeaksGR),
                    end(chrPeaksGR),
                    wt_norm_peakCov,
                    kss_norm_peakCov,
                    wt_norm_peakCov-kss_norm_peakCov)
  allCovPeaks <- rbind(allCovPeaks, chrPeaks)
}
colnames(allCovPeaks) <- c("Chr", "Start", "End",
                           "wt_norm_peakCov", "kss_norm_peakCov", "diff_wt_kss_norm_peakCov")

write.table(allCovPeaks,
            file = paste0(outDir,
                          "REC8_HA_Rep1_peaks_mean_log2_REC8_HA_Rep1_wt_kss_diff_cov.txt"),
            sep = "\t", row.names = F, quote = F)

#allCovPeaks <- read.table(paste0(outDir,
#                                 "REC8_HA_Rep1_peaks_mean_log2_MNase_cov.txt"),
#                          header = T)

## Select peaks with MNase coverage <= 0
#lowNucPeaks <- allCovPeaks[allCovPeaks$wt_norm_peakCov <= 0,]                      
#write.table(lowNucPeaks,
#            file = paste0(outDir,
#                          "REC8_HA_Rep1_peaks_mean_log2_MNase_cov_lessThanOrEqualToZero.txt"),
#            sep = "\t", row.names = F, quote = F)


