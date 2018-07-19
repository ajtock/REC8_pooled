
# Create histograms of width distributions of REC8-HA Rep1 peaks in each differential set

# Usage:
# Rscript peak_set_width_dists.R

library(GenomicRanges)

inDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/"
plotDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/genome_wide/peak_set_overlaps_and_width_dists/plots/"

load(paste0(inDir,
            "REC8_HA_Rep1_armrangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData"))
load(paste0(inDir,
            "REC8_HA_Rep1_perirangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData"))
peaksGR <- sort(c(armrangerPeaksGRmergedOverlaps, perirangerPeaksGRmergedOverlaps))
strand(peaksGR) <- "*"
load("/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/genome_wide/REC8_HA_Rep1_wt_kss_diff/REC8_HA_Rep1_peaks_mean_log2_REC8_HA_Rep1_wt_kss_diff_cov_top10percentGR.RData")
REC8_GR <- allCovPeaks_top10percentGR
allCovPeaks_top10percentGR <- NULL
load("/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/genome_wide/RNAseq_Rep1_kss_wt_diff/REC8_HA_Rep1_peaks_mean_RNAseq_Rep1_kss_wt_diff_cov_top10percentGR.RData")
RNAseq_GR <- allCovPeaks_top10percentGR
allCovPeaks_top10percentGR <- NULL
load("/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/genome_wide/SPO11oligo_Rep1_kss_wt_diff/REC8_HA_Rep1_peaks_mean_log2_SPO11oligo_Rep1_kss_wt_diff_cov_top10percentGR.RData")
SPO11oligo_GR <- allCovPeaks_top10percentGR
allCovPeaks_top10percentGR <- NULL
load("/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/genome_wide/H3K9me2_wt_kss_diff/REC8_HA_Rep1_peaks_mean_log2_H3K9me2_wt_kss_diff_cov_top10percentGR.RData")
H3K9me2_GR <- allCovPeaks_top10percentGR
allCovPeaks_top10percentGR <- NULL

peakSetList <- list(peaksGR,
                    REC8_GR,
                    RNAseq_GR,
                    SPO11oligo_GR,
                    H3K9me2_GR)

names <- c("REC8-HA Rep1",
           "REC8-HA Rep1 wt-kss",
           "RNA-seq kss-wt",
           "SPO11-1-oligos kss-wt",
           "H3K9me2 wt-kss")

mycols <- c("black", "red", "blue", "green", "darkmagenta")

pdf(paste0(plotDir, "REC8_HA_Rep1_peaks_sets_differential_wt_kss_width_dists.pdf"),
    height = 5, width = 10)
par(mfrow = c(1,1), mar = c(6, 6, 2, 2), mgp = c(4, 1.5, 0))
maxDensity <- max(max(density(width(peakSetList[[1]]))$y), max(density(width(peakSetList[[2]]))$y),
                  max(density(width(peakSetList[[3]]))$y), max(density(width(peakSetList[[4]]))$y),
                  max(density(width(peakSetList[[5]]))$y))
plot(density(width(peakSetList[[1]])), col = mycols[1], lwd = 1.5,
     xlim = c(min(width(peakSetList[[1]]), width(peakSetList[[2]]), width(peakSetList[[3]]), width(peakSetList[[4]]), width(peakSetList[[5]])),
              max(width(peakSetList[[1]]), width(peakSetList[[2]]), width(peakSetList[[3]]), width(peakSetList[[4]]), width(peakSetList[[5]]))),
     ylim = c(0, maxDensity),
     xlab = "", ylab = "", main = "", cex.lab = 1, cex.axis = 1)
abline(v = mean(width(peakSetList[[1]])), col = mycols[1], lty = 2, lwd = 1)
for(i in 2:length(peakSetList)) {
  lines(density(width(peakSetList[[i]])), col = mycols[i], lwd = 1.5) 
  abline(v = mean(width(peakSetList[[i]])), col = mycols[i], lty = 2, lwd = 1)
}
mtext(side = 1, line = 3, cex = 1, text = "Peak width (bp)")
mtext(side = 2, line = 3, cex = 1, text = "Density of REC8-HA Rep1 peaks")
legend("right",
       legend = names,
       col = mycols,
       text.col = mycols,
       ncol = 1, cex = 0.5, lwd = 1.5, bty = "n")
dev.off()



