####################################################################################
# Plot peak density genome profiles, with cumulative chromosomal coordinates       #
####################################################################################

# Usage:
# Rscript peak_density_genomeProfilesPlot.R 10kb

args <- commandArgs(trailingOnly = TRUE)
winName <- as.character(args[1])

inDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/genome_wide/peak_density_genomeProfiles/"
plotDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/genome_wide/peak_density_genomeProfiles/plots/"

names <- c("REC8_HA_Rep1_peaks",
           "REC8_HA_Rep1_peaks_mean_log2_REC8_HA_Rep1_wt_kss_diff_cov_top10percent",
           "REC8_HA_Rep1_peaks_mean_RNAseq_Rep1_kss_wt_diff_cov_top10percent",
           "REC8_HA_Rep1_peaks_mean_log2_SPO11oligo_Rep1_kss_wt_diff_cov_top10percent",
           "REC8_HA_Rep1_peaks_mean_log2_H3K9me2_wt_kss_diff_cov_top10percent")
filt_peak_density_genomeProfiles <- lapply(seq_along(names), function(x) {
  read.table(paste0(inDir, "filt_", names[x], "_density_genome_", winName, ".txt"),
             header = T)
})
filt_noNA_peak_density_genomeProfiles <- lapply(seq_along(names), function(x) {
  read.table(paste0(inDir, "filt_noNA_", names[x], "_density_genome_", winName, ".txt"),
             header = T)
})

# Redefine names for plotting
names <- c("REC8-HA Rep1",
           "REC8-HA Rep1 wt-kss",
           "RNA-seq kss-wt",
           "SPO11-1-oligos kss-wt",
           "H3K9me2 wt-kss")
mycols <- c("black", "red", "blue", "green", "darkmagenta") 

# genomic definitions
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
# pericentromeric regions are as defined in Supplemental Table S26 of Ziolkowski et al. (2017) Genes Dev. 31
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)

#################################
# make cumulative genomes       #
#################################

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

# Function to plot 5 overlaid peak density genome profiles
peakDensityGenomePlot <- function(xplot, datList) {
  plot(xplot, datList[[1]][,2], type = "l", lwd = 1.5, col = mycols[1],
       ylim = c(min(datList[[1]][,2], datList[[2]][,2], datList[[3]][,2], datList[[4]][,2], datList[[5]][,2]),
                max(datList[[1]][,2], datList[[2]][,2], datList[[3]][,2], datList[[4]][,2], datList[[5]][,2])),
       xlab = "", ylab = "", yaxt = "n", main = "")
  axis(side = 2, at = pretty(datList[[1]][,2], datList[[2]][,2], datList[[3]][,2], datList[[4]][,2], datList[[5]][,2]), lwd.tick = 1.5)
  mtext(side = 2, line = 2.25, cex = 1, text = "REC8-HA Rep1 peaks", col = mycols[1])
  for(i in 2:5) {
     lines(xplot, datList[[i]][,2], lwd = 1.5, col = mycols[i])
  }
  abline(v = sumchr, lty = 1, lwd = 0.75)
  abline(v = centromeres, lty = 5, lwd = 0.75)
  box(lwd = 1.5)
  legend("left",
         legend = names,
         col = mycols,
         text.col = mycols,
         ncol = 1, cex = 0.37, lwd = 1.5, bty = "n")
}

# Function to plot 5 overlaid peak density genome profiles
peakDensityGenomePlot2 <- function(xplot, datList) {
  plot(xplot, datList[[1]][,2], type = "l", lwd = 1.5, col = mycols[1],
       ylim = c(min(datList[[1]][,2]),
                max(datList[[1]][,2])),
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", main = "")
  axis(side = 4, at = pretty(datList[[1]][,2]), lwd.tick = 1.5)
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 1, adj = c(0.5, -4), labels = "REC8-HA Rep1 peaks", xpd = NA, srt = -90, col = mycols[1])
  par(new = T)
  plot(xplot, datList[[2]][,2], type = "l", lwd = 1.5, col = mycols[2],
       ylim = c(min(datList[[2]][,2], datList[[3]][,2], datList[[4]][,2], datList[[5]][,2]),
                max(datList[[2]][,2], datList[[3]][,2], datList[[4]][,2], datList[[5]][,2])),
       xlab = "", ylab = "", main = "")
  for(i in 3:5) {
     lines(xplot, datList[[i]][,2], lwd = 1.5, col = mycols[i])
  }
  #axis(side = 2, at = pretty(datList[[2]][,2], datList[[3]][,2], datList[[4]][,2], datList[[5]][,2]), lwd.tick = 1.5) 
  mtext(side = 2, line = 2.25, cex = 1, text = "Differential REC8-HA Rep1 peaks", col = mycols[1])
  abline(v = sumchr, lty = 1, lwd = 0.75)
  abline(v = centromeres, lty = 5, lwd = 0.75)
  box(lwd = 1.5)
  legend("left",
         legend = names,
         col = mycols,
         text.col = mycols,
         ncol = 1, cex = 0.37, lwd = 1.5, bty = "n")
}

# Plot
pdf(file = paste0(plotDir, "REC8_HA_Rep1_peak_density_genomeProfiles_wt_kss_diff_v100718.pdf"),
    height = 5, width = 7.5)
par(mfcol = c(2, 1))
par(mar = c(2.1, 4.1, 2.1, 4.1))
par(mgp = c(3, 1, 0))
peakDensityGenomePlot(xplot = filt_peak_density_genomeProfiles[[1]][,1],
                      datList = filt_peak_density_genomeProfiles)
peakDensityGenomePlot2(xplot = filt_peak_density_genomeProfiles[[1]][,1],
                       datList = filt_peak_density_genomeProfiles)
dev.off()
 
