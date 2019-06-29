#!/applications/R/R-3.5.0/bin/Rscript

# Generate genome-wide plots of chromosome profiles

# REC8-HA Rep2 (red) vs:
# Nucleosomes (purple4)

# REC8-HA Rep2 (red) vs:
# H3K9me2 (green2)
# TEs (darkgreen)

# REC8-HA Rep2 (red) vs:
# DNA methylation ("navy", "blue", "deepskyblue1")

# REC8-HA Rep2 (red) vs:
# H3K4me3 (goldenrod1)
# genes (goldenrod4)

# SPO11-1-oligos (dodgerblue2)
# Crossovers (darkblue)

# Usage:
# Rscript ./log2_ChIPinput_genomewideZscore_perwin_PLOTONLY_v190619.R 10kb

args <- commandArgs(trailingOnly = T)
winName <- as.character(args[1])

library(parallel)

# Genomic definitions
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
# Pericentromeric regions are as defined in Supplemental Table S26
# of Ziolkowski et al. (2017) Genes Dev. 31
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)

# Make cumulative genomes
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

plotDir <- "./plots/"
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

REC8HADir <- "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/genomeProfiles/genomewideZscore/"
REC8HAprefix <- "REC8_HA_Rep2_ChIP_REC8_MYC_Rep1_input"
REC8HAname <- "REC8-HA"
REC8HAcolour <- "red"

inDir1 <- "/projects/ajt200/BAM_masters/nucleosomes/WT/coverage/nakedDNA_untrimmed_input/log2ChIPinput/genomeProfiles/genomewideZscore/"
prefix1 <- "WT_nuc_nakedDNAuntrimmed"
name1 <- "Nucleosomes"
colour1 <- "purple4"

inDir2 <- "/projects/ajt200/BAM_masters/H3K9me2/WT/coverage/log2ChIPinput/genomeProfiles/genomewideZscore/"
prefix2 <- "WT_H3K9me2_ChIP_WT_H3K9me2_input"
name2 <- "H3K9me2"
colour2 <- "green2"

inDir3 <- "/projects/ajt200/BAM_masters/H3K4me3/replicates/coverage/log2ChIPinput/genomeProfiles/genomewideZscore/"
prefix3 <- "WT_H3K4me3_ChIP14_WT_H3K9me2_input"
name3 <- "H3K4me3"
colour3 <- "goldenrod1"

inDir4 <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/genomeProfiles/genomewideZscore/"
prefix4 <- "WT_SPO11oligo_RPI1_WT_nakedDNA_R1"
name4 <- "SPO11-1"
colour4 <- "dodgerblue2"

# REC8
inDirREC8 <- c(REC8HADir)
prefixREC8 <- c(REC8HAprefix)
nameREC8 <- c(REC8HAname)
colourREC8 <- c(REC8HAcolour)

filt_REC8log2trans_list <- mclapply(seq_along(nameREC8), function(x) {
  read.table(paste0(inDirREC8[x], "filt_log2_", prefixREC8[x],
                    "_genome_norm_coverage_Zscore_", winName, ".txt"))
}, mc.cores = length(nameREC8))

# All
inDirAll <- c(inDir1, inDir2, inDir3, inDir4)
prefixAll <- c(prefix1, prefix2, prefix3, prefix4)
nameAll <- c(name1, name2, name3, name4)
colourAll <- c(colour1, colour2, colour3, colour4)

filt_log2trans_list <- mclapply(seq_along(nameAll), function(x) {
  read.table(paste0(inDirAll[x], "filt_log2_", prefixAll[x],
                    "_genome_norm_coverage_Zscore_", winName, ".txt"))
}, mc.cores = length(nameAll))

# Feature frequency
inDirFeatures <- "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/genomeProfiles/genomewideZscore/"

# TEs
filt_TEs <- read.table(paste0(inDirFeatures, "filt_TE_frequency_genome_",
                              winName, ".txt"))
# genes
filt_genes <- read.table(paste0(inDirFeatures, "filt_gene_frequency_genome_",
                                winName, ".txt"))
# COs
filt_COs <- read.table(paste0(inDirFeatures, "filt_CO_frequency_genome_",
                              winName, ".txt"))
# DNA methylation in 200-kb windows
filt_meth <- read.table(paste0(inDirFeatures, "filt_DNAmeth_GSM980986_WT_rep2_genome_200kb.txt"))


# Function to plot REC8 genome-scale coverage overlaid with other datasets
dat1to2diffY <- function(xplot,
                         dat1, dat2,
                         dat1Lab, dat2Lab,
                         dat1Colour, dat2Colour) {
  plot(xplot, dat1, type = "l", lwd = 2, col = dat1Colour,
       ylim = c(min(dat1),
                max(dat1)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  axis(side = 2, lwd.tick = 2, cex.axis = 2.00)
  mtext(side = 2, line = 3.50, cex = 2,
        text = dat1Lab, col = dat1Colour)
  par(new = T)
  plot(xplot, dat2, type = "l", lwd = 2, col = dat2Colour,
       ylim = c(min(dat2),
                max(dat2)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 3.00, adj = c(0.5, -2.0), labels = dat2Lab, xpd = NA, srt = -90,
       col = dat2Colour)
  axis(side = 4, at = pretty(dat2), lwd.tick = 2, cex.axis = 2.00)
  axis(side = 1, lwd.tick = 2,
       at = c(0, 2e7, 4e7, 6e7, 8e7, 1e8, 1.2e8),
       labels = c("0", "20", "40", "60", "80", "100", "120"), cex.axis = 2.00)
  mtext(side = 1, line = 3.50, cex = 2,
        text = "Coordinates (Mb)")
  abline(v = sumchr, lty = 1, lwd = 2)
  abline(v = centromeres, lty = 5, lwd = 2)
  rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 1, lwd = 1, col = "slategray1")
  rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 3, lwd = 1, col = "slategray1")
  box(lwd = 2)
}

# Function to plot REC8 genome-scale coverage overlaid with other datasets 
dat1to2GenomePlot <- function(xplot,
                              dat1, dat2,
                              legendNames, plotColours) {
  plot(xplot, dat2, type = "l", lwd = 2, col = plotColours[2],
       ylim = c(min(dat1, dat2),
                max(dat1, dat2)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  lines(xplot, dat1, type = "l", lwd = 2, col = plotColours[1])
  axis(side = 1, lwd.tick = 2,
       at = c(0, 2e7, 4e7, 6e7, 8e7, 1e8, 1.2e8),
       labels = c("0", "20", "40", "60", "80", "100", "120"), cex.axis = 2.00)
  mtext(side = 1, line = 3.50, cex = 2,
        text = "Coordinates (Mb)")
  axis(side = 2, lwd.tick = 2, cex.axis = 2.00)
  mtext(side = 2, line = 3.00, cex = 2,
        text = expression("Log"[2]*"(ChIP/input)"))
  abline(v = sumchr, lty = 1, lwd = 2)
  abline(v = centromeres, lty = 5, lwd = 2)
  rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 1, lwd = 1, col = "slategray1")
  rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 3, lwd = 1, col = "slategray1")
  box(lwd = 2)
  legend("topleft",
         legend = legendNames,
         col = "white",
         text.col = plotColours,
         ncol = 1, cex = 1.2, lwd = 1.5, bty = "n")
}

# Function to plot REC8 genome-scale coverage overlaid with other datasets
dat1to2vFeatureGenomePlot <- function(xplot,
                                      dat1, dat2, feature1,
                                      legendNames, legendPos, feature1Lab,
                                      plotColours, feature1Colour) {
  plot(xplot, feature1, type = "l", lwd = 2, col = feature1Colour,
       ylim = c(0, max(feature1)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 3.00, adj = c(0.5, -2.0), labels = feature1Lab, xpd = NA, srt = -90,
       col = feature1Colour)
  axis(side = 4, at = pretty(feature1), lwd.tick = 2, cex.axis = 2.00)
  par(new = T)
  plot(xplot, dat2, type = "l", lwd = 2, col = plotColours[2],
       ylim = c(min(dat1, dat2),
                max(dat1, dat2)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  lines(xplot, dat1, type = "l", lwd = 2, col = plotColours[1])
  axis(side = 1, lwd.tick = 2,
       at = c(0, 2e7, 4e7, 6e7, 8e7, 1e8, 1.2e8),
       labels = c("0", "20", "40", "60", "80", "100", "120"), cex.axis = 2.00)
  mtext(side = 1, line = 3.50, cex = 2,
        text = "Coordinates (Mb)")
  axis(side = 2, lwd.tick = 2, cex.axis = 2.00)
  mtext(side = 2, line = 3.00, cex = 2,
        text = expression("Log"[2]*"(ChIP/input)"))
  abline(v = sumchr, lty = 1, lwd = 2)
  abline(v = centromeres, lty = 5, lwd = 2)
  rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 1, lwd = 1, col = "slategray1")
  rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 3, lwd = 1, col = "slategray1")
  box(lwd = 2)
  legend(legendPos,
         legend = legendNames,
         col = "white",
         text.col = plotColours,
         ncol = 1, cex = 1.2, lwd = 1.5, bty = "n")
}

# Function to plot REC8 genome-scale coverage overlaid with other datasets
dat1vFeatureGenomePlot <- function(xplot,
                                   dat1, feature1,
                                   legendNames, feature1Lab,
                                   plotColours, feature1Colour) {
  plot(xplot, feature1, type = "l", lwd = 2, col = feature1Colour,
       ylim = c(0, max(feature1)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 3.00, adj = c(0.5, -2.0), labels = feature1Lab, xpd = NA, srt = -90,
       col = feature1Colour)
  axis(side = 4, at = pretty(feature1), lwd.tick = 2, cex.axis = 2.00)
  par(new = T)
  plot(xplot, dat1, type = "l", lwd = 2, col = plotColours[1],
       ylim = c(min(dat1),
                max(dat1)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  axis(side = 1, lwd.tick = 2,
       at = c(0, 2e7, 4e7, 6e7, 8e7, 1e8, 1.2e8),
       labels = c("0", "20", "40", "60", "80", "100", "120"), cex.axis = 2.00)
  mtext(side = 1, line = 3.50, cex = 2,
        text = "Coordinates (Mb)")
  axis(side = 2, lwd.tick = 2, cex.axis = 2.00)
  mtext(side = 2, line = 3.00, cex = 2,
        text = expression("Log"[2]*"(ChIP/input)"))
  abline(v = sumchr, lty = 1, lwd = 2)
  abline(v = centromeres, lty = 5, lwd = 2)
  rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 1, lwd = 1, col = "slategray1")
  rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 3, lwd = 1, col = "slategray1")
  box(lwd = 2)
  legend("topleft",
         legend = legendNames,
         col = "white",
         text.col = plotColours,
         ncol = 1, cex = 1.2, lwd = 1.5, bty = "n")
}

# Function to plot REC8 genome-scale coverage overlaid with other datasets 
dat1to4GenomePlot <- function(xplot,
                              dat1, dat2, dat3, dat4,
                              legendNames, plotColours) {
  plot(xplot, dat4, type = "l", lwd = 2, col = plotColours[4],
       ylim = c(min(dat1, dat2, dat3, dat4),
                max(dat1, dat2, dat3, dat4)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  lines(xplot, dat3, type = "l", lwd = 2, col = plotColours[3])
  lines(xplot, dat2, type = "l", lwd = 2, col = plotColours[2])
  lines(xplot, dat1, type = "l", lwd = 2, col = plotColours[1])
  axis(side = 1, lwd.tick = 2,
       at = c(0, 2e7, 4e7, 6e7, 8e7, 1e8, 1.2e8),
       labels = c("0", "20", "40", "60", "80", "100", "120"), cex.axis = 2.00)
  mtext(side = 1, line = 3.50, cex = 2,
        text = "Coordinates (Mb)")
  axis(side = 2, lwd.tick = 2, cex.axis = 2.00)
  mtext(side = 2, line = 3.00, cex = 2,
        text = expression("Log"[2]*"(ChIP/input)"))
  abline(v = sumchr, lty = 1, lwd = 2)
  abline(v = centromeres, lty = 5, lwd = 2)
  rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 1, lwd = 1, col = "slategray1")
  rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 3, lwd = 1, col = "slategray1")
  box(lwd = 2)
  legend("topleft",
         legend = legendNames,
         col = "white",
         text.col = plotColours,
         ncol = 1, cex = 1.2, lwd = 1.5, bty = "n")
}

# Function to plot REC8 genome-scale coverage overlaid with DNA methylation (each context)
REC8vsDNAmethSepGenomePlot <- function(xplot1,
                                       xplot2,
                                       dat1A,
                                       dat2,
                                       dat1ALab,
                                       legendNames, plotColours) {
  plot(xplot1, dat1A, type = "l", lwd = 2, col = plotColours[1],
       ylim = c(min(dat1A), max(dat1A)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  axis(side = 2, lwd.tick = 2, cex.axis = 2.00)
  mtext(side = 2, line = 3.50, cex = 2,
        #text = expression("Log"[2]*"(ChIP/input)"))
        text = dat1ALab, col = plotColours[1])
  par(new = T)
  plot(xplot2, dat2[,3], type = "l", lwd = 2, col = plotColours[2],
       ylim = c(min(dat2[,3], dat2[,4], dat2[,5]), max(dat2[,3], dat2[,4], dat2[,5])),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 3.00, adj = c(0.5, -2.0), labels = "DNA methylation", xpd = NA, srt = -90,
       col = "blue")
  lines(xplot2, dat2[,4], type = "l", lwd = 2, col = plotColours[3])
  lines(xplot2, dat2[,5], type = "l", lwd = 2, col = plotColours[4])
  axis(side = 4, at = pretty(c(dat2[,3], dat2[,4], dat2[,5])), lwd.tick = 2, cex.axis = 2.00)
  axis(side = 1, lwd.tick = 2,
       at = c(0, 2e7, 4e7, 6e7, 8e7, 1e8, 1.2e8),
       labels = c("0", "20", "40", "60", "80", "100", "120"), cex.axis = 2.00)
  mtext(side = 1, line = 3.50, cex = 2,
        text = "Coordinates (Mb)")
  abline(v = sumchr, lty = 1, lwd = 2)
  abline(v = centromeres, lty = 5, lwd = 2)
  rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 1, lwd = 1, col = "slategray1")
  rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 3, lwd = 1, col = "slategray1")
  box(lwd = 2)
  legend("topleft",
         legend = legendNames,
         col = "white",
         text.col = plotColours[2:4],
         ncol = 1, cex = 1.2, lwd = 1.5, bty = "n")
}


#
pdf(paste0(plotDir, "genomewideZscore_log2_", REC8HAprefix, "_", winName, "_genomeProfiles_v190619.pdf"),
    height = 40, width = 12)
par(mfcol = c(10, 1))
par(mar = c(5.1, 6.1, 2.1, 6.1))
par(mgp = c(3, 1.5, 0))
dat1to2diffY(xplot = filt_REC8log2trans_list[[1]]$cumWindows,
             dat1 = filt_REC8log2trans_list[[1]]$filt_ZscoreLog2cov,
             dat2 = filt_log2trans_list[[1]]$filt_ZscoreLog2cov,
             dat1Lab = nameREC8,
             dat2Lab = nameAll[1],
             dat1Colour = colourREC8,
             dat2Colour = colourAll[1])
dat1to2diffY(xplot = filt_REC8log2trans_list[[1]]$cumWindows,
             dat1 = filt_REC8log2trans_list[[1]]$filt_ZscoreLog2cov,
             dat2 = filt_log2trans_list[[2]]$filt_ZscoreLog2cov,
             dat1Lab = nameREC8,
             dat2Lab = nameAll[2],
             dat1Colour = colourREC8,
             dat2Colour = colourAll[2])
dat1to2diffY(xplot = filt_REC8log2trans_list[[1]]$cumWindows,
             dat1 = filt_REC8log2trans_list[[1]]$filt_ZscoreLog2cov,
             dat2 = filt_TEs$filt_TEs,
             dat1Lab = nameREC8,
             dat2Lab = "TEs",
             dat1Colour = colourREC8,
             dat2Colour = "darkgreen")
REC8vsDNAmethSepGenomePlot(xplot1 = filt_REC8log2trans_list[[1]]$cumWindows,
                           xplot2 = filt_meth$cumWindows,
                           dat1A = filt_REC8log2trans_list[[1]]$filt_ZscoreLog2cov,
                           dat2 = filt_meth,
                           dat1ALab = nameREC8,
                           legendNames = c("mCG", "mCHG", "mCHH"),
                           plotColours = c(colourREC8[1], "navy", "blue", "deepskyblue1"))
dat1to2diffY(xplot = filt_REC8log2trans_list[[1]]$cumWindows,
             dat1 = filt_REC8log2trans_list[[1]]$filt_ZscoreLog2cov,
             dat2 = filt_log2trans_list[[3]]$filt_ZscoreLog2cov,
             dat1Lab = nameREC8,
             dat2Lab = nameAll[3],
             dat1Colour = colourREC8,
             dat2Colour = colourAll[3])
dat1to2diffY(xplot = filt_REC8log2trans_list[[1]]$cumWindows,
             dat1 = filt_REC8log2trans_list[[1]]$filt_ZscoreLog2cov,
             dat2 = filt_genes$filt_genes,
             dat1Lab = nameREC8,
             dat2Lab = "Genes",
             dat1Colour = colourREC8,
             dat2Colour = "goldenrod4")
dat1to2diffY(xplot = filt_REC8log2trans_list[[1]]$cumWindows,
             dat1 = filt_REC8log2trans_list[[1]]$filt_ZscoreLog2cov,
             dat2 = filt_log2trans_list[[4]]$filt_ZscoreLog2cov,
             dat1Lab = nameREC8,
             dat2Lab = nameAll[4],
             dat1Colour = colourREC8,
             dat2Colour = colourAll[4])
dat1to2diffY(xplot = filt_REC8log2trans_list[[1]]$cumWindows,
             dat1 = filt_REC8log2trans_list[[1]]$filt_ZscoreLog2cov,
             dat2 = filt_COs$filt_COs,
             dat1Lab = nameREC8,
             dat2Lab = "Crossovers",
             dat1Colour = colourREC8,
             dat2Colour = "darkblue")
dat1to2diffY(xplot = filt_log2trans_list[[3]]$cumWindows,
             dat1 = filt_log2trans_list[[3]]$filt_ZscoreLog2cov,
             dat2 = filt_log2trans_list[[4]]$filt_ZscoreLog2cov,
             dat1Lab = nameAll[3],
             dat2Lab = nameAll[4],
             dat1Colour = colourAll[3],
             dat2Colour = colourAll[4])
dat1to2diffY(xplot = filt_REC8log2trans_list[[1]]$cumWindows,
             dat1 = filt_genes$filt_genes,
             dat2 = filt_COs$filt_COs,
             dat1Lab = "Genes",
             dat2Lab = "Crossovers",
             dat1Colour = "goldenrod4",
             dat2Colour = "darkblue")
dev.off()
