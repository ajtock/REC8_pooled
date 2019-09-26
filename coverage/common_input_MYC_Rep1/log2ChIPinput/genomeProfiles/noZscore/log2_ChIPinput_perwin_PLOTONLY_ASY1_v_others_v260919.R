#!/applications/R/R-3.5.0/bin/Rscript

#########################################
#########################################
# Generate plots of chromosome profiles
#########################################
#########################################

# Usage:
# ./log2_ChIPinput_perwin_PLOTONLY_ASY1_v_others_v260919.R 10kb

#winName <- "10kb"

args <- commandArgs(trailingOnly = T)
winName <- as.character(args[1])

libNames <- c(
              "ASY1_Rep2_ChIP_REC8_MYC_Rep1_input",
              "REC8_HA_Rep2_ChIP_REC8_MYC_Rep1_input",
              "WT_nuc_WT_nakedDNA_untrimmed",
              "WT_H3K9me2_ChIP_WT_H3K9me2_input",
              "WT_H3K4me3_ChIP14_WT_H3K9me2_input",
              "WT_SPO11oligo_RPI1_WT_nakedDNA_R1",
              "MTOPVIB_HA_Rep1_ChIP_REC8_MYC_Rep1_input",
              "MTOPVIB_HA_Rep2_ChIP_REC8_MYC_Rep1_input",
              "WT_SPO11_ChIP4_REC8_MYC_Rep1_input",
              "WT_SPO11_ChIP13_REC8_MYC_Rep1_input",
              "WT_H3K4me1_Rep1_ChIP_WT_H3K9me2_input",
              "WT_H3K4me2_Rep1_ChIP_WT_H3K9me2_input",
              "WT_H3K27me1_Rep1_ChIP_WT_H3K9me2_input",
              "H3K27me3_ChIP_SRR1509478_WT_H3K9me2_input",
              "H2A_input",
              "H2AW_input",
              "H2AX_input",
              "H2AZ_input"
             )
libNamesPlot <- c(
                  "ASY1",
                  "REC8-HA",
                  "Nucleosomes",
                  "H3K9me2",
                  "H3K4me3",
                  "SPO11-1-oligos",
                  "MTOPVIB Rep1",
                  "MTOPVIB Rep2",
                  "SPO11-1 ChIP Rep1",
                  "SPO11-1 ChIP Rep2",
                  "H3K4me1",
                  "H3K4me2",
                  "H3K27me1",
                  "H3K27me3",
                  "H2A",
                  "H2A.W",
                  "H2A.X",
                  "H2A.Z"
                 )
libColours <- c(
                "red",
                "green2",
                "purple4",
                "magenta3",
                "goldenrod1",
                "dodgerblue2",
                "darkcyan",
                "cyan",
                "grey40",
                "grey40",
                "goldenrod2",
                "goldenrod3",
                "purple",
                "navy",
                "grey40",
                "green",
                "blue",
                "orange"
               )

library(parallel)

# Genomic definitions
chrs <- c("Chr2", "Chr2", "Chr3", "Chr4", "Chr5")
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
# Pericentromeric regions are as defined in Supplemental Table S26
# of Ziolkowski et al. (2017) Genes Dev. 31
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)

# Make chromosomal coordinates cumulative
# such that the first coordinate of Chr2 is
# equal to the last coordinate of Chr1 + 1
sumchr <- cumsum(c(0, chrLens))
print(sumchr)
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

plotDir <- "plots/"
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

lib_filt_winDF_list <- mclapply(seq_along(libNames), function(x) {
  read.table(paste0("filt_", libNames[x],
                    "_genome_norm_coverage_",
                    winName, "_noZscore.tsv"),
             header = T)
}, mc.cores = length(libNames))

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

# Function to plot genome-scale coverage overlaid with other datasets
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
  text(p[2], mean(p[3:4]), cex = 3.00, adj = c(0.5, -1.75), labels = dat2Lab, xpd = NA, srt = -90,
       col = dat2Colour)
  axis(side = 4, at = pretty(dat2), lwd.tick = 2, cex.axis = 2.00)
  axis(side = 1, lwd.tick = 2,
       at = c(0, 2e7, 4e7, 6e7, 8e7, 1e8, 1.2e8),
       labels = c("0", "20", "40", "60", "80", "100", "120"), cex.axis = 2.00)
  mtext(side = 1, line = 3.75, cex = 2,
        text = "Coordinates (Mb)")
  abline(v = sumchr, lty = 1, lwd = 2)
  abline(v = centromeres, lty = 5, lwd = 2)
  rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 1, lwd = 1, col = "slategray1")
  rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 3, lwd = 1, col = "slategray1")
  box(lwd = 2)
}

# Function to plot genome-scale coverage overlaid with DNA methylation (each context)
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
  text(p[2], mean(p[3:4]), cex = 3.00, adj = c(0.5, -1.75), labels = "DNA methylation", xpd = NA, srt = -90,
       col = "blue")
  lines(xplot2, dat2[,4], type = "l", lwd = 2, col = plotColours[3])
  lines(xplot2, dat2[,5], type = "l", lwd = 2, col = plotColours[4])
  axis(side = 4, at = pretty(c(dat2[,3], dat2[,4], dat2[,5])), lwd.tick = 2, cex.axis = 2.00)
  axis(side = 1, lwd.tick = 2,
       at = c(0, 2e7, 4e7, 6e7, 8e7, 1e8, 1.2e8),
       labels = c("0", "20", "40", "60", "80", "100", "120"), cex.axis = 2.00)
  mtext(side = 1, line = 3.75, cex = 2,
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

# Function to plot genome-scale coverage overlaid with other datasets
dat1to2vFeatureGenomePlot <- function(xplot,
                                      dat1, dat2, feature1,
                                      legendNames, legendPos, feature1Lab,
                                      plotColours, feature1Colour) {
  plot(xplot, feature1, type = "l", lwd = 2, col = feature1Colour,
       ylim = c(0, max(feature1)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 3.00, adj = c(0.5, -1.75), labels = feature1Lab, xpd = NA, srt = -90,
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
  mtext(side = 1, line = 3.75, cex = 2,
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

# Function to plot genome-scale coverage overlaid with other datasets
dat1vFeatureGenomePlot <- function(xplot,
                                   dat1, feature1,
                                   legendNames, feature1Lab,
                                   plotColours, feature1Colour) {
  plot(xplot, feature1, type = "l", lwd = 2, col = feature1Colour,
       ylim = c(0, max(feature1)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 3.00, adj = c(0.5, -1.75), labels = feature1Lab, xpd = NA, srt = -90,
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
  mtext(side = 1, line = 3.75, cex = 2,
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

# Function to plot genome-scale coverage overlaid with other datasets 
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
  mtext(side = 1, line = 3.75, cex = 2,
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

# Function to plot genome-scale coverage overlaid with other datasets 
dat1to3GenomePlot <- function(xplot,
                              dat1, dat2, dat3,
                              legendNames, plotColours) {
  plot(xplot, dat3, type = "l", lwd = 2, col = plotColours[3],
       ylim = c(min(dat1, dat2, dat3),
                max(dat1, dat2, dat3)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  lines(xplot, dat2, type = "l", lwd = 2, col = plotColours[2])
  lines(xplot, dat1, type = "l", lwd = 2, col = plotColours[1])
  axis(side = 1, lwd.tick = 2,
       at = c(0, 2e7, 4e7, 6e7, 8e7, 1e8, 1.2e8),
       labels = c("0", "20", "40", "60", "80", "100", "120"), cex.axis = 2.00)
  mtext(side = 1, line = 3.75, cex = 2,
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

# Function to plot genome-scale coverage overlaid with other datasets 
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
  mtext(side = 1, line = 3.75, cex = 2,
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

# Plot
pdf(paste0(plotDir, "ASY1_Rep2_vs_others_genomeProfiles_noZscore_v260919.pdf"),
    height = 116, width = 12)
par(mfcol = c(29, 1))
par(mar = c(5.1, 6.1, 2.1, 6.1))
par(mgp = c(3, 1.5, 0))
dat1to3GenomePlot(xplot = lib_filt_winDF_list[[1]]$cumwindow,
                  dat1 = lib_filt_winDF_list[[1]]$filt_log2,
                  dat2 = lib_filt_winDF_list[[2]]$filt_log2,
                  dat3 = lib_filt_winDF_list[[6]]$filt_log2,
                  legendNames = c(libNamesPlot[1], libNamesPlot[2], libNamesPlot[6]),
                  plotColours = c(libColours[1], libColours[2], libColours[6]))
# MTOBVIB Rep1 v Rep2
dat1to2GenomePlot(xplot = lib_filt_winDF_list[[7]]$cumwindow,
                  dat1 = lib_filt_winDF_list[[7]]$filt_log2,
                  dat2 = lib_filt_winDF_list[[8]]$filt_log2,
                  legendNames = c(libNamesPlot[7], libNamesPlot[8]),
                  plotColours = c(libColours[7], libColours[8]))
# MTOPVIB Rep1 vs SPO11-1-oligos
dat1to2diffY(xplot = lib_filt_winDF_list[[6]]$cumwindow,
             dat1 = lib_filt_winDF_list[[7]]$filt_log2,
             dat2 = lib_filt_winDF_list[[6]]$filt_log2,
             dat1Lab = libNamesPlot[7],
             dat2Lab = libNamesPlot[6],
             dat1Colour = libColours[7],
             dat2Colour = libColours[6])
# MTOPVIB Rep2 vs SPO11-1-oligos
dat1to2diffY(xplot = lib_filt_winDF_list[[6]]$cumwindow,
             dat1 = lib_filt_winDF_list[[8]]$filt_log2,
             dat2 = lib_filt_winDF_list[[6]]$filt_log2,
             dat1Lab = libNamesPlot[8],
             dat2Lab = libNamesPlot[6],
             dat1Colour = libColours[8],
             dat2Colour = libColours[6])
# MTOPVIB Rep1 vs ASY1
dat1to2diffY(xplot = lib_filt_winDF_list[[1]]$cumwindow,
             dat1 = lib_filt_winDF_list[[7]]$filt_log2,
             dat2 = lib_filt_winDF_list[[1]]$filt_log2,
             dat1Lab = libNamesPlot[7],
             dat2Lab = libNamesPlot[1],
             dat1Colour = libColours[7],
             dat2Colour = libColours[1])
# MTOPVIB Rep2 vs ASY1
dat1to2diffY(xplot = lib_filt_winDF_list[[1]]$cumwindow,
             dat1 = lib_filt_winDF_list[[8]]$filt_log2,
             dat2 = lib_filt_winDF_list[[1]]$filt_log2,
             dat1Lab = libNamesPlot[8],
             dat2Lab = libNamesPlot[1],
             dat1Colour = libColours[8],
             dat2Colour = libColours[1])
# MTOPVIB Rep1 vs REC8-HA
dat1to2diffY(xplot = lib_filt_winDF_list[[2]]$cumwindow,
             dat1 = lib_filt_winDF_list[[7]]$filt_log2,
             dat2 = lib_filt_winDF_list[[2]]$filt_log2,
             dat1Lab = libNamesPlot[7],
             dat2Lab = libNamesPlot[2],
             dat1Colour = libColours[7],
             dat2Colour = libColours[2])
# MTOPVIB Rep2 vs REC8-HA
dat1to2diffY(xplot = lib_filt_winDF_list[[2]]$cumwindow,
             dat1 = lib_filt_winDF_list[[8]]$filt_log2,
             dat2 = lib_filt_winDF_list[[2]]$filt_log2,
             dat1Lab = libNamesPlot[8],
             dat2Lab = libNamesPlot[2],
             dat1Colour = libColours[8],
             dat2Colour = libColours[2])
# REC8
dat1to2diffY(xplot = lib_filt_winDF_list[[1]]$cumwindow,
             dat1 = lib_filt_winDF_list[[1]]$filt_log2,
             dat2 = lib_filt_winDF_list[[2]]$filt_log2,
             dat1Lab = libNamesPlot[1],
             dat2Lab = libNamesPlot[2],
             dat1Colour = libColours[1],
             dat2Colour = libColours[2])
# nucleosomes
dat1to2diffY(xplot = lib_filt_winDF_list[[1]]$cumwindow,
             dat1 = lib_filt_winDF_list[[1]]$filt_log2,
             dat2 = lib_filt_winDF_list[[3]]$filt_log2,
             dat1Lab = libNamesPlot[1],
             dat2Lab = libNamesPlot[3],
             dat1Colour = libColours[1],
             dat2Colour = libColours[3])
# H3K9me2
dat1to2diffY(xplot = lib_filt_winDF_list[[1]]$cumwindow,
             dat1 = lib_filt_winDF_list[[1]]$filt_log2,
             dat2 = lib_filt_winDF_list[[4]]$filt_log2,
             dat1Lab = libNamesPlot[1],
             dat2Lab = libNamesPlot[4],
             dat1Colour = libColours[1],
             dat2Colour = libColours[4])
# TEs
dat1to2diffY(xplot = lib_filt_winDF_list[[1]]$cumwindow,
             dat1 = lib_filt_winDF_list[[1]]$filt_log2,
             dat2 = filt_TEs$filt_TEs,
             dat1Lab = libNamesPlot[1],
             dat2Lab = "TEs",
             dat1Colour = libColours[1],
             dat2Colour = "darkgreen")
# DNA methylation
REC8vsDNAmethSepGenomePlot(xplot1 = lib_filt_winDF_list[[1]]$cumwindow,
                           xplot2 = filt_meth$cumWindows,
                           dat1A = lib_filt_winDF_list[[1]]$filt_log2,
                           dat2 = filt_meth,
                           dat1ALab = libNamesPlot[1],
                           legendNames = c("mCG", "mCHG", "mCHH"),
                           plotColours = c(libColours[1], "navy", "blue", "deepskyblue1"))
# H3K4me3
dat1to2diffY(xplot = lib_filt_winDF_list[[1]]$cumwindow,
             dat1 = lib_filt_winDF_list[[1]]$filt_log2,
             dat2 = lib_filt_winDF_list[[5]]$filt_log2,
             dat1Lab = libNamesPlot[1],
             dat2Lab = libNamesPlot[5],
             dat1Colour = libColours[1],
             dat2Colour = libColours[5])
# Genes
dat1to2diffY(xplot = lib_filt_winDF_list[[1]]$cumwindow,
             dat1 = lib_filt_winDF_list[[1]]$filt_log2,
             dat2 = filt_genes$filt_genes,
             dat1Lab = libNamesPlot[1],
             dat2Lab = "Genes",
             dat1Colour = libColours[1],
             dat2Colour = "goldenrod4")
# SPO11-1-oligos
dat1to2diffY(xplot = lib_filt_winDF_list[[1]]$cumwindow,
             dat1 = lib_filt_winDF_list[[1]]$filt_log2,
             dat2 = lib_filt_winDF_list[[6]]$filt_log2,
             dat1Lab = libNamesPlot[1],
             dat2Lab = libNamesPlot[6],
             dat1Colour = libColours[1],
             dat2Colour = libColours[6])
# Crossovers
dat1to2diffY(xplot = lib_filt_winDF_list[[1]]$cumwindow,
             dat1 = lib_filt_winDF_list[[1]]$filt_log2,
             dat2 = filt_COs$filt_COs,
             dat1Lab = libNamesPlot[1],
             dat2Lab = "Crossovers",
             dat1Colour = libColours[1],
             dat2Colour = "darkblue")
# MTOPVIB, SPO11-1 ChIP, and the rest
for(x in 7:length(libNames)) {
  dat1to2diffY(xplot = lib_filt_winDF_list[[1]]$cumwindow,
               dat1 = lib_filt_winDF_list[[1]]$filt_log2,
               dat2 = lib_filt_winDF_list[[x]]$filt_log2,
               dat1Lab = libNamesPlot[1],
               dat2Lab = libNamesPlot[x],
               dat1Colour = libColours[1],
               dat2Colour = libColours[x])
}
dev.off()
