#!/applications/R/R-3.5.0/bin/Rscript

######################################
######################################
# Generate plots chromosome profiles
######################################
######################################

# REC8_HA_Rep2 (red)
# kss_REC8_HA_Rep1 (red4)

# Usage:
# Rscript ./log2_ChIPinput_perwin_PLOTONLY_v040919.R 10kb 'REC8_HA_Rep2,REC8_HA_Rep1,REC8_MYC_Rep1' 'kss_REC8_HA_Rep1,kss_REC8_HA_Rep2' 'wt REC8-HA_Rep2,wt REC8-HA_Rep1,wt REC8-Myc_Rep1' 'kss REC8-HA_Rep1,kss REC8-HA_Rep2' '

winName <- "10kb"
geno1Names <- unlist(strsplit("REC8_HA_Rep2,REC8_HA_Rep1,REC8_MYC_Rep1",
                              split = ","))
geno2Names <- unlist(strsplit("kss_REC8_HA_Rep1,kss_REC8_HA_Rep2",
                              split = ","))
geno1NamesPlot <- unlist(strsplit("wt REC8-HA_Rep2,wt REC8-HA_Rep1,wt REC8-Myc_Rep1",
                                  split = ","))
geno2NamesPlot <- unlist(strsplit("kss REC8-HA_Rep1,kss REC8-HA_Rep2",
                                  split = ","))
genoColours <- unlist(strsplit("red,red4",
                               split = ","))

args <- commandArgs(trailingOnly = T)
winName <- as.character(args[1])
geno1Names <- unlist(strsplit(args[2],
                              split = ","))
geno2Names <- unlist(strsplit(args[3],
                              split = ","))
geno1NamesPlot <- unlist(strsplit(args[4],
                                  split = ","))
geno2NamesPlot <- unlist(strsplit(args[5],
                                  split = ","))
genoColours <- unlist(strsplit(args[6],
                               split = ","))

library(parallel)

# Genomic definitions
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
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

geno1_filt_winDF_list <- mclapply(seq_along(geno1Names), function(x) {
  read.table(paste0("filt_", geno1Names[x],
                    "_ChIP_REC8_MYC_Rep1_input_genome_norm_coverage_",
                    winName, "_noZscore.tsv"),
             header = T)
}, mc.cores = length(geno1Names))

geno2_filt_winDF_list <- mclapply(seq_along(geno2Names), function(x) {
  read.table(paste0("filt_", geno2Names[x],
                    "_ChIP_kss_REC8_HA_Rep1_input_genome_norm_coverage_",
                    winName, "_noZscore.tsv"),
             header = T)
}, mc.cores = length(geno2Names))

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

# Function to plot genome-scale coverage overlaid with other datasets 
dat1to2GenomePlot <- function(xplot,
                              dat1, dat2,
                              Ylab,
                              legendNames,
                              plotColours) {
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
        text = Ylab)
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


# Plot
pdf(paste0(plotDir, "wt_v_kss_REC8_ChIP_and_log2ChIPinput_genomeProfiles_noZscore_v040919.pdf"),
    height = 24, width = 24)
par(mfcol = c(6, 2))
par(mar = c(5.1, 6.1, 2.1, 6.1))
par(mgp = c(3, 1.5, 0))
# ChIP
dat1to2GenomePlot(xplot = geno1_filt_winDF_list[[1]]$cumwindow,
                  dat1 = geno1_filt_winDF_list[[1]][,4],
                  dat2 = geno2_filt_winDF_list[[1]][,4],
                  Ylab = expression("Normalized ChIP"),
                  legendNames = c(bquote(.(geno1NamesPlot[1])),
                                  as.expression(bquote(italic(.(unlist(strsplit(geno2NamesPlot[1],
                                                                                split = " "))[1])) ~
                                                       .(unlist(strsplit(geno2NamesPlot[1],
                                                                         split = " "))[2])))),
                  plotColours = genoColours)
dat1to2GenomePlot(xplot = geno1_filt_winDF_list[[1]]$cumwindow,
                  dat1 = geno1_filt_winDF_list[[1]][,4],
                  dat2 = geno2_filt_winDF_list[[2]][,4],
                  Ylab = expression("Normalized ChIP"),
                  legendNames = c(bquote(.(geno1NamesPlot[1])),
                                  as.expression(bquote(italic(.(unlist(strsplit(geno2NamesPlot[2],
                                                                                split = " "))[1])) ~
                                                       .(unlist(strsplit(geno2NamesPlot[2],
                                                                         split = " "))[2])))),
                  plotColours = genoColours)
dat1to2GenomePlot(xplot = geno1_filt_winDF_list[[2]]$cumwindow,
                  dat1 = geno1_filt_winDF_list[[2]][,4],
                  dat2 = geno2_filt_winDF_list[[1]][,4],
                  Ylab = expression("Normalized ChIP"),
                  legendNames = c(bquote(.(geno1NamesPlot[2])),
                                  as.expression(bquote(italic(.(unlist(strsplit(geno2NamesPlot[1],
                                                                                split = " "))[1])) ~
                                                       .(unlist(strsplit(geno2NamesPlot[1],
                                                                         split = " "))[2])))),
                  plotColours = genoColours)
dat1to2GenomePlot(xplot = geno1_filt_winDF_list[[2]]$cumwindow,
                  dat1 = geno1_filt_winDF_list[[2]][,4],
                  dat2 = geno2_filt_winDF_list[[2]][,4],
                  Ylab = expression("Normalized ChIP"),
                  legendNames = c(bquote(.(geno1NamesPlot[2])),
                                  as.expression(bquote(italic(.(unlist(strsplit(geno2NamesPlot[2],
                                                                                split = " "))[1])) ~
                                                       .(unlist(strsplit(geno2NamesPlot[2],
                                                                         split = " "))[2])))),
                  plotColours = genoColours)
dat1to2GenomePlot(xplot = geno1_filt_winDF_list[[3]]$cumwindow,
                  dat1 = geno1_filt_winDF_list[[3]][,4],
                  dat2 = geno2_filt_winDF_list[[1]][,4],
                  Ylab = expression("Normalized ChIP"),
                  legendNames = c(bquote(.(geno1NamesPlot[3])),
                                  as.expression(bquote(italic(.(unlist(strsplit(geno2NamesPlot[1],
                                                                                split = " "))[1])) ~
                                                       .(unlist(strsplit(geno2NamesPlot[1],
                                                                         split = " "))[2])))),
                  plotColours = genoColours)
dat1to2GenomePlot(xplot = geno1_filt_winDF_list[[3]]$cumwindow,
                  dat1 = geno1_filt_winDF_list[[3]][,4],
                  dat2 = geno2_filt_winDF_list[[2]][,4],
                  Ylab = expression("Normalized ChIP"),
                  legendNames = c(bquote(.(geno1NamesPlot[3])),
                                  as.expression(bquote(italic(.(unlist(strsplit(geno2NamesPlot[2],
                                                                                split = " "))[1])) ~
                                                       .(unlist(strsplit(geno2NamesPlot[2],
                                                                         split = " "))[2])))),
                  plotColours = genoColours)
# log2(ChIP/input)
dat1to2GenomePlot(xplot = geno1_filt_winDF_list[[1]]$cumwindow,
                  dat1 = geno1_filt_winDF_list[[1]][,6],
                  dat2 = geno2_filt_winDF_list[[1]][,6],
                  Ylab = expression("Log"[2]*"(ChIP/input)"),
                  legendNames = c(bquote(.(geno1NamesPlot[1])),
                                  as.expression(bquote(italic(.(unlist(strsplit(geno2NamesPlot[1],
                                                                                split = " "))[1])) ~
                                                       .(unlist(strsplit(geno2NamesPlot[1],
                                                                         split = " "))[2])))),
                  plotColours = genoColours)
dat1to2GenomePlot(xplot = geno1_filt_winDF_list[[1]]$cumwindow,
                  dat1 = geno1_filt_winDF_list[[1]][,6],
                  dat2 = geno2_filt_winDF_list[[2]][,6],
                  Ylab = expression("Log"[2]*"(ChIP/input)"),
                  legendNames = c(bquote(.(geno1NamesPlot[1])),
                                  as.expression(bquote(italic(.(unlist(strsplit(geno2NamesPlot[2],
                                                                                split = " "))[1])) ~
                                                       .(unlist(strsplit(geno2NamesPlot[2],
                                                                         split = " "))[2])))),
                  plotColours = genoColours)
dat1to2GenomePlot(xplot = geno1_filt_winDF_list[[2]]$cumwindow,
                  dat1 = geno1_filt_winDF_list[[2]][,6],
                  dat2 = geno2_filt_winDF_list[[1]][,6],
                  Ylab = expression("Log"[2]*"(ChIP/input)"),
                  legendNames = c(bquote(.(geno1NamesPlot[2])),
                                  as.expression(bquote(italic(.(unlist(strsplit(geno2NamesPlot[1],
                                                                                split = " "))[1])) ~
                                                       .(unlist(strsplit(geno2NamesPlot[1],
                                                                         split = " "))[2])))),
                  plotColours = genoColours)
dat1to2GenomePlot(xplot = geno1_filt_winDF_list[[2]]$cumwindow,
                  dat1 = geno1_filt_winDF_list[[2]][,6],
                  dat2 = geno2_filt_winDF_list[[2]][,6],
                  Ylab = expression("Log"[2]*"(ChIP/input)"),
                  legendNames = c(bquote(.(geno1NamesPlot[2])),
                                  as.expression(bquote(italic(.(unlist(strsplit(geno2NamesPlot[2],
                                                                                split = " "))[1])) ~
                                                       .(unlist(strsplit(geno2NamesPlot[2],
                                                                         split = " "))[2])))),
                  plotColours = genoColours)
dat1to2GenomePlot(xplot = geno1_filt_winDF_list[[3]]$cumwindow,
                  dat1 = geno1_filt_winDF_list[[3]][,6],
                  dat2 = geno2_filt_winDF_list[[1]][,6],
                  Ylab = expression("Log"[2]*"(ChIP/input)"),
                  legendNames = c(bquote(.(geno1NamesPlot[3])),
                                  as.expression(bquote(italic(.(unlist(strsplit(geno2NamesPlot[1],
                                                                                split = " "))[1])) ~
                                                       .(unlist(strsplit(geno2NamesPlot[1],
                                                                         split = " "))[2])))),
                  plotColours = genoColours)
dat1to2GenomePlot(xplot = geno1_filt_winDF_list[[3]]$cumwindow,
                  dat1 = geno1_filt_winDF_list[[3]][,6],
                  dat2 = geno2_filt_winDF_list[[2]][,6],
                  Ylab = expression("Log"[2]*"(ChIP/input)"),
                  legendNames = c(bquote(.(geno1NamesPlot[3])),
                                  as.expression(bquote(italic(.(unlist(strsplit(geno2NamesPlot[2],
                                                                                split = " "))[1])) ~
                                                       .(unlist(strsplit(geno2NamesPlot[2],
                                                                         split = " "))[2])))),
                  plotColours = genoColours)
dev.off()

