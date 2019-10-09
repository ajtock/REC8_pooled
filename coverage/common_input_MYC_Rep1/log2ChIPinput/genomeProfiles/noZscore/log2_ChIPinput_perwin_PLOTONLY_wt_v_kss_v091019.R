#!/applications/R/R-3.5.0/bin/Rscript

######################################
######################################
# Generate plots chromosome profiles
######################################
######################################

# REC8_HA_Rep2 (red)
# kss_REC8_HA_Rep1 (red4)

# Usage:
# Rscript ./log2_ChIPinput_perwin_PLOTONLY_wt_v_kss_v091019.R 10kb 'REC8_HA_Rep2_ChIP_REC8_MYC_Rep1_input,WT_SPO11oligo_RPI1_WT_nakedDNA_R1,WT_H3K9me2_ChIP_WT_H3K9me2_input' 'kss_REC8_HA_Rep1_ChIP_kss_REC8_HA_Rep1_input,suvh456_SPO11oligo_RPI34_WT_nakedDNA_R1,kss_H3K9me2_ChIP_kss_H3K9me2_input' 'wt REC8-HA_Rep2,wt SPO11-1-oligos_Rep1,wt H3K9me2,wt mCHG' 'kss REC8-HA_Rep1,kss SPO11-1-oligos_Rep2,kss H3K9me2,kss mCHG' 'red,dodgerblue2,green2,orange' 'red4,navy,darkgreen,orange4'

#winName <- "10kb"
#geno1Names <- unlist(strsplit("REC8_HA_Rep2_ChIP_REC8_MYC_Rep1_input,WT_SPO11oligo_RPI1_WT_nakedDNA_R1,WT_H3K9me2_ChIP_WT_H3K9me2_input",
#                              split = ","))
#geno2Names <- unlist(strsplit("kss_REC8_HA_Rep1_ChIP_kss_REC8_HA_Rep1_input,suvh456_SPO11oligo_RPI34_WT_nakedDNA_R1,kss_H3K9me2_ChIP_kss_H3K9me2_input",
#                              split = ","))
#geno1NamesPlot <- unlist(strsplit("wt REC8-HA_Rep2,wt SPO11-1-oligos_Rep1,wt H3K9me2,wt mCHG",
#                                  split = ","))
#geno2NamesPlot <- unlist(strsplit("kss REC8-HA_Rep1,kss SPO11-1-oligos_Rep2,kss H3K9me2,kss mCHG",
#                                  split = ","))
#geno1Colours <- unlist(strsplit("red,dodgerblue2,green2,orange",
#                                split = ","))
#geno2Colours <- unlist(strsplit("red4,navy,darkgreen,orange4",
#                                split = ","))

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
geno1Colours <- unlist(strsplit(args[6],
                                split = ","))
geno2Colours <- unlist(strsplit(args[7],
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
                    "_genome_norm_coverage_",
                    winName, "_noZscore.tsv"),
             header = T)
}, mc.cores = length(geno1Names))

geno2_filt_winDF_list <- mclapply(seq_along(geno2Names), function(x) {
  read.table(paste0("filt_", geno2Names[x],
                    "_genome_norm_coverage_",
                    winName, "_noZscore.tsv"),
             header = T)
}, mc.cores = length(geno2Names))

# Feature frequency
inDirFeatures <- "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/genomeProfiles/genomewideZscore/"

# DNA methylation in 200-kb windows
geno1_filt_meth <- read.table(paste0(inDirFeatures, "filt_DNAmeth_GSM980986_WT_rep2_genome_200kb.txt"))
geno2_filt_meth <- read.table(paste0(inDirFeatures, "filt_DNAmeth_GSM981060_suvh456_genome_200kb.txt"))

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

# Plot
pdf(paste0(plotDir, "wt_v_kss_ChIP_and_mCHG_genomeProfiles_noZscore_v091019.pdf"),
    height = 16, width = 12)
par(mfcol = c(4, 1))
par(mar = c(5.1, 6.1, 2.1, 6.1))
par(mgp = c(3, 1.5, 0))
# REC8 ChIP
dat1to2GenomePlot(xplot = geno1_filt_winDF_list[[1]]$cumwindow,
                  dat1 = geno1_filt_winDF_list[[1]][,4],
                  dat2 = geno2_filt_winDF_list[[1]][,4],
                  Ylab = expression("REC8-HA"),
                  legendNames = c("wt",
                                  expression(italic("kss"))),
                  plotColours = c(geno1Colours[1], geno2Colours[1]))
# SPO11-1-oligos
dat1to2GenomePlot(xplot = geno1_filt_winDF_list[[2]]$cumwindow,
                  dat1 = geno1_filt_winDF_list[[2]][,4],
                  dat2 = geno2_filt_winDF_list[[2]][,4],
                  Ylab = expression("SPO11-1-oligos"),
                  legendNames = c("wt",
                                  expression(italic("kss"))),
                  plotColours = c(geno1Colours[2], geno2Colours[2]))
# H3K9me2 ChIP 
dat1to2GenomePlot(xplot = geno1_filt_winDF_list[[3]]$cumwindow,
                  dat1 = geno1_filt_winDF_list[[3]][,4],
                  dat2 = geno2_filt_winDF_list[[3]][,4],
                  Ylab = expression("H3K9me2"),
                  legendNames = c("wt",
                                  expression(italic("kss"))),
                  plotColours = c(geno1Colours[3], geno2Colours[3]))
# mCHG
dat1to2GenomePlot(xplot = geno1_filt_meth$cumWindows,
                  dat1 = geno1_filt_meth[,4],
                  dat2 = geno2_filt_meth[,4],
                  Ylab = expression("CHG methylation"),
                  legendNames = c("wt",
                                  expression(italic("kss"))),
                  plotColours = c(geno1Colours[4], geno2Colours[4]))
dev.off()
## REC8 ChIP
#dat1to2GenomePlot(xplot = geno1_filt_winDF_list[[1]]$cumwindow,
#                  dat1 = geno1_filt_winDF_list[[1]][,4],
#                  dat2 = geno2_filt_winDF_list[[1]][,4],
#                  Ylab = expression("Normalized ChIP"),
#                  legendNames = c(bquote(.(geno1NamesPlot[1])),
#                                  as.expression(bquote(italic(.(unlist(strsplit(geno2NamesPlot[1],
#                                                                                split = " "))[1])) ~
#                                                       .(unlist(strsplit(geno2NamesPlot[1],
#                                                                         split = " "))[2])))),
#                  plotColours = c(geno1Colours[1], geno2Colours[1]))
## SPO11-1-oligos
#dat1to2GenomePlot(xplot = geno1_filt_winDF_list[[2]]$cumwindow,
#                  dat1 = geno1_filt_winDF_list[[2]][,4],
#                  dat2 = geno2_filt_winDF_list[[2]][,4],
#                  Ylab = expression("Normalized ChIP"),
#                  legendNames = c(bquote(.(geno1NamesPlot[2])),
#                                  as.expression(bquote(italic(.(unlist(strsplit(geno2NamesPlot[2],
#                                                                                split = " "))[1])) ~
#                                                       .(unlist(strsplit(geno2NamesPlot[2],
#                                                                         split = " "))[2])))),
#                  plotColours = c(geno1Colours[2], geno2Colours[2]))
## H3K9me2 ChIP 
#dat1to2GenomePlot(xplot = geno1_filt_winDF_list[[3]]$cumwindow,
#                  dat1 = geno1_filt_winDF_list[[3]][,4],
#                  dat2 = geno2_filt_winDF_list[[3]][,4],
#                  Ylab = expression("Normalized ChIP"),
#                  legendNames = c(bquote(.(geno1NamesPlot[3])),
#                                  as.expression(bquote(italic(.(unlist(strsplit(geno2NamesPlot[3],
#                                                                                split = " "))[1])) ~
#                                                       .(unlist(strsplit(geno2NamesPlot[3],
#                                                                         split = " "))[2])))),
#                  plotColours = c(geno1Colours[3], geno2Colours[3]))
## mCHG
#dat1to2GenomePlot(xplot = geno1_filt_meth$cumWindows,
#                  dat1 = geno1_filt_meth[,4],
#                  dat2 = geno2_filt_meth[,4],
#                  Ylab = expression("DNA methylation"),
#                  legendNames = c(bquote(.(geno1NamesPlot[4])),
#                                  as.expression(bquote(italic(.(unlist(strsplit(geno2NamesPlot[4],
#                                                                                split = " "))[1])) ~
#                                                       .(unlist(strsplit(geno2NamesPlot[4],
#                                                                         split = " "))[2])))),
#                  plotColours = c(geno1Colours[4], geno2Colours[4]))
#dev.off()
