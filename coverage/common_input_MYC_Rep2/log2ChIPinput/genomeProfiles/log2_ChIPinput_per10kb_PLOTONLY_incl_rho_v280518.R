#############################################################################################
# Generate chromosome-scale plots                                                           #
#############################################################################################

library(corrplot)
library(parallel)

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

winNames <- c("10kb")

plotDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/genomeProfiles/plots/"

outDir1 <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/genomeProfiles/"
ChIPnames1 <- "REC8_HA_Rep1_ChIP"
names1 <- "REC8-HA Rep1"

outDir2 <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/genomeProfiles/"
ChIPnames2 <- "REC8_MYC_Rep1_ChIP"
names2 <- "REC8-MYC Rep1"

outDir3 <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/genomeProfiles/"
ChIPnames3 <- "REC8_MYC_Rep2_ChIP"
names3 <- "REC8-MYC Rep2"

outDir4 <- "/projects/ajt200/BAM_masters/H3K9me2/WT/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames4 <- "WT_H3K9me2_ChIP"
names4 <- "wt H3K9me2"

outDir5 <- "/projects/ajt200/BAM_masters/H3K9me2/kss/coverage/Xiaohui/log2ChIPinput/genomeProfiles/"
ChIPnames5 <- "kss_H3K9me2_ChIP"
names5 <- "kss H3K9me2"

outDir6 <- "/projects/ajt200/BAM_masters/H3K9me2/cmt3/coverage/Xiaohui/log2ChIPinput/genomeProfiles/"
ChIPnames6 <- "cmt3_H3K9me2_ChIP"
names6 <- "cmt3 H3K9me2"

outDir7 <- "/home/meiosis/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/WT/coverage/genomeProfiles/"
ChIPnames7 <- "WT_RNAseq_Chris_Rep1"
names7 <- "wt RNA-seq (Chris) Rep1"

outDir8 <- "/home/meiosis/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/WT/coverage/genomeProfiles/"
ChIPnames8 <- "WT_RNAseq_Chris_Rep2"
names8 <- "wt RNA-seq (Chris) Rep2"

outDir9 <- "/home/meiosis/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/kss/coverage/genomeProfiles/"
ChIPnames9 <- "kss_RNAseq_Chris_Rep1"
names9 <- "kss RNA-seq (Chris) Rep1"

outDir10 <- "/home/meiosis/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/kss/coverage/genomeProfiles/"
ChIPnames10 <- "kss_RNAseq_Chris_Rep2"
names10 <- "kss RNA-seq (Chris) Rep2"

outDirCombined <- c(outDir1, outDir2, outDir3, outDir4, outDir5, outDir6, outDir7, outDir8, outDir9, outDir10)
ChIPnamesCombined <- c(ChIPnames1, ChIPnames2, ChIPnames3, ChIPnames4, ChIPnames5, ChIPnames6, ChIPnames7, ChIPnames8, ChIPnames9, ChIPnames10)
namesCombined <- c(names1, names2, names3, names4, names5, names6, names7, names8, names9, names10)

# Load concatenated filtered, log2-transformed coverage datasets
filt_log2trans_list <- mclapply(seq_along(ChIPnamesCombined), function(k) {
  read.table(file = paste0(outDirCombined[k], "filt_log2_", ChIPnamesCombined[k], "_genome_norm_coverage_", winNames[1], ".txt"))
}, mc.cores = length(ChIPnamesCombined))
filt_noNA_log2trans_list <- mclapply(seq_along(ChIPnamesCombined), function(k) {
  read.table(file = paste0(outDirCombined[k], "filt_noNA_log2_", ChIPnamesCombined[k], "_genome_norm_coverage_", winNames[1], ".txt"))
}, mc.cores = length(ChIPnamesCombined))

# Function to plot REC8 genome-scale coverage overlaid with other datasets 
REC8vsOthersGenomePlot <- function(xplot, filt_REC8_norm, filt_REC8_norm_noNA, otherFiltNorm, otherFiltNormNoNA, Ylabel, otherYlabel) {
  plot(xplot, otherFiltNorm, type = "l", lwd = 1.5, col = "blue",
       ylim = c(min(otherFiltNormNoNA), max(otherFiltNormNoNA)),
       xlab = "",
       ylab = "", xaxt = "n", yaxt = "n",
       main = bquote(italic("r"[s]) ~ " = " ~ .(round(cor(filt_REC8_norm, otherFiltNorm, method = "spearman"), digits = 2))))
  axis(side = 4, at = pretty(otherFiltNorm), lwd.tick = 1.5)
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 1.5, adj = c(0.5, -2.5), labels = otherYlabel, xpd = NA, srt = -90, col = "blue")
#  mtext(side = 4, line = 2.5, cex = 1, text = otherYlabel, col = "blue")
  par(new = T)
  plot(xplot, filt_REC8_norm, type = "l", lwd = 1.5, col = "red",
       ylim = c(min(filt_REC8_norm_noNA), max(filt_REC8_norm_noNA)),
       xlab = "",
       ylab = "")
  axis(side = 2, lwd.tick = 1.5)
  mtext(side = 2, line = 2.25, cex = 1, text = Ylabel, col = "red")
  abline(v = sumchr, lty = 1, lwd = 0.75)
  abline(v = centromeres, lty = 5, lwd = 0.75, col = "black")
  box(lwd = 1.5)
}

# Function to plot REC8 genome-scale coverage overlaid with DNA methylation (each context)
REC8vsDNAmethSepGenomePlot <- function(xplot, xplotMeth, filt_REC8_norm, filt_REC8_norm_noNA, otherFiltNorm, otherFiltNormNoNA, Ylabel, otherYlabel) {
  plot(xplotMeth, otherFiltNorm[,2], type = "l", lwd = 1.5, col = "navy",
       ylim = c(min(otherFiltNormNoNA[,2], otherFiltNormNoNA[,3], otherFiltNormNoNA[,4]), max(otherFiltNormNoNA[,2], otherFiltNormNoNA[,3], otherFiltNormNoNA[,4])),
       xlab = "",
       ylab = "", xaxt = "n", yaxt = "n")
  axis(side = 4, at = pretty(otherFiltNorm[,2]), lwd.tick = 1.5)
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 1.5, adj = c(0.5, -2.5), labels = otherYlabel, xpd = NA, srt = -90, col = "blue")
#  mtext(side = 4, line = 2.5, cex = 1, text = otherYlabel, col = "blue")
  par(new = T)
  plot(xplotMeth, otherFiltNorm[,3], type = "l", lwd = 1.5, col = "blue",
       ylim = c(min(otherFiltNormNoNA[,2], otherFiltNormNoNA[,3], otherFiltNormNoNA[,4]), max(otherFiltNormNoNA[,2], otherFiltNormNoNA[,3], otherFiltNormNoNA[,4])),
       xlab = "",
       ylab = "", xaxt = "n", yaxt = "n")
  par(new = T)
  plot(xplotMeth, otherFiltNorm[,4], type = "l", lwd = 1.5, col = "deepskyblue1",
       ylim = c(min(otherFiltNormNoNA[,2], otherFiltNormNoNA[,3], otherFiltNormNoNA[,4]), max(otherFiltNormNoNA[,2], otherFiltNormNoNA[,3], otherFiltNormNoNA[,4])),
       xlab = "",
       ylab = "", xaxt = "n", yaxt = "n")
  par(new = T)
  plot(xplot, filt_REC8_norm, type = "l", lwd = 1.5, col = "red",
       ylim = c(min(filt_REC8_norm_noNA), max(filt_REC8_norm_noNA)),
       xlab = "",
       ylab = "")
  axis(side = 2, lwd.tick = 1.5)
  axis(side = 1, lwd.tick = 1.5)
  mtext(side = 2, line = 2.25, cex = 1, text = Ylabel, col = "red")
  abline(v = sumchr, lty = 1, lwd = 0.75)
  abline(v = centromeres, lty = 5, lwd = 0.75, col = "black")
#  abline(h = mean(filt_REC8_norm, na.rm = T), lty = 5, col = "red")
  box(lwd = 1.5)
  legend("topleft",
         legend = c("CG", "CHG", "CHH"),
         col = c("navy", "blue", "deepskyblue1"),
         text.col = c("navy", "blue", "deepskyblue1"),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
}

# Function to plot genome-scale coverage overlaid with alternate genotype
genoVsOtherGenoGenomePlot <- function(xplot, filt_norm, filt_norm_noNA, filt_other_norm, filt_other_norm_noNA, Ylabel, otherYlabel) {
  plot(xplot, filt_norm, type = "l", lwd = 1.5, col = "blue",
       ylim = c(min(filt_norm_noNA, filt_other_norm_noNA), max(filt_norm_noNA, filt_other_norm_noNA)),
       xlab = "",
       ylab = "", xaxt = "n", yaxt = "n",
       main = bquote(italic("r"[s]) ~ " = " ~ .(round(cor(filt_norm, filt_other_norm, method = "spearman"), digits = 2))))
  mtext(side = 2, line = 2.25, cex = 1, text = Ylabel, col = "blue")
  axis(side = 2, lwd.tick = 1.5)
  lines(xplot, filt_other_norm, lwd = 1.5, col = "green")
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 1.5, adj = c(0.5, -2.5), labels = otherYlabel, xpd = NA, srt = -90, col = "green")
  axis(side = 1, lwd.tick = 1.5)
  abline(v = sumchr, lty = 1, lwd = 0.75)
  abline(v = centromeres, lty = 5, lwd = 0.75, col = "black")
  box(lwd = 1.5)
}


# Function to plot REC8 genome-scale coverage (all replicates) 
REC8GenomePlot <- function(xplot, dat1, dat2, dat3, dat1noNA, dat2noNA, dat3noNA) {
  plot(xplot, dat1, type = "l", lwd = 1.5, col = "red",
       ylim = c(min(dat1noNA, dat2noNA, dat3noNA), max(dat1noNA, dat2noNA, dat3noNA)),
       xlab = "",
       ylab = "")
  axis(side = 1, lwd.tick = 1.5)
  axis(side = 2, lwd.tick = 1.5)
  lines(xplot, dat2, lwd = 1.5, col = "green")
  lines(xplot, dat3, lwd = 1.5, col = "blue")
  mtext(side = 2, line = 2.25, cex = 1, text = "log2(ChIP/input)")
  abline(v = sumchr, lty = 1, lwd = 0.75)
  abline(v = centromeres, lty = 5, lwd = 0.75)
  #abline(h = mean(dat1, na.rm = T), lty = 5, col = "red")
  #abline(h = mean(dat2, na.rm = T), lty = 5, col = "green")
  #abline(h = mean(dat3, na.rm = T), lty = 5, col = "blue")
  box(lwd = 1.5)
  legend("topleft",
         legend = c("REC8_HA_Rep1", "REC8_MYC_Rep2", "REC8_MYC_Rep1"),
         col = c("red", "green", "blue"),
         text.col = c("red", "green", "blue"),
         ncol = 1, cex = 0.6, lwd = 1.5, bty = "n")
}


##################
## REC8-HA Rep1 ##
##################
pdf(file = paste0(plotDir, "REC8_HA_Rep1_vs_H3K9me2_RNAseq_log2_zscore_genomeplot_", winNames[1], "_rho.pdf"), height = 15, width = 15)
par(mfcol = c(6, 2))
par(mar = c(2.1, 4.1, 2.1, 4.1))
par(mgp = c(3, 1, 0))

for(i in 4:length(filt_log2trans_list)) {
  REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[1]][,1],
                         filt_REC8_norm = filt_log2trans_list[[1]][,2],
                         filt_REC8_norm_noNA = filt_noNA_log2trans_list[[1]][,1],
                         otherFiltNorm = filt_log2trans_list[[i]][,2],
                         otherFiltNormNoNA = filt_noNA_log2trans_list[[i]][,1],
                         Ylabel = namesCombined[1],
                         otherYlabel = namesCombined[i])
}
genoVsOtherGenoGenomePlot(xplot = filt_log2trans_list[[4]][,1],
                          filt_norm = filt_log2trans_list[[4]][,2],
                          filt_norm_noNA = filt_noNA_log2trans_list[[4]][,1],
                          filt_other_norm = filt_log2trans_list[[5]][,2],
                          filt_other_norm_noNA = filt_noNA_log2trans_list[[5]][,1],
                          Ylabel = namesCombined[4],
                          otherYlabel = namesCombined[5])
genoVsOtherGenoGenomePlot(xplot = filt_log2trans_list[[4]][,1],
                          filt_norm = filt_log2trans_list[[4]][,2],
                          filt_norm_noNA = filt_noNA_log2trans_list[[4]][,1],
                          filt_other_norm = filt_log2trans_list[[6]][,2],
                          filt_other_norm_noNA = filt_noNA_log2trans_list[[6]][,1],
                          Ylabel = namesCombined[4],
                          otherYlabel = namesCombined[6])
genoVsOtherGenoGenomePlot(xplot = filt_log2trans_list[[5]][,1],
                          filt_norm = filt_log2trans_list[[5]][,2],
                          filt_norm_noNA = filt_noNA_log2trans_list[[5]][,1],
                          filt_other_norm = filt_log2trans_list[[6]][,2],
                          filt_other_norm_noNA = filt_noNA_log2trans_list[[6]][,1],
                          Ylabel = namesCombined[5],
                          otherYlabel = namesCombined[6])
genoVsOtherGenoGenomePlot(xplot = filt_log2trans_list[[7]][,1],
                          filt_norm = filt_log2trans_list[[7]][,2],
                          filt_norm_noNA = filt_noNA_log2trans_list[[7]][,1],
                          filt_other_norm = filt_log2trans_list[[9]][,2],
                          filt_other_norm_noNA = filt_noNA_log2trans_list[[9]][,1],
                          Ylabel = namesCombined[7],
                          otherYlabel = namesCombined[9])
genoVsOtherGenoGenomePlot(xplot = filt_log2trans_list[[8]][,1],
                          filt_norm = filt_log2trans_list[[8]][,2],
                          filt_norm_noNA = filt_noNA_log2trans_list[[8]][,1],
                          filt_other_norm = filt_log2trans_list[[10]][,2],
                          filt_other_norm_noNA = filt_noNA_log2trans_list[[10]][,1],
                          Ylabel = namesCombined[8],
                          otherYlabel = namesCombined[10])
dev.off()

##################
## REC8-MYC Rep1 ##
##################
pdf(file = paste0(plotDir, "REC8_MYC_Rep1_vs_H3K9me2_RNAseq_log2_zscore_genomeplot_", winNames[1], "_rho.pdf"), height = 15, width = 15)
par(mfcol = c(6, 2))
par(mar = c(2.1, 4.1, 2.1, 4.1))
par(mgp = c(3, 1, 0))

for(i in 4:length(filt_log2trans_list)) {
  REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[2]][,1],
                         filt_REC8_norm = filt_log2trans_list[[2]][,2],
                         filt_REC8_norm_noNA = filt_noNA_log2trans_list[[2]][,1],
                         otherFiltNorm = filt_log2trans_list[[i]][,2],
                         otherFiltNormNoNA = filt_noNA_log2trans_list[[i]][,1],
                         Ylabel = namesCombined[2],
                         otherYlabel = namesCombined[i])
}
genoVsOtherGenoGenomePlot(xplot = filt_log2trans_list[[4]][,1],
                          filt_norm = filt_log2trans_list[[4]][,2],
                          filt_norm_noNA = filt_noNA_log2trans_list[[4]][,1],
                          filt_other_norm = filt_log2trans_list[[5]][,2],
                          filt_other_norm_noNA = filt_noNA_log2trans_list[[5]][,1],
                          Ylabel = namesCombined[4],
                          otherYlabel = namesCombined[5])
genoVsOtherGenoGenomePlot(xplot = filt_log2trans_list[[4]][,1],
                          filt_norm = filt_log2trans_list[[4]][,2],
                          filt_norm_noNA = filt_noNA_log2trans_list[[4]][,1],
                          filt_other_norm = filt_log2trans_list[[6]][,2],
                          filt_other_norm_noNA = filt_noNA_log2trans_list[[6]][,1],
                          Ylabel = namesCombined[4],
                          otherYlabel = namesCombined[6])
genoVsOtherGenoGenomePlot(xplot = filt_log2trans_list[[5]][,1],
                          filt_norm = filt_log2trans_list[[5]][,2],
                          filt_norm_noNA = filt_noNA_log2trans_list[[5]][,1],
                          filt_other_norm = filt_log2trans_list[[6]][,2],
                          filt_other_norm_noNA = filt_noNA_log2trans_list[[6]][,1],
                          Ylabel = namesCombined[5],
                          otherYlabel = namesCombined[6])
genoVsOtherGenoGenomePlot(xplot = filt_log2trans_list[[7]][,1],
                          filt_norm = filt_log2trans_list[[7]][,2],
                          filt_norm_noNA = filt_noNA_log2trans_list[[7]][,1],
                          filt_other_norm = filt_log2trans_list[[9]][,2],
                          filt_other_norm_noNA = filt_noNA_log2trans_list[[9]][,1],
                          Ylabel = namesCombined[7],
                          otherYlabel = namesCombined[9])
genoVsOtherGenoGenomePlot(xplot = filt_log2trans_list[[8]][,1],
                          filt_norm = filt_log2trans_list[[8]][,2],
                          filt_norm_noNA = filt_noNA_log2trans_list[[8]][,1],
                          filt_other_norm = filt_log2trans_list[[10]][,2],
                          filt_other_norm_noNA = filt_noNA_log2trans_list[[10]][,1],
                          Ylabel = namesCombined[8],
                          otherYlabel = namesCombined[10])
dev.off()

##################
## REC8-MYC Rep2 ##
##################
pdf(file = paste0(plotDir, "REC8_MYC_Rep2_vs_H3K9me2_RNAseq_log2_zscore_genomeplot_", winNames[1], "_rho.pdf"), height = 15, width = 15)
par(mfcol = c(6, 2))
par(mar = c(2.1, 4.1, 2.1, 4.1))
par(mgp = c(3, 1, 0))

for(i in 4:length(filt_log2trans_list)) {
  REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[3]][,1],
                         filt_REC8_norm = filt_log2trans_list[[3]][,2],
                         filt_REC8_norm_noNA = filt_noNA_log2trans_list[[3]][,1],
                         otherFiltNorm = filt_log2trans_list[[i]][,2],
                         otherFiltNormNoNA = filt_noNA_log2trans_list[[i]][,1],
                         Ylabel = namesCombined[3],
                         otherYlabel = namesCombined[i])
}
genoVsOtherGenoGenomePlot(xplot = filt_log2trans_list[[4]][,1],
                          filt_norm = filt_log2trans_list[[4]][,2],
                          filt_norm_noNA = filt_noNA_log2trans_list[[4]][,1],
                          filt_other_norm = filt_log2trans_list[[5]][,2],
                          filt_other_norm_noNA = filt_noNA_log2trans_list[[5]][,1],
                          Ylabel = namesCombined[4],
                          otherYlabel = namesCombined[5])
genoVsOtherGenoGenomePlot(xplot = filt_log2trans_list[[4]][,1],
                          filt_norm = filt_log2trans_list[[4]][,2],
                          filt_norm_noNA = filt_noNA_log2trans_list[[4]][,1],
                          filt_other_norm = filt_log2trans_list[[6]][,2],
                          filt_other_norm_noNA = filt_noNA_log2trans_list[[6]][,1],
                          Ylabel = namesCombined[4],
                          otherYlabel = namesCombined[6])
genoVsOtherGenoGenomePlot(xplot = filt_log2trans_list[[5]][,1],
                          filt_norm = filt_log2trans_list[[5]][,2],
                          filt_norm_noNA = filt_noNA_log2trans_list[[5]][,1],
                          filt_other_norm = filt_log2trans_list[[6]][,2],
                          filt_other_norm_noNA = filt_noNA_log2trans_list[[6]][,1],
                          Ylabel = namesCombined[5],
                          otherYlabel = namesCombined[6])
genoVsOtherGenoGenomePlot(xplot = filt_log2trans_list[[7]][,1],
                          filt_norm = filt_log2trans_list[[7]][,2],
                          filt_norm_noNA = filt_noNA_log2trans_list[[7]][,1],
                          filt_other_norm = filt_log2trans_list[[9]][,2],
                          filt_other_norm_noNA = filt_noNA_log2trans_list[[9]][,1],
                          Ylabel = namesCombined[7],
                          otherYlabel = namesCombined[9])
genoVsOtherGenoGenomePlot(xplot = filt_log2trans_list[[8]][,1],
                          filt_norm = filt_log2trans_list[[8]][,2],
                          filt_norm_noNA = filt_noNA_log2trans_list[[8]][,1],
                          filt_other_norm = filt_log2trans_list[[10]][,2],
                          filt_other_norm_noNA = filt_noNA_log2trans_list[[10]][,1],
                          Ylabel = namesCombined[8],
                          otherYlabel = namesCombined[10])
dev.off()

