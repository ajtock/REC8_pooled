#############################################################################################
# Generate chromosome-scale plots                                                           #
#############################################################################################

#MNase
#SPO11-1 ChIP
#H2A.Z
#H3K4me3
#Genes
#DNA methylation CG CHG CHH
#H2A.W
#H3K9me2
#TEs
#H3K4me1
#H3K4me2
#H3K27me1
#H3K27me3
#AT:GC
#SPO11-1-oligos
#Crossovers
#MSH4
#REC8 reps

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

outDir2 <- "/home/meiosis/ajt200/analysis/180622_Chris_lambing_ChIP_REC8_HA_Col_kss/WT/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames2 <- "REC8_HA_Rep2_ChIP"
names2 <- "REC8-HA Rep2"

outDir3 <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/genomeProfiles/"
ChIPnames3 <- "REC8_MYC_Rep1_ChIP"
names3 <- "REC8-MYC Rep1"

outDir4 <- "/projects/ajt200/BAM_masters/nucleosomes/WT/coverage/nakedDNA_untrimmed_input/log2ChIPinput/genomeProfiles/"
ChIPnames4 <- "WT_nuc"
names4 <- "MNase"

outDir5 <- "/projects/ajt200/BAM_masters/SPO11_ChIP/WT/coverage/REC8_MYC_Rep2_input/log2ChIPinput/genomeProfiles/"
ChIPnames5 <- "WT_SPO11_ChIP4"
names5 <- "SPO11-1 ChIP"

outDir6 <- "/projects/ajt200/BAM_masters/H2A/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames6 <- "H2AZ_ChIP"
names6 <- "H2A.Z"

outDir7 <- "/projects/ajt200/BAM_masters/H3K4me3/replicates/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames7 <- "WT_H3K4me3_ChIP14"
names7 <- "H3K4me3"

outDir8 <- "/projects/ajt200/BAM_masters/H2A/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames8 <- "H2AW_ChIP"
names8 <- "H2A.W"

outDir9 <- "/projects/ajt200/BAM_masters/H3K9me2/WT/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames9 <- "WT_H3K9me2_ChIP"
names9 <- "H3K9me2"

outDir10 <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K4me1/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames10 <- "WT_H3K4me1_Rep1_ChIP"
names10 <- "H3K4me1"

outDir11 <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K4me2/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames11 <- "WT_H3K4me2_Rep1_ChIP"
names11 <- "H3K4me2"

outDir12 <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K27me1/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames12 <- "WT_H3K27me1_Rep1_ChIP"
names12 <- "H3K27me1"

outDir13 <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K27me3/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames13 <- "WT_H3K27me3_Rep1_ChIP"
names13 <- "H3K27me3"

outDir14 <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames14 <- "WT_SPO11_oligo_RPI1"
names14 <- "SPO11-1-oligos"

outDir15 <- "/projects/ajt200/BAM_masters/MSH4/WT/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames15 <- "WT_MSH4_ChIP"
names15 <- "MSH4"

outDir16 <- "/projects/ajt200/BAM_masters/PolIV_Law_Jacobsen_2013_Nature/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames16 <- "PolIV_Rep2_ChIP"
names16 <- "Pol IV"

outDir17 <- "/projects/ajt200/BAM_masters/PolV_Liu_Jacobsen_2018_NatPlants/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames17 <- "PolV_ChIP"
names17 <- "Pol V"

outDir18 <- "/projects/ajt200/BAM_masters/RNAseq_meiocyte_Walker_Feng_2018_NatGenet/WT/meiocyte/coverage/genomeProfiles/"
ChIPnames18 <- "WT_RNAseq_meiocyte_Rep1"
names18 <- "wt RNA-seq meiocyte Rep1"

outDir19 <- "/projects/ajt200/BAM_masters/RNAseq_meiocyte_Walker_Feng_2018_NatGenet/WT/meiocyte/coverage/genomeProfiles/"
ChIPnames19 <- "WT_RNAseq_meiocyte_Rep2"
names19 <- "wt RNA-seq meiocyte Rep2"

outDir20 <- "/projects/ajt200/BAM_masters/RNAseq_meiocyte_Walker_Feng_2018_NatGenet/WT/meiocyte/coverage/genomeProfiles/"
ChIPnames20 <- "WT_RNAseq_meiocyte_Rep3"
names20 <- "wt RNA-seq meiocyte Rep3"

outDir21 <- "/projects/ajt200/BAM_masters/RNAseq/WT/coverage/genomeProfiles/"
ChIPnames21 <- "WT_RNAseq_Rep1"
names21 <- "wt RNA-seq (Kyuha) Rep1"

outDir22 <- "/projects/ajt200/BAM_masters/RNAseq/WT/coverage/genomeProfiles/"
ChIPnames22 <- "WT_RNAseq_Rep2"
names22 <- "wt RNA-seq (Kyuha) Rep2"

outDir23 <- "/projects/ajt200/BAM_masters/RNAseq/WT/coverage/genomeProfiles/"
ChIPnames23 <- "WT_RNAseq_Rep3"
names23 <- "wt RNA-seq (Kyuha) Rep3"

outDir24 <- "/home/meiosis/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/WT/coverage/genomeProfiles/"
ChIPnames24 <- "WT_RNAseq_Chris_Rep1"
names24 <- "wt RNA-seq (Chris) Rep1"

outDir25 <- "/home/meiosis/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/WT/coverage/genomeProfiles/"
ChIPnames25 <- "WT_RNAseq_Chris_Rep2"
names25 <- "wt RNA-seq (Chris) Rep2"

outDirCombined <- c(outDir1, outDir2, outDir3, outDir4, outDir5, outDir6, outDir7, outDir8, outDir9,
                    outDir10, outDir11, outDir12, outDir13, outDir14, outDir15, outDir16, outDir17,
                    outDir18, outDir19, outDir20, outDir21, outDir22, outDir23, outDir24, outDir25)
ChIPnamesCombined <- c(ChIPnames1, ChIPnames2, ChIPnames3, ChIPnames4, ChIPnames5, ChIPnames6,
                       ChIPnames7, ChIPnames8, ChIPnames9, ChIPnames10, ChIPnames11, ChIPnames12,
                       ChIPnames13, ChIPnames14, ChIPnames15, ChIPnames16, ChIPnames17,
                       ChIPnames18, ChIPnames19, ChIPnames20, ChIPnames21, ChIPnames22, ChIPnames23, ChIPnames24, ChIPnames25)
namesCombined <- c(names1, names2, names3, names4, names5, names6,
                   names7, names8, names9, names10, names11, names12,
                   names13, names14, names15, names16, names17,
                   names18, names19, names20, names21, names22, names23, names24, names25)

# Load concatenated filtered, log2-transformed coverage datasets
filt_log2trans_list <- mclapply(seq_along(ChIPnamesCombined), function(k) {
  read.table(file = paste0(outDirCombined[k], "filt_log2_", ChIPnamesCombined[k], "_genome_norm_coverage_", winNames[1], ".txt"))
}, mc.cores = length(ChIPnamesCombined))
filt_noNA_log2trans_list <- mclapply(seq_along(ChIPnamesCombined), function(k) {
  read.table(file = paste0(outDirCombined[k], "filt_noNA_log2_", ChIPnamesCombined[k], "_genome_norm_coverage_", winNames[1], ".txt"))
}, mc.cores = length(ChIPnamesCombined))

inDirFeatures <- "/projects/ajt200/REC8_MSH4/data_merged_fastq/coverage/log2ChIPinput/genomeProfiles/"
inDirATGC <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/REC8/coverage/common_input_MYC_Rep2/log2ChIPinput/genomeProfiles/"

# Load concatenated filtered gene dataset
filt_GeneDat <- read.table(file = paste0(inDirFeatures, "filt_gene_density_genome_", winNames[1], ".txt"))
filt_GeneDat_noNA <- read.table(file = paste0(inDirFeatures, "filt_noNA_gene_density_genome_", winNames[1], ".txt"))

# Load concatenated filtered DNA meth datasets
filt_methDat <- read.table(file = paste0(inDirFeatures, "filt_wtMeth_genome_200kb.txt"))
filt_noNA_methDat <- read.table(file = paste0(inDirFeatures, "filt_noNA_wtMeth_genome_200kb.txt"))

# Load concatenated filtered TE dataset
filt_TEDat <- read.table(file = paste0(inDirFeatures, "filt_TE_density_genome_", winNames[1], ".txt"))
filt_TEDat_noNA <- read.table(file = paste0(inDirFeatures, "filt_noNA_TE_density_genome_", winNames[1], ".txt"))

# Load concatenated filtered AT and GC content datasets
filt_ATcontent <- read.table(file = paste0(inDirATGC, "filt_AT_content_genome_", winNames[1], ".txt"))
filt_ATcontent_noNA <- read.table(file = paste0(inDirATGC, "filt_noNA_AT_content_genome_", winNames[1], ".txt"))
filt_GCcontent <- read.table(file = paste0(inDirATGC, "filt_GC_content_genome_", winNames[1], ".txt"))
filt_GCcontent_noNA <- read.table(file = paste0(inDirATGC, "filt_noNA_GC_content_genome_", winNames[1], ".txt"))

# Load concatenated filtered CO density datasets
filt_CODat <- read.table(file = paste0(inDirFeatures, "filt_WTCO_density_genome_", winNames[1], ".txt"))
filt_CODat_noNA <- read.table(file = paste0(inDirFeatures, "filt_noNA_WTCO_density_genome_", winNames[1], ".txt"))


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

# Function to plot REC8 genome-scale coverage overlaid with feature density
REC8vsFeaturesGenomePlot <- function(xplot, filt_REC8_norm, filt_REC8_norm_noNA, otherFiltNorm, otherFiltNormNoNA, Ylabel, otherYlabel) {
  plot(xplot, otherFiltNorm, type = "l", lwd = 1.5, col = "blue",
       ylim = c(min(otherFiltNormNoNA), max(otherFiltNormNoNA)),
       xlab = "",
       ylab = "", xaxt = "n", yaxt = "n",
       main = bquote(italic("r"[s]) ~ " = " ~ .(round(cor(filt_REC8_norm, otherFiltNorm, method = "spearman"), digits = 2))))
  axis(side = 4, at = pretty(otherFiltNorm), lwd.tick = 1.5)
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 1.5, adj = c(0.5, -2.5), labels = otherYlabel, xpd = NA, srt = -90, col = "blue")
  #mtext(side = 4, line = 2.5, cex = 1, text = otherYlabel, col = "blue")
#  abline(h = mean(otherFiltNorm, na.rm = T), lty = 5, col = "blue")
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
#  abline(h = mean(filt_REC8_norm, na.rm = T), lty = 5, col = "black")
  box(lwd = 1.5)
}

# Function to plot REC8 genome-scale coverage overlaid with % AT and GC content
REC8vsATGCcontentGenomePlot <- function(xplot, filt_REC8_norm, filt_REC8_norm_noNA, ATcontentFiltNorm, ATcontentFiltNormNoNA, GCcontentFiltNorm, GCcontentFiltNormNoNA, Ylabel, otherYlabel) {
  plot(xplot, ATcontentFiltNorm, type = "l", lwd = 1.5, col = "blue",
       ylim = c(min(ATcontentFiltNormNoNA, GCcontentFiltNormNoNA), max(ATcontentFiltNormNoNA, GCcontentFiltNormNoNA)),
       xlab = "",
       ylab = "", xaxt = "n", yaxt = "n")
  axis(side = 4, lwd.tick = 1.5)
  lines(xplot, GCcontentFiltNorm, lwd = 1.5, col = "green")
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 1.5, adj = c(0.5, -2.5), labels = otherYlabel, xpd = NA, srt = -90, col = "blue")
#  mtext(side = 4, line = 2.5, cex = 1, text = otherYlabel, col = "black")
  #par(new = T)
  #plot(xplot, GCcontentFiltNorm, type = "l", lwd = 1.5, col = "green",
  #     ylim = c(min(ATcontentFiltNormNoNA, GCcontentFiltNormNoNA), max(ATcontentFiltNormNoNA, GCcontentFiltNormNoNA)),
  #     xlab = "",
  #     ylab = "", xaxt = "n", yaxt = "n")
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
         legend = c("AT", "GC"),
         col = c("blue", "green"),
         text.col = c("blue", "green"),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
}

# Function to plot REC8 genome-scale coverage overlaid with % AT or GC content
REC8vsbaseContentGenomePlot <- function(xplot, filt_REC8_norm, filt_REC8_norm_noNA, baseContentFiltNorm, baseContentFiltNormNoNA, Ylabel, otherYlabel) {
  plot(xplot, baseContentFiltNorm, type = "l", lwd = 1.5, col = "blue",
       ylim = c(min(baseContentFiltNormNoNA), max(baseContentFiltNormNoNA)),
       xlab = "",
       ylab = "", xaxt = "n", yaxt = "n",
       main = bquote(italic("r"[s]) ~ " = " ~ .(round(cor(filt_REC8_norm, baseContentFiltNorm, method = "spearman"), digits = 2))))
  axis(side = 4, lwd.tick = 1.5)
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 1.5, adj = c(0.5, -2.5), labels = otherYlabel, xpd = NA, srt = -90, col = "blue")
#  mtext(side = 4, line = 2.5, cex = 1, text = otherYlabel, col = "blue")
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
         legend = c("REC8-HA Rep1", "REC8-HA Rep2", "REC8-MYC Rep1"),
         col = c("red", "green", "blue"),
         text.col = c("red", "green", "blue"),
         ncol = 1, cex = 0.6, lwd = 1.5, bty = "n")
}

###################
## REC8 all reps ##
###################
pdf(file = paste0(plotDir, "REC8_HA_Rep1_Rep2_MYC_Rep1_log2_zscore_genomeplot_",
                  winNames[1], "_h3w9_rho_MYC_Rep2_input_v290618.pdf"), height = 3, width = 9)
par(mfcol = c(1, 1))
par(mar = c(2.1, 4.1, 2.1, 2.1))
par(mgp = c(3, 1, 0))

REC8GenomePlot(xplot = filt_log2trans_list[[1]][,1],
               dat1 = filt_log2trans_list[[1]][,2],
               dat2 = filt_log2trans_list[[2]][,2],
               dat3 = filt_log2trans_list[[3]][,2],
               dat1noNA = filt_noNA_log2trans_list[[1]][,1],
               dat2noNA = filt_noNA_log2trans_list[[2]][,1],
               dat3noNA = filt_noNA_log2trans_list[[3]][,1])
dev.off()



##################
## REC8-HA Rep1 ##
##################
pdf(file = paste0(plotDir, "REC8_HA_Rep1_vs_others_log2_zscore_genomeplot_", winNames[1], "_h35w15_rho_MYC_Rep2_input_v290618.pdf"), height = 35, width = 15)
par(mfcol = c(14, 2))
par(mar = c(2.1, 4.1, 2.1, 4.1))
par(mgp = c(3, 1, 0))

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[1]][,1], filt_REC8_norm = filt_log2trans_list[[1]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[1]][,1], otherFiltNorm = filt_log2trans_list[[4]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[4]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[4])
REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[1]][,1], filt_REC8_norm = filt_log2trans_list[[1]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[1]][,1], otherFiltNorm = filt_log2trans_list[[5]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[5]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[5])
REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[1]][,1], filt_REC8_norm = filt_log2trans_list[[1]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[1]][,1], otherFiltNorm = filt_log2trans_list[[6]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[6]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[6])
REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[1]][,1], filt_REC8_norm = filt_log2trans_list[[1]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[1]][,1], otherFiltNorm = filt_log2trans_list[[7]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[7]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[7])

# Genes
REC8vsFeaturesGenomePlot(xplot = filt_log2trans_list[[1]][,1], filt_REC8_norm = filt_log2trans_list[[1]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[1]][,1], otherFiltNorm = filt_GeneDat[,2], otherFiltNormNoNA = filt_GeneDat_noNA[,1], Ylabel = namesCombined[1], otherYlabel = "Genes")
# DNA methylation
REC8vsDNAmethSepGenomePlot(xplot = filt_log2trans_list[[1]][,1], xplotMeth = filt_methDat[,1], filt_REC8_norm = filt_log2trans_list[[1]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[1]][,1], otherFiltNorm = filt_methDat, otherFiltNormNoNA = filt_noNA_methDat, Ylabel = namesCombined[1], otherYlabel = "DNA methylation") 

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[1]][,1], filt_REC8_norm = filt_log2trans_list[[1]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[1]][,1], otherFiltNorm = filt_log2trans_list[[8]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[8]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[8])
REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[1]][,1], filt_REC8_norm = filt_log2trans_list[[1]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[1]][,1], otherFiltNorm = filt_log2trans_list[[9]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[9]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[9])

# TEs
REC8vsFeaturesGenomePlot(xplot = filt_log2trans_list[[1]][,1], filt_REC8_norm = filt_log2trans_list[[1]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[1]][,1], otherFiltNorm = filt_TEDat[,2], otherFiltNormNoNA = filt_TEDat_noNA[,1], Ylabel = namesCombined[1], otherYlabel = "TEs")

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[1]][,1], filt_REC8_norm = filt_log2trans_list[[1]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[1]][,1], otherFiltNorm = filt_log2trans_list[[10]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[10]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[10])
REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[1]][,1], filt_REC8_norm = filt_log2trans_list[[1]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[1]][,1], otherFiltNorm = filt_log2trans_list[[11]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[11]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[11])
REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[1]][,1], filt_REC8_norm = filt_log2trans_list[[1]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[1]][,1], otherFiltNorm = filt_log2trans_list[[12]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[12]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[12])
REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[1]][,1], filt_REC8_norm = filt_log2trans_list[[1]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[1]][,1], otherFiltNorm = filt_log2trans_list[[13]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[13]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[13])

# AT content
REC8vsbaseContentGenomePlot(xplot = filt_log2trans_list[[1]][,1], filt_REC8_norm = filt_log2trans_list[[1]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[1]][,1], baseContentFiltNorm = filt_ATcontent[,2], baseContentFiltNormNoNA = filt_ATcontent_noNA[,1],Ylabel = namesCombined[1], otherYlabel = "AT content (%)")

# GC content
REC8vsbaseContentGenomePlot(xplot = filt_log2trans_list[[1]][,1], filt_REC8_norm = filt_log2trans_list[[1]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[1]][,1], baseContentFiltNorm = filt_GCcontent[,2], baseContentFiltNormNoNA = filt_GCcontent_noNA[,1],Ylabel = namesCombined[1], otherYlabel = "GC content (%)")

#REC8vsATGCcontentGenomePlot(xplot = filt_log2trans_list[[1]][,1], filt_REC8_norm = filt_log2trans_list[[1]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[1]][,1], ATcontentFiltNorm = filt_ATcontent[,2], ATcontentFiltNormNoNA = filt_ATcontent_noNA[,1], GCcontentFiltNorm = filt_GCcontent[,2], GCcontentFiltNormNoNA = filt_GCcontent_noNA[,1], Ylabel = namesCombined[1], otherYlabel = "AT:GC content (%)")

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[1]][,1], filt_REC8_norm = filt_log2trans_list[[1]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[1]][,1], otherFiltNorm = filt_log2trans_list[[14]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[14]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[14])

# COs
REC8vsFeaturesGenomePlot(xplot = filt_log2trans_list[[1]][,1], filt_REC8_norm = filt_log2trans_list[[1]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[1]][,1], otherFiltNorm = filt_CODat[,2], otherFiltNormNoNA = filt_CODat_noNA[,1], Ylabel = namesCombined[1], otherYlabel = "Crossovers")

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[1]][,1], filt_REC8_norm = filt_log2trans_list[[1]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[1]][,1], otherFiltNorm = filt_log2trans_list[[15]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[15]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[15])

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[1]][,1], filt_REC8_norm = filt_log2trans_list[[1]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[1]][,1], otherFiltNorm = filt_log2trans_list[[16]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[16]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[16])

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[1]][,1], filt_REC8_norm = filt_log2trans_list[[1]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[1]][,1], otherFiltNorm = filt_log2trans_list[[17]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[17]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[17])

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[1]][,1], filt_REC8_norm = filt_log2trans_list[[1]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[1]][,1], otherFiltNorm = filt_log2trans_list[[18]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[18]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[18])

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[1]][,1], filt_REC8_norm = filt_log2trans_list[[1]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[1]][,1], otherFiltNorm = filt_log2trans_list[[19]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[19]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[19])

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[1]][,1], filt_REC8_norm = filt_log2trans_list[[1]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[1]][,1], otherFiltNorm = filt_log2trans_list[[20]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[20]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[20])

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[1]][,1], filt_REC8_norm = filt_log2trans_list[[1]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[1]][,1], otherFiltNorm = filt_log2trans_list[[21]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[21]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[21])

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[1]][,1], filt_REC8_norm = filt_log2trans_list[[1]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[1]][,1], otherFiltNorm = filt_log2trans_list[[22]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[22]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[22])

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[1]][,1], filt_REC8_norm = filt_log2trans_list[[1]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[1]][,1], otherFiltNorm = filt_log2trans_list[[23]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[23]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[23])

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[1]][,1], filt_REC8_norm = filt_log2trans_list[[1]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[1]][,1], otherFiltNorm = filt_log2trans_list[[24]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[24]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[24])

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[1]][,1], filt_REC8_norm = filt_log2trans_list[[1]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[1]][,1], otherFiltNorm = filt_log2trans_list[[25]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[25]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[25])

dev.off()


###################
## REC8-HA Rep2 ##
###################
pdf(file = paste0(plotDir, "REC8_HA_Rep2_vs_others_log2_zscore_genomeplot_", winNames[1], "_h35w15_rho_MYC_Rep2_input_v290618.pdf"), height = 35, width = 15)
par(mfcol = c(14, 2))
par(mar = c(2.1, 4.1, 2.1, 4.1))
par(mgp = c(3, 1, 0))

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[2]][,1], filt_REC8_norm = filt_log2trans_list[[2]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[2]][,1], otherFiltNorm = filt_log2trans_list[[4]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[4]][,1], Ylabel = namesCombined[2], otherYlabel = namesCombined[4])
REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[2]][,1], filt_REC8_norm = filt_log2trans_list[[2]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[2]][,1], otherFiltNorm = filt_log2trans_list[[5]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[5]][,1], Ylabel = namesCombined[2], otherYlabel = namesCombined[5])
REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[2]][,1], filt_REC8_norm = filt_log2trans_list[[2]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[2]][,1], otherFiltNorm = filt_log2trans_list[[6]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[6]][,1], Ylabel = namesCombined[2], otherYlabel = namesCombined[6])
REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[2]][,1], filt_REC8_norm = filt_log2trans_list[[2]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[2]][,1], otherFiltNorm = filt_log2trans_list[[7]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[7]][,1], Ylabel = namesCombined[2], otherYlabel = namesCombined[7])

# Genes
REC8vsFeaturesGenomePlot(xplot = filt_log2trans_list[[2]][,1], filt_REC8_norm = filt_log2trans_list[[2]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[2]][,1], otherFiltNorm = filt_GeneDat[,2], otherFiltNormNoNA = filt_GeneDat_noNA[,1], Ylabel = namesCombined[2], otherYlabel = "Genes")
# DNA methylation
REC8vsDNAmethSepGenomePlot(xplot = filt_log2trans_list[[2]][,1], xplotMeth = filt_methDat[,1], filt_REC8_norm = filt_log2trans_list[[2]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[2]][,1], otherFiltNorm = filt_methDat, otherFiltNormNoNA = filt_noNA_methDat, Ylabel = namesCombined[2], otherYlabel = "DNA methylation")

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[2]][,1], filt_REC8_norm = filt_log2trans_list[[2]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[2]][,1], otherFiltNorm = filt_log2trans_list[[8]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[8]][,1], Ylabel = namesCombined[2], otherYlabel = namesCombined[8])
REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[2]][,1], filt_REC8_norm = filt_log2trans_list[[2]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[2]][,1], otherFiltNorm = filt_log2trans_list[[9]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[9]][,1], Ylabel = namesCombined[2], otherYlabel = namesCombined[9])

# TEs
REC8vsFeaturesGenomePlot(xplot = filt_log2trans_list[[2]][,1], filt_REC8_norm = filt_log2trans_list[[2]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[2]][,1], otherFiltNorm = filt_TEDat[,2], otherFiltNormNoNA = filt_TEDat_noNA[,1], Ylabel = namesCombined[2], otherYlabel = "TEs")

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[2]][,1], filt_REC8_norm = filt_log2trans_list[[2]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[2]][,1], otherFiltNorm = filt_log2trans_list[[10]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[10]][,1], Ylabel = namesCombined[2], otherYlabel = namesCombined[10])
REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[2]][,1], filt_REC8_norm = filt_log2trans_list[[2]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[2]][,1], otherFiltNorm = filt_log2trans_list[[11]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[11]][,1], Ylabel = namesCombined[2], otherYlabel = namesCombined[11])
REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[2]][,1], filt_REC8_norm = filt_log2trans_list[[2]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[2]][,1], otherFiltNorm = filt_log2trans_list[[12]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[12]][,1], Ylabel = namesCombined[2], otherYlabel = namesCombined[12])
REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[2]][,1], filt_REC8_norm = filt_log2trans_list[[2]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[2]][,1], otherFiltNorm = filt_log2trans_list[[13]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[13]][,1], Ylabel = namesCombined[2], otherYlabel = namesCombined[13])

# AT content
REC8vsbaseContentGenomePlot(xplot = filt_log2trans_list[[2]][,1], filt_REC8_norm = filt_log2trans_list[[2]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[2]][,1], baseContentFiltNorm = filt_ATcontent[,2], baseContentFiltNormNoNA = filt_ATcontent_noNA[,1],Ylabel = namesCombined[2], otherYlabel = "AT content (%)")

# GC content
REC8vsbaseContentGenomePlot(xplot = filt_log2trans_list[[2]][,1], filt_REC8_norm = filt_log2trans_list[[2]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[2]][,1], baseContentFiltNorm = filt_GCcontent[,2], baseContentFiltNormNoNA = filt_GCcontent_noNA[,1],Ylabel = namesCombined[2], otherYlabel = "GC content (%)")

#REC8vsATGCcontentGenomePlot(xplot = filt_log2trans_list[[2]][,1], filt_REC8_norm = filt_log2trans_list[[2]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[2]][,1], ATcontentFiltNorm = filt_ATcontent[,2], ATcontentFiltNormNoNA = filt_ATcontent_noNA[,1], GCcontentFiltNorm = filt_GCcontent[,2], GCcontentFiltNormNoNA = filt_GCcontent_noNA[,1], Ylabel = namesCombined[2], otherYlabel = "AT:GC content (%)")

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[2]][,1], filt_REC8_norm = filt_log2trans_list[[2]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[2]][,1], otherFiltNorm = filt_log2trans_list[[14]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[14]][,1], Ylabel = namesCombined[2], otherYlabel = namesCombined[14])

# COs
REC8vsFeaturesGenomePlot(xplot = filt_log2trans_list[[2]][,1], filt_REC8_norm = filt_log2trans_list[[2]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[2]][,1], otherFiltNorm = filt_CODat[,2], otherFiltNormNoNA = filt_CODat_noNA[,1], Ylabel = namesCombined[2], otherYlabel = "Crossovers")

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[2]][,1], filt_REC8_norm = filt_log2trans_list[[2]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[2]][,1], otherFiltNorm = filt_log2trans_list[[15]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[15]][,1], Ylabel = namesCombined[2], otherYlabel = namesCombined[15])

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[2]][,1], filt_REC8_norm = filt_log2trans_list[[2]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[2]][,1], otherFiltNorm = filt_log2trans_list[[16]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[16]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[16])

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[2]][,1], filt_REC8_norm = filt_log2trans_list[[2]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[2]][,1], otherFiltNorm = filt_log2trans_list[[17]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[17]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[17])

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[2]][,1], filt_REC8_norm = filt_log2trans_list[[2]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[2]][,1], otherFiltNorm = filt_log2trans_list[[18]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[18]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[18])

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[2]][,1], filt_REC8_norm = filt_log2trans_list[[2]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[2]][,1], otherFiltNorm = filt_log2trans_list[[19]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[19]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[19])

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[2]][,1], filt_REC8_norm = filt_log2trans_list[[2]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[2]][,1], otherFiltNorm = filt_log2trans_list[[20]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[20]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[20])

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[2]][,1], filt_REC8_norm = filt_log2trans_list[[2]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[2]][,1], otherFiltNorm = filt_log2trans_list[[21]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[21]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[21])

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[2]][,1], filt_REC8_norm = filt_log2trans_list[[2]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[2]][,1], otherFiltNorm = filt_log2trans_list[[22]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[22]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[22])

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[2]][,1], filt_REC8_norm = filt_log2trans_list[[2]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[2]][,1], otherFiltNorm = filt_log2trans_list[[23]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[23]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[23])

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[2]][,1], filt_REC8_norm = filt_log2trans_list[[2]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[2]][,1], otherFiltNorm = filt_log2trans_list[[24]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[24]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[24])

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[2]][,1], filt_REC8_norm = filt_log2trans_list[[2]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[2]][,1], otherFiltNorm = filt_log2trans_list[[25]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[25]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[25])

dev.off()


###################
## REC8-MYC Rep1 ##
###################
pdf(file = paste0(plotDir, "REC8_MYC_Rep1_vs_others_log2_zscore_genomeplot_", winNames[1], "_h35w15_rho_MYC_Rep2_input_v290618.pdf"), height = 35, width = 15)
par(mfcol = c(14, 2))
par(mar = c(2.1, 4.1, 2.1, 4.1))
par(mgp = c(3, 1, 0))

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[3]][,1], filt_REC8_norm = filt_log2trans_list[[3]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[3]][,1], otherFiltNorm = filt_log2trans_list[[4]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[4]][,1], Ylabel = namesCombined[3], otherYlabel = namesCombined[4])
REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[3]][,1], filt_REC8_norm = filt_log2trans_list[[3]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[3]][,1], otherFiltNorm = filt_log2trans_list[[5]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[5]][,1], Ylabel = namesCombined[3], otherYlabel = namesCombined[5])
REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[3]][,1], filt_REC8_norm = filt_log2trans_list[[3]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[3]][,1], otherFiltNorm = filt_log2trans_list[[6]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[6]][,1], Ylabel = namesCombined[3], otherYlabel = namesCombined[6])
REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[3]][,1], filt_REC8_norm = filt_log2trans_list[[3]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[3]][,1], otherFiltNorm = filt_log2trans_list[[7]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[7]][,1], Ylabel = namesCombined[3], otherYlabel = namesCombined[7])

# Genes
REC8vsFeaturesGenomePlot(xplot = filt_log2trans_list[[3]][,1], filt_REC8_norm = filt_log2trans_list[[3]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[3]][,1], otherFiltNorm = filt_GeneDat[,2], otherFiltNormNoNA = filt_GeneDat_noNA[,1], Ylabel = namesCombined[3], otherYlabel = "Genes")
# DNA methylation
REC8vsDNAmethSepGenomePlot(xplot = filt_log2trans_list[[3]][,1], xplotMeth = filt_methDat[,1], filt_REC8_norm = filt_log2trans_list[[3]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[3]][,1], otherFiltNorm = filt_methDat, otherFiltNormNoNA = filt_noNA_methDat, Ylabel = namesCombined[3], otherYlabel = "DNA methylation")

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[3]][,1], filt_REC8_norm = filt_log2trans_list[[3]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[3]][,1], otherFiltNorm = filt_log2trans_list[[8]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[8]][,1], Ylabel = namesCombined[3], otherYlabel = namesCombined[8])
REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[3]][,1], filt_REC8_norm = filt_log2trans_list[[3]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[3]][,1], otherFiltNorm = filt_log2trans_list[[9]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[9]][,1], Ylabel = namesCombined[3], otherYlabel = namesCombined[9])

# TEs
REC8vsFeaturesGenomePlot(xplot = filt_log2trans_list[[3]][,1], filt_REC8_norm = filt_log2trans_list[[3]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[3]][,1], otherFiltNorm = filt_TEDat[,2], otherFiltNormNoNA = filt_TEDat_noNA[,1], Ylabel = namesCombined[3], otherYlabel = "TEs")

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[3]][,1], filt_REC8_norm = filt_log2trans_list[[3]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[3]][,1], otherFiltNorm = filt_log2trans_list[[10]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[10]][,1], Ylabel = namesCombined[3], otherYlabel = namesCombined[10])
REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[3]][,1], filt_REC8_norm = filt_log2trans_list[[3]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[3]][,1], otherFiltNorm = filt_log2trans_list[[11]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[11]][,1], Ylabel = namesCombined[3], otherYlabel = namesCombined[11])
REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[3]][,1], filt_REC8_norm = filt_log2trans_list[[3]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[3]][,1], otherFiltNorm = filt_log2trans_list[[12]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[12]][,1], Ylabel = namesCombined[3], otherYlabel = namesCombined[12])
REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[3]][,1], filt_REC8_norm = filt_log2trans_list[[3]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[3]][,1], otherFiltNorm = filt_log2trans_list[[13]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[13]][,1], Ylabel = namesCombined[3], otherYlabel = namesCombined[13])

# AT content
REC8vsbaseContentGenomePlot(xplot = filt_log2trans_list[[3]][,1], filt_REC8_norm = filt_log2trans_list[[3]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[3]][,1], baseContentFiltNorm = filt_ATcontent[,2], baseContentFiltNormNoNA = filt_ATcontent_noNA[,1],Ylabel = namesCombined[3], otherYlabel = "AT content (%)")

# GC content
REC8vsbaseContentGenomePlot(xplot = filt_log2trans_list[[3]][,1], filt_REC8_norm = filt_log2trans_list[[3]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[3]][,1], baseContentFiltNorm = filt_GCcontent[,2], baseContentFiltNormNoNA = filt_GCcontent_noNA[,1],Ylabel = namesCombined[3], otherYlabel = "GC content (%)")

#REC8vsATGCcontentGenomePlot(xplot = filt_log2trans_list[[3]][,1], filt_REC8_norm = filt_log2trans_list[[3]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[3]][,1], ATcontentFiltNorm = filt_ATcontent[,2], ATcontentFiltNormNoNA = filt_ATcontent_noNA[,1], GCcontentFiltNorm = filt_GCcontent[,2], GCcontentFiltNormNoNA = filt_GCcontent_noNA[,1], Ylabel = namesCombined[3], otherYlabel = "AT:GC content (%)")

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[3]][,1], filt_REC8_norm = filt_log2trans_list[[3]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[3]][,1], otherFiltNorm = filt_log2trans_list[[14]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[14]][,1], Ylabel = namesCombined[3], otherYlabel = namesCombined[14])

# COs
REC8vsFeaturesGenomePlot(xplot = filt_log2trans_list[[3]][,1], filt_REC8_norm = filt_log2trans_list[[3]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[3]][,1], otherFiltNorm = filt_CODat[,2], otherFiltNormNoNA = filt_CODat_noNA[,1], Ylabel = namesCombined[3], otherYlabel = "Crossovers")

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[3]][,1], filt_REC8_norm = filt_log2trans_list[[3]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[3]][,1], otherFiltNorm = filt_log2trans_list[[15]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[15]][,1], Ylabel = namesCombined[3], otherYlabel = namesCombined[15])

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[3]][,1], filt_REC8_norm = filt_log2trans_list[[3]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[3]][,1], otherFiltNorm = filt_log2trans_list[[16]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[16]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[16])

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[3]][,1], filt_REC8_norm = filt_log2trans_list[[3]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[3]][,1], otherFiltNorm = filt_log2trans_list[[17]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[17]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[17])

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[3]][,1], filt_REC8_norm = filt_log2trans_list[[3]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[3]][,1], otherFiltNorm = filt_log2trans_list[[18]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[18]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[18])

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[3]][,1], filt_REC8_norm = filt_log2trans_list[[3]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[3]][,1], otherFiltNorm = filt_log2trans_list[[19]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[19]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[19])

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[3]][,1], filt_REC8_norm = filt_log2trans_list[[3]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[3]][,1], otherFiltNorm = filt_log2trans_list[[20]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[20]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[20])

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[3]][,1], filt_REC8_norm = filt_log2trans_list[[3]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[3]][,1], otherFiltNorm = filt_log2trans_list[[21]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[21]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[21])

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[3]][,1], filt_REC8_norm = filt_log2trans_list[[3]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[3]][,1], otherFiltNorm = filt_log2trans_list[[22]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[22]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[22])

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[3]][,1], filt_REC8_norm = filt_log2trans_list[[3]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[3]][,1], otherFiltNorm = filt_log2trans_list[[23]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[23]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[23])

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[3]][,1], filt_REC8_norm = filt_log2trans_list[[3]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[3]][,1], otherFiltNorm = filt_log2trans_list[[24]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[24]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[24])

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[3]][,1], filt_REC8_norm = filt_log2trans_list[[3]][,2], filt_REC8_norm_noNA = filt_noNA_log2trans_list[[3]][,1], otherFiltNorm = filt_log2trans_list[[25]][,2], otherFiltNormNoNA = filt_noNA_log2trans_list[[25]][,1], Ylabel = namesCombined[1], otherYlabel = namesCombined[25])

dev.off()



library(corrplot)

# Create FILTERED, log2-transformed genome-wide Spearman's rho correlation matrix
allDF <- data.frame(filt_log2trans_list[[1]][,2], filt_log2trans_list[[2]][,2], filt_log2trans_list[[3]][,2])
colnames(allDF) <- c("REC8_HA_Rep1", "REC8_HA_Rep2", "REC8_MYC_Rep1")
allDF_corMat <- cor(allDF, method = "spearman", use = "pairwise.complete.obs")
col1 <- colorRampPalette(c("red", "white", "blue"))
pdf(file = paste0(plotDir, "REC8_HA_Rep1_Rep2_MYC_Rep1_filtered_genome-wide_correlation_matrix_", winNames[1], "_colouronly_MYC_Rep2_input_v290618.pdf"))
corrplot(allDF_corMat, method = "color", type = "upper", col = col1(20), tl.col = "black",
         addgrid.col = "white", addCoef.col = "grey90", mar = c(0,0,1,0), tl.srt = 0, tl.cex = 1, cl.cex = 0.8, number.cex = 1,
         title = paste0("Filtered genome-wide Spearman correlation matrix (10-kb windows)"))
dev.off()

