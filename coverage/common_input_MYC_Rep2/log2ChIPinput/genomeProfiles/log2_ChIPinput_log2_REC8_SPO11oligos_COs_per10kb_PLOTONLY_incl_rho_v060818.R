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

plotDir <- "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/genomeProfiles/plots/"

outDir1 <- "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/genomeProfiles/"
ChIPnames1 <- "REC8_HA_Rep1_ChIP"
names1 <- "wt REC8-HA Rep1"

outDir2 <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames2 <- "WT_SPO11_oligo_RPI1"
names2 <- "wt SPO11-1-oligos Rep1"

outDirCombined <- c(outDir1, outDir2)
ChIPnamesCombined <- c(ChIPnames1, ChIPnames2)
namesCombined <- c(names1, names2)

inDirFeatures <- "/projects/ajt200/REC8_MSH4/data_merged_fastq/coverage/log2ChIPinput/genomeProfiles/"

# Load concatenated filtered, log2-transformed coverage datasets
filt_log2trans_list <- mclapply(seq_along(ChIPnamesCombined), function(k) {
  read.table(file = paste0(outDirCombined[k], "filt_log2_", ChIPnamesCombined[k], "_genome_norm_coverage_", winNames[1], ".txt"))
}, mc.cores = length(ChIPnamesCombined))
filt_noNA_log2trans_list <- mclapply(seq_along(ChIPnamesCombined), function(k) {
  read.table(file = paste0(outDirCombined[k], "filt_noNA_log2_", ChIPnamesCombined[k], "_genome_norm_coverage_", winNames[1], ".txt"))
}, mc.cores = length(ChIPnamesCombined))

filt_COsDat <- read.table(paste0(inDirFeatures, "filt_WTCO_density_genome_10kb.txt"))
filt_noNA_COsDat <- read.table(paste0(inDirFeatures, "filt_noNA_WTCO_density_genome_10kb.txt"))

filt_log2_REC8_SPO11oligos <- log2((filt_log2trans_list[[1]][,2]+2)/(filt_log2trans_list[[2]][,2]+2))
filt_noNA_log2_REC8_SPO11oligos <- log2((filt_noNA_log2trans_list[[1]][,1]+2)/(filt_noNA_log2trans_list[[2]][,1]+2))

filt_log2_REC8_COs <- log2((filt_log2trans_list[[1]][,2]+2)/(filt_COsDat[,2]+2))
filt_noNA_log2_REC8_COs <- log2((filt_noNA_log2trans_list[[1]][,1]+2)/(filt_noNA_COsDat[,1]+2))

# Function to plot REC8 genome-scale coverage overlaid with other datasets 
REC8vsOthersGenomePlot <- function(xplot, filt_REC8_norm, filt_REC8_norm_noNA, otherFiltNorm, otherFiltNormNoNA, Ylabel, otherYlabel) {
  plot(xplot, otherFiltNorm, type = "l", lwd = 1.5, col = "blue",
       ylim = c(min(otherFiltNormNoNA), max(otherFiltNormNoNA)),
       xlab = "",
       ylab = "", xaxt = "n", yaxt = "n",
       main = bquote(italic("r"[s]) ~ " = " ~ .(round(cor(filt_REC8_norm, otherFiltNorm, method = "spearman"), digits = 2))))
  axis(side = 4, at = pretty(otherFiltNorm), lwd.tick = 1.5)
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 1, adj = c(0.5, -4.0), labels = otherYlabel, xpd = NA, srt = -90, col = "blue")
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

##################
## REC8-HA Rep1 ##
##################
pdf(file = paste0(plotDir, "log2_REC8_HA_Rep1_SPO11oligos_log2_REC8_HA_Rep1_COs_genomeplot_", winNames[1], "_rho_v060818.pdf"), height = 3.125, width = 9.375)
par(mfcol = c(1, 1))
par(mar = c(2.1, 4.1, 2.1, 4.1))
par(mgp = c(3, 1, 0))

REC8vsOthersGenomePlot(xplot = filt_log2trans_list[[1]][,1],
                       filt_REC8_norm = filt_log2_REC8_SPO11oligos,
                       filt_REC8_norm_noNA = filt_noNA_log2_REC8_SPO11oligos,
                       otherFiltNorm = filt_log2_REC8_COs,
                       otherFiltNormNoNA = filt_noNA_log2_REC8_COs,
                       Ylabel = "log2(REC8:SPO11-1-oligos)",
                       otherYlabel = "log2(REC8:Crossovers)")
dev.off()
