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

outDir2 <- "/home/ajt200/analysis/180622_Chris_lambing_ChIP_REC8_HA_Col_kss/WT/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames2 <- "REC8_HA_Rep2_ChIP" 
names2 <- "wt REC8-HA Rep2"

outDir3 <- "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/genomeProfiles/"
ChIPnames3 <- "REC8_MYC_Rep1_ChIP"
names3 <- "wt REC8-MYC Rep1"

outDir4 <- "/home/ajt200/analysis/180622_Chris_lambing_ChIP_REC8_HA_Col_kss/kss/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames4 <- "kss_REC8_HA_Rep1_ChIP"
names4 <- "kss REC8-HA Rep1"

outDir5 <- "/projects/ajt200/BAM_masters/H3K9me2/WT/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames5 <- "WT_H3K9me2_ChIP"
names5 <- "wt H3K9me2"

outDir6 <- "/projects/ajt200/BAM_masters/H3K9me2/kss/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames6 <- "kss_H3K9me2_ChIP"
names6 <- "kss H3K9me2"

outDir7 <- "/projects/ajt200/BAM_masters/H3K9me2/cmt3/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames7 <- "cmt3_H3K9me2_ChIP"
names7 <- "cmt3 H3K9me2"

outDir8 <- "/home/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/WT/coverage/genomeProfiles/"
ChIPnames8 <- "WT_RNAseq_Chris_Rep1"
names8 <- "wt RNA-seq (Chris) Rep1"

outDir9 <- "/home/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/WT/coverage/genomeProfiles/"
ChIPnames9 <- "WT_RNAseq_Chris_Rep2"
names9 <- "wt RNA-seq (Chris) Rep2"

outDir10 <- "/home/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/kss/coverage/genomeProfiles/"
ChIPnames10 <- "kss_RNAseq_Chris_Rep1"
names10 <- "kss RNA-seq (Chris) Rep1"

outDir11 <- "/home/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/kss/coverage/genomeProfiles/"
ChIPnames11 <- "kss_RNAseq_Chris_Rep2"
names11 <- "kss RNA-seq (Chris) Rep2"

outDir12 <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames12 <- "WT_SPO11_oligo_RPI1"
names12 <- "wt SPO11-1-oligos Rep1"

outDir13 <- "/projects/ajt200/BAM_masters/SPO11-oligo/suvh456/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames13 <- "suvh456_SPO11_oligo_RPI34"
names13 <- "kss SPO11-1-oligos Rep1"

outDir14 <- "/projects/ajt200/BAM_masters/SPO11-oligo/suvh456/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames14 <- "suvh456_SPO11_oligo_RPI35"
names14 <- "kss SPO11-1-oligos Rep2"

outDirCombined <- c(outDir1, outDir2, outDir3, outDir4, outDir5, outDir6, outDir7, outDir8, outDir9, outDir10, outDir11, outDir12, outDir13, outDir14)
ChIPnamesCombined <- c(ChIPnames1, ChIPnames2, ChIPnames3, ChIPnames4, ChIPnames5, ChIPnames6, ChIPnames7, ChIPnames8, ChIPnames9, ChIPnames10, ChIPnames11, ChIPnames12, ChIPnames13, ChIPnames14)
namesCombined <- c(names1, names2, names3, names4, names5, names6, names7, names8, names9, names10, names11, names12, names13, names14)

inDirFeatures <- "/projects/ajt200/REC8_MSH4/data_merged_fastq/coverage/log2ChIPinput/genomeProfiles/"

# Load concatenated filtered, log2-transformed coverage datasets
filt_log2trans_list <- mclapply(seq_along(ChIPnamesCombined), function(k) {
  read.table(file = paste0(outDirCombined[k], "filt_log2_", ChIPnamesCombined[k], "_genome_norm_coverage_", winNames[1], ".txt"))
}, mc.cores = length(ChIPnamesCombined))
filt_noNA_log2trans_list <- mclapply(seq_along(ChIPnamesCombined), function(k) {
  read.table(file = paste0(outDirCombined[k], "filt_noNA_log2_", ChIPnamesCombined[k], "_genome_norm_coverage_", winNames[1], ".txt"))
}, mc.cores = length(ChIPnamesCombined))

# Load concatenated filtered DNA meth datasets
filt_wtMethDat <- read.table(file = paste0(inDirFeatures, "filt_wtMeth_genome_200kb.txt"))
filt_noNA_wtMethDat <- read.table(file = paste0(inDirFeatures, "filt_noNA_wtMeth_genome_200kb.txt"))
# Load concatenated filtered DNA meth datasets
filt_kssMethDat <- read.table(file = paste0(inDirFeatures, "filt_kssMeth_genome_200kb.txt"))
filt_noNA_kssMethDat <- read.table(file = paste0(inDirFeatures, "filt_noNA_kssMeth_genome_200kb.txt"))

# Function to plot genome-scale coverage overlaid with alternate genotype
genoVsOtherGenoGenomePlot <- function(xplot, filt_norm, filt_norm_noNA, filt_other_norm, filt_other_norm_noNA, Ylabel, otherYlabel, wtCol, mutCol) {
  plot(xplot, filt_norm, type = "l", lwd = 1.5, col = wtCol,
       ylim = c(min(filt_norm_noNA, filt_other_norm_noNA), max(filt_norm_noNA, filt_other_norm_noNA)),
       xlab = "",
       ylab = "", xaxt = "n", yaxt = "n",
       main = bquote(italic("r"[s]) ~ " = " ~ .(round(cor(filt_norm, filt_other_norm, method = "spearman"), digits = 2))))
  mtext(side = 2, line = 2.25, cex = 1, text = Ylabel, col = wtCol)
  axis(side = 2, lwd.tick = 1.5)
  lines(xplot, filt_other_norm, lwd = 1.5, col = mutCol)
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 1.5, adj = c(0.5, -2.5), labels = otherYlabel, xpd = NA, srt = -90, col = mutCol)
  axis(side = 1, lwd.tick = 1.5)
  abline(v = sumchr, lty = 1, lwd = 0.75)
  abline(v = centromeres, lty = 5, lwd = 0.75, col = "black")
  box(lwd = 1.5)
}

##################
## REC8-HA Rep1 ##
##################
pdf(file = paste0(plotDir, "REC8_HA_Rep1_H3K9me2_RNAseq_SPO11oligos_DNAmeth_wt_vs_kss_log2_zscore_genomeplot_", winNames[1], "_rho_MYC_Rep2_input_v310718_talk.pdf"), height = 17.5, width = 7.5)
par(mfcol = c(7, 1))
par(mar = c(2.1, 4.1, 2.1, 4.1))
par(mgp = c(3, 1, 0))

genoVsOtherGenoGenomePlot(xplot = filt_log2trans_list[[1]][,1],
                          filt_norm = filt_log2trans_list[[1]][,2],
                          filt_norm_noNA = filt_noNA_log2trans_list[[1]][,1],
                          filt_other_norm = filt_log2trans_list[[4]][,2],
                          filt_other_norm_noNA = filt_noNA_log2trans_list[[4]][,1],
                          Ylabel = namesCombined[1],
                          otherYlabel = namesCombined[4],
                          wtCol = "red", mutCol = "red4")
genoVsOtherGenoGenomePlot(xplot = filt_log2trans_list[[5]][,1],
                          filt_norm = filt_log2trans_list[[5]][,2],
                          filt_norm_noNA = filt_noNA_log2trans_list[[5]][,1],
                          filt_other_norm = filt_log2trans_list[[6]][,2],
                          filt_other_norm_noNA = filt_noNA_log2trans_list[[6]][,1],
                          Ylabel = namesCombined[5],
                          otherYlabel = namesCombined[6],
                          wtCol = "deepskyblue", mutCol = "deepskyblue4")
genoVsOtherGenoGenomePlot(xplot = filt_wtMethDat[,1],
                          filt_norm = filt_wtMethDat[,2],
                          filt_norm_noNA = filt_noNA_wtMethDat[,2],
                          filt_other_norm = filt_kssMethDat[,2],
                          filt_other_norm_noNA = filt_noNA_kssMethDat[,2],
                          Ylabel = "wt CG methylation",
                          otherYlabel = "kss CG methylation",
                          wtCol = "yellow3", mutCol = "yellow4")
genoVsOtherGenoGenomePlot(xplot = filt_wtMethDat[,1],
                          filt_norm = filt_wtMethDat[,3],
                          filt_norm_noNA = filt_noNA_wtMethDat[,3],
                          filt_other_norm = filt_kssMethDat[,3],
                          filt_other_norm_noNA = filt_noNA_kssMethDat[,3],
                          Ylabel = "wt CHG methylation",
                          otherYlabel = "kss CHG methylation",
                          wtCol = "green", mutCol = "darkgreen")
genoVsOtherGenoGenomePlot(xplot = filt_wtMethDat[,1],
                          filt_norm = filt_wtMethDat[,4],
                          filt_norm_noNA = filt_noNA_wtMethDat[,4],
                          filt_other_norm = filt_kssMethDat[,4],
                          filt_other_norm_noNA = filt_noNA_kssMethDat[,4],
                          Ylabel = "wt CHH methylation",
                          otherYlabel = "kss CHH methylation",
                          wtCol = "orange", mutCol = "orange4")
genoVsOtherGenoGenomePlot(xplot = filt_log2trans_list[[12]][,1],
                          filt_norm = filt_log2trans_list[[12]][,2],
                          filt_norm_noNA = filt_noNA_log2trans_list[[12]][,1],
                          filt_other_norm = filt_log2trans_list[[13]][,2],
                          filt_other_norm_noNA = filt_noNA_log2trans_list[[13]][,1],
                          Ylabel = namesCombined[12],
                          otherYlabel = namesCombined[13],
                          wtCol = "dodgerblue", mutCol = "navy")
genoVsOtherGenoGenomePlot(xplot = filt_log2trans_list[[8]][,1],
                          filt_norm = filt_log2trans_list[[8]][,2],
                          filt_norm_noNA = filt_noNA_log2trans_list[[8]][,1],
                          filt_other_norm = filt_log2trans_list[[10]][,2],
                          filt_other_norm_noNA = filt_noNA_log2trans_list[[10]][,1],
                          Ylabel = namesCombined[8],
                          otherYlabel = namesCombined[10],
                          wtCol = "magenta", mutCol = "magenta4")
dev.off()
