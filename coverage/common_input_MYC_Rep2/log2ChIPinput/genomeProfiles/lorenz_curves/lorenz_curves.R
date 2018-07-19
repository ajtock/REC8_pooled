# Plot Lorenz curves for each ChIP-seq and MNase-seq dataset

library(ineq)
library(RColorBrewer)

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

windows <- c(10000)
winNames <- c("10kb")

plotDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/genomeProfiles/lorenz_curves/"

outDir1 <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/genomeProfiles/"
ChIPnames1 <- "REC8_HA_Rep1_ChIP"
names1 <- "REC8-HA"

# Euchromatic marks
outDir2 <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K4me1/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames2 <- "WT_H3K4me1_Rep1_ChIP"
names2 <- "H3K4me1"

outDir3 <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K4me2/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames3 <- "WT_H3K4me2_Rep1_ChIP"
names3 <- "H3K4me2"

outDir4 <- "/projects/ajt200/BAM_masters/H3K4me3/replicates/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames4 <- "WT_H3K4me3_ChIP14"
names4 <- "H3K4me3"

outDir5 <- "/projects/ajt200/BAM_masters/H2A/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames5 <- "H2AZ_ChIP"
names5 <- "H2A.Z"

outDir6 <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K27me3/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames6 <- "WT_H3K27me3_Rep1_ChIP"
names6 <- "H3K27me3"

# Heterochromatic marks
outDir7 <- "/projects/ajt200/BAM_masters/H2A/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames7 <- "H2AW_ChIP"
names7 <- "H2A.W"

outDir8 <- "/projects/ajt200/BAM_masters/H3K9me2/WT/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames8 <- "WT_H3K9me2_ChIP"
names8 <- "H3K9me2"

outDir9 <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K27me1/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames9 <- "WT_H3K27me1_Rep1_ChIP"
names9 <- "H3K27me1"

outDir10 <- "/projects/ajt200/BAM_masters/nucleosomes/WT/coverage/nakedDNA_untrimmed_input/log2ChIPinput/genomeProfiles/"
ChIPnames10 <- "WT_nuc"
names10 <- "Nucleosomes"

outDirCombined <- c(outDir1, outDir2, outDir3, outDir4, outDir5,
                    outDir6, outDir7, outDir8, outDir9, outDir10)
ChIPnamesCombined <- c(ChIPnames1, ChIPnames2, ChIPnames3, ChIPnames4, ChIPnames5,
                       ChIPnames6, ChIPnames7, ChIPnames8, ChIPnames9, ChIPnames10)
namesCombined <- c(names1, names2, names3, names4, names5,
                   names6, names7, names8, names9, names10)

# Load concatenated filtered, log2-transformed coverage datasets
# and add 2.18 for representation as Lorenz curves
filt_log2trans_list <- lapply(seq_along(ChIPnamesCombined), function(k) {
  read.table(file = paste0(outDirCombined[k],
                           "filt_log2_",
                           ChIPnamesCombined[k],
                           "_genome_norm_coverage_",
                           winNames[1],
                           ".txt"))$filt_chrlog2trans+2.18
})
#filt_noNA_log2trans_list <- lapply(seq_along(ChIPnamesCombined), function(k) {
#  read.table(file = paste0(outDirCombined[k], "filt_noNA_log2_", ChIPnamesCombined[k], "_genome_norm_coverage_", winNames[1], ".txt"))
#})

mycols <- c(
            # REC8
            "green",
            # Euchromatin
            "red4", "red3", "red", "lightcoral", "lightpink",
            # Heterochromatin
            "navy", "blue", "lightskyblue",
            # MNase
            "darkorange"
           )
#mycols <- c("green",
#            colorRampPalette(c("red4", "lightcoral"))(5),
#            colorRampPalette(c("navy", "lightskyblue1"))(3),
#            "darkorange") 
#mycols <- brewer.pal(n = 10, name = "Paired")
#mycols[1:2] <- c("navy", "#A6CEE3")
pdf(paste0(plotDir, "Lorenz_curves_log2_coverage_winSize_", winNames[1], ".pdf"))
plot(Lc(filt_log2trans_list[[1]]), col = mycols[1], lwd = 2, main = "10-kb windows", xlab = "Cumulative proportion of log2-transformed library-size-normalized coverage + 2.18", ylab = "Cumulative proportion of physical map")
for(i in 2:length(filt_log2trans_list)) {
  lines(Lc(filt_log2trans_list[[i]]), col = mycols[i], lwd = 2)
}
legend("topleft",
       legend = namesCombined,
       col = mycols,
       text.col = mycols,
       ncol = 1, cex = 1, lwd = 1.5, bty = "n")
box(lwd = 1.5)
dev.off()


ChIP_10kb_list <- lapply(seq_along(ChIPnamesCombined), function(k) {
  read.table(file = paste0(outDirCombined[k],
                           ChIPnamesCombined[k],
                           "_genome_norm_coverage_",
                           winNames[1],
                           ".txt"))$ChIPcovWinVals
})

pdf(paste0(plotDir, "Lorenz_curves_coverage_winSize_", winNames[1], ".pdf"))
plot(Lc(ChIP_10kb_list[[1]]), col = mycols[1], lwd = 2, main = "10-kb windows", xlab = "Cumulative proportion of library-size-normalized coverage", ylab = "Cumulative proportion of physical map")
for(i in 2:length(ChIP_10kb_list)) {
  lines(Lc(ChIP_10kb_list[[i]]), col = mycols[i], lwd = 2)
}
legend("topleft",
       legend = namesCombined,
       col = mycols,
       text.col = mycols,
       ncol = 1, cex = 1, lwd = 1.5, bty = "n")
box(lwd = 1.5)
dev.off()


# Load concatenated GBS crossover density dataset
# HOWEVER, GBS CROSSOVER DATA TOO SPARSE

inDirFeatures <- "/projects/ajt200/REC8_MSH4/data_merged_fastq/coverage/log2ChIPinput/genomeProfiles/"

CODat <- read.table(file = paste0(inDirFeatures, "WTCO_density_genome_", winNames[1], ".txt"))
filt_CODat <- read.table(file = paste0(inDirFeatures, "filt_WTCO_density_genome_", winNames[1], ".txt"))

mycols <- c("red")
pdf(paste0(plotDir, "Lorenz_curve_GBS_crossovers_wins_10kb.pdf"))
plot(Lc(CODat[,2]), col = mycols[1], lwd = 2, main = "10-kb windows", xlab = "Cumulative proportion of GBS crossovers", ylab = "Cumulative proportion of physical map")
box(lwd = 1.5)
dev.off()

mycols <- c("red")
pdf(paste0(plotDir, "Lorenz_curve_filtered_GBS_crossovers_wins_10kb.pdf"))
plot(Lc(filt_CODat[,2]), col = mycols[1], lwd = 2, main = "10-kb windows", xlab = "Cumulative proportion of filtered GBS crossovers", ylab = "Cumulative proportion of physical map")
box(lwd = 1.5)
dev.off()



## DOESN'T WORK:
#mycols <- brewer.pal(n = 10, name = "Paired")
#pdf(paste0(plotDir, "test.pdf"))
#windows <- filt_log2trans_list[[1]][,1]
#cov <- filt_log2trans_list[[1]][,2]
#cumul_wins <- cumsum(windows*cov)/max(cumsum(windows*cov))
#cumul_cov <- cumsum(cov)/max(cumsum(cov))
#plot(cumul_cov, cumul_wins, type = "l", col = mycols[1], lwd = 2, main = "10-kb windows", xlab = "Cumulative proportion of coverage", ylab = "Cumulative proportion of physical map")
#for(i in 2:length(filt_log2trans_list)) {
#  windows <- filt_log2trans_list[[i]][,1]
#  cov <- filt_log2trans_list[[i]][,2]
#  cumul_wins <- cumsum(windows*cov)/max(cumsum(windows*cov))
#  cumul_cov <- cumsum(cov)/max(cumsum(cov))
#  lines(cumul_cov, cumul_wins, type = "l", col = mycols[i], lwd = 2)
#}
#legend("topleft",
#       legend = namesCombined,
#       col = mycols,
#       text.col = mycols,
#       ncol = 1, cex = 1, lwd = 1.5, bty = "n")
#box(lwd = 1.5)
#dev.off()

