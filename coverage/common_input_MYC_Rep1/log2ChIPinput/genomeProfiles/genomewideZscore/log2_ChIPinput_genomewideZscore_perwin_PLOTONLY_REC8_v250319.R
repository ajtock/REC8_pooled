#!/applications/R/R-3.5.0/bin/Rscript

# Generate genome-wide plots of chromosome profiles

# Usage:
# Rscript ./log2_ChIPinput_genomewideZscore_perwin_PLOTONLY_REC8_v250319.R 10kb

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
REC8HAprefix <- "REC8_HA_Rep1_ChIP_REC8_MYC_Rep1_input"
REC8HAname <- "REC8-HA"
 
REC8MYCDir <- "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/genomeProfiles/genomewideZscore/"
REC8MYCprefix <- "REC8_MYC_Rep1_ChIP_REC8_MYC_Rep1_input"
REC8MYCname <- "REC8-Myc"

controlHADir <- "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/genomeProfiles/genomewideZscore/"
controlHAprefix <- "control_HA_Rep1_ChIP_REC8_MYC_Rep1_input"
controlHAname <- "control-HA"

controlMYCDir <- "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/genomeProfiles/genomewideZscore/"
controlMYCprefix <- "control_MYC_Rep1_ChIP_REC8_MYC_Rep1_input"
controlMYCname <- "control-Myc"

# REC8
inDirREC8 <- c(REC8HADir, REC8MYCDir, controlHADir, controlMYCDir)
prefixREC8 <- c(REC8HAprefix, REC8MYCprefix, controlHAprefix, controlMYCprefix)
nameREC8 <- c(REC8HAname, REC8MYCname, controlHAname, controlMYCname)

filt_REC8log2trans_list <- mclapply(seq_along(nameREC8), function(x) {
  read.table(paste0(inDirREC8[x], "filt_log2_", prefixREC8[x],
                    "_genome_norm_coverage_Zscore_", winName, ".txt"))
}, mc.cores = length(nameREC8))

# Function to plot REC8 genome-scale coverage (all replicates) 
REC8GenomePlot <- function(xplot,
                           dat1,
                           dat2,
                           dat3,
                           dat4) {
  plot(xplot, dat1, type = "l", lwd = 1.5, col = "red",
       ylim = c(min(dat1, dat2, dat3, dat4), max(dat1, dat2, dat3, dat4)),
       xlab = "",
       ylab = "",
       xaxt = "n")
#       main = bquote(italic("r"[s]) ~ " = " ~ .(round(cor(dat1, dat2, method = "spearman"), digits = 2))))
  axis(side = 1, lwd.tick = 1.5,
       at = c(0, 2e7, 4e7, 6e7, 8e7, 1e8, 1.2e8),
       labels = c("0", "20", "40", "60", "80", "100", "120"), cex.axis = 1)
  axis(side = 2, lwd.tick = 1.5)
  lines(xplot, dat2, lwd = 1.5, col = "blue")
  lines(xplot, dat3, lwd = 1.5, col = "magenta3")
  lines(xplot, dat4, lwd = 1.5, col = "green2")
  mtext(side = 1, line = 2.25, cex = 1.3, text = "Coordinates (Mb)")
  mtext(side = 2, line = 2.25, cex = 1.3, text = expression("Log"[2]*"(ChIP/REC8-Myc Rep1 input)"))
  abline(v = sumchr, lty = 1, lwd = 0.75)
  abline(v = centromeres, lty = 5, lwd = 0.75)
  box(lwd = 1.5)
  legend("topleft",
         legend = c("REC8-HA", "REC8-Myc", "Control-HA", "Control-Myc"),
         col = c("red", "blue", "magenta3", "green2"),
         text.col = c("red", "blue", "magenta3", "green2"),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
}


## Fig. S5A
pdf(paste0(plotDir, "genomewideZscore_log2_REC8_control_ChIP_REC8_MYC_Rep1_input_", winName, "_genomeProfiles_FigS5A.pdf"),
    height = 4, width = 10)
par(mfcol = c(1, 1))
par(mar = c(5.1, 5.1, 2.1, 2.1))
par(mgp = c(3, 1, 0))
# Fig. S5A (Ylab: Log[2](ChIP/input))
# REC8-HA Rep1
# REC8-Myc Rep1 
# Control-Myc Rep1
REC8GenomePlot(xplot = filt_REC8log2trans_list[[1]]$cumWindows,
               dat1  = filt_REC8log2trans_list[[1]]$filt_ZscoreLog2cov,
               dat2  = filt_REC8log2trans_list[[2]]$filt_ZscoreLog2cov,
               dat3  = filt_REC8log2trans_list[[3]]$filt_ZscoreLog2cov,
               dat4  = filt_REC8log2trans_list[[4]]$filt_ZscoreLog2cov)
dev.off()

