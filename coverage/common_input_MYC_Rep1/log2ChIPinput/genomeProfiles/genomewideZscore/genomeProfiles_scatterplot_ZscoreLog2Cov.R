#!/applications/R/R-3.5.0/bin/Rscript

# Usage:
#./genomeProfiles_scatterplot_ZscoreLog2Cov.R filt_log2_REC8_HA_Rep2_ChIP_REC8_MYC_Rep1_input filt_log2_REC8_MYC_Rep1_ChIP_REC8_MYC_Rep1_input 'REC8-HA' 'REC8-Myc' 10kb

#lib1Name <- "filt_log2_REC8_HA_Rep2_ChIP_REC8_MYC_Rep1_input"
#lib2Name <- "filt_log2_REC8_MYC_Rep1_ChIP_REC8_MYC_Rep1_input"
#lib1NamePlot <- "REC8-HA"
#lib2NamePlot <- "REC8-Myc"
#winName <- "10kb"

args <- commandArgs(trailingOnly = T)
lib1Name <- args[1]
lib2Name <- args[2]
lib1NamePlot <- args[3]
lib2NamePlot <- args[4]
winName <- args[5]

lib1 <- read.table(paste0(lib1Name,
                          "_genome_norm_coverage_Zscore_",
                          winName, ".txt"))
lib2 <- read.table(paste0(lib2Name,
                          "_genome_norm_coverage_Zscore_",
                          winName, ".txt"))

pdf(paste0("./plots/", lib1Name, "_vs_", lib2Name,
           "_genomeProfiles_scatterplot_", winName, ".pdf"),
    height = 5, width = 5)
par(mfcol = c(1, 1))
par(mar = c(4.1, 4.1, 4.1, 4.1))
par(mgp = c(3, 1, 0))
plot(lib1$filt_ZscoreLog2cov, lib2$filt_ZscoreLog2cov,
     type = "p", pch = 16, cex = 0.25, col = "grey30",
     xlim = c(min(lib1$filt_ZscoreLog2cov),
              max(lib1$filt_ZscoreLog2cov)),
     ylim = c(min(lib2$filt_ZscoreLog2cov),
              max(lib2$filt_ZscoreLog2cov)),
     xlab = "", ylab = "",
     main = bquote(italic("r"[s]) ~ " = " ~
                   .(round(cor(lib1$filt_ZscoreLog2cov,
                               lib2$filt_ZscoreLog2cov,
                               method = "spearman",
                               use = "pairwise.complete.obs"),
                           digits = 2))))
mtext(side = 1, line = 2.30, cex = 1.3,
      text = lib1NamePlot)
mtext(side = 2, line = 2.30, cex = 1.3,
      text = lib2NamePlot)
#lines(x = c(min(lib1$filt_ZscoreLog2cov),
#            max(lib1$filt_ZscoreLog2cov)),
#      y = c(min(lib2$filt_ZscoreLog2cov), 
#            max(lib2$filt_ZscoreLog2cov)),
#      type = "l",
#      lty = 5,
#      lwd = 2,
#      col = "red")
abline(lm(lib2$filt_ZscoreLog2cov~lib1$filt_ZscoreLog2cov),
       lty = 5, lwd = 2, col = "red")
abline(v = mean(lib1$filt_ZscoreLog2cov),
       lty = 5, lwd = 2, col = "blue")
abline(h = mean(lib2$filt_ZscoreLog2cov),
       lty = 5, lwd = 2, col = "blue")
box(lwd = 2)
dev.off()
