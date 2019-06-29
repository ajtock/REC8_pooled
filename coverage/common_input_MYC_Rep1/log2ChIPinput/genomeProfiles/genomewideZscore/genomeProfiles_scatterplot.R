#!/applications/R/R-3.5.0/bin/Rscript

# Usage:
#./genomeProfiles_scatterplot.R REC8_HA_Rep2_ChIP REC8_MYC_Rep1_ChIP 'REC8-HA' 'REC8-Myc' 10kb 20

#lib1Name <- "REC8_HA_Rep2_ChIP"
#lib2Name <- "REC8_MYC_Rep1_ChIP"
#lib1NamePlot <- "REC8-HA"
#lib2NamePlot <- "REC8-Myc"
#winName <- "10kb"
#covMax <- 20

args <- commandArgs(trailingOnly = T)
lib1Name <- args[1]
lib2Name <- args[2]
lib1NamePlot <- args[3]
lib2NamePlot <- args[4]
winName <- args[5]
covMax <- as.numeric(args[6])

lib1 <- read.table(paste0(lib1Name,
                          "_genome_norm_coverage_",
                          winName, ".txt"))
lib2 <- read.table(paste0(lib2Name,
                          "_genome_norm_coverage_",
                          winName, ".txt"))
lib1[lib1$coverage > covMax,]$coverage <- NA
lib2[lib2$coverage > covMax,]$coverage <- NA

pdf(paste0("./plots/", lib1Name, "_vs_", lib2Name,
           "_genomeProfiles_scatterplot_", winName,
           "_covMax", as.character(covMax), ".pdf"),
    height = 5, width = 5)
par(mfcol = c(1, 1))
par(mar = c(4.1, 4.1, 4.1, 4.1))
par(mgp = c(3, 1, 0))
plot(lib1$coverage, lib2$coverage,
     type = "p", pch = 16, cex = 0.25, col = "grey30",
     xlim = c(0, covMax),
     ylim = c(0, covMax),
     xlab = "", ylab = "",
     main = bquote(italic("r"[s]) ~ " = " ~
                   .(round(cor(lib1$coverage,
                               lib2$coverage,
                               method = "spearman",
                               use = "pairwise.complete.obs"),
                           digits = 2))))
mtext(side = 1, line = 2.30, cex = 1.3,
      text = lib1NamePlot)
mtext(side = 2, line = 2.30, cex = 1.3,
      text = lib2NamePlot)
#lines(x = c(min(lib1$coverage),
#            max(lib1$coverage)),
#      y = c(min(lib2$coverage), 
#            max(lib2$coverage)),
#      type = "l",
#      lty = 5,
#      lwd = 2,
#      col = "red")
abline(lm(lib2$coverage~lib1$coverage),
       lty = 5, lwd = 2, col = "red")
abline(v = mean(lib1$coverage, na.rm = T),
       lty = 5, lwd = 2, col = "blue")
abline(h = mean(lib2$coverage, na.rm = T),
       lty = 5, lwd = 2, col = "blue")
box(lwd = 2)
dev.off()
