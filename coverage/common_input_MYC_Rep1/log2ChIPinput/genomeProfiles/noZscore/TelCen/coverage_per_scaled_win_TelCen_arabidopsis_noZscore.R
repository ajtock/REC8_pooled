#!/applications/R/R-3.5.0/bin/Rscript

# Usage:
# ./coverage_per_scaled_win_TelCen_arabidopsis_noZscore.R ASY1_Rep2_ChIP_REC8_MYC_Rep1_input ASY1 10kb 10000 100 100ths red 

#libName <- "ASY1_Rep2_ChIP_REC8_MYC_Rep1_input"
#libNamePlot <- "ASY1"
#winName <- "10kb"
#winSize <- 10000
#prop <- 100
#propName <- "100ths" 
#profileColour <- "red"

args <- commandArgs(trailingOnly = T)
libName <- args[1]
libNamePlot <- args[2]
winName <- args[3]
winSize <- as.numeric(args[4])
prop <- as.numeric(args[5])
propName <- args[6]
profileColour <- args[7]

library(parallel)
library(data.table)
library(plyr)
library(GenomicRanges)
library(segmentSeq)

makeTransparent <- function(thisColour, alpha = 150)
{
  newColour <- col2rgb(thisColour)
  apply(newColour, 2, function(x) {
    rgb(red = x[1], green = x[2], blue = x[3],
        alpha = alpha, maxColorValue = 255)
  })
}
profileColour <- sapply(seq_along(profileColour), function(x) {
  makeTransparent(profileColour[x])
})

# Genomic definitions
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromere <- c(15086045, 3607929, 13587786, 3956021, 11725024)

# Windowed coverage profile
profileLog2 <- read.table(paste0("../filt_", libName, "_genome_norm_coverage_",
                                 winName, "_noZscore.tsv"),
                          header = TRUE)
colnames(profileLog2) <- c("chr", "window", "cumwindow", "ChIP", "input", "log2ChIPinput")

chrProfiles <- mclapply(seq_along(chrs), function(x) {
  profileLog2[profileLog2$chr == chrs[x],]
}, mc.cores = length(chrs))

# Calculate average coverage in proportionally scaled
# windows along left and right chromosome arms
TelCenMatrix <- NULL
for(x in seq_along(chrs)) {
  print(chrs[x])
  # Left arm
  seqWindowStartsL <- as.integer(seq(from = 1,
                                     to = centromere[x],
                                     by = centromere[x]/prop))
  # Check that last window start coordinate is
  # < centromere[x]-1
  stopifnot(seqWindowStartsL[length(seqWindowStartsL)] < centromere[x]-1)
  seqWindowEndsL <- c(seqWindowStartsL[2:length(seqWindowStartsL)]-1,
                      centromere[x]-1)
  windowGRangesL <- GRanges(seqnames = chrs[x],
                            ranges = IRanges(start = seqWindowStartsL,
                                             end = seqWindowEndsL),
                            strand = "*")
  print(windowGRangesL)
  # Check that last window start coordinate == centromere[x]-1
  # (i.e., most proximal coordinate)
  stopifnot(end(windowGRangesL[length(windowGRangesL)]) == centromere[x]-1)
  seqWindowStartsR <- as.integer(seq(from = centromere[x],
                                     to = chrLens[x],
                                     by = ((chrLens[x]-centromere[x])+1)/prop))
  # Check that last window start coordinate is
  # < chrLens[x]
  stopifnot(seqWindowStartsR[length(seqWindowStartsR)] < chrLens[x])

  # Right arm
  seqWindowEndsR <- c(seqWindowStartsR[2:length(seqWindowStartsR)]-1,
                      chrLens[x])
  windowGRangesR <- rev(GRanges(seqnames = chrs[x],
                                 ranges = IRanges(start = seqWindowStartsR,
                                                  end = seqWindowEndsR),
                                 strand = "*"))
  print(windowGRangesR)
  # Check that first window start coordinate == chrLens[x]
  # (i.e., most distal coordinate)
  stopifnot(end(windowGRangesR[1]) == chrLens[x])

  # Create GRanges object of winName windows
  # corresponding to winName windowed coverage values
  covWindowStarts <- chrProfiles[[x]]$window
  covWindowEnds <- c((chrProfiles[[x]]$window[2:length(chrProfiles[[x]]$window)])-1,
                      chrLens[x])
  covWindowGRanges <- GRanges(seqnames = chrs[x],
                              ranges = IRanges(start = covWindowStarts,
                                               end = covWindowEnds),
                              strand = "*")

  # Calculate mean log2(markChIPA/markControlA) in each
  # scaled window using coverage values in winName windows
  overlapsL <- getOverlaps(coordinates = windowGRangesL,
                           segments = covWindowGRanges,
                           overlapType = "overlapping",
                           whichOverlaps = TRUE)
  overlapsR <- getOverlaps(coordinates = windowGRangesR,
                           segments = covWindowGRanges,
                           overlapType = "overlapping",
                           whichOverlaps = TRUE)
  scaledWinAvgCovL <- sapply(overlapsL, function(y) {
                        mean(chrProfiles[[x]]$log2ChIPinput[y])
                      })
  scaledWinAvgCovR <- sapply(overlapsR, function(y) {
                        mean(chrProfiles[[x]]$log2ChIPinput[y])
                      })
  scaledWinAvgCovLR <- cbind(scaledWinAvgCovL,
                             scaledWinAvgCovR)
  TelCenMatrix <- cbind(TelCenMatrix, scaledWinAvgCovLR)
  #scaledWinAvgCovMeanLR <- sapply(seq_along(scaledWinAvgCovL), function(y) {
  #                           mean(c(scaledWinAvgCovL[y], scaledWinAvgCovR[y]))
  #                         })
  #scaledWinAvgCov <- scaledWinAvgCovL+scaledWinAvgCovR
  #TelCenProfile <- TelCenProfile+scaledWinAvgCov
}
write.table(TelCenMatrix,
            file = paste0("log2_",
                          libName, "_",
                          winName, "_",
                          propName, "_TelCenMatrix_noZscore.txt"))

# Load TelCenMatrix
TelCenDF <- read.table(paste0("log2_",
                              libName, "_",
                              winName, "_",
                              propName, "_TelCenMatrix_noZscore.txt"))
TelCenProfile <- as.vector(rowMeans(TelCenDF))

# Function to plot telomere to centromere (Tel-Cen)
# profile of log2(markChIPA/markControlA) 
TelCenPlot <- function(xplot,
                       profile,
                       proportions,
                       proportionsName,
                       profileColour,
                       Ylabel,
                       Ylim,
                       legendLabs,
                       legendLoc) {
  plot(xplot, profile, type = "h", lwd = 4.9, col = profileColour,
       ylim = Ylim,
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = "")
  axis(side = 1, cex.axis = 1, lwd.tick = 1.5,
       at = c(1, seq(10, proportions, by = 10)),
       labels = c(expression(italic("TEL")),
                  seq(10, proportions-10, by = 10),
                  expression(italic("CEN"))))
  mtext(side = 1, line = 2.1, cex = 1,
        text = paste0("Scaled windows (", proportionsName, ")"))
  axis(side = 2, at = pretty(c(profile)), cex.axis = 1, lwd.tick = 1.5)
  mtext(side = 2, line = 2.1, cex = 1, text = Ylabel, col = "black")
  abline(h = 0, lwd = 1.5, lty = 1)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c(profileColour),
         text.col = c(profileColour),
         text.font = c(1),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
}

pdf(paste0("./log2_",
           libName, "_",
           winName, "_",
           propName, "_TelCenProfile_noZscore.pdf"),
    height = 3.5, width = 7)
par(mfrow = c(1, 1))
par(mar = c(3.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
TelCenPlot(xplot = 1:length(TelCenProfile),
           profile = TelCenProfile,
           proportions = prop,
           proportionsName = propName,
           profileColour = profileColour,
           Ylabel = bquote("Log"[2]*"(ChIP/input)"),
           Ylim = c(min(TelCenProfile),
                    max(TelCenProfile)),
           legendLabs = libNamePlot,
           legendLoc = "top")
dev.off()
