############################################################################################################################
# select peaks for hexile formation based on mean log2-normalised coverage levels between start and end coordinates        #   
############################################################################################################################

library(GenomicAlignments)
library(segmentSeq)

inDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/"

load(paste0(inDir,
            "REC8_HA_Rep1_armrangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData"))
peaksGR <- armrangerPeaksGRmergedOverlaps
allPeaks <- data.frame(peaksGR)[,-4]
print(dim(allPeaks))

##############################################################
# mean log2-normalised coverage levels between start and end #
##############################################################

chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)

hexDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/arms/peakHexiles/log2_REC8_HA_Rep1_ChIP_MYC_Rep2_input/"
covDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/"
REC8_norm <- read.table(file = paste0(covDir,
                                      "log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab.bed"))
allCovPeaks <- NULL
for(i in 1:5) {
  print(i)
  chrCov <- REC8_norm[REC8_norm[,1] == chrs[i],]
  chrCov <- chrCov[,4]
  print(i)
  covCoords <- seq(1, length(chrCov), by = 1)
  covIRcoords <- IRanges(start = covCoords, width = 1)
  covGRcoords <- GRanges(seqnames = chrs[i],
                         strand = "+",
                         ranges = covIRcoords)
  print(i)
  chrPeaks <- allPeaks[which(allPeaks[,1] == chrs[i]),]
  peaksIRcoords <- IRanges(start = chrPeaks[,2],
                           end = chrPeaks[,3])
  peaksGRcoords <- GRanges(seqnames = chrs[i],
                           strand = "+",
                           ranges = peaksIRcoords)
  print(i)
  peaksOverlaps <- getOverlaps(peaksGRcoords,
                               covGRcoords,
                               whichOverlaps = T)
  REC8_norm_peakCov <- sapply(peaksOverlaps,
                              function(x) mean(chrCov[x]))
  print(i)
  chrPeaks <- cbind(chrPeaks, REC8_norm_peakCov)
  allCovPeaks <- rbind(allCovPeaks, chrPeaks)
}
write.table(allCovPeaks,
            file = paste0(hexDir,
                          "log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_peakMeanCov.txt"))

#############################################################
# Rank peaks into hexiles based on log2-normalised coverage #
#############################################################

print(head(allCovPeaks))
allCovPeaks <- allCovPeaks[order(allCovPeaks$REC8_norm_peakCov,
                                 decreasing = T),]
print(head(allCovPeaks))
print(tail(allCovPeaks))
hex <- round(length(allCovPeaks[,1])/6)
hex <- cumsum(c(1, rep(hex, times = 6)))
hex <- c(hex, length(allCovPeaks[,1]))
hexDat <- NULL
for(j in 1:length(hex)) {
  print(j)
  if(j == 1) {
    hexjPeaks <- allCovPeaks[hex[j]:hex[j+1]-1,]
    print(paste0("condition 1: hexile ", j))
    print(dim(hexjPeaks))
  }
  if(j > 1 & j < 6) {
    hexjPeaks <- allCovPeaks[(hex[j]):(hex[j+1]-1),]
    print(paste0("condition 2: hexile ", j))
    print(dim(hexjPeaks))
  }
  if(j == 6) {
    hexjPeaks <- allCovPeaks[(hex[j]):(hex[length(hex)]),]
    print(paste0("condition 3: hexile ", j))
    print(dim(hexjPeaks))
  } 
  if(j <= 6) {    
    dat <- c(j, length(hexjPeaks[,1]),
             round(mean(hexjPeaks[,3]-hexjPeaks[,2])),
             sum(hexjPeaks[,3]-hexjPeaks[,2]),
             mean(hexjPeaks[,5]))
    hexDat <- rbind(hexDat, dat)
    write.table(hexjPeaks,
                file = paste0(hexDir,
                              "hexile_", j,
                               "_log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_peakCov.txt"))
  }
}
colnames(hexDat) <- c("hexile",
                      "n",
                      "mean_bp",
                      "total_bp",
                      "REC8_norm_mean_cov")
write.table(hexDat,
            file = paste0(hexDir,
                          "log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_peakCov_hexile_stats.txt"),
            row.names = F)


print(sessionInfo())

