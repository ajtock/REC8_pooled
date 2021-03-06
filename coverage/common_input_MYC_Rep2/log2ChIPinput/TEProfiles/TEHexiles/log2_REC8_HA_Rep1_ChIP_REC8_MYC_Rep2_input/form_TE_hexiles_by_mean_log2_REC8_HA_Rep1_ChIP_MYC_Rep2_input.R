#################################################################
# select TEs in TAIR10 annotation for hexile formation based    #
# on mean log2-normalised coverage levels between start and end #   
#################################################################

library(GenomicAlignments)
library(segmentSeq)

allTEs <- read.table("/projects/ajt200/TAIR10/TAIR10_Buisine_TEs_strand_tab_ann.txt",
                     header = T)
print(dim(allTEs))
#[1] 31189     7

##############################################################
# mean log2-normalised coverage levels between start and end #
##############################################################

chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)

hexDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/TEProfiles/TEHexiles/log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input/"
covDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/"

REC8_norm <- read.table(file = paste0(covDir, "log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab.bed"))
allCovTEs <- NULL
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
  chrTEs <- allTEs[which(allTEs[,1] == chrs[i]),]
  TEsIRcoords <- IRanges(start = chrTEs[,2],
                         end = chrTEs[,3])
  TEsGRcoords <- GRanges(seqnames = chrs[i],
                         strand = "+",
                         ranges = TEsIRcoords)
  print(i)
  TEsOverlaps <- getOverlaps(TEsGRcoords,
                             covGRcoords,
                             whichOverlaps = T)
  REC8_norm_TECov <- sapply(TEsOverlaps,
                            function(x) mean(chrCov[x]))
  print(i)
  chrTEs <- cbind(chrTEs, REC8_norm_TECov)
  allCovTEs <- rbind(allCovTEs, chrTEs)
}
write.table(allCovTEs,
            file = paste0(hexDir,
                          "log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_TEMeanCov.txt"))

###########################################################
# Rank TEs into hexiles based on log2-normalised coverage #
###########################################################

print(head(allCovTEs))
allCovTEs <- allCovTEs[order(allCovTEs$REC8_norm_TECov,
                             decreasing = T),]
print(head(allCovTEs))
print(tail(allCovTEs))
hex <- round(length(allCovTEs[,1])/6)
hex <- cumsum(c(1, rep(hex, times = 6)))
hex <- c(hex, length(allCovTEs[,1]))
hexDat <- NULL
for(j in 1:length(hex)) {
  print(j)
  if(j == 1) {
    hexjTEs <- allCovTEs[hex[j]:hex[j+1]-1,]
    print(paste0("condition 1: hexile ", j))
    print(dim(hexjTEs))
  }
  if(j > 1 & j < 6) {
    hexjTEs <- allCovTEs[(hex[j]):(hex[j+1]-1),]
    print(paste0("condition 2: hexile ", j))
    print(dim(hexjTEs))
  }
  if(j == 6) {
    hexjTEs <- allCovTEs[(hex[j]):(hex[length(hex)]),]
    print(paste0("condition 3: hexile ", j))
    print(dim(hexjTEs))
  } 
  if(j <= 6) {    
    dat <- c(j, length(hexjTEs[,1]),
             round(mean(hexjTEs[,3]-hexjTEs[,2])),
             sum(hexjTEs[,3]-hexjTEs[,2]),
             mean(hexjTEs[,8]))
    hexDat <- rbind(hexDat, dat)
    write.table(hexjTEs,
                file = paste0(hexDir,
                              "hexile_", j,
                              "_log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_TECov.txt"))
  }
}
colnames(hexDat) <- c("hexile",
                      "n",
                      "mean_bp",
                      "total_bp", "REC8_norm_mean_cov")
write.table(hexDat,
            file = paste0(hexDir,
                          "log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_TECov_hexile_stats.txt"),
            row.names = F)


