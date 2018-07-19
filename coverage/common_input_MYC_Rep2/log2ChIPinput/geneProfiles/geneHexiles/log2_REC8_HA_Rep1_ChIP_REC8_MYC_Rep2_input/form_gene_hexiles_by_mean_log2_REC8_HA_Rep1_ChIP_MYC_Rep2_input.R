################################################################
# Select genes in TAIR10 annotation for hexile formation based #
# on mean log2-normalised coverage levels between TSS and TTS  #
################################################################

library(GenomicAlignments)
library(segmentSeq)

allGenes <- read.table("/projects/ajt200/TAIR10/representative_genes/representative_genes_uniq_fmt_strand.txt",
                       header = T)
allGenes <- cbind(chr = paste0("Chr",
                  allGenes[,1]),
                  allGenes[,-1])
print(dim(allGenes))
#[1] 27204     5

############################################################
# Mean log2-normalised coverage levels between TSS and TTS #
############################################################

chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)

hexDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/geneProfiles/geneHexiles/log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input/"
covDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/"

REC8_norm <- read.table(file = paste0(covDir, "log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab.bed"))
allCovGenes <- NULL
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
  chrGenes <- allGenes[which(allGenes[,1] == chrs[i]),]
  genesIRcoords <- IRanges(start = chrGenes[,2],
                           end = chrGenes[,3])
  genesGRcoords <- GRanges(seqnames = chrs[i],
                           strand = "+",
                           ranges = genesIRcoords)
  print(i)
  genesOverlaps <- getOverlaps(genesGRcoords,
                               covGRcoords,
                               whichOverlaps = T)
  REC8_norm_geneCov <- sapply(genesOverlaps,
                              function(x) mean(chrCov[x]))
  print(i)
  chrGenes <- cbind(chrGenes, REC8_norm_geneCov)
  allCovGenes <- rbind(allCovGenes, chrGenes)
}
write.table(allCovGenes,
            file = paste0(hexDir,
                          "log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_geneMeanCov.txt"))

#############################################################
# Rank genes into hexiles based on log2-normalised coverage #
#############################################################

print(head(allCovGenes))
allCovGenes <- allCovGenes[order(allCovGenes$REC8_norm_geneCov,
                                 decreasing = T),]
print(head(allCovGenes))
print(tail(allCovGenes))
hex <- round(length(allCovGenes[,1])/6)
hex <- cumsum(c(1, rep(hex, times = 6)))
hex <- c(hex, length(allCovGenes[,1]))
hexDat <- NULL
for(j in 1:length(hex)) {
  print(j)
  if(j == 1) {
    hexjGenes <- allCovGenes[hex[j]:hex[j+1]-1,]
    print(paste0("condition 1: hexile ", j))
    print(dim(hexjGenes))
  }
  if(j > 1 & j < 6) {
    hexjGenes <- allCovGenes[(hex[j]):(hex[j+1]-1),]
    print(paste0("condition 2: hexile ", j))
    print(dim(hexjGenes))
  }
  if(j == 6) {
    hexjGenes <- allCovGenes[(hex[j]):(hex[length(hex)]),]
    print(paste0("condition 3: hexile ", j))
    print(dim(hexjGenes))
  } 
  if(j <= 6) {    
    dat <- c(j, length(hexjGenes[,1]),
             round(mean(hexjGenes[,3]-hexjGenes[,2])),
             sum(hexjGenes[,3]-hexjGenes[,2]),
             mean(hexjGenes[,6]))
    hexDat <- rbind(hexDat, dat)
    write.table(hexjGenes,
                file = paste0(hexDir,
                              "hexile_", j,
                              "_log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_geneCov.txt"))
  }
}
colnames(hexDat) <- c("hexile",
                      "n",
                      "mean_bp",
                      "total_bp",
                      "REC8_norm_mean_cov")
write.table(hexDat,
            file = paste0(hexDir,
                          "log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_geneCov_hexile_stats.txt"),
            row.names = F)

