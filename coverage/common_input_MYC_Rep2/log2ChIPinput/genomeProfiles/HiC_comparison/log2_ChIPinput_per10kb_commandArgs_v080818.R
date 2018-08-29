#!/applications/R/R-3.3.2/bin/Rscript

###########################################################################
# Calculate log2 ratio of ChIP and input coverage values in 10-kb windows #
# for generating chromosome-scale plots                                   #
###########################################################################

# Usage on hydrogen node7
# csmit -m 20G -c 1 "./log2_ChIPinput_per10kb_commandArgs_v080818.R REC8_HA_Rep1_ChIP REC8_MYC_Rep2_input 10000 10kb"

library(segmentSeq)

args <- commandArgs(trailingOnly = TRUE)
ChIPname <- args[1]
inputname <- args[2]
winSize <- as.numeric(args[3])
winName <- as.character(args[4])

ChIPcovFile <- paste0("../../../../", ChIPname,
                      "_norm_allchrs_coverage_coord_tab.bed")
inputcovFile <- paste0("../../../../", inputname,
                       "_norm_allchrs_coverage_coord_tab.bed")
inputnameChIP <- paste0(inputname, "_", ChIPname)
outDir <- "../../../log2ChIPinput/genomeProfiles/HiC_comparison/" 

ChIPcovDat <- read.table(ChIPcovFile,
                         colClasses = c(NA, rep("NULL", 2), NA))
inputcovDat <- read.table(inputcovFile, 
                          colClasses = c(NA, rep("NULL", 2), NA))

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

# Calculate mean coverage values within windows and log2(ChIP/input) ratios
print(winName)
cumChIPwinDat <- NULL
cuminputwinDat <- NULL
for(i in 1:5) {
  print(i)
  chrChIP <- ChIPcovDat[ChIPcovDat[,1] == chrs[i],]
  chrinput <- inputcovDat[inputcovDat[,1] == chrs[i],]
  ChIPcov <- chrChIP[,2]
  inputcov <- chrinput[,2]
  covCoords <- seq(1, length(ChIPcov), by = 1)
  covIRcoords <- IRanges(start = covCoords, width = 1)
  covGRcoords <- GRanges(seqnames = chrs[i], strand = "*",
                         ranges = covIRcoords)
  seqWindows <- seq(1, chrLens[i], by = winSize)
#  seqWindows <- c(seqWindows, chrLens[i])
  cumWindows <- seqWindows + sumchr[i]
  windowsIRanges <- IRanges(start = seqWindows,
                            width = winSize)
  windowsIRanges <- windowsIRanges[-length(windowsIRanges)]
  windowsIRanges <- append(windowsIRanges,
                           IRanges(start = seqWindows[length(seqWindows)],
                                   end = chrLens[i]))
  windowsGRanges <- GRanges(seqnames = chrs[i], strand = "*",
                            ranges = windowsIRanges)
  print(windowsGRanges)
  overlaps <- getOverlaps(windowsGRanges, covGRcoords,
                          whichOverlaps = TRUE)
  #+1 is offset to avoid infinite values
  ChIPcovWinVals <- sapply(overlaps, function(x) mean(ChIPcov[x])+1)
  inputcovWinVals <- sapply(overlaps, function(x) mean(inputcov[x])+1)
  ChIPwinDat <- cbind(cumWindows, ChIPcovWinVals)
  inputwinDat <- cbind(cumWindows, inputcovWinVals)
  cumChIPwinDat <- rbind(cumChIPwinDat, ChIPwinDat)
  cuminputwinDat <- rbind(cuminputwinDat, inputwinDat)
  write.table(ChIPwinDat,
              file = paste0(outDir, ChIPname, "_chr", i, "_norm_coverage_",
                            winName, ".txt"))
  write.table(inputwinDat,
              file = paste0(outDir, inputnameChIP, "_chr", i, "_norm_coverage_",
                            winName, ".txt"))
  log2trans <- log2(ChIPwinDat[,2]/inputwinDat[,2])
  log2trans <- (log2trans-mean(log2trans, na.rm = T))/sd(log2trans, na.rm = T)
  log2trans <- cbind(ChIPwinDat[,1], log2trans)
  colnames(log2trans) <- c("cumWindows", "log2covWinVals")
  write.table(log2trans,
              file = paste0(outDir, "log2_", ChIPname, "_chr", i, "_norm_coverage_",
                            winName, ".txt"))
}
write.table(cumChIPwinDat,
            file = paste0(outDir, ChIPname, "_genome_norm_coverage_",
                          winName, ".txt"))
write.table(cuminputwinDat,
            file = paste0(outDir, inputnameChIP, "_genome_norm_coverage_",
                          winName, ".txt"))
cumlog2trans <- log2(cumChIPwinDat[,2]/cuminputwinDat[,2])
cumlog2trans <- (cumlog2trans-mean(cumlog2trans, na.rm = T))/sd(cumlog2trans, na.rm = T)
cumlog2trans <- cbind(cumChIPwinDat[,1], cumlog2trans)
colnames(cumlog2trans) <- c("cumWindows", "log2covWinVals")
write.table(cumlog2trans,
            file = paste0(outDir, "log2_", ChIPname, "_genome_norm_coverage_",
                          winName, ".txt"))

chrlog2trans <- lapply(1:5, function(i) {
  read.table(file = paste0(outDir, "log2_", ChIPname, "_chr", i,
                           "_norm_coverage_", winName, ".txt"))
})

# smooth coverage values with MA filter
test <- seq(1, 1000, by = 1)
j = 100
ma <- rep(1, test[j])/test[j]

filt_log2trans <- NULL
filt_log2trans_noNA <- NULL  
for(i in 1:5) {
  filt_chrlog2trans <- stats::filter(chrlog2trans[i][[1]][,2], ma)
  which_na <- which(is.na(filt_chrlog2trans) == TRUE)
  left_na <- which_na[which(which_na < 100)]
  left_val <- filt_chrlog2trans[left_na[length(left_na)]+1]
  filt_chrlog2trans[left_na] <- left_val
  right_na <- which_na[which(which_na > 100)]
  right_val <- filt_chrlog2trans[right_na[1]-1]
  filt_chrlog2trans[right_na] <- right_val
  filt_chrlog2trans_noNA <- filt_chrlog2trans[!is.na(filt_chrlog2trans)]
  filt_chrlog2trans <- cbind(chrlog2trans[i][[1]][,1], filt_chrlog2trans)
  write.table(filt_chrlog2trans,
              file = paste0(outDir, "filt_log2_", ChIPname, "_chr", i,
                            "_norm_coverage_", winName, ".txt"))
  write.table(filt_chrlog2trans_noNA,
              file = paste0(outDir, "filt_noNA_log2_", ChIPname, "_chr", i,
                            "_norm_coverage_", winName, ".txt"))
  filt_log2trans <- rbind(filt_log2trans, filt_chrlog2trans)
  filt_log2trans_noNA <- c(filt_log2trans_noNA, filt_chrlog2trans_noNA)
}
write.table(filt_log2trans,
            file = paste0(outDir, "filt_log2_", ChIPname,
                          "_genome_norm_coverage_", winName, ".txt"))
write.table(filt_log2trans_noNA,
            file = paste0(outDir, "filt_noNA_log2_", ChIPname,
                          "_genome_norm_coverage_", winName, ".txt"))

