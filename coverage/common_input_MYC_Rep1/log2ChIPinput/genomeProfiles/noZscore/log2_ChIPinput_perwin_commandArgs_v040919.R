#!/applications/R/R-3.5.0/bin/Rscript

###########################################################
###########################################################
# Calculate log2 ratio of ChIP to input coverage
# values in windows of defined size for generating
# chromosome-scale plots and for calculating correlations
###########################################################
###########################################################

# Usage on hydrogen node7:
# csmit -m 20G -c 1 "Rscript ./log2_ChIPinput_perwin_commandArgs_v040919.R '/home/ajt200/analysis/REC8_pooled/coverage/' '/home/ajt200/analysis/REC8_pooled/coverage/' REC8_HA_Rep2_ChIP REC8_MYC_Rep1_input 10000 10kb 101"

library(segmentSeq)

#ChIPDir <- "/home/ajt200/analysis/REC8_pooled/coverage/"
#inputDir <- "/home/ajt200/analysis/REC8_pooled/coverage/"
#ChIPname <- "REC8_HA_Rep2_ChIP"
#inputname <- "REC8_MYC_Rep1_input"
#winSize <- 10000
#winName <- "10kb"
#smoothN <- 101

args <- commandArgs(trailingOnly = TRUE)
ChIPDir <- args[1]
inputDir <- args[2]
ChIPname <- args[3]
inputname <- args[4]
winSize <- as.numeric(args[5])
winName <- as.character(args[6])
smoothN <- as.numeric(args[7])
 
ChIPfile <- paste0(ChIPDir, ChIPname, "_norm_allchrs_coverage_coord_tab.bed")
inputfile <- paste0(inputDir, inputname, "_norm_allchrs_coverage_coord_tab.bed")

ChIPcov <- read.table(ChIPfile,
                      colClasses = c(NA, rep("NULL", 2), NA))
inputcov <- read.table(inputfile,
                       colClasses = c(NA, rep("NULL", 2), NA))

# Genomic definitions
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
# Pericentromeric regions are as defined in Supplemental Table S26
# of Ziolkowski et al. (2017) Genes Dev. 31
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)

# Make chromosomal coordinates cumulative
# such that the first coordinate of Chr2 is
# equal to the last coordinate of Chr1 + 1
sumchr <- cumsum(c(0, chrLens))
print(sumchr)
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

# For each winSize-bp window, calculate the mean of the winSize
# per-base coverage values within that window, and
# calculate windowed log2(ChIP/input) ratios
print(winName)
ChIPwinDF <- NULL
inputwinDF <- NULL



