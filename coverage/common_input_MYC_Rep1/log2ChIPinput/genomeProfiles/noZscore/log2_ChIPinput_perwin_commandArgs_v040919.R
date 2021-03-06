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
library(GenomicRanges)
library(parallel)

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
# smoothN must be an odd number
# see profile smoothing below
smoothN <- as.numeric(args[7])
 
ChIPfile <- paste0(ChIPDir, ChIPname, "_norm_allchrs_coverage_coord_tab.bed")
inputfile <- paste0(inputDir, inputname, "_norm_allchrs_coverage_coord_tab.bed")

ChIP <- read.table(ChIPfile,
                   colClasses = c(NA, rep("NULL", 2), NA))
input <- read.table(inputfile,
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
winDF <- NULL
for(i in seq_along(chrs)) {
  print(chrs[i])
  # Subset coverage to chrs[i]
  ChIPchr <- ChIP[ChIP[,1] == chrs[i],]
  inputchr <- input[input[,1] == chrs[i],]
  ChIPcov <- ChIPchr[,2]
  inputcov <- inputchr[,2]
  # Define coverage coordinates
  covLocGR <- GRanges(seqnames = chrs[i],
                      ranges = IRanges(start = seq(from = 1,
                                                   to = length(ChIPcov),
                                                   by = 1),
                                       width = 1),
                      strand = "*")
  print(covLocGR)
  # Define adjacent windows
  winSeq <- seq(from = 1, to = chrLens[i], by = winSize)
  winCum <- winSeq + sumchr[i]
  winIR <- IRanges(start = winSeq,
                   width = winSize)
  winIR <- winIR[-length(winIR)]
  winIR <- append(winIR,
                  IRanges(start = winSeq[length(winSeq)],
                          end = chrLens[i]))
  winGR <- GRanges(seqnames = chrs[i],
                   ranges = winIR,
                   strand = "*")
  print(winGR)
  # Identify overlapping windows and coverage coordinates
  #fOverlaps <- findOverlaps(query = winGR,
  #                          subject = covLocGR,
  #                          type = "any",
  #                          select = "all",
  #                          ignore.strand = TRUE)
  ## Convert fOverlaps into list object equivalent to that
  ## generated by segmentSeq::getOverlaps(), in which list
  ## element corresponds to a genomic window of
  ## sequentially numbered coverage coordinates
  #fOverlapsList <- mclapply(seq_along(unique(queryHits(fOverlaps))),
  #                 function(x) {
  #                   subjectHits(fOverlaps[queryHits(fOverlaps) == x])
  #                 }, mc.cores = detectCores())
  gOverlaps <- getOverlaps(coordinates = winGR,
                           segments = covLocGR,
                           overlapType = "overlapping",
                           whichOverlaps = TRUE,
                           ignoreStrand = TRUE)
  ChIPwinCov <- sapply(gOverlaps, function(x) {
                  mean(ChIPcov[x])
                })
  inputwinCov <- sapply(gOverlaps, function(x) {
                  mean(inputcov[x])
                })
  # Calculate log2(ChIP+1/input+1) ratios;
  # +1 is offset to avoid infinite values
  log2ChIPinput <- log2((ChIPwinCov+1)/(inputwinCov+1)) 
  # Combine in a dataframe
  winDFchr <- data.frame(chr = as.character(chrs[i]),
                         window = as.integer(start(winGR)),
                         cumwindow = as.integer(start(winGR) + sumchr[i]),
                         ChIP = as.numeric(ChIPwinCov),
                         input = as.numeric(inputwinCov),
                         log2ChIPinput = as.numeric(log2ChIPinput))
  winDF <- rbind(winDF, winDFchr)
}
colnames(winDF) <- c("chr", "window", "cumwindow",
                     ChIPname, inputname,
                     paste0("log2_", ChIPname, "_", inputname))
write.table(winDF,
            file = paste0(ChIPname, "_", inputname,
                          "_genome_norm_coverage_",
                          winName, "_noZscore.tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T) 

# Create list in which each list element is a dataframe
# of windowed coverage values for a chromosome to to
# enable per-chromosome smoothing
winDFchrList <- lapply(seq_along(chrs), function(i) {
  winDF[winDF$chr == chrs[i],]
})

# Calculate moving average (MA) of current window,
# (smoothN/2)-0.5 previous windows (where smoothN is odd),
# and
# (smoothN/2)-0.5 subsequent windows (where smoothN is odd)
# The higher smoothN is, the greater the smoothing
# Use modulo operation (%%) to confirm that the remainder of
# integer division smoothN %/% 2 is 1 (i.e., smoothN is odd)
stopifnot(smoothN %% 2 == 1)
flank <- smoothN %/% 2
# Define MA filter coefficients
f <- rep(x = 1/smoothN, times = smoothN)

filt_winDF <- NULL
for(i in seq_along(chrs)) {
  print(chrs[i])
  # ChIP
  filt_ChIP <- stats::filter(winDFchrList[[i]][,4],
                             filter = f,
                             sides = 2)
  filt_ChIP[1:flank] <- filt_ChIP[flank+1]
  filt_ChIP[ (length(filt_ChIP)-flank+1) :
              length(filt_ChIP) ] <- filt_ChIP[ (length(filt_ChIP)-flank) ]
  # input
  filt_input <- stats::filter(winDFchrList[[i]][,5],
                              filter = f,
                              sides = 2)
  filt_input[1:flank] <- filt_input[flank+1]
  filt_input[ (length(filt_input)-flank+1) :
               length(filt_input) ] <- filt_input[ (length(filt_input)-flank) ]
  # log2
  filt_log2 <- stats::filter(winDFchrList[[i]][,6],
                             filter = f,
                             sides = 2)
  filt_log2[1:flank] <- filt_log2[flank+1]
  filt_log2[ (length(filt_log2)-flank+1) :
              length(filt_log2) ] <- filt_log2[ (length(filt_log2)-flank) ]
  # Combine in dataframe
  filt_winDFchr <- data.frame(chr = as.character(chrs[i]),
                              window = as.integer(winDFchrList[[i]]$window),
                              cumwindow = as.integer(winDFchrList[[i]]$cumwindow),
                              filt_ChIP = as.numeric(filt_ChIP),
                              filt_input = as.numeric(filt_input),
                              filt_log2ChIPinput = as.numeric(filt_log2))
  filt_winDF <- rbind(filt_winDF, filt_winDFchr)
}
colnames(filt_winDF) <- c("chr", "window", "cumwindow",
                          paste0("filt_", ChIPname),
                          paste0("filt_", inputname),
                          paste0("filt_log2_", ChIPname, "_", inputname))
write.table(filt_winDF,
            file = paste0("filt_", ChIPname, "_", inputname,
                          "_genome_norm_coverage_",
                          winName, "_noZscore.tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
