#!/applications/R/R-3.5.0/bin/Rscript

###################################################################################
# Generate genome profiles of feature frequency in adjacent windows of given size # 
###################################################################################

# Usage:
# ./CO_frequency_perwin.R 10000 10kb 201

args <- commandArgs(trailingOnly = T)
winSize <- as.numeric(args[1])
winName <- as.character(args[2])
N <- as.numeric(args[3])

library(GenomicRanges)
library(parallel)

outDir <- "./"
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
# pericentromeric regions are as defined in Supplemental Table S26 of Ziolkowski et al. (2017) Genes Dev. 31
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)

# Genomic definitions with cumulative coordinates
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

# Count number of features in each adjacent genomic windows using countOverlaps()
#COs <- read.csv("/projects/ajt200/GBS_CO/BethRowan_190219/BR_wt_ColLer_F2/CO.all.with.420.flagged.csv")
#COsGR <- GRanges(seqnames = paste0("Chr", COs$chr),
#                 ranges = IRanges(start = COs$block1.end.pos,
#                                  end = COs$block2.start.pos),
#                 strand = "*")
#print(length(COsGR))
load("/projects/ajt200/GBS_CO/HS_CU_080617/wt/COsGRcoords.RData")
COsGR <- COsGRcoords
print(length(COsGR))

cumWinCOsFreq <- NULL
for(i in 1:5) {
  # Define adjacent windows
  seqWindows <- seq(1, chrLens[i], by = winSize)
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
  # Count features within windows
  chrCOs <- COsGR[seqnames(COsGR) == chrs[i]]    
  winCOs <- countOverlaps(query = windowsGRanges,
                          subject = chrCOs,
                          ignore.strand = T)
  COsDF <- data.frame(chr = seqnames(windowsGRanges),
                      cumWindows = as.integer(cumWindows),
                      COs = as.numeric(winCOs))
  cumWinCOsFreq <- rbind(cumWinCOsFreq, COsDF)
}
write.table(cumWinCOsFreq,
            file = paste0(outDir, "CO_frequency_genome_", winName, ".txt"))

chrCOsDat <- lapply(seq_along(chrs), function(x) {
  cumWinCOsFreq[cumWinCOsFreq$chr == chrs[x],]
})

# Smooth feature frequency values with moving-average filter
# Calculate moving average of current window, ((N/2)-0.5) previous windows,
# and ((N/2)-0.5) subsequent windows
# (the higher N is, the greater the smoothing)
stopifnot(N %% 2 != 0)
flank <- (N/2)-0.5
# Define MA filter coefficients
f <- rep(1/N, N)

filt_COsDat <- NULL
for(i in seq_along(chrs)) {
  filt_chrCOsDat <- stats::filter(x = chrCOsDat[[i]]$COs,
                                  filter = f,
                                  sides = 2)
  filt_chrCOsDat[1:flank] <- filt_chrCOsDat[flank+1]
  filt_chrCOsDat[(length(filt_chrCOsDat)-flank+1):length(filt_chrCOsDat)] <- filt_chrCOsDat[(length(filt_chrCOsDat)-flank)]
  filt_chrCOsDatDF <- data.frame(chr = chrCOsDat[[i]]$chr,
                                 cumWindows = as.integer(chrCOsDat[[i]]$cumWindows),
                                 filt_COs = as.numeric(filt_chrCOsDat))
  filt_COsDat <- rbind(filt_COsDat, filt_chrCOsDatDF)
}
write.table(filt_COsDat, file = paste0(outDir, "filt_CO_frequency_genome_", winName, ".txt"))

