#!/applications/R/R-3.5.0/bin/Rscript

###################################################################################
# Generate genome profiles of feature frequency in adjacent windows of given size # 
###################################################################################

# Usage:
# ./TE_frequency_perwin.R 10000 10kb 101

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
TEs <- read.table("/projects/ajt200/TAIR10/TAIR10_Buisine_TEs_strand_tab_ann.txt",
                  header = T)
TEsGR <- GRanges(seqnames = TEs$Chr,
                 ranges = IRanges(start = TEs$start,
                                  end = TEs$end),
                 strand = TEs$strand,
                 TE_model = TEs$Transposon_Name)
print("******TEsGR******")
print(TEsGR)

cumWinTEsFreq <- NULL
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
  chrTEs <- TEsGR[seqnames(TEsGR) == chrs[i]]    
  winTEs <- countOverlaps(query = windowsGRanges,
                          subject = chrTEs,
                          ignore.strand = T)
  TEsDF <- data.frame(chr = seqnames(windowsGRanges),
                      cumWindows = as.integer(cumWindows),
                      TEs = as.numeric(winTEs))
  cumWinTEsFreq <- rbind(cumWinTEsFreq, TEsDF)
}
write.table(cumWinTEsFreq,
            file = paste0(outDir, "TE_frequency_genome_", winName, ".txt"))

chrTEsDat <- lapply(seq_along(chrs), function(x) {
  cumWinTEsFreq[cumWinTEsFreq$chr == chrs[x],]
})

# Smooth feature frequency values with moving-average filter
# Calculate moving average of current window, ((N/2)-0.5) previous windows,
# and ((N/2)-0.5) subsequent windows
# (the higher N is, the greater the smoothing)
stopifnot(N %% 2 != 0)
flank <- (N/2)-0.5
# Define MA filter coefficients
f <- rep(1/N, N)

filt_TEsDat <- NULL
for(i in seq_along(chrs)) {
  filt_chrTEsDat <- stats::filter(x = chrTEsDat[[i]]$TEs,
                                  filter = f,
                                  sides = 2)
  filt_chrTEsDat[1:flank] <- filt_chrTEsDat[flank+1]
  filt_chrTEsDat[(length(filt_chrTEsDat)-flank+1):length(filt_chrTEsDat)] <- filt_chrTEsDat[(length(filt_chrTEsDat)-flank)]
  filt_chrTEsDatDF <- data.frame(chr = chrTEsDat[[i]]$chr,
                                 cumWindows = as.integer(chrTEsDat[[i]]$cumWindows),
                                 filt_TEs = as.numeric(filt_chrTEsDat))
  filt_TEsDat <- rbind(filt_TEsDat, filt_chrTEsDatDF)
}
write.table(filt_TEsDat, file = paste0(outDir, "filt_TE_frequency_genome_", winName, ".txt"))

