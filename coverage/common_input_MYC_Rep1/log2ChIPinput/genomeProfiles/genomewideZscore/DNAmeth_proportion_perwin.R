#!/applications/R/R-3.5.0/bin/Rscript

# Calculate mean DNA methylation proportions in adjacent windows
# for generating chromosome-scale plots

# Usage:
# Rscript ./DNAmeth_proportion_perwin.R '/home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/' GSM980986_WT_rep2 200000 200kb 5

args <- commandArgs(trailingOnly = T)
bedDir <- args[1]
libName <- args[2]
winSize <- as.numeric(args[3])
winName <- args[4]
N <- as.numeric(args[5])

outDir <- "./"

library(segmentSeq)

# DNA methylation
CG  <- read.table(paste0(bedDir, libName, "_CG.wig.bed.gr.tab.bed"),
                  colClasses = c(NA, NA, "NULL", NA))
CHG <- read.table(paste0(bedDir, libName, "_CHG.wig.bed.gr.tab.bed"),
                  colClasses = c(NA, NA, "NULL", NA))
CHH <- read.table(paste0(bedDir, libName, "_CHH.wig.bed.gr.tab.bed"),
                  colClasses = c(NA, NA, "NULL", NA))

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

# Calculate mean DNA methylation proportions min adjacent windows
print(winName)
cumMethDat <- NULL
for(i in 1:5) {
  # Define adjacent windows as GRanges object
  seqWindows <- seq(1, chrLens[i], by = winSize)
  cumWindows <- seqWindows + sumchr[i]
  windowsIRanges <- IRanges(start = seqWindows,
                            width = winSize)
  windowsIRanges <- windowsIRanges[-length(windowsIRanges)]
  windowsIRanges <- append(windowsIRanges,
                           IRanges(start = seqWindows[length(seqWindows)],
                                   end = chrLens[i]))
  windowsGRanges <- GRanges(seqnames = chrs[i],
                            ranges = windowsIRanges,
                            strand = "*")
  print(windowsGRanges)

  # Define DNA methylation coordinates as GRanges objects
  # and calculate mean methylation proportions in each window
  # CG
  chrCG <- CG[CG[,1] == paste0("chr", i),]
  chrCG_GRanges <- GRanges(seqnames = chrs[i],
                           ranges = IRanges(start = chrCG[,2],
                                            width = 1),
                           strand = "*")
  CGoverlaps <- getOverlaps(coordinates = windowsGRanges,
                            segments = chrCG_GRanges,
                            overlapType = "overlapping",
                            whichOverlaps = TRUE)
  CGwinVals <- sapply(CGoverlaps, function(x) mean(as.numeric(chrCG[,3][x])))
  # CHG
  chrCHG <- CHG[CHG[,1] == paste0("chr", i),]
  chrCHG_GRanges <- GRanges(seqnames = chrs[i],
                            ranges = IRanges(start = chrCHG[,2],
                                             width = 1),
                            strand = "*")
  CHGoverlaps <- getOverlaps(coordinates = windowsGRanges,
                             segments = chrCHG_GRanges,
                             overlapType = "overlapping",
                             whichOverlaps = TRUE)
  CHGwinVals <- sapply(CHGoverlaps, function(x) mean(as.numeric(chrCHG[,3][x])))
  # CHH
  chrCHH <- CHH[CHH[,1] == paste0("chr", i),]
  chrCHH_GRanges <- GRanges(seqnames = chrs[i],
                            ranges = IRanges(start = chrCHH[,2],
                                             width = 1),
                            strand = "*")
  CHHoverlaps <- getOverlaps(coordinates = windowsGRanges,
                             segments = chrCHH_GRanges,
                             overlapType = "overlapping",
                             whichOverlaps = TRUE)
  CHHwinVals <- sapply(CHHoverlaps, function(x) mean(as.numeric(chrCHH[,3][x])))
  # Mean of all 3 contexts
  CwinVals <- sapply(seq_along(CGwinVals), function(x) {
                      mean(c(CGwinVals[x], CHGwinVals[x], CHHwinVals[x]))
                    })
  # Combine in data frame
  methDat <- data.frame(chr = rep(chrs[i], times = length(CGwinVals)),
                        cumWindows = as.integer(cumWindows),
                        mCG  = CGwinVals,
                        mCHG = CHGwinVals,
                        mCHH = CHHwinVals,
                        mC   = CwinVals)
  cumMethDat <- rbind(cumMethDat, methDat)
}
write.table(cumMethDat,
            file = paste0(outDir, "DNAmeth_",
                          libName, "_genome_", winName, ".txt"))

chrMethDat <- lapply(seq_along(chrs), function(x) {
  cumMethDat[cumMethDat$chr == chrs[x],]
})

# Calculate moving average of current window,
#### (N/2) previous windows (where N is even) OR
# (N/2)-0.5 previous windows (where N is odd),
# and
#### (N/2) subsequent windows (where N is even) OR
# (N/2)-0.5 subsequent windows (where N is odd)
# (the higher N is, the greater the smoothing)
stopifnot(N %% 2 != 0)
flank <- (N/2)-0.5
# Define MA filter coefficients
f <- rep(1/N, N)

filt_methDat <- NULL
for(i in seq_along(chrs)) {
  # mCG
  filt_chr_mCG <- stats::filter(chrMethDat[[i]]$mCG,
                                filter = f,
                                sides = 2)
  filt_chr_mCG[1:flank] <- filt_chr_mCG[flank+1]
  filt_chr_mCG[(length(filt_chr_mCG)-flank+1):length(filt_chr_mCG)] <- filt_chr_mCG[(length(filt_chr_mCG)-flank)]
  # mCHG
  filt_chr_mCHG <- stats::filter(chrMethDat[[i]]$mCHG,
                                 filter = f,
                                 sides = 2)
  filt_chr_mCHG[1:flank] <- filt_chr_mCHG[flank+1]
  filt_chr_mCHG[(length(filt_chr_mCHG)-flank+1):length(filt_chr_mCHG)] <- filt_chr_mCHG[(length(filt_chr_mCHG)-flank)]
  # mCHH
  filt_chr_mCHH <- stats::filter(chrMethDat[[i]]$mCHH,
                                 filter = f,
                                 sides = 2)
  filt_chr_mCHH[1:flank] <- filt_chr_mCHH[flank+1]
  filt_chr_mCHH[(length(filt_chr_mCHH)-flank+1):length(filt_chr_mCHH)] <- filt_chr_mCHH[(length(filt_chr_mCHH)-flank)]
  # mC
  filt_chr_mC <- stats::filter(chrMethDat[[i]]$mC,
                               filter = f,
                               sides = 2)
  filt_chr_mC[1:flank] <- filt_chr_mC[flank+1]
  filt_chr_mC[(length(filt_chr_mC)-flank+1):length(filt_chr_mC)] <- filt_chr_mC[(length(filt_chr_mC)-flank)]
  # Combine in data frame
  filt_chrMethDat <- data.frame(chr = chrMethDat[[i]]$chr,
                                cumWindows = as.integer(chrMethDat[[i]]$cumWindows),
                                filt_mCG  = as.numeric(filt_chr_mCG),
                                filt_mCHG = as.numeric(filt_chr_mCHG),
                                filt_mCHH = as.numeric(filt_chr_mCHH),
                                filt_mC   = as.numeric(filt_chr_mC))
  filt_methDat <- rbind(filt_methDat, filt_chrMethDat)
}
write.table(filt_methDat,
            file = paste0(outDir, "filt_DNAmeth_",
                          libName, "_genome_", winName, ".txt"))
