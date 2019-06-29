#!/applications/R/R-3.5.0/bin/Rscript

###################################################################################
# Generate genome profiles of feature frequency in adjacent windows of given size # 
###################################################################################

# Usage:
# ./gene_frequency_perwin.R 10000 10kb 101

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
genes <- read.table("/projects/ajt200/TAIR10/representative_genes/representative_genes_uniq_fmt_strand.txt",
                    header = T)
genesGR <- GRanges(seqnames = paste0("Chr", genes$chr),
                   ranges = IRanges(start = genes$start,
                                    end = genes$end),
                   strand = genes$strand,
                   gene_model = genes$gene_model)
print("******genesGR******")
print(genesGR)

cumWinGenesFreq <- NULL
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
  chrGenes <- genesGR[seqnames(genesGR) == chrs[i]]    
  winGenes <- countOverlaps(query = windowsGRanges,
                            subject = chrGenes,
                            ignore.strand = T)
  genesDF <- data.frame(chr = seqnames(windowsGRanges),
                        cumWindows = as.integer(cumWindows),
                        genes = as.numeric(winGenes))
  cumWinGenesFreq <- rbind(cumWinGenesFreq, genesDF)
}
write.table(cumWinGenesFreq,
            file = paste0(outDir, "gene_frequency_genome_", winName, ".txt"))

chrGenesDat <- lapply(seq_along(chrs), function(x) {
  cumWinGenesFreq[cumWinGenesFreq$chr == chrs[x],]
})

# Smooth feature frequency values with moving-average filter
# Calculate moving average of current window, ((N/2)-0.5) previous windows,
# and ((N/2)-0.5) subsequent windows
# (the higher N is, the greater the smoothing)
stopifnot(N %% 2 != 0)
flank <- (N/2)-0.5
# Define MA filter coefficients
f <- rep(1/N, N)

filt_GenesDat <- NULL
for(i in seq_along(chrs)) {
  filt_chrGenesDat <- stats::filter(x = chrGenesDat[[i]]$genes,
                                    filter = f,
                                    sides = 2)
  filt_chrGenesDat[1:flank] <- filt_chrGenesDat[flank+1]
  filt_chrGenesDat[(length(filt_chrGenesDat)-flank+1):length(filt_chrGenesDat)] <- filt_chrGenesDat[(length(filt_chrGenesDat)-flank)]
  filt_chrGenesDatDF <- data.frame(chr = chrGenesDat[[i]]$chr,
                                   cumWindows = as.integer(chrGenesDat[[i]]$cumWindows),
                                   filt_genes = as.numeric(filt_chrGenesDat))
  filt_GenesDat <- rbind(filt_GenesDat, filt_chrGenesDatDF)
}
write.table(filt_GenesDat, file = paste0(outDir, "filt_gene_frequency_genome_", winName, ".txt"))

