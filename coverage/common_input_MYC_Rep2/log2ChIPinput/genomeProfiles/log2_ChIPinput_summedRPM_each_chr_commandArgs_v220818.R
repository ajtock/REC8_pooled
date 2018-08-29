#!/applications/R/R-3.4.0/bin/Rscript

# Calculate summed RPM for each chromosome and take log2(ChIP/input) ratio

# Usage on cluster node7:
# csmit -m 20G -c 1 "./log2_ChIPinput_summedRPM_each_chr_commandArgs_v220818.R /home/ajt200/analysis/REC8_pooled/REC8_HA_Rep1_ChIP_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam /home/ajt200/analysis/REC8_pooled/REC8_MYC_Rep2_input_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam REC8_HA_Rep1_ChIP REC8_MYC_Rep2_input"

library(GenomicAlignments)

chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")

args <- commandArgs(trailingOnly = TRUE)
ChIP_bamFile <- args[1]
input_bamFile <- args[2]
ChIP_libName <- args[3]
input_libName <- args[4]

outDir <- "./" 

# ChIP
# Load BAM and create RangedData object
ChIP_readGAlignmentPairs <- readGAlignmentPairs(ChIP_bamFile)
ChIP_ranges <- ranges(ChIP_readGAlignmentPairs)
ChIP_chrs <- as.character(as.data.frame(seqnames(ChIP_readGAlignmentPairs))[,1])
ChIP_strand <- as.character(strand(ChIP_readGAlignmentPairs))
ChIP_ranged <- RangedData(space = ChIP_chrs,
                          ranges = ChIP_ranges,
                          strand = "*")
ChIP_ranged <- ChIP_ranged[space(ChIP_ranged) != "chloroplast" & space(ChIP_ranged) != "mitochondria",]
names(ChIP_ranged) <- sub("", "Chr", names(ChIP_ranged))

# Convert library ranged data into GRanges object
ChIP_rangedGR <- as(ChIP_ranged, "GRanges")
save(ChIP_rangedGR, file = paste0(outDir, ChIP_libName, "_GRanges.RData"))

# Calculate library size
ChIP_libSize <- length(ChIP_rangedGR)

# Calculate "per million" scaling factor
ChIP_RPM_scaling_factor <- ChIP_libSize/1e+06

ChIP_rangedGR_chr_RPM_allchrs <- NULL
for(i in 1:5) {
  ChIP_rangedGR_chr <- ChIP_rangedGR[seqnames(ChIP_rangedGR) == chrs[i]]
  ChIP_rangedGR_chr_reads <- length(ChIP_rangedGR_chr)
  ChIP_rangedGR_chr_RPM <- ChIP_rangedGR_chr_reads/ChIP_RPM_scaling_factor
  ChIP_rangedGR_chr_RPM_allchrs <- c(ChIP_rangedGR_chr_RPM_allchrs, ChIP_rangedGR_chr_RPM)
}
write.table(ChIP_rangedGR_chr_RPM_allchrs,
            file = paste0(outDir, ChIP_libName, "_summedRPM_each_chr.txt"))

# input
# Load BAM and create RangedData object
input_readGAlignmentPairs <- readGAlignmentPairs(input_bamFile)
input_ranges <- ranges(input_readGAlignmentPairs)
input_chrs <- as.character(as.data.frame(seqnames(input_readGAlignmentPairs))[,1])
input_strand <- as.character(strand(input_readGAlignmentPairs))
input_ranged <- RangedData(space = input_chrs,
                          ranges = input_ranges,
                          strand = "*")
input_ranged <- input_ranged[space(input_ranged) != "chloroplast" & space(input_ranged) != "mitochondria",]
names(input_ranged) <- sub("", "Chr", names(input_ranged))

# Convert library ranged data into GRanges object
input_rangedGR <- as(input_ranged, "GRanges")
save(input_rangedGR, file = paste0(outDir, input_libName, "_GRanges.RData"))

# Calculate library size
input_libSize <- length(input_rangedGR)

# Calculate "per million" scaling factor
input_RPM_scaling_factor <- input_libSize/1e+06

input_rangedGR_chr_RPM_allchrs <- NULL
for(i in 1:5) {
  input_rangedGR_chr <- input_rangedGR[seqnames(input_rangedGR) == chrs[i]]
  input_rangedGR_chr_reads <- length(input_rangedGR_chr)
  input_rangedGR_chr_RPM <- input_rangedGR_chr_reads/input_RPM_scaling_factor
  input_rangedGR_chr_RPM_allchrs <- c(input_rangedGR_chr_RPM_allchrs, input_rangedGR_chr_RPM)
}
write.table(input_rangedGR_chr_RPM_allchrs,
            file = paste0(outDir, input_libName, "_summedRPM_each_chr.txt"))

# Calculate log2(ChIP/input) ratio
log2_ChIPinput_rangerGR_chr_RPM_allchrs <- log2(ChIP_rangedGR_chr_RPM_allchrs/input_rangedGR_chr_RPM_allchrs)
write.table(log2_ChIPinput_rangerGR_chr_RPM_allchrs,
            file = paste0(outDir, "log2_", ChIP_libName, "_", input_libName, "_summedRPM_each_chr.txt"))

