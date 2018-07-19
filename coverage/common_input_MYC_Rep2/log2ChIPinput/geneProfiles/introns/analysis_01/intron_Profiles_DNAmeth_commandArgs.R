#!/applications/R/R-3.3.2/bin/Rscript

# Profile mean DNA methylation around gene introns and random loci

# Usage via Condor submission system on node7:
# csmit -m 50G -c 1 "Rscript intron_Profiles_DNAmeth_commandArgs.R /home/meiosis/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM980986_WT_rep2_CG.wig.bed.gr.tab.bed CGmeth"

# Source functions to be used in this script
source("/projects/ajt200/Rfunctions/covMatrix_DNAmethMatrix_target_ranLoc.R")

library(EnrichedHeatmap)
library(genomation)
library(regioneR)

args <- commandArgs(trailingOnly = T)
covDatPath <- as.character(args[1])
libName <- as.character(args[2])

matDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/geneProfiles/introns/analysis_01/matrices/"

# Chromosome definitions
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chrStart <- c(1, 1, 1, 1, 1)
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))

# Import introns tables and convert to GRanges object
intronsPlus <- read.table("/projects/ajt200/TAIR10/all_plus_introns.txt", header = T)
intronsMinus <- read.table("/projects/ajt200/TAIR10/all_minus_introns.txt", header = T)
intronsPlusGR <- GRanges(seqnames = paste0("Chr", intronsPlus$chr),
                         ranges = (IRanges(start = intronsPlus$all.intron.starts+1,
                                           end = intronsPlus$all.intron.stops-1)),
                         strand = "+", gene_model = intronsPlus$all.atgs)
intronsMinusGR <- GRanges(seqnames = paste0("Chr", intronsMinus$chr),
                          ranges = (IRanges(start = intronsMinus$all.intron.stops+1,
                                            end = intronsMinus$all.intron.starts-1)),
                          strand = "-", gene_model = intronsMinus$all.atgs)
intronsGR <- sort(append(intronsPlusGR, intronsMinusGR), by = ~ seqnames + start + end)
print("***********introns***********")
print(intronsGR)
# Generate GRanges object containing random loci of same number and
# size distribution as intronsGR
ranLocGR <- randomizeRegions(intronsGR,
                             genome = genome,
                             per.chromosome = TRUE,
                             allow.overlaps = TRUE)

# Specify locations of normalised per base coverage files
libPath <- system(paste0("ls ", covDatPath), intern = T)

# Import DNA methylation proportion files as GRanges objects
# and assign to library names
gr <- readGeneric(libPath, meta.col = list(coverage = 4))
seqlevels(gr) <- sub("chr", "Chr", seqlevels(gr))
assign(paste0(libName), gr)

# Define column mean DNA methylation outfile (mean profiles)
outDFCM <- list(list(paste0(matDir, libName,
                            "_norm_cov_introns_mat1_smoothed_target_and_flank_dataframe_colMeans.txt"),
                     paste0(matDir, libName,
                            "_norm_cov_ranLoc_mat2_smoothed_target_and_flank_dataframe_colMeans.txt")))

# Run DNAmethMatrixHexile() function on each DNA methylation GRanges object to obtain matrices
## containing normalised DNA methylation values around target and random loci
DNAmethMatrix(signal = gr, target = intronsGR, ranLoc = ranLocGR,
              targetSize = mean(width(intronsGR)), flankSize = 1000, winSize = 20,
              DNAmethOutDFCM = outDFCM, x = 1)
print(paste0(libName, " profile calculation complete"))


