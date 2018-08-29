#!/applications/R/R-3.4.0/bin/Rscript

# Profile mean DNA methylation around peaks and random loci

# Usage via Condor submission system on node7:
# csmit -m 10G -c 1 "Rscript peak_Profiles_DNAmeth_commandArgs.R /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM980986_WT_rep2_CG.wig.bed.gr.tab.bed CGmeth"

# Source functions to be used in this script
source("/projects/ajt200/Rfunctions/covMatrix_DNAmethMatrix_target_ranLoc_R3.4.0.R")

library(EnrichedHeatmap)
library(genomation)
library(regioneR)

args <- commandArgs(trailingOnly = T)
covDatPath <- as.character(args[1])
libName <- as.character(args[2])

inDir <- "/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/DiffBind/"
matDir <- "/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/DiffBind/analysis_01_DB_sites_kss_exceeds_wt/matrices/"

# Chromosome definitions
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chrStart <- c(1, 1, 1, 1, 1)
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))

# Import peaks as GRanges object
load(paste0(inDir,
            "DB_sites_kss_exceeds_wt_wt_kss_REC8_HA_peaks_0.05.RData"))
peaksGR <- dbaObjAnalyzeDB_kss_exceeds_wt
dbaObjAnalyzeDB_kss_exceeds_wt <- NULL
seqlevels(peaksGR) <- sub("chr", "Chr", seqlevels(peaksGR))
strand(peaksGR) <- "*"
print(length(peaksGR))
#[1] 2781
# Generate GRanges object containing random loci of same number
# and size distribution as peaksGR
set.seed(937234)
ranLocGR <- randomizeRegions(peaksGR,
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
                            "_norm_cov_peaks_mat1_smoothed_target_and_flank_dataframe_colMeans.txt"),
                     paste0(matDir, libName,
                            "_norm_cov_ranLoc_mat2_smoothed_target_and_flank_dataframe_colMeans.txt")))

# Run DNAmethMatrixHexile() function on each DNA methylation GRanges object to obtain matrices
## containing normalised DNA methylation values around target and random loci
DNAmethMatrix(signal = gr, target = peaksGR, ranLoc = ranLocGR,
              targetSize = mean(width(peaksGR)), flankSize = 2000, winSize = 20,
              DNAmethOutDFCM = outDFCM, x = 1)
print(paste0(libName, " profile calculation complete"))


