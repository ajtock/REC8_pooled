#!/applications/R/R-3.3.2/bin/Rscript

# Profile mean coverage around TEs and random loci

# Usage via Condor submission system on node7:
# csmit -m 50G -c 1 "Rscript TE_Profiles_commandArgs.R /home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab.bed REC8_HA_Rep1"

# Source functions to be used in this script
source("/projects/ajt200/Rfunctions/covMatrix_DNAmethMatrix_target_ranLoc.R")

library(EnrichedHeatmap)
library(genomation)
library(regioneR)

args <- commandArgs(trailingOnly = T)
covDatPath <- as.character(args[1])
libName <- as.character(args[2])

matDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/TEProfiles/analysis_01/matrices/"

# Chromosome definitions
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chrStart <- c(1, 1, 1, 1, 1)
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))

# Import TEs as GRanges object
TEs <- system("ls /projects/ajt200/TAIR10/TAIR10_Buisine_TEs_strand_tab_ann.txt",
              intern = T)
TEsGR <- readGeneric(TEs,
                     header = T,
                     strand = 4,
                     meta.col = list(Transposon_Name = 5))
print("***********TEs***********")
print(TEsGR)
# Generate GRanges object containing random loci of same number and
# size distribution as TEsGR
ranLocGR <- randomizeRegions(TEsGR,
                             genome = genome,
                             per.chromosome = TRUE,
                             allow.overlaps = TRUE)

# Specify locations of normalised per base coverage files
libPath <- system(paste0("ls ", covDatPath), intern = T)

# Import coverage files as GRanges objects and assign to library names
gr <- readGeneric(libPath, meta.col = list(coverage = 4))
assign(paste0(libName), gr)

# Define column mean coverage outfile (mean profiles)
outDFCM <- list(list(paste0(matDir, libName,
                            "_norm_cov_TEs_mat1_smoothed_target_and_flank_dataframe_colMeans.txt"),
                     paste0(matDir, libName,
                            "_norm_cov_ranLoc_mat2_smoothed_target_and_flank_dataframe_colMeans.txt")))

# Run covMatrix() function on each coverage GRanges object to obtain matrices
## containing normalised coverage values around target and random loci
covMatrix(signal = gr, target = TEsGR, ranLoc = ranLocGR,
          targetSize = mean(width(TEsGR)), flankSize = 2000, winSize = 20,
          outDFCM = outDFCM, x = 1)
print(paste0(libName, " profile calculation complete"))


