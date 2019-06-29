#!/applications/R/R-3.4.0/bin/Rscript

# Profile mean coverage around RNA TEs and random loci

# Usage via Condor submission system on node7:
# csmit -m 20G -c 1 "./TE_RNAfams_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep1_input_norm_allchrs_coverage_coord_tab.bed REC8_HA_Rep1"

# Source functions to be used in this script
source("/projects/ajt200/Rfunctions/covMatrix_DNAmethMatrix_target_ranLoc_R3.4.0.R")

library(EnrichedHeatmap)
library(genomation)
library(regioneR)

args <- commandArgs(trailingOnly = T)
flankSize <- as.numeric(args[1])
flankName <- as.character(args[2])
winSize <- as.numeric(args[3])
covDatPath <- as.character(args[4])
libName <- as.character(args[5])

matDir <- "./matrices/"
plotDir <- "./plots/"
system(paste0("[ -d ", matDir, " ] || mkdir ", matDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

DNAfamNames <- c("dna", "heli", "ptmari", "mudr", "enspm", "hat", "harbinger")
RNAfamNames <- c("rna", "gypsy", "copia", "linel1", "sine")
DNAdir <- "/projects/ajt200/TAIR10/TE_classes/DNA/"
RNAdir <- "/projects/ajt200/TAIR10/TE_classes/RNA/"

# Chromosome definitions
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chrStart <- c(1, 1, 1, 1, 1)
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))

# Import RNA TEs as GRanges object
for(h in seq_along(RNAfamNames)) {
  # Import TEs as GRanges object
  TEs <- system(paste0("ls ", RNAdir, "TAIR10_Buisine_TEs_strand_tab_ann_", RNAfamNames[h], ".txt"),
                intern = T)
  TEsGR <- readGeneric(TEs,
                       header = T,
                       strand = 4,
                       meta.col = list(Transposon_Name = 5))
  print("***********TEs***********")
  print(TEsGR)
  # Generate GRanges object containing random loci of same number and
  # size distribution as TEsGR
  set.seed(374592)
  ranLocGR <- randomizeRegions(TEsGR,
                               genome = genome,
                               per.chromosome = TRUE,
                               allow.overlaps = TRUE)
  
  # Specify locations of normalised per base coverage files
  libPath <- system(paste0("ls ", covDatPath), intern = T)
  
  # Import coverage files as GRanges objects and assign to library names
  covGR <- readGeneric(libPath, meta.col = list(coverage = 4))
  assign(paste0(libName), covGR)
  
  # Define matrix and column mean coverage outfile (mean profiles)
  outDF <- list(paste0(matDir, libName,
                       "_norm_cov_", RNAfamNames[h], "_TEs_smoothed_target_and_", flankName, "_flank_dataframe.txt"),
                paste0(matDir, libName,
                       "_norm_cov_", RNAfamNames[h], "_ranLoc_smoothed_target_and_", flankName, "_flank_dataframe.txt"))
  outDFcolMeans <- list(paste0(matDir, libName,
                               "_norm_cov_", RNAfamNames[h], "_TEs_smoothed_target_and_", flankName, "_flank_dataframe_colMeans.txt"),
                        paste0(matDir, libName,
                               "_norm_cov_", RNAfamNames[h], "_ranLoc_smoothed_target_and_", flankName, "_flank_dataframe_colMeans.txt"))
  
  # Run covMatrix() function on each coverage GRanges object to obtain matrices
  ## containing normalised coverage values around target and random loci
  covMatrix(signal = covGR,
            feature = TEsGR,
            ranLoc = ranLocGR,
            featureSize = mean(width(TEsGR)),
            flankSize = flankSize,
            winSize = winSize,
            outDF = outDF,
            outDFcolMeans = outDFcolMeans)
  print(paste0(libName, " ", RNAfamNames[h], " profile calculation complete"))
}
