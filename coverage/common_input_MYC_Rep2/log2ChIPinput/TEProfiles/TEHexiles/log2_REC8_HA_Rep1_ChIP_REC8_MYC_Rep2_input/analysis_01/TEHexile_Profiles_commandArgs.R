#!/applications/R/R-3.3.2/bin/Rscript

# Profile mean coverage around TEs and random loci
# separated into hexiles

# Usage via Condor submission system on node7:
# csmit -m 20G -c 1 "Rscript TEHexile_Profiles_commandArgs.R /home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab.bed REC8_HA_Rep1" 

# Source functions to be used in this script
source("/projects/ajt200/Rfunctions/covMatrix_DNAmethMatrix_target_ranLoc_hexiles.R")

library(EnrichedHeatmap)
library(genomation)
library(regioneR)

args <- commandArgs(trailingOnly = T)
covDatPath <- as.character(args[1])
libName <- as.character(args[2])

hexDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/TEProfiles/TEHexiles/log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input/"
matDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/TEProfiles/TEHexiles/log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input/analysis_01/matrices/"

# Define hexile locations
targetsHex <- paste0(hexDir, "hexile_", c(1:6), "_log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_TECov.txt")

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
# Generate GRanges object containing random loci of same number and size distribution as TEsGR
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
                            "_norm_cov_TE_hexile",
                            c(1:6),
                            "_by_log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_mat1_smoothed_target_and_flank_dataframe_colMeans.txt"),
                     paste0(matDir, libName,
                            "_norm_cov_ranLoc_hexile",
                            c(1:6),
                            "_mat2_smoothed_target_and_flank_dataframe_colMeans.txt")))

# Run covMatrixHexile() function on each coverage GRanges object to obtain matrices
## containing normalised coverage values around target and random loci
for(y in 1:6) {
  covMatrixHexile(signal = gr, targetsHex = targetsHex,
                  targetSize = mean(width(TEsGR)), flankSize = 2000, winSize = 20,
		  outDFCM = outDFCM, x = 1, y = y)
  print(paste(libName, " hexile", y, " profile calculation complete"))
}
 
