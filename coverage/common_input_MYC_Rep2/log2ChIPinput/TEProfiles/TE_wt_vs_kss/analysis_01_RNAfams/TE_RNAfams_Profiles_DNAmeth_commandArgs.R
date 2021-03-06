#!/applications/R/R-3.3.2/bin/Rscript

# Profile mean DNA methylation around TEs and random loci

# Usage via Condor submission system on node7:
# csmit -m 5G -c 1 "Rscript TE_RNAfams_Profiles_DNAmeth_commandArgs.R /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM980986_WT_rep2_CG.wig.bed.gr.tab.bed CGmeth"

# Source functions to be used in this script
source("/projects/ajt200/Rfunctions/covMatrix_DNAmethMatrix_target_ranLoc.R")

library(EnrichedHeatmap)
library(genomation)
library(regioneR)

args <- commandArgs(trailingOnly = T)
covDatPath <- as.character(args[1])
libName <- as.character(args[2])

matDir <- "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/TEProfiles/TE_wt_vs_kss/analysis_01_RNAfams/matrices/"

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
  ranLocGR <- randomizeRegions(TEsGR, 
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
                              "_norm_cov_", RNAfamNames[h], "_TEs_mat1_smoothed_target_and_flank_dataframe_colMeans.txt"),
                       paste0(matDir, libName,
                              "_norm_cov_", RNAfamNames[h], "_ranLoc_mat2_smoothed_target_and_flank_dataframe_colMeans.txt")))
  
  # Run DNAmethMatrixHexile() function on each DNA methylation GRanges object to obtain matrices
  ## containing normalised DNA methylation values around target and random loci
  DNAmethMatrix(signal = gr, target = TEsGR, ranLoc = ranLocGR,
                targetSize = mean(width(TEsGR)), flankSize = 2000, winSize = 20,
                DNAmethOutDFCM = outDFCM, x = 1)
  print(paste0(libName, " ", RNAfamNames[h], " profile calculation complete"))
}
