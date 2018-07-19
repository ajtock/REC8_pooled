#!/applications/R/R-3.3.2/bin/Rscript

# Calculate mean z-score standardised log2-transformed library-size-normalised
# wt and kss coverage levels at TEs (mean over TE width) within each family
# Calculate correlation coefficients for levels in wt vs kss and for
# dataset vs dataset

# Usage via Condor submission system on node7:
# csmit -m 50G -c 1 "Rscript crossover_Profiles_commandArgs.R /home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab.bed REC8_HA_Rep1"


library(regioneR)
library(GenomicRanges)
library(genomation)
library(doParallel)
registerDoParallel(cores=7)
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())


args <- commandArgs(trailingOnly = T)
covDatPath <- as.character(args[1])
libName <- as.character(args[2])

chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chrStart <- c(1, 1, 1, 1, 1)
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))
#seqlevels(genome) <- sub("Chr", "", seqlevels(genome))

DNAfamNames <- c("dna", "heli", "ptmari", "mudr", "enspm", "hat", "harbinger")
RNAfamNames <- c("rna", "gypsy", "copia", "linel1", "sine")
DNAdir <- "/projects/ajt200/TAIR10/TE_classes/DNA/"
RNAdir <- "/projects/ajt200/TAIR10/TE_classes/RNA/"


### DNA transposons

TEsDNAGR <- lapply(seq_along(DNAfamNames), function(x) {
  TEsDNA <- read.table(file = paste0(DNAdir, "TAIR10_Buisine_TEs_strand_tab_ann_", DNAfamNames[x], ".txt"), header = T)
  GRanges(seqnames = TEsDNA$chr, ranges = IRanges(start = TEsDNA$start, end = TEsDNA$end), strand = "*")
})

r

############################################################
# mean log2-normalised coverage levels at TEs              #
############################################################

inDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/crossover_hotspots/"
outDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/crossover_hotspots/regioneR/"
plotDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/crossover_hotspots/regioneR/plots/"

load(paste0(inDir, "all_CO_hotspots_log2wtSPO11oligoMeanAllRepsNakedDNA_levels_GRanges.RData"))
allCovHotspots <- c(allCovHotspotsGR[[1]], allCovHotspotsGR[[2]], allCovHotspotsGR[[3]], allCovHotspotsGR[[4]], allCovHotspotsGR[[5]], allCovHotspotsGR[[6]], allCovHotspotsGR[[7]], allCovHotspotsGR[[8]], allCovHotspotsGR[[9]], allCovHotspotsGR[[10]], allCovHotspotsGR[[11]], allCovHotspotsGR[[12]], allCovHotspotsGR[[13]])
RAC1 <- allCovHotspots[allCovHotspots$hotspot == "RAC1"]
a3 <- allCovHotspots[allCovHotspots$hotspot == "3a"]
b3 <- allCovHotspots[allCovHotspots$hotspot == "3b"]
x130 <- allCovHotspots[allCovHotspots$hotspot == "130x"]
I4a <- allCovHotspots[allCovHotspots$hotspot == "I4a"]

hotspotsReducedGR <- c(reduce(RAC1, min.gapwidth = 2L), reduce(a3, min.gapwidth = 2L), reduce(b3, min.gapwidth = 2L),
                       reduce(x130, min.gapwidth = 2L), reduce(I4a, min.gapwidth = 2L))

SPO11oligos <- system("ls /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/log2wtSPO11oligoRPI1NakedDNA_norm_allchrs_coverage_coord_tab.bed", intern = T)

SPO11oligosGR <- readGeneric(SPO11oligos, meta.col = list(value = 4))


# Perform permutation test with randomized regions generated on a per chromosome basis
set.seed(123)
ptHotspotsSPO11cov <- permTest(A = hotspotsReducedGR, genome = genome,
                               randomize.function = randomizeRegions,
                               allow.overlaps = FALSE, per.chromosome = TRUE,
                               evaluate.function = meanInRegions, x = SPO11oligosGR,
                               ntimes = 1000, mc.set.seed = FALSE, mc.cores = 40)
save(ptHotspotsSPO11cov, file = paste0(outDir, "pt_mean_SPO11_coverage_at_CO_hotspots_RAC1_3a_3b_130x_I4a_nperm1000.RData"))

pdf(file = paste0(plotDir, "permTest_mean_SPO11_coverage_at_CO_hotspots_RAC1_3a_3b_130x_I4a_nperm1000.pdf"), width = 10, height = 7)
plot(ptHotspotsSPO11cov, main = "", xlab = "Mean log2(SPO11-1-oligos/naked DNA) coverage", ylab = "Relative frequency")
dev.off()

set.seed(123)
ptHotspotsSPO11cov <- permTest(A = hotspotsReducedGR, genome = genome,
                               randomize.function = randomizeRegions,
                               allow.overlaps = FALSE, per.chromosome = TRUE,
                               evaluate.function = meanInRegions, x = SPO11oligosGR,
                               ntimes = 10000, mc.set.seed = FALSE, mc.cores = 40)
save(ptHotspotsSPO11cov, file = paste0(outDir, "pt_mean_SPO11_coverage_at_CO_hotspots_RAC1_3a_3b_130x_I4a_nperm10000.RData"))

pdf(file = paste0(plotDir, "permTest_mean_SPO11_coverage_at_CO_hotspots_RAC1_3a_3b_130x_I4a_nperm10000.pdf"), width = 10, height = 7)
plot(ptHotspotsSPO11cov, main = "", xlab = "Mean log2(SPO11-1-oligos/naked DNA) coverage", ylab = "Relative frequency")
dev.off()


 


