#!/applications/R/R-3.3.2/bin/Rscript

# Calculate mean coverage between target start and end coordinates of targets for
# each wt and mutant library and tabulate
# This will enable Mann-Whitney U tests comparing each wt library with each kss library

# Usage on node7:
# /scripts/csmit -m 20G -c 1 "./targets_meanCov.R /home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab.bed wt_REC8_HA_Rep1" 

args <- commandArgs(trailingOnly = T)
libPath <- args[1]
libName <- args[2]

library(segmentSeq)
library(genomation)

outDir <- "/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/genome_wide/H3K9me2peakProfiles/U_test/"

inDir <- "/projects/ajt200/BAM_masters/H3K9me2/WT/peaks/PeakRanger1.18/ranger/p0.05_q0.05/"
# Import peaks as GRanges object
load(paste0(inDir,
            "WT_H3K9me2_ChIP_armrangerPeaksGRmergedOverlaps_minuslog10_p0.05_q0.05_noMinWidth.RData"))
load(paste0(inDir,
            "WT_H3K9me2_ChIP_perirangerPeaksGRmergedOverlaps_minuslog10_p0.05_q0.05_noMinWidth.RData"))
peaksGR <- sort(c(armrangerPeaksGRmergedOverlaps, perirangerPeaksGRmergedOverlaps))
armrangerPeaksGRmergedOverlaps <- NULL
perirangerPeaksGRmergedOverlaps <- NULL
strand(peaksGR) <- "*"
print(peaksGR)

gr <- readGeneric(libPath, meta.col = list(coverage = 4))

targetsMeanCov <- function(targetsName, targetsGR, signalGR, libName) {
  targetsOverlaps <- getOverlaps(targetsGR,
                                 signalGR,
                                 whichOverlaps = T,
                                 ignoreStrand = T)
  targetsCov <- sapply(targetsOverlaps,
                       function(x) mean(signalGR$coverage[x]))
  write.table(targetsCov,
              file = paste0(outDir,
                            libName,
                            "_",
                            targetsName,
                            "_meanCov.txt"),
              col.names = libName)
}

targetsMeanCov(targetsName = "H3K9me2_peaks",
               targetsGR = peaksGR,
               signalGR = gr,
               libName = libName)

