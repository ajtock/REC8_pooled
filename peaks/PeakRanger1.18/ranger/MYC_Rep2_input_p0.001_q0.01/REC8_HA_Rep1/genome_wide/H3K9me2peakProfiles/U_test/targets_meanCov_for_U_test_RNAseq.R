#!/applications/R/R-3.3.2/bin/Rscript

# Calculate mean coverage between target start and end coordinates of targets for
# each wt and mutant library and tabulate
# This will enable Mann-Whitney U tests comparing each wt library with each kss library

library(segmentSeq)
library(parallel)
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

wtRNAseqDir <- "/home/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/WT/coverage/"
kssRNAseqDir <- "/home/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/kss/coverage/"

wt_RNAseq_Rep1 <- system(paste0("ls ", wtRNAseqDir, "WT_RNAseq_Chris_Rep1_norm_allchrs_coverage_coord_tab.bed"), intern = T)
wt_RNAseq_Rep2 <- system(paste0("ls ", wtRNAseqDir, "WT_RNAseq_Chris_Rep2_norm_allchrs_coverage_coord_tab.bed"), intern = T)
kss_RNAseq_Rep1 <- system(paste0("ls ", kssRNAseqDir, "kss_RNAseq_Chris_Rep1_norm_allchrs_coverage_coord_tab.bed"), intern = T)
kss_RNAseq_Rep2 <- system(paste0("ls ", kssRNAseqDir, "kss_RNAseq_Chris_Rep2_norm_allchrs_coverage_coord_tab.bed"), intern = T)

libPaths <- list(
                 wt_RNAseq_Rep1,
                 wt_RNAseq_Rep2,
                 kss_RNAseq_Rep1, 
                 kss_RNAseq_Rep2
                ) 

libNames <- list(
                 "wt_RNAseq_Rep1",
                 "wt_RNAseq_Rep2",
                 "kss_RNAseq_Rep1",
                 "kss_RNAseq_Rep2"
                )

grTmp <- mclapply(seq_along(libPaths), function(x) {
  readGeneric(libPaths[[x]], meta.col = list(coverage = 4))
}, mc.cores = length(libPaths), mc.preschedule = F)
for(i in 1:length(grTmp)) {
  #seqlevels(grTmp[[i]]) <- sub("Chr", "", seqlevels(grTmp[[i]]))
  assign(paste0(libNames[i]), grTmp[[i]])
}

print("grTmp loaded")

grl <- GRangesList(
                   "wt_RNAseq_Rep1" = wt_RNAseq_Rep1,
                   "wt_RNAseq_Rep2" = wt_RNAseq_Rep2,
                   "kss_RNAseq_Rep1" = kss_RNAseq_Rep1,
                   "kss_RNAseq_Rep2" = kss_RNAseq_Rep2
                  )

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

mclapply(seq_along(grl), function(i) {
  targetsMeanCov(targetName = "H3K9me2_peaks",
                 targetsGR = peaksGR,
                 signalGR = grl[[i]],
                 libName = libNames[i])
}, mc.cores = 4)

covTabs <- mclapply(seq_along(libNames), function(i) {
  read.table(file = paste0(outDir,
                           libNames[i],
                           "_H3K9me2_peaks_meanCov.txt"))
}, mc.cores = length(libNames))

targetsCovBind <- cbind(as.vector(seqnames(peaksGR)),
                        start(peaksGR),
                        end(peaksGR),
                        as.vector(strand(peaksGR)))
colnames(targetsCovBind) <- c("chr", "start", "end", "strand")
for(i in 1:length(covTabs)) {
  targetsCovBind <- cbind(targetsCovBind, covTabs[[i]])
}

write.table(targetsCovBind,
            file = paste0(outDir,
                          "H3K9me2_peaks_meanCov_wt_kss_RNAseq.txt"),
            sep = "\t",
            quote = F,
            row.names = F)

 
