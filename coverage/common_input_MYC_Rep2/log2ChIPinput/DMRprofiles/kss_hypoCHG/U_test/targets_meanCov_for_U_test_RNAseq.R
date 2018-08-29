#!/applications/R/R-3.3.2/bin/Rscript

# Calculate mean coverage between target start and end coordinates of targets for
# each wt and mutant library and tabulate
# This will enable Mann-Whitney U tests comparing each wt library with each kss library

library(segmentSeq)
library(parallel)
library(genomation)

outDir <- "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/DMRprofiles/kss_hypoCHG/U_test/"

# Import DMRs as GRanges object
DMRs <- read.table("/home/ajt200/BS_Seq/Stroud_2013/DMRs/suvh456_hypoCHG_DMR_vs3reps_min4filter_mg200.bed", header = F)
DMRsGR <- GRanges(seqnames = DMRs[,1],
                  ranges = IRanges(start = DMRs[,2], end = DMRs[,3]),
                  strand = "*")
seqlevels(DMRsGR) <- sub("chr", "Chr", seqlevels(DMRsGR))
print("***********DMRs***********")
print(DMRsGR)

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
  targetsMeanCov(targetsName = "kss_hypoCHG_DMRs",
                 targetsGR = DMRsGR,
                 signalGR = grl[[i]],
                 libName = libNames[i])
}, mc.cores = 4)

covTabs <- mclapply(seq_along(libNames), function(i) {
  read.table(file = paste0(outDir,
                           libNames[i],
                           "_kss_hypoCHG_DMRs_meanCov.txt"))
}, mc.cores = length(libNames))

targetsCovBind <- cbind(as.vector(seqnames(DMRsGR)),
                        start(DMRsGR),
                        end(DMRsGR),
                        as.vector(strand(DMRsGR)))
colnames(targetsCovBind) <- c("chr", "start", "end", "strand")
for(i in 1:length(covTabs)) {
  targetsCovBind <- cbind(targetsCovBind, covTabs[[i]])
}

write.table(targetsCovBind,
            file = paste0(outDir,
                          "kss_hypoCHG_DMRs_meanCov_wt_kss_RNAseq.txt"),
            sep = "\t",
            quote = F,
            row.names = F)

 
