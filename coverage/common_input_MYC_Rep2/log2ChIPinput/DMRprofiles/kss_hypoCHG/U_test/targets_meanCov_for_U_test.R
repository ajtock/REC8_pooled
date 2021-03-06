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

REC8Dir <- "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/"
wtH3K9me2Dir <- "/projects/ajt200/BAM_masters/H3K9me2/WT/coverage/log2ChIPinput/"
kssH3K9me2Dir <- "/projects/ajt200/BAM_masters/H3K9me2/kss/coverage/log2ChIPinput/"
wtSPO11Dir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/"
kssSPO11Dir <- "/projects/ajt200/BAM_masters/SPO11-oligo/suvh456/coverage/log2ChIPinput/"

wt_REC8_HA_Rep1 <- system(paste0("ls ", 
                                 REC8Dir,
                                 "log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab.bed"),
                          intern = T)
wt_REC8_HA_Rep2 <- system(paste0("ls ", 
                                 REC8Dir,
                                 "log2_REC8_HA_Rep2_ChIP_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab.bed"),
                          intern = T)
wt_REC8_MYC_Rep1 <- system(paste0("ls ", 
                                 REC8Dir,
                                 "log2_REC8_MYC_Rep1_ChIP_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab.bed"),
                          intern = T)
kss_REC8_HA_Rep1 <- system(paste0("ls ", 
                                  REC8Dir,
                                  "log2_kss_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab.bed"),
                           intern = T)
wt_H3K9me2 <- system(paste0("ls ", wtH3K9me2Dir, "log2_WT_H3K9me2_ChIP_WT_H3K9me2_input_norm_allchrs_coverage_coord_tab.bed"), intern = T)
kss_H3K9me2 <- system(paste0("ls ", kssH3K9me2Dir, "log2_kss_H3K9me2_ChIP_kss_H3K9me2_input_norm_allchrs_coverage_coord_tab.bed"), intern = T)
wt_SPO11oligos_RPI1 <- system(paste0("ls ", wtSPO11Dir, "log2wtSPO11oligoRPI1NakedDNA_norm_allchrs_coverage_coord_tab.bed"), intern = T)
wt_SPO11oligos_RPI3 <- system(paste0("ls ", wtSPO11Dir, "log2wtSPO11oligoRPI3NakedDNA_norm_allchrs_coverage_coord_tab.bed"), intern = T)
wt_SPO11oligos_RPI8 <- system(paste0("ls ", wtSPO11Dir, "log2wtSPO11oligoRPI8NakedDNA_norm_allchrs_coverage_coord_tab.bed"), intern = T)
kss_SPO11oligos_RPI34 <- system(paste0("ls ", kssSPO11Dir, "log2suvh456SPO11oligoRPI34NakedDNA_norm_allchrs_coverage_coord_tab.bed"), intern = T)
kss_SPO11oligos_RPI35 <- system(paste0("ls ", kssSPO11Dir, "log2suvh456SPO11oligoRPI35NakedDNA_norm_allchrs_coverage_coord_tab.bed"), intern = T)

libPaths <- list(
                 wt_REC8_HA_Rep1,
                 wt_REC8_HA_Rep2,
                 wt_REC8_MYC_Rep1,
                 kss_REC8_HA_Rep1,
                 wt_H3K9me2,
                 kss_H3K9me2,
                 wt_SPO11oligos_RPI1,
                 wt_SPO11oligos_RPI3,
                 wt_SPO11oligos_RPI8,
                 kss_SPO11oligos_RPI34,
                 kss_SPO11oligos_RPI35
                )

libNames <- list(
                 "wt_REC8_HA_Rep1",
                 "wt_REC8_HA_Rep2",
                 "wt_REC8_MYC_Rep1",
                 "kss_REC8_HA_Rep1",
                 "wt_H3K9me2",
                 "kss_H3K9me2",
                 "wt_SPO11oligos_RPI1",
                 "wt_SPO11oligos_RPI3",
                 "wt_SPO11oligos_RPI8",
                 "kss_SPO11oligos_RPI34",
                 "kss_SPO11oligos_RPI35"
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
                   "wt_REC8_HA_Rep1" = wt_REC8_HA_Rep1,
                   "wt_REC8_HA_Rep2" = wt_REC8_HA_Rep2,
                   "wt_REC8_MYC_Rep1" = wt_REC8_MYC_Rep1,
                   "kss_REC8_HA_Rep1" = kss_REC8_HA_Rep1,
                   "wt_H3K9me2" = wt_H3K9me2,
                   "kss_H3K9me2" = kss_H3K9me2,
                   "wt_SPO11oligos_RPI1" = wt_SPO11oligos_RPI1,
                   "wt_SPO11oligos_RPI3" = wt_SPO11oligos_RPI3,
                   "wt_SPO11oligos_RPI8" = wt_SPO11oligos_RPI8,
                   "kss_SPO11oligos_RPI34" = kss_SPO11oligos_RPI34,
                   "kss_SPO11oligos_RPI35" = kss_SPO11oligos_RPI35
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
                          "kss_hypoCHG_DMRs_meanCov_wt_kss_REC8_H3K9me2_SPO11oligos.txt"),
            sep = "\t",
            quote = F,
            row.names = F)

 
