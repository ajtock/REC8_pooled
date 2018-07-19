## Calculate mean RNA-seq, H3K9me2 and SPO11-1-oligo coverage
## in wt and kss for each TE between start and end

library(GenomicAlignments)
library(segmentSeq)
library(parallel)
library(genomation)

TEs <- read.table("/projects/ajt200/TAIR10/TAIR10_Buisine_TEs_strand_tab_ann.txt",
                  header = T)
print(dim(TEs))
#[1] 31189     7

TEsGR <- GRanges(seqnames = TEs$Chr,
                 ranges = IRanges(start = TEs$start,
                                  end = TEs$end),
                 strand = TEs$strand,
                 Transposon_Name = TEs$Transposon_Name,
                 Transposon_Family = TEs$Transposon_Family,
                 Transposon_Superfamily = TEs$Transposon_Super_Family)

chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)

outDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/TEProfiles/TE_wt_vs_kss/"
wtRNAseqDir <- "/home/meiosis/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/WT/coverage/"
kssRNAseqDir <- "/home/meiosis/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/kss/coverage/"
wtH3K9me2Dir <- "/projects/ajt200/BAM_masters/H3K9me2/WT/coverage/log2ChIPinput/"
kssH3K9me2Dir <- "/projects/ajt200/BAM_masters/H3K9me2/kss/coverage/Xiaohui/log2ChIPinput/"
wtSPO11Dir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/"
kssSPO11Dir <- "/projects/ajt200/BAM_masters/SPO11-oligo/suvh456/coverage/log2ChIPinput/"

wt_RNAseq_Rep1 <- system(paste0("ls ", wtRNAseqDir, "WT_RNAseq_Chris_Rep1_norm_allchrs_coverage_coord_tab.bed"), intern = T)
wt_RNAseq_Rep2 <- system(paste0("ls ", wtRNAseqDir, "WT_RNAseq_Chris_Rep2_norm_allchrs_coverage_coord_tab.bed"), intern = T)
kss_RNAseq_Rep1 <- system(paste0("ls ", kssRNAseqDir, "kss_RNAseq_Chris_Rep1_norm_allchrs_coverage_coord_tab.bed"), intern = T)
kss_RNAseq_Rep2 <- system(paste0("ls ", kssRNAseqDir, "kss_RNAseq_Chris_Rep2_norm_allchrs_coverage_coord_tab.bed"), intern = T)
wt_H3K9me2 <- system(paste0("ls ", wtH3K9me2Dir, "WT_H3K9me2_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed"), intern = T)
kss_H3K9me2 <- system(paste0("ls ", kssH3K9me2Dir, "log2_kss_H3K9me2_ChIP_kss_H3K9me2_input_norm_allchrs_coverage_coord_tab.bed"), intern = T)
wt_SPO11oligos_RPI1 <- system(paste0("ls ", wtSPO11Dir, "log2wtSPO11oligoRPI1NakedDNA_norm_allchrs_coverage_coord_tab.bed"), intern = T)
wt_SPO11oligos_RPI3 <- system(paste0("ls ", wtSPO11Dir, "log2wtSPO11oligoRPI3NakedDNA_norm_allchrs_coverage_coord_tab.bed"), intern = T)
wt_SPO11oligos_RPI8 <- system(paste0("ls ", wtSPO11Dir, "log2wtSPO11oligoRPI8NakedDNA_norm_allchrs_coverage_coord_tab.bed"), intern = T)
kss_SPO11oligos_RPI34 <- system(paste0("ls ", kssSPO11Dir, "log2suvh456SPO11oligoRPI34NakedDNA_norm_allchrs_coverage_coord_tab.bed"), intern = T)
kss_SPO11oligos_RPI35 <- system(paste0("ls ", kssSPO11Dir, "log2suvh456SPO11oligoRPI35NakedDNA_norm_allchrs_coverage_coord_tab.bed"), intern = T)

libPaths <- list(
                 wt_RNAseq_Rep1,
                 wt_RNAseq_Rep2,
                 kss_RNAseq_Rep1,
                 kss_RNAseq_Rep2,
                 wt_H3K9me2,
                 kss_H3K9me2,
                 wt_SPO11oligos_RPI1,
                 wt_SPO11oligos_RPI3,
                 wt_SPO11oligos_RPI8,
                 kss_SPO11oligos_RPI34,
                 kss_SPO11oligos_RPI35
                )

libNames <- c(
              "wt_RNAseq_Rep1",
              "wt_RNAseq_Rep2",
              "kss_RNAseq_Rep1",
              "kss_RNAseq_Rep2",
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
                   "wt_RNAseq_Rep1" = wt_RNAseq_Rep1,
                   "wt_RNAseq_Rep2" = wt_RNAseq_Rep2,
                   "kss_RNAseq_Rep1" = kss_RNAseq_Rep1,
                   "kss_RNAseq_Rep2" = kss_RNAseq_Rep2,
                   "wt_H3K9me2" = wt_H3K9me2,
                   "kss_H3K9me2" = kss_H3K9me2,
                   "wt_SPO11oligos_RPI1" = wt_SPO11oligos_RPI1,
                   "wt_SPO11oligos_RPI3" = wt_SPO11oligos_RPI3,
                   "wt_SPO11oligos_RPI8" = wt_SPO11oligos_RPI8,
                   "kss_SPO11oligos_RPI34" = kss_SPO11oligos_RPI34,
                   "kss_SPO11oligos_RPI35" = kss_SPO11oligos_RPI35
                  )

targetMeanCov <- function(targetName, targetsGR, signalGR, libName) {
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
                            targetName,
                            "MeanCov.txt"),
              col.names = libName)
}

mclapply(seq_along(grl), function(i) {
  targetMeanCov(targetName = "TEs",
                targetsGR = TEsGR,
                signalGR = grl[[i]],
                libName = libNames[i])
}, mc.cores = 4)

covTabs <- mclapply(seq_along(libNames), function(i) {
  read.table(file = paste0(outDir,
                           libNames[i],
                           "_TEsMeanCov.txt"))
}, mc.cores = length(libNames))

TEsCovBind <- TEs
for(i in 1:length(covTabs)) {
  TEsCovBind <- cbind(TEsCovBind, covTabs[[i]])
}

write.table(TEsCovBind,
            file = paste0(outDir,
                          "TEs_mean_coverage_wt_kss_RNAseq_H3K9me2_SPO11oligos.txt"),
            sep = "\t",
            quote = F,
            row.names = F)

