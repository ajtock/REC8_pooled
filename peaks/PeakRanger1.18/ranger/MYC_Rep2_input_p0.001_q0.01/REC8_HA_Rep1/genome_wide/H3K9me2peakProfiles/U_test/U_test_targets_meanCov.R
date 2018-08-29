#!/applications/R/R-3.3.2/bin/Rscript

# Calculate mean coverage between target start and end coordinates of targets for
# each wt and mutant library and tabulate
# This will enable Mann-Whitney U tests comparing each wt library with each kss library

library(parallel)
library(GenomicRanges)

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

libNames <- list(
                 "wt_REC8_HA_Rep1",
                 "wt_REC8_HA_Rep2",
                 "wt_REC8_MYC_Rep1",
                 "kss_REC8_HA_Rep1",
                 "wt_H3K9me2",
                 "kss_H3K9me2",
                 "wt_SPO11_1_oligos_RPI1",
                 "wt_SPO11_1_oligos_RPI3",
                 "wt_SPO11_1_oligos_RPI8",
                 "kss_SPO11_1_oligos_RPI34",
                 "kss_SPO11_1_oligos_RPI35",
                 "wt_RNAseq_Rep1",
                 "wt_RNAseq_Rep2",
                 "kss_RNAseq_Rep1",
                 "kss_RNAseq_Rep2"
                )

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
                          "H3K9me2_peaks_meanCov_wt_kss_REC8_H3K9me2_SPO11oligos_RNAseq.txt"),
            sep = "\t",
            quote = F,
            row.names = F)

# U tests
# REC8
REC8_HA_Rep1_u_test <- wilcox.test(x = c(targetsCovBind$wt_REC8_HA_Rep1),
                                   y = c(targetsCovBind$kss_REC8_HA_Rep1),
                                   alternative = "greater")
print("REC8_HA_Rep1_u_test")
print(REC8_HA_Rep1_u_test)
#W = 222040000, p-value < 2.2e-16
#alternative hypothesis: true location shift is greater than 0
REC8_HA_Rep2_u_test <- wilcox.test(x = c(targetsCovBind$wt_REC8_HA_Rep2),
                                   y = c(targetsCovBind$kss_REC8_HA_Rep1),
                                   alternative = "greater")
print("REC8_HA_Rep2_u_test")
print(REC8_HA_Rep1_u_test)
#W = 222040000, p-value < 2.2e-16
#alternative hypothesis: true location shift is greater than 0
REC8_MYC_Rep1_u_test <- wilcox.test(x = c(targetsCovBind$wt_REC8_MYC_Rep1),
                                    y = c(targetsCovBind$kss_REC8_HA_Rep1),
                                    alternative = "greater")
print("REC8_MYC_Rep1_u_test")
print(REC8_MYC_Rep1_u_test)
#W = 204780000, p-value = 0.8107
#alternative hypothesis: true location shift is greater than 0

# H3K9me2
H3K9me2_u_test <- wilcox.test(x = c(targetsCovBind$wt_H3K9me2),
                              y = c(targetsCovBind$kss_H3K9me2),
                              alternative = "greater")
print("H3K9me2_u_test")
print(H3K9me2_u_test)
#W = 240020000, p-value < 2.2e-16
#alternative hypothesis: true location shift is greater than 0

# SPO11-1-oligos
SPO11_1_oligos_Rep1_u_test <- wilcox.test(x = c(targetsCovBind$wt_SPO11_1_oligos_RPI1),
                                          y = c(targetsCovBind$kss_SPO11_1_oligos_RPI34),
                                          alternative = "less")
print("SPO11_1_oligos_Rep1_u_test")
print(SPO11_1_oligos_Rep1_u_test)
#W = 182860000, p-value < 2.2e-16
#alternative hypothesis: true location shift is less than 0
SPO11_1_oligos_Rep2_u_test <- wilcox.test(x = c(targetsCovBind$wt_SPO11_1_oligos_RPI3),
                                          y = c(targetsCovBind$kss_SPO11_1_oligos_RPI35),
                                          alternative = "less")
print("SPO11_1_oligos_Rep2_u_test")
print(SPO11_1_oligos_Rep2_u_test)
#W = 193440000, p-value < 2.2e-16
#alternative hypothesis: true location shift is less than 0

# RNA-seq
RNAseq_Rep1_u_test <- wilcox.test(x = c(targetsCovBind$wt_RNAseq_Rep1),
                                  y = c(targetsCovBind$kss_RNAseq_Rep1),
                                  alternative = "less")
print("RNAseq_Rep1_u_test")
print(RNAseq_Rep1_u_test)
#W = 197570000, p-value = 5.305e-14
#alternative hypothesis: true location shift is less than 0
RNAseq_Rep2_u_test <- wilcox.test(x = c(targetsCovBind$wt_RNAseq_Rep2),
                                  y = c(targetsCovBind$kss_RNAseq_Rep2),
                                  alternative = "less")
print("RNAseq_Rep2_u_test")
print(RNAseq_Rep2_u_test)
#W = 200120000, p-value = 1.098e-07
#alternative hypothesis: true location shift is less than 0

