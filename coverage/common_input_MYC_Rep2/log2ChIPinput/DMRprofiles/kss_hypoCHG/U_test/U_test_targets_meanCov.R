#!/applications/R/R-3.3.2/bin/Rscript

# Calculate mean coverage between target start and end coordinates of targets for
# each wt and mutant library and tabulate
# This will enable Mann-Whitney U tests comparing each wt library with each kss library

library(parallel)
library(GenomicRanges)

outDir <- "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/DMRprofiles/kss_hypoCHG/U_test/"

# Import DMRs as GRanges object
DMRs <- read.table("/home/ajt200/BS_Seq/Stroud_2013/DMRs/suvh456_hypoCHG_DMR_vs3reps_min4filter_mg200.bed", header = F)
DMRsGR <- GRanges(seqnames = DMRs[,1],
                  ranges = IRanges(start = DMRs[,2], end = DMRs[,3]),
                  strand = "*")
seqlevels(DMRsGR) <- sub("chr", "Chr", seqlevels(DMRsGR))
print("***********DMRs***********")
print(DMRsGR)

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
                          "kss_hypoCHG_DMRs_meanCov_wt_kss_REC8_H3K9me2_SPO11oligos_RNAseq.txt"),
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
#W = 144750000, p-value < 2.2e-16
#alternative hypothesis: true location shift is greater than 0
REC8_HA_Rep2_u_test <- wilcox.test(x = c(targetsCovBind$wt_REC8_HA_Rep2),
                                   y = c(targetsCovBind$kss_REC8_HA_Rep1),
                                   alternative = "greater")
print("REC8_HA_Rep2_u_test")
print(REC8_HA_Rep1_u_test)
#W = 1.4e+08, p-value < 2.2e-16
#alternative hypothesis: true location shift is greater than 0
REC8_MYC_Rep1_u_test <- wilcox.test(x = c(targetsCovBind$wt_REC8_MYC_Rep1),
                                    y = c(targetsCovBind$kss_REC8_HA_Rep1),
                                    alternative = "greater")
print("REC8_MYC_Rep1_u_test")
print(REC8_MYC_Rep1_u_test)
#W = 124520000, p-value = 7.42e-06
#alternative hypothesis: true location shift is greater than 0

# H3K9me2
H3K9me2_u_test <- wilcox.test(x = c(targetsCovBind$wt_H3K9me2),
                              y = c(targetsCovBind$kss_H3K9me2),
                              alternative = "greater")
print("H3K9me2_u_test")
print(H3K9me2_u_test)
#W = 2.07e+08, p-value < 2.2e-16
#alternative hypothesis: true location shift is greater than 0

# SPO11-1-oligos
SPO11_1_oligos_Rep1_u_test <- wilcox.test(x = c(targetsCovBind$wt_SPO11_1_oligos_RPI1),
                                          y = c(targetsCovBind$kss_SPO11_1_oligos_RPI34),
                                          alternative = "less")
print("SPO11_1_oligos_Rep1_u_test")
print(SPO11_1_oligos_Rep1_u_test)
#W = 80197000, p-value < 2.2e-16
#alternative hypothesis: true location shift is less than 0
SPO11_1_oligos_Rep2_u_test <- wilcox.test(x = c(targetsCovBind$wt_SPO11_1_oligos_RPI3),
                                          y = c(targetsCovBind$kss_SPO11_1_oligos_RPI35),
                                          alternative = "less")
print("SPO11_1_oligos_Rep2_u_test")
print(SPO11_1_oligos_Rep2_u_test)
#W = 89540000, p-value < 2.2e-16
#alternative hypothesis: true location shift is less than 0

# RNA-seq
RNAseq_Rep1_u_test <- wilcox.test(x = c(targetsCovBind$wt_RNAseq_Rep1),
                                  y = c(targetsCovBind$kss_RNAseq_Rep1),
                                  alternative = "less")
print("RNAseq_Rep1_u_test")
print(RNAseq_Rep1_u_test)
#W = 102740000, p-value < 2.2e-16
#alternative hypothesis: true location shift is less than 0
RNAseq_Rep2_u_test <- wilcox.test(x = c(targetsCovBind$wt_RNAseq_Rep2),
                                  y = c(targetsCovBind$kss_RNAseq_Rep2),
                                  alternative = "less")
print("RNAseq_Rep2_u_test")
print(RNAseq_Rep2_u_test)
#W = 105910000, p-value < 2.2e-16
#alternative hypothesis: true location shift is less than 0

