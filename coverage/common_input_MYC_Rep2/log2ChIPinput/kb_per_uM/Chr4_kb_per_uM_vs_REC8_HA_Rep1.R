# Calculate factor coverage within intervals on chromosome 4 and evaluate relationship
# with kb/uM values

library(segmentSeq)

chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")

inDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/kb_per_uM/"
plotDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/kb_per_uM/plots/" 
targets <- read.table(paste0(inDir,
                             "kb_per_uM.txt"), header = T)

covDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/"
factor_norm <- read.table(paste0(covDir,
                                 "log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab_Chr4.bed"),
                          colClasses = c(NA, NA, "NULL", NA))

i <- 4
print(i)
chrCov <- factor_norm[factor_norm[,1] == chrs[i],]
chrCov <- chrCov[,3]
print(i)
covCoords <- seq(1, length(chrCov), by = 1)
covGR <- GRanges(seqnames = chrs[i],
                 ranges = IRanges(start = covCoords, width = 1),
                 strand = "*")
chrTargets <- targets[targets[,1] == chrs[i],]
targetsGR <- GRanges(seqnames = chrs[i],
                     ranges = IRanges(start = chrTargets$start,
                                      end = chrTargets$end),
                     strand = "*")
targetsOverlaps <- getOverlaps(targetsGR,
                               covGR,
                               whichOverlaps = T)
factor_norm_targetMeanCov <- sapply(targetsOverlaps,
                                    function(x) mean(chrCov[x]))
factor_norm_targetSumCov <- sapply(targetsOverlaps,
                                   function(x) sum(chrCov[x]))
chrTargets <- cbind(chrTargets, factor_norm_targetMeanCov, factor_norm_targetSumCov)

pdf(paste0(plotDir,
           "Chr4_kb_per_uM_vs_mean_REC8_HA_Rep1_sum_REC8_HA_Rep1.pdf"), height = 4, width = 8)
par(mfrow = c(1, 2), mar = c(6, 6, 2, 2), mgp = c(4, 1.5, 0))
plot(x = chrTargets$kb.uM,
     y = chrTargets$factor_norm_targetMeanCov,
     pch = 19,
     main = bquote(italic("r"[s]) ~ " = " ~ .(round(cor(chrTargets$kb.uM, chrTargets$factor_norm_targetMeanCov, method = "spearman"), digits = 2))),
     xlab = "kb/uM",
     ylab = "Mean log2-transformed REC8-HA Rep1")
plot(x = chrTargets$kb.uM,
     y = chrTargets$factor_norm_targetSumCov,
     pch = 19,
     main = bquote(italic("r"[s]) ~ " = " ~ .(round(cor(chrTargets$kb.uM, chrTargets$factor_norm_targetSumCov, method = "spearman"), digits = 2))),
     xlab = "kb/uM",
     ylab = "Total log2-transformed REC8-HA Rep1")
dev.off()



