#!/applications/R/R-3.3.2/bin/Rscript

# Calculate mean library-size-normalised wt and kss
# coverage between the start and end of each peak
# Select peaks with greatest change in kss (top 10%)

library(dplyr)
library(GenomicRanges)

outDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/genome_wide/RNAseq_Rep1_kss_wt_diff/"

# Chromosome definitions
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)

allCovPeaks <- read.table(paste0(outDir,
                                 "REC8_HA_Rep1_peaks_mean_RNAseq_Rep1_kss_wt_diff_cov.txt"),
                          header = T)
allCovPeaks_percent_rank <- mutate(allCovPeaks,
                                   diff_kss_wt_norm_peakRank = percent_rank(diff_kss_wt_norm_peakCov))
allCovPeaks_top10percent <- allCovPeaks_percent_rank[allCovPeaks_percent_rank$diff_kss_wt_norm_peakRank >= 0.9,]

write.table(allCovPeaks_top10percent,
            paste0(outDir,
                   "REC8_HA_Rep1_peaks_mean_RNAseq_Rep1_kss_wt_diff_cov_top10percent.txt"),
            sep = "\t", quote = F)

allCovPeaks_top10percent_arm <- NULL
allCovPeaks_top10percent_peri <- NULL
for(i in 1:5) {
  chr_allCovPeaks_top10percent <- allCovPeaks_top10percent[allCovPeaks_top10percent$Chr == chrs[i],]
  chr_allCovPeaks_top10percent_arm <- chr_allCovPeaks_top10percent %>%
    filter(End < pericenStart[i] | Start > pericenEnd[i])
  chr_allCovPeaks_top10percent_peri <- chr_allCovPeaks_top10percent %>%
    filter(End > pericenStart[i] & Start < pericenEnd[i])
  allCovPeaks_top10percent_arm <- rbind(allCovPeaks_top10percent_arm,
                                        chr_allCovPeaks_top10percent_arm)
  allCovPeaks_top10percent_peri <- rbind(allCovPeaks_top10percent_peri,
                                         chr_allCovPeaks_top10percent_peri)
}

allCovPeaks_top10percentGR <- GRanges(seqnames = allCovPeaks_top10percent$Chr,
                                      ranges = IRanges(start = allCovPeaks_top10percent$Start,
                                                       end = allCovPeaks_top10percent$End),
                                      strand = "*",
                                      diff_kss_wt_norm_peakCov = allCovPeaks_top10percent$diff_kss_wt_norm_peakCov, 
                                      diff_kss_wt_norm_peakRank = allCovPeaks_top10percent$diff_kss_wt_norm_peakRank)
save(allCovPeaks_top10percentGR,
     file = paste0(outDir,
                   "REC8_HA_Rep1_peaks_mean_RNAseq_Rep1_kss_wt_diff_cov_top10percentGR.RData"))

allCovPeaks_top10percent_armGR <- GRanges(seqnames = allCovPeaks_top10percent_arm$Chr,
                                          ranges = IRanges(start = allCovPeaks_top10percent_arm$Start,
                                                           end = allCovPeaks_top10percent_arm$End),
                                          strand = "*",
                                          diff_kss_wt_norm_peakCov = allCovPeaks_top10percent_arm$diff_kss_wt_norm_peakCov,
                                          diff_kss_wt_norm_peakRank = allCovPeaks_top10percent_arm$diff_kss_wt_norm_peakRank)
save(allCovPeaks_top10percent_armGR,
     file = paste0(outDir,
                   "REC8_HA_Rep1_peaks_mean_RNAseq_Rep1_kss_wt_diff_cov_top10percent_armGR.RData"))

allCovPeaks_top10percent_periGR <- GRanges(seqnames = allCovPeaks_top10percent_peri$Chr,
                                          ranges = IRanges(start = allCovPeaks_top10percent_peri$Start,
                                                           end = allCovPeaks_top10percent_peri$End),
                                          strand = "*",
                                          diff_kss_wt_norm_peakCov = allCovPeaks_top10percent_peri$diff_kss_wt_norm_peakCov,
                                          diff_kss_wt_norm_peakRank = allCovPeaks_top10percent_peri$diff_kss_wt_norm_peakRank)
save(allCovPeaks_top10percent_periGR,
     file = paste0(outDir,
                   "REC8_HA_Rep1_peaks_mean_RNAseq_Rep1_kss_wt_diff_cov_top10percent_periGR.RData"))

