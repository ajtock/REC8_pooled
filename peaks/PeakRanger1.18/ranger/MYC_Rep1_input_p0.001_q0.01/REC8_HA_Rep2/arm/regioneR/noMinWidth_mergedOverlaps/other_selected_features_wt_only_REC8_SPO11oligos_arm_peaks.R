#!/applications/R/R-3.5.0/bin/Rscript

# Plot bar graph of log2(observed:expected) peaks overlapping other features

# Selected other features:
# Well-positioned nucleosomes [3 for REC8, 2 for SPO11]
# SPO11-1-oligo hotspots/REC8-HA peaks [4, for REC8, 3 for SPO11]
# H3K9me2 peaks [13 for REC8, 12 for SPO11]
# Crossovers [14 for REC8, 13 for SPO11]

# Usage:
# ./other_selected_features_wt_only_REC8_SPO11oligos_arm_peaks.R "Arm REC8-HA peaks and SPO11-1-oligo hotspots" "REC8-HA" "SPO11-1" REC8_HA_Rep2 wt_SPO11oligo REC8_HA_Rep2 arm 10000

library(ggplot2)
library(ggthemes)

#plotTitle <- "Arm REC8-HA peaks and SPO11-1-oligo hotspots"
#pt1Name <- "wt REC8-HA"
#pt3Name <- "wt SPO11-1"
#pt1LibName <- "REC8_HA_Rep2" 
#pt3LibName <- "wt_SPO11oligo"
#ptOrderLibName <- "REC8_HA_Rep2"
#region <- "arm" 
## Number of permutations (randomisations) performed
#perms <- "10000"

args <- commandArgs(trailingOnly = T)
plotTitle <- args[1]
pt1Name <- args[2]
pt3Name <- args[3]
pt1LibName <- args[4]
pt3LibName <- args[5]
ptOrderLibName <- args[6]
region <- args[7]
# Number of permutations (randomisations) performed
perms <- as.character(args[8])

ptOrderDir <- paste0("/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep1_input_p0.001_q0.01/",
                     ptOrderLibName,
                     "/peri/regioneR/noMinWidth_mergedOverlaps/")
pt1Dir <- paste0("/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep1_input_p0.001_q0.01/",
                 pt1LibName, "/",
                 region, "/regioneR/noMinWidth_mergedOverlaps/")
pt3Dir <- paste0("/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/RPI1_RPI8_profiles_rangerPeaks_idr0.05/",
                 region, "/regioneR/noMinWidth_mergedOverlaps/")

plotDir <- "./bar_graphs/"
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Selected other features:
# Well-positioned nucleosomes [3 for REC8, 2 for SPO11]
# SPO11-1-oligo hotspots/REC8-HA peaks [4 for REC8, 3 for SPO11]
# H3K9me2 peaks [13 for REC8, 12 for SPO11]
# Crossovers [14 for REC8, 13 for SPO11]

otherNumbersREC8 <- c(3, 4, 13, 14)
otherNumbersSPO11 <- c(2, 3, 12, 13)
otherNames <- c("REC8_MYC_Rep1GR",
                "REC8_HA_Rep1GR",
                "nucleRnucsGR",
                "SPO11GR",
                "SPO11_ChIP4GR",
                "SPO11_ChIP13GR",
                "H3K4me3_ChIP14GR",
                "H3K4me3_ChIP15GR",
                "H3K4me1GR",
                "H3K4me2GR",
                "H3K27me1GR",
                "H3K27me3GR",
                "H3K9me2GR",
                "COsGR",
                "genesGR",
                "promotersGR",
                "terminatorsGR",
                "TSSdownstream500GR",
                "TTSupstream500GR",
                "exonsGR",
                "intronsGR",
                "TEsGR",
                "kss_hypoCHG_DMRsGR",
                "kss_hypoCHH_DMRsGR")
otherNames <- otherNames[otherNumbersREC8]
otherNamesPlot <- c("REC8-Myc",
                    "REC8-HA Rep1",
                    "Nucleosomes",
                    "SPO11-1/REC8-HA",
                    "SPO11-1 ChIP Rep1",
                    "SPO11-1 ChIP Rep2",
                    "H3K4me3 Rep2",
                    "H3K4me3 Rep3",
                    "H3K4me1",
                    "H3K4me2",
                    "H3K27me1",
                    "H3K27me3",
                    "H3K9me2",
                    "Crossovers",
                    "Genes",
                    "Gene promoters",
                    "Gene terminators",
                    "Gene 5' ends",
                    "Gene 3' ends",
                    "Gene exons",
                    "Gene introns",
                    "Transposons",
                    "Hypo-mCHG",
                    "Hypo-mCHH")
otherNamesPlot <- otherNamesPlot[otherNumbersREC8]

# Load permutation test results for peak set to be used for ordering
# of other features in bar graph
load(paste0(ptOrderDir, "permTest_", perms, "perms_",
            ptOrderLibName, "_peri_peaks_vs_others.RData"))
ptOrder <- ptPeaksOtherPerChrom
ptPeaksOtherPerChrom <- NULL

# Load permutation test results to be used for plotting
load(paste0(pt1Dir, "permTest_", perms, "perms_",
            pt1LibName, "_", region, "_peaks_vs_others.RData"))
pt1 <- ptPeaksOtherPerChrom
ptPeaksOtherPerChrom <- NULL

# Re-order pt3 and pt4 to make consistent with pt1 and pt2
load(paste0(pt3Dir, "permTest_", perms, "perms_",
            pt3LibName, "_hotspots_", region, "_vs_others.RData"))
pt3 <- ptPeaksOtherPerChrom
pt3 <- pt3[c(3, 2, 1, 4:length(pt3))]
ptPeaksOtherPerChrom <- NULL

# Permutation test results for selected features
ptOrder <- ptOrder[otherNumbersREC8]
pt1 <- pt1[otherNumbersREC8]
pt3 <- pt3[otherNumbersSPO11]

# ptOrder
ptOrder_Pval <- lapply(seq_along(ptOrder), function(x) {
  ptOrder[[x]]$numOverlaps$pval
})
ptOrder_Obs <- lapply(seq_along(ptOrder), function(x) {
  ptOrder[[x]]$numOverlaps$observed
})
ptOrder_Perm <- lapply(seq_along(ptOrder), function(x) {
  ptOrder[[x]]$numOverlaps$permuted
})
ptOrder_Exp <- lapply(seq_along(ptOrder), function(x) {
  mean(ptOrder[[x]]$numOverlaps$permuted)
})
ptOrder_log2ObsExp <- lapply(seq_along(ptOrder_Obs), function(x) {
  log2(ptOrder_Obs[[x]]/ptOrder_Exp[[x]])
})
ptOrder_Zscore <- lapply(seq_along(ptOrder), function(x) {
  ptOrder[[x]]$numOverlaps$zscore
})
ptOrder_AltHyp <- lapply(seq_along(ptOrder), function(x) {
  ptOrder[[x]]$numOverlaps$alternative
})
ptOrder_alpha0.05 <- lapply(seq_along(ptOrder_Perm), function(x) {
  if(ptOrder_AltHyp[[x]] == "greater") {
    quantile(ptOrder_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(ptOrder_Perm[[x]], 0.05)[[1]]
  }
})
ptOrder_log2alpha0.05 <- lapply(seq_along(ptOrder_alpha0.05), function(x) {
  log2(ptOrder_alpha0.05[[x]]/ptOrder_Exp[[x]])
})

# Order permutation test statistics and feature names by
# decreasing log2(observed/expected) overlaps of wt treatment peaks
ptOrder_log2ObsExp_sorted <- sort.int(unlist(ptOrder_log2ObsExp),
                                      decreasing = T)
ptOrder_log2alpha0.05_sorted <- unlist(ptOrder_log2alpha0.05[sort.int(unlist(ptOrder_log2ObsExp),
                                                                      decreasing = T,
                                                                      index.return = T)$ix])
ptOrder_otherNames_sorted <- otherNames[sort.int(unlist(ptOrder_log2ObsExp),
                                                 decreasing = T,
                                                 index.return = T)$ix]
ptOrder_otherNamesPlot_sorted <- otherNamesPlot[sort.int(unlist(ptOrder_log2ObsExp),
                                                         decreasing = T,
                                                         index.return = T)$ix]

# pt1
pt1_Pval <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$pval
})
pt1_Obs <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$observed
})
pt1_Perm <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$permuted
})
pt1_Exp <- lapply(seq_along(pt1), function(x) {
  mean(pt1[[x]]$numOverlaps$permuted)
})
pt1_log2ObsExp <- lapply(seq_along(pt1_Obs), function(x) {
  log2(pt1_Obs[[x]]/pt1_Exp[[x]])
})
pt1_Zscore <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$zscore
})
pt1_AltHyp <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$alternative
})
pt1_alpha0.05 <- lapply(seq_along(pt1_Perm), function(x) {
  if(pt1_AltHyp[[x]] == "greater") {
    quantile(pt1_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(pt1_Perm[[x]], 0.05)[[1]]
  }
})
pt1_log2alpha0.05 <- lapply(seq_along(pt1_alpha0.05), function(x) {
  log2(pt1_alpha0.05[[x]]/pt1_Exp[[x]])
})

# Order permutation test statistics and feature names by
# decreasing log2(observed/expected) overlaps of wt treatment peaks
pt1_log2ObsExp_sorted <- unlist(pt1_log2ObsExp[sort.int(unlist(ptOrder_log2ObsExp),
                                                        decreasing = T,
                                                        index.return = T)$ix])
pt1_log2alpha0.05_sorted <- unlist(pt1_log2alpha0.05[sort.int(unlist(ptOrder_log2ObsExp),
                                                              decreasing = T,
                                                              index.return = T)$ix])
pt1_otherNames_sorted <- otherNames[sort.int(unlist(ptOrder_log2ObsExp),
                                             decreasing = T,
                                             index.return = T)$ix]
pt1_otherNamesPlot_sorted <- otherNamesPlot[sort.int(unlist(ptOrder_log2ObsExp),
                                                     decreasing = T,
                                                     index.return = T)$ix]

# pt3
pt3_Pval <- lapply(seq_along(pt3), function(x) {
  pt3[[x]]$numOverlaps$pval
})
pt3_Obs <- lapply(seq_along(pt3), function(x) {
  pt3[[x]]$numOverlaps$observed
})
pt3_Perm <- lapply(seq_along(pt3), function(x) {
  pt3[[x]]$numOverlaps$permuted
})
pt3_Exp <- lapply(seq_along(pt3), function(x) {
  mean(pt3[[x]]$numOverlaps$permuted)
})
pt3_log2ObsExp <- lapply(seq_along(pt3_Obs), function(x) {
  log2(pt3_Obs[[x]]/pt3_Exp[[x]])
})
pt3_Zscore <- lapply(seq_along(pt3), function(x) {
  pt3[[x]]$numOverlaps$zscore
})
pt3_AltHyp <- lapply(seq_along(pt3), function(x) {
  pt3[[x]]$numOverlaps$alternative
})
pt3_alpha0.05 <- lapply(seq_along(pt3_Perm), function(x) {
  if(pt3_AltHyp[[x]] == "greater") {
    quantile(pt3_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(pt3_Perm[[x]], 0.05)[[1]]
  }
})
pt3_log2alpha0.05 <- lapply(seq_along(pt3_alpha0.05), function(x) {
  log2(pt3_alpha0.05[[x]]/pt3_Exp[[x]])
})

# Order permutation test statistics and feature names by
# decreasing log2(observed/expected) overlaps of wt treatment peaks
pt3_log2ObsExp_sorted <- unlist(pt3_log2ObsExp[sort.int(unlist(ptOrder_log2ObsExp),
                                                        decreasing = T,
                                                        index.return = T)$ix])
pt3_log2alpha0.05_sorted <- unlist(pt3_log2alpha0.05[sort.int(unlist(ptOrder_log2ObsExp),
                                                              decreasing = T,
                                                              index.return = T)$ix])
pt3_otherNames_sorted <- otherNames[sort.int(unlist(ptOrder_log2ObsExp),
                                             decreasing = T,
                                             index.return = T)$ix]
pt3_otherNamesPlot_sorted <- otherNamesPlot[sort.int(unlist(ptOrder_log2ObsExp),
                                                     decreasing = T,
                                                     index.return = T)$ix]

# Combine in data.frame
df <- data.frame(Sample = rep(c(pt1Name, pt3Name),
                              each = length(ptOrder_log2ObsExp_sorted)),
                 Feature = rep(pt1_otherNamesPlot_sorted, 2),
                 log2ObsExp = c(pt1_log2ObsExp_sorted,
                                pt3_log2ObsExp_sorted),
                 log2alpha0.05 = c(pt1_log2alpha0.05_sorted,
                                   pt3_log2alpha0.05_sorted))

df$Feature <- factor(df$Feature,
                     levels = pt1_otherNamesPlot_sorted)
df$Sample <- factor(df$Sample,
                    levels = c(pt1Name, pt3Name))

bp <- ggplot(data = df,
             mapping = aes(x = Feature,
                           y = log2ObsExp,
                           fill = Sample)) +
  geom_bar(stat = "identity",
           position = position_dodge()) +
  scale_fill_manual(name = "",
                    values = c("red",
                               "dodgerblue2"),
                    labels = c(pt1Name,
                               pt3Name)) + 
  geom_point(mapping = aes(Feature, log2alpha0.05),
             position = position_dodge(0.90),
             shape = "-", colour  = "grey80", size = 8) +
  labs(y = expression("Log"[2]*"(observed/expected) peak overlap")) +
  scale_y_continuous(limits = c(-3.1, 1)) +
  scale_x_discrete(position = "top") +
  guides(fill = guide_legend(direction = "horizontal",
                             label.position = "top",
                             label.theme = element_text(size = 10, hjust = 0, vjust = 0.5, angle = 90),
                             nrow = 1,
                             byrow = TRUE)) +
  theme_bw() +
  theme(axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.ticks.y = element_line(size = 0.75, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black", hjust = 0.5, vjust = 0.5, angle = 90),
        axis.title.y = element_text(size = 10, colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black", hjust = 0, vjust = 0.5, angle = 45),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        #legend.position = c(0.05, 0.30),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(colour = "transparent",
                                  fill = "transparent"),
        plot.margin = unit(c(5.5, 5.5, 6.0, 5.5), "pt"),
        plot.title = element_text(size = 6, colour = "black", hjust = 0.25)) +
  ggtitle(paste0(plotTitle, " (", prettyNum(as.character(perms),
                                            big.mark = ",", trim = "T"),
                 " permutations)"))
ggsave(paste0(plotDir, "barplot_selected_other_features_permTestResults_",
              as.character(perms), "perms_",
              "log2_Observed_Expected_",
              pt1LibName, "_",
              pt3LibName, "_",
              region, "_peaks_v250619.pdf"),
       plot = bp,
       height = 5, width = 3.4)
save(bp,
     file = paste0(plotDir, "barplot_selected_other_features_permTestResults_",
                   as.character(perms), "perms_",
                   "log2_Observed_Expected_",
                   pt1LibName, "_",
                   pt3LibName, "_",
                   region, "_peaks.RData"))
