#!/applications/R/R-3.5.0/bin/Rscript

# Plot bar graph of log2(observed:expected) peaks overlapping other features

# Usage:
# ./other_features_REC8_vs_control_arm_peaks.R "REC8-Myc Rep1 and Control-Myc Rep1 arm peaks" "REC8-Myc Rep1" "Control-Myc Rep1" REC8_MYC_Rep1 control_MYC_Rep1 REC8_MYC_Rep1 arm 10000

library(ggplot2)
library(ggthemes)

#plotTitle <- "REC8-Myc Rep1 and control-Myc Rep1 arm peaks"
#wt_pt1Name <- "REC8-Myc Rep1"
#wt_pt2Name <- "Control-Myc Rep1"
#wt_pt1LibName <- "REC8_MYC_Rep1" 
#wt_pt2LibName <- "control_MYC_Rep1"
#wt_ptOrderLibName <- "REC8_MYC_Rep1"
#region <- "arm" 
## Number of permutations (randomisations) performed
#perms <- "10000"

args <- commandArgs(trailingOnly = T)
plotTitle <- args[1]
wt_pt1Name <- args[2]
wt_pt2Name <- args[3]
wt_pt1LibName <- args[4]
wt_pt2LibName <- args[5]
wt_ptOrderLibName <- args[6]
region <- args[7]
# Number of permutations (randomisations) performed
perms <- as.character(args[8])

wt_ptOrderDir <- paste0("/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep1_input_p0.001_q0.01/",
                        wt_ptOrderLibName, "/",
                        region, "/regioneR/noMinWidth_mergedOverlaps/")
wt_pt1Dir <- paste0("/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep1_input_p0.001_q0.01/",
                    wt_pt1LibName, "/",
                    region, "/regioneR/noMinWidth_mergedOverlaps/")
wt_pt2Dir <- paste0("/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep1_input_p0.001_q0.01/",
                    wt_pt2LibName, "/",
                    region, "/regioneR/noMinWidth_mergedOverlaps/")

plotDir <- "./bar_graphs/"
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

otherNames <- c(paste0(wt_pt2LibName, "GR"),
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
                "TEsGR")

otherNamesPlot <- c(paste0(wt_pt1Name, "/", wt_pt2Name),
                    "Nucleosomes",
                    "SPO11-1-oligo hotspots",
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
                    "Transposons")

# Load permutation test results for peak set to be used for ordering
# of other features in bar graph
load(paste0(wt_ptOrderDir, "permTest_", perms, "perms_",
            wt_ptOrderLibName, "_arm_peaks_vs_others.RData"))
wt_ptOrder <- ptPeaksOtherPerChrom
ptPeaksOtherPerChrom <- NULL

# Load permutation test results to be used for plotting
load(paste0(wt_pt1Dir, "permTest_", perms, "perms_",
            wt_pt1LibName, "_arm_peaks_vs_others.RData"))
wt_pt1 <- ptPeaksOtherPerChrom
ptPeaksOtherPerChrom <- NULL

load(paste0(wt_pt2Dir, "permTest_", perms, "perms_",
            wt_pt2LibName, "_arm_peaks_vs_others.RData"))
wt_pt2 <- ptPeaksOtherPerChrom
ptPeaksOtherPerChrom <- NULL

# wt_ptOrder
wt_ptOrder_Pval <- lapply(seq_along(wt_ptOrder), function(x) {
  wt_ptOrder[[x]]$numOverlaps$pval
})
wt_ptOrder_Obs <- lapply(seq_along(wt_ptOrder), function(x) {
  wt_ptOrder[[x]]$numOverlaps$observed
})
wt_ptOrder_Perm <- lapply(seq_along(wt_ptOrder), function(x) {
  wt_ptOrder[[x]]$numOverlaps$permuted
})
wt_ptOrder_Exp <- lapply(seq_along(wt_ptOrder), function(x) {
  mean(wt_ptOrder[[x]]$numOverlaps$permuted)
})
wt_ptOrder_log2ObsExp <- lapply(seq_along(wt_ptOrder_Obs), function(x) {
  log2(wt_ptOrder_Obs[[x]]/wt_ptOrder_Exp[[x]])
})
wt_ptOrder_Zscore <- lapply(seq_along(wt_ptOrder), function(x) {
  wt_ptOrder[[x]]$numOverlaps$zscore
})
wt_ptOrder_AltHyp <- lapply(seq_along(wt_ptOrder), function(x) {
  wt_ptOrder[[x]]$numOverlaps$alternative
})
wt_ptOrder_alpha0.05 <- lapply(seq_along(wt_ptOrder_Perm), function(x) {
  if(wt_ptOrder_AltHyp[[x]] == "greater") {
    quantile(wt_ptOrder_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(wt_ptOrder_Perm[[x]], 0.05)[[1]]
  }
})
wt_ptOrder_log2alpha0.05 <- lapply(seq_along(wt_ptOrder_alpha0.05), function(x) {
  log2(wt_ptOrder_alpha0.05[[x]]/wt_ptOrder_Exp[[x]])
})

# Order permutation test statistics and feature names by
# decreasing log2(observed/expected) overlaps of wt treatment peaks
wt_ptOrder_log2ObsExp_sorted <- sort.int(unlist(wt_ptOrder_log2ObsExp),
                                         decreasing = T)
wt_ptOrder_log2alpha0.05_sorted <- unlist(wt_ptOrder_log2alpha0.05[sort.int(unlist(wt_ptOrder_log2ObsExp),
                                                                            decreasing = T,
                                                                            index.return = T)$ix])
wt_ptOrder_otherNames_sorted <- otherNames[sort.int(unlist(wt_ptOrder_log2ObsExp),
                                                    decreasing = T,
                                                    index.return = T)$ix]
wt_ptOrder_otherNamesPlot_sorted <- otherNamesPlot[sort.int(unlist(wt_ptOrder_log2ObsExp),
                                                            decreasing = T,
                                                            index.return = T)$ix]

# wt_pt1
wt_pt1_Pval <- lapply(seq_along(wt_pt1), function(x) {
  wt_pt1[[x]]$numOverlaps$pval
})
wt_pt1_Obs <- lapply(seq_along(wt_pt1), function(x) {
  wt_pt1[[x]]$numOverlaps$observed
})
wt_pt1_Perm <- lapply(seq_along(wt_pt1), function(x) {
  wt_pt1[[x]]$numOverlaps$permuted
})
wt_pt1_Exp <- lapply(seq_along(wt_pt1), function(x) {
  mean(wt_pt1[[x]]$numOverlaps$permuted)
})
wt_pt1_log2ObsExp <- lapply(seq_along(wt_pt1_Obs), function(x) {
  log2(wt_pt1_Obs[[x]]/wt_pt1_Exp[[x]])
})
wt_pt1_Zscore <- lapply(seq_along(wt_pt1), function(x) {
  wt_pt1[[x]]$numOverlaps$zscore
})
wt_pt1_AltHyp <- lapply(seq_along(wt_pt1), function(x) {
  wt_pt1[[x]]$numOverlaps$alternative
})
wt_pt1_alpha0.05 <- lapply(seq_along(wt_pt1_Perm), function(x) {
  if(wt_pt1_AltHyp[[x]] == "greater") {
    quantile(wt_pt1_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(wt_pt1_Perm[[x]], 0.05)[[1]]
  }
})
wt_pt1_log2alpha0.05 <- lapply(seq_along(wt_pt1_alpha0.05), function(x) {
  log2(wt_pt1_alpha0.05[[x]]/wt_pt1_Exp[[x]])
})

# Order permutation test statistics and feature names by
# decreasing log2(observed/expected) overlaps of wt treatment peaks
wt_pt1_log2ObsExp_sorted <- unlist(wt_pt1_log2ObsExp[sort.int(unlist(wt_ptOrder_log2ObsExp),
                                                              decreasing = T,
                                                              index.return = T)$ix])
wt_pt1_log2alpha0.05_sorted <- unlist(wt_pt1_log2alpha0.05[sort.int(unlist(wt_ptOrder_log2ObsExp),
                                                                    decreasing = T,
                                                                    index.return = T)$ix])
wt_pt1_otherNames_sorted <- otherNames[sort.int(unlist(wt_ptOrder_log2ObsExp),
                                                decreasing = T,
                                                index.return = T)$ix]
wt_pt1_otherNamesPlot_sorted <- otherNamesPlot[sort.int(unlist(wt_ptOrder_log2ObsExp),
                                                        decreasing = T,
                                                        index.return = T)$ix]

# wt_pt2
wt_pt2_Pval <- lapply(seq_along(wt_pt2), function(x) {
  wt_pt2[[x]]$numOverlaps$pval
})
wt_pt2_Obs <- lapply(seq_along(wt_pt2), function(x) {
  wt_pt2[[x]]$numOverlaps$observed
})
wt_pt2_Perm <- lapply(seq_along(wt_pt2), function(x) {
  wt_pt2[[x]]$numOverlaps$permuted
})
wt_pt2_Exp <- lapply(seq_along(wt_pt2), function(x) {
  mean(wt_pt2[[x]]$numOverlaps$permuted)
})
wt_pt2_log2ObsExp <- lapply(seq_along(wt_pt2_Obs), function(x) {
  log2(wt_pt2_Obs[[x]]/wt_pt2_Exp[[x]])
})
wt_pt2_Zscore <- lapply(seq_along(wt_pt2), function(x) {
  wt_pt2[[x]]$numOverlaps$zscore
})
wt_pt2_AltHyp <- lapply(seq_along(wt_pt2), function(x) {
  wt_pt2[[x]]$numOverlaps$alternative
})
wt_pt2_alpha0.05 <- lapply(seq_along(wt_pt2_Perm), function(x) {
  if(wt_pt2_AltHyp[[x]] == "greater") {
    quantile(wt_pt2_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(wt_pt2_Perm[[x]], 0.05)[[1]]
  }
})
wt_pt2_log2alpha0.05 <- lapply(seq_along(wt_pt2_alpha0.05), function(x) {
  log2(wt_pt2_alpha0.05[[x]]/wt_pt2_Exp[[x]])
})

# Order permutation test statistics and feature names by
# decreasing log2(observed/expected) overlaps of wt treatment peaks
wt_pt2_log2ObsExp_sorted <- unlist(wt_pt2_log2ObsExp[sort.int(unlist(wt_ptOrder_log2ObsExp),
                                                              decreasing = T,
                                                              index.return = T)$ix])
wt_pt2_log2alpha0.05_sorted <- unlist(wt_pt2_log2alpha0.05[sort.int(unlist(wt_ptOrder_log2ObsExp),
                                                                    decreasing = T,
                                                                    index.return = T)$ix])
wt_pt2_otherNames_sorted <- otherNames[sort.int(unlist(wt_ptOrder_log2ObsExp),
                                                decreasing = T,
                                                index.return = T)$ix]
wt_pt2_otherNamesPlot_sorted <- otherNamesPlot[sort.int(unlist(wt_ptOrder_log2ObsExp),
                                                        decreasing = T,
                                                        index.return = T)$ix]

df <- data.frame(Sample = rep(c(wt_pt1Name, wt_pt2Name),
                              each = length(wt_ptOrder_log2ObsExp_sorted)),
                 Feature = rep(wt_pt1_otherNamesPlot_sorted, 2),
                 log2ObsExp = c(wt_pt1_log2ObsExp_sorted,
                                wt_pt2_log2ObsExp_sorted),
                 log2alpha0.05 = c(wt_pt1_log2alpha0.05_sorted,
                                   wt_pt2_log2alpha0.05_sorted))

df$Feature <- factor(df$Feature,
                               levels = wt_pt1_otherNamesPlot_sorted)
df$Sample <- factor(df$Sample,
                    levels = c(wt_pt1Name, wt_pt2Name))

bp <- ggplot(data = df,
             mapping = aes(x = Feature,
                           y = log2ObsExp,
                           fill = Sample)) +
      geom_bar(stat = "identity",
               position = position_dodge()) +
      scale_fill_manual(name = "Sample",
                        values = c("red",
                                   "grey30"),
                        labels = c(wt_pt1Name,
                                   wt_pt2Name)) +
      geom_point(mapping = aes(Feature, log2alpha0.05),
                 position = position_dodge(0.9),
                 shape = "-", colour  = "grey70", size = 3.5) +
      labs(x = "Feature",
           y = expression("Log"[2]*"(observed:expected) peak overlap")) +
      theme_bw() +
      theme(axis.line.y = element_line(size = 0.5, colour = "black"),
            axis.ticks.y = element_line(size = 0.25, colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.ticks.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 7),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      ggtitle(paste0(plotTitle, " (", as.character(perms), " permutations)"))
ggsave(paste0(plotDir, "barplot_other_features_permTestResults_",
              as.character(perms), "perms_",
              "log2_Observed_Expected_",
               wt_pt1LibName, "_", wt_pt2LibName, "_", region, "_peaks.pdf"),
       plot = bp,
       height = 5, width = 9)
save(bp,
     file = paste0(plotDir, "barplot_other_features_permTestResults_",
                   as.character(perms), "perms_",
                  "log2_Observed_Expected_",
                   wt_pt1LibName, "_", wt_pt2LibName, "_", region, "_peaks.RData"))
