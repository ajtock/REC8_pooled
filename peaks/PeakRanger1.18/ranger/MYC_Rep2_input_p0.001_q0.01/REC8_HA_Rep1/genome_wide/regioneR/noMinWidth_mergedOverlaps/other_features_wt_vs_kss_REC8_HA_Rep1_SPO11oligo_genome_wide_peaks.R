#!/applications/R/R-3.3.2/bin/Rscript

# Plot bar chart of log2(observed:expected) peaks overlapping other features

# Usage:
# Rscript other_features_wt_vs_kss_REC8_HA_Rep1_SPO11oligo_genome_wide_peaks.R "REC8 peaks and SPO11-1-oligo hotspots" "wt REC8-HA Rep1" "kss REC8-HA Rep1" "wt SPO11-1-oligos" "kss SPO11-1-oligos" wt_REC8_HA_Rep1 kss_REC8_HA_Rep1 wt_SPO11oligo kss_SPO11oligo 10000 

library(ggplot2)
library(ggthemes)

dataName <- "REC8 peaks and SPO11-1-oligo hotspots"
pt1REC8Name <- "wt REC8-HA Rep1"
pt2REC8Name <- "kss REC8-HA Rep1"
pt1Name <- "wt SPO11-1-oligos" 
pt2Name <-  "kss SPO11-1-oligos" 
pt1REC8LibName <- "wt_REC8_HA_Rep1"
pt2REC8LibName <- "kss_REC8_HA_Rep1"
pt1LibName <- "wt_SPO11oligo"
pt2LibName <- "kss_SPO11oligo"
# Number of permutations (randomisations) performed
perms <- 10000

args <- commandArgs(trailingOnly = T)
dataName <- as.character(args[1])
pt1REC8Name <- as.character(args[2])
pt2REC8Name <- as.character(args[3])
pt1Name <- as.character(args[4])
pt2Name <- as.character(args[5])
pt1REC8LibName <- as.character(args[6])
pt2REC8LibName <- as.character(args[7])
pt1LibName <- as.character(args[8])
pt2LibName <- as.character(args[9])
# Number of permutations (randomisations) performed
perms <- as.numeric(args[10])

REC8Dir <- "/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/genome_wide/regioneR/noMinWidth_mergedOverlaps/"
REC8Dir1 <- "/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/genome_wide/regioneR/noMinWidth_mergedOverlaps/"
REC8Dir2 <- "/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/kss_REC8_HA_Rep1/genome_wide/regioneR/noMinWidth_mergedOverlaps/"
inDir1 <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/RPI1_RPI8_profiles_rangerPeaks_idr0.05/genome_wide/regioneR/"
inDir2 <- "/projects/ajt200/BAM_masters/SPO11-oligo/suvh456/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/RPI34_RPI35_profiles_rangerPeaks_idr0.05/genome_wide/regioneR/"
plotDir <- "/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/genome_wide/regioneR/noMinWidth_mergedOverlaps/plots/"

otherNames <- c("REC8_HA_Rep2GR", "REC8_MYC_Rep1GR", "kss_REC8_HA_Rep1GR",
                "nucleRnucsGR", "rangernucsGR",
                "SPO11GR", "SPO11_ChIP4GR", "SPO11_ChIP13GR",
                "H3K4me3GR", "H3K9me2GR", "H3K9me2GRbcp", "COsGR",
                "genesGR", "promotersGR", "terminatorsGR",
                "TSSdownstream500GR", "TTSupstream500GR",
                "exonsGR", "intronsGR", "TEsGR")
otherNamesPlot <- c("REC8-HA Rep2 peaks", "REC8-MYC Rep1 peaks", "kss REC8-HA Rep1 peaks",
                    "Nucleosomes", "Nucleosomes (ranger)",
                    "SPO11-1-oligo hotspots/REC8 peaks", "SPO11-1 ChIP peaks", "SPO11-1 ChIP13 peaks",
                    "H3K4me3 peaks", "H3K9me2 peaks", "H3K9me2 peaks (bcp)", "Crossovers",
                    "Genes", "Gene promoters", "Gene terminators",
                    "Gene 5' ends", "Gene 3' ends",
                    "Gene exons", "Gene introns", "Transposons")

otherNames <- otherNames[c(-1:-3, -5, -8, -11)]
otherNamesPlot <- otherNamesPlot[c(-1:-3, -5, -8, -11)]

load(paste0(REC8Dir, "permTest_REC8_HA_Rep1_rangerPeaks_vs_others.RData"))
pt_REC8 <- ptPeaksOtherPerChrom
ptPeaksOtherPerChrom <- NULL
pt_REC8 <- pt_REC8[c(-1:-3, -5, -8, -11)]

load(paste0(REC8Dir1, "permTest_REC8_HA_Rep1_rangerPeaks_vs_others.RData"))
pt1_REC8 <- ptPeaksOtherPerChrom
ptPeaksOtherPerChrom <- NULL
pt1_REC8 <- pt1_REC8[c(-1:-3, -5, -8, -11)]

load(paste0(REC8Dir2, "permTest_kss_REC8_HA_Rep1_rangerPeaks_vs_others.RData"))
pt2_REC8 <- ptPeaksOtherPerChrom
ptPeaksOtherPerChrom <- NULL
pt2_REC8 <- pt2_REC8[c(-1:-3, -5, -8, -11)]

load(paste0(inDir1, "permTest_SPO11_RPI1_RPI8_rangerPeaks_vs_others.RData"))
pt1 <- ptPeaksOtherPerChrom
ptPeaksOtherPerChrom <- NULL
pt1 <- pt1[c(-1:-3, -5, -8, -11)]

load(paste0(inDir2, "permTest_kss_SPO11_RPI34_RPI35_rangerPeaks_vs_others.RData"))
pt2 <- ptPeaksOtherPerChrom
ptPeaksOtherPerChrom <- NULL
pt2 <- pt2[c(-1:-3, -5, -8, -11)]

# pt_REC8
pt_REC8_Pval <- lapply(seq_along(pt_REC8), function(x) {
  pt_REC8[[x]]$numOverlaps$pval
})
pt_REC8_Obs <- lapply(seq_along(pt_REC8), function(x) {
  pt_REC8[[x]]$numOverlaps$observed
})
pt_REC8_Perm <- lapply(seq_along(pt_REC8), function(x) {
  pt_REC8[[x]]$numOverlaps$permuted
})
pt_REC8_Exp <- lapply(seq_along(pt_REC8), function(x) {
  mean(pt_REC8[[x]]$numOverlaps$permuted)
})
pt_REC8_log2ObsExp <- lapply(seq_along(pt_REC8_Obs), function(x) {
  log2(pt_REC8_Obs[[x]]/pt_REC8_Exp[[x]])
})
pt_REC8_Zscore <- lapply(seq_along(pt_REC8), function(x) {
  pt_REC8[[x]]$numOverlaps$zscore
})
pt_REC8_AltHyp <- lapply(seq_along(pt_REC8), function(x) {
  pt_REC8[[x]]$numOverlaps$alternative
})
pt_REC8_alpha0.05 <- lapply(seq_along(pt_REC8_Perm), function(x) {
  if(pt_REC8_AltHyp[[x]] == "greater") {
    quantile(pt_REC8_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(pt_REC8_Perm[[x]], 0.05)[[1]]
  }
})
pt_REC8_log2alpha0.05 <- lapply(seq_along(pt_REC8_alpha0.05), function(x) {
  log2(pt_REC8_alpha0.05[[x]]/pt_REC8_Exp[[x]])
})

pt_REC8_log2ObsExp_sorted <- sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T)
pt_REC8_log2alpha0.05_sorted <- unlist(pt_REC8_log2alpha0.05[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix])
pt_REC8_otherNames_sorted <- otherNames[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix]
pt_REC8_otherNamesPlot_sorted <- otherNamesPlot[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix]

# pt1_REC8
pt1_REC8_Pval <- lapply(seq_along(pt1_REC8), function(x) {
  pt1_REC8[[x]]$numOverlaps$pval
})
pt1_REC8_Obs <- lapply(seq_along(pt1_REC8), function(x) {
  pt1_REC8[[x]]$numOverlaps$observed
})
pt1_REC8_Perm <- lapply(seq_along(pt1_REC8), function(x) {
  pt1_REC8[[x]]$numOverlaps$permuted
})
pt1_REC8_Exp <- lapply(seq_along(pt1_REC8), function(x) {
  mean(pt1_REC8[[x]]$numOverlaps$permuted)
})
pt1_REC8_log2ObsExp <- lapply(seq_along(pt1_REC8_Obs), function(x) {
  log2(pt1_REC8_Obs[[x]]/pt1_REC8_Exp[[x]])
})
pt1_REC8_Zscore <- lapply(seq_along(pt1_REC8), function(x) {
  pt1_REC8[[x]]$numOverlaps$zscore
})
pt1_REC8_AltHyp <- lapply(seq_along(pt1_REC8), function(x) {
  pt1_REC8[[x]]$numOverlaps$alternative
})
pt1_REC8_alpha0.05 <- lapply(seq_along(pt1_REC8_Perm), function(x) {
  if(pt1_REC8_AltHyp[[x]] == "greater") {
    quantile(pt1_REC8_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(pt1_REC8_Perm[[x]], 0.05)[[1]]
  }
})
pt1_REC8_log2alpha0.05 <- lapply(seq_along(pt1_REC8_alpha0.05), function(x) {
  log2(pt1_REC8_alpha0.05[[x]]/pt1_REC8_Exp[[x]])
})

pt1_REC8_log2ObsExp_sorted <- unlist(pt1_REC8_log2ObsExp[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix])
pt1_REC8_log2alpha0.05_sorted <- unlist(pt1_REC8_log2alpha0.05[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix])
pt1_REC8_otherNames_sorted <- otherNames[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix]
pt1_REC8_otherNamesPlot_sorted <- otherNamesPlot[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix]

# pt2_REC8
pt2_REC8_Pval <- lapply(seq_along(pt2_REC8), function(x) {
  pt2_REC8[[x]]$numOverlaps$pval
})
pt2_REC8_Obs <- lapply(seq_along(pt2_REC8), function(x) {
  pt2_REC8[[x]]$numOverlaps$observed
})
pt2_REC8_Perm <- lapply(seq_along(pt2_REC8), function(x) {
  pt2_REC8[[x]]$numOverlaps$permuted
})
pt2_REC8_Exp <- lapply(seq_along(pt2_REC8), function(x) {
  mean(pt2_REC8[[x]]$numOverlaps$permuted)
})
pt2_REC8_log2ObsExp <- lapply(seq_along(pt2_REC8_Obs), function(x) {
  log2(pt2_REC8_Obs[[x]]/pt2_REC8_Exp[[x]])
})
pt2_REC8_Zscore <- lapply(seq_along(pt2_REC8), function(x) {
  pt2_REC8[[x]]$numOverlaps$zscore
})
pt2_REC8_AltHyp <- lapply(seq_along(pt2_REC8), function(x) {
  pt2_REC8[[x]]$numOverlaps$alternative
})
pt2_REC8_alpha0.05 <- lapply(seq_along(pt2_REC8_Perm), function(x) {
  if(pt2_REC8_AltHyp[[x]] == "greater") {
    quantile(pt2_REC8_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(pt2_REC8_Perm[[x]], 0.05)[[1]]
  }
})
pt2_REC8_log2alpha0.05 <- lapply(seq_along(pt2_REC8_alpha0.05), function(x) {
  log2(pt2_REC8_alpha0.05[[x]]/pt2_REC8_Exp[[x]])
})

pt2_REC8_log2ObsExp_sorted <- unlist(pt2_REC8_log2ObsExp[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix])
pt2_REC8_log2alpha0.05_sorted <- unlist(pt2_REC8_log2alpha0.05[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix])
pt2_REC8_otherNames_sorted <- otherNames[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix]
pt2_REC8_otherNamesPlot_sorted <- otherNamesPlot[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix]

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

pt1_log2ObsExp_sorted <- unlist(pt1_log2ObsExp[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix])
pt1_log2alpha0.05_sorted <- unlist(pt1_log2alpha0.05[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix])
pt1_otherNames_sorted <- otherNames[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix]
pt1_otherNamesPlot_sorted <- otherNamesPlot[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix]

# pt2
pt2_Pval <- lapply(seq_along(pt2), function(x) {
  pt2[[x]]$numOverlaps$pval
})
pt2_Obs <- lapply(seq_along(pt2), function(x) {
  pt2[[x]]$numOverlaps$observed
})
pt2_Perm <- lapply(seq_along(pt2), function(x) {
  pt2[[x]]$numOverlaps$permuted
})
pt2_Exp <- lapply(seq_along(pt2), function(x) {
  mean(pt2[[x]]$numOverlaps$permuted)
})
pt2_log2ObsExp <- lapply(seq_along(pt2_Obs), function(x) {
  log2(pt2_Obs[[x]]/pt2_Exp[[x]])
})
pt2_Zscore <- lapply(seq_along(pt2), function(x) {
  pt2[[x]]$numOverlaps$zscore
})
pt2_AltHyp <- lapply(seq_along(pt2), function(x) {
  pt2[[x]]$numOverlaps$alternative
})
pt2_alpha0.05 <- lapply(seq_along(pt2_Perm), function(x) {
  if(pt2_AltHyp[[x]] == "greater") {
    quantile(pt2_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(pt2_Perm[[x]], 0.05)[[1]]
  }
})
pt2_log2alpha0.05 <- lapply(seq_along(pt2_alpha0.05), function(x) {
  log2(pt2_alpha0.05[[x]]/pt2_Exp[[x]])
})

pt2_log2ObsExp_sorted <- unlist(pt2_log2ObsExp[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix])
pt2_log2alpha0.05_sorted <- unlist(pt2_log2alpha0.05[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix])
pt2_otherNames_sorted <- otherNames[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix]
pt2_otherNamesPlot_sorted <- otherNamesPlot[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix]

df <- data.frame(Sample = rep(c(pt1REC8Name, pt2REC8Name, pt1Name, pt2Name),
                              each = length(pt1_log2ObsExp_sorted)),
                 Annotation_feature = rep(pt1_otherNamesPlot_sorted, 2),
                 log2ObsExp = c(pt1_REC8_log2ObsExp_sorted,
                                pt2_REC8_log2ObsExp_sorted,
                                pt1_log2ObsExp_sorted,
                                pt2_log2ObsExp_sorted),
                 log2alpha0.05 = c(pt1_REC8_log2alpha0.05_sorted,
                                   pt2_REC8_log2alpha0.05_sorted,
                                   pt1_log2alpha0.05_sorted,
                                   pt2_log2alpha0.05_sorted))

df$Annotation_feature <- factor(df$Annotation_feature,
                                levels = pt1_otherNamesPlot_sorted)
df$Sample <- factor(df$Sample,
                    levels = c(pt1REC8Name, pt2REC8Name,
                               pt1Name, pt2Name))

bp <- ggplot(data = df,
             mapping = aes(x = Annotation_feature,
                           y = log2ObsExp,
                           fill = Sample)) +
      geom_bar(stat = "identity",
               position = position_dodge()) +
      scale_fill_manual(name = "Genotype",
                        values = c("red",
                                   "red4",
                                   "blue",
                                   "navy"),
                        labels = c(pt1REC8Name,
                                   pt2REC8Name,
                                   pt1Name,
                                   pt2Name)) +
      geom_point(mapping = aes(Annotation_feature, log2alpha0.05),
                 position = position_dodge(0.9),
                 shape = "-", colour  = "grey70", size = 3.5) +
      labs(x = "Annotation feature",
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
      ggtitle(paste0(dataName, " (", as.character(perms), " permutations)"))
ggsave(paste0(plotDir, "barplot_other_features_permTestResults_",
              as.character(perms), "perms_",
              "log2_Observed_Expected_",
              pt1REC8LibName, "_", pt2REC8LibName, "_peaks_",
              pt1LibName, "_", pt2LibName, "_peaks.pdf"),
       plot = bp,
       height = 4.5, width = 8)
save(bp,
     file = paste0(plotDir, "barplot_other_features_permTestResults_",
                   as.character(perms), "perms_",
                   "log2_Observed_Expected_",
                   pt1REC8LibName, "_", pt2REC8LibName, "_peaks_",
                   pt1LibName, "_", pt2LibName, "_peaks.RData"))
