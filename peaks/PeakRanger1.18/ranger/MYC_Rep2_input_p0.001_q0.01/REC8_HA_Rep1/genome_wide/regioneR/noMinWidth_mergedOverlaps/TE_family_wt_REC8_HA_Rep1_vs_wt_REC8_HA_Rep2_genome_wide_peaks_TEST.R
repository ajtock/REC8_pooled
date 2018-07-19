#!/applications/R/R-3.3.2/bin/Rscript

# Plot bar chart of log2(observed:expected) peaks overlapping TEs within each family

# Usage:
# Rscript TE_family_wt_REC8_HA_Rep1_vs_wt_REC8_HA_Rep2_genome_wide_peaks_TEST.R "REC8 peaks" "wt REC8-HA Rep1" "wt REC8-HA Rep2" wt_REC8_HA_Rep1 wt_REC8_HA_Rep2 10000

library(ggplot2)
library(ggthemes)

dataName <- "REC8 peaks" 
pt1Name <- "wt REC8-HA Rep1" 
pt2Name <-  "wt REC8-HA Rep2" 
pt1LibName <- "wt_REC8_HA_Rep1"
pt2LibName <- "wt_REC8_HA_Rep2"
# Number of permutations (randomisations) performed
perms <- 10000

args <- commandArgs(trailingOnly = T)
dataName <- as.character(args[1])
pt1Name <- as.character(args[2])
pt2Name <- as.character(args[3])
pt1LibName <- as.character(args[4])
pt2LibName <- as.character(args[5])
# Number of permutations (randomisations) performed
perms <- as.numeric(args[6])


inDir1 <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/genome_wide/regioneR/noMinWidth_mergedOverlaps/"
inDir2 <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep2/genome_wide/regioneR/noMinWidth_mergedOverlaps/"
plotDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/genome_wide/regioneR/noMinWidth_mergedOverlaps/plots/"

famNames <- c("dna", "heli", "ptmari", "mudr", "enspm", "hat", "harbinger",
              "rna", "gypsy", "copia", "linel1", "sine")
famNamesPlot <- c("DNA", "Helitrons", "Pogo/Tc1/Mariner", "MuDR", "EnSpm", "hAT", "Harbinger",
                  "RNA", "Gypsy LTR", "Copia LTR", "LINE-1", "SINE")

load(paste0(inDir1, "permTest_REC8_HA_Rep1_rangerPeaks_vs_TEsDNA.RData"))
pt1_DNA <- ptPeaksTEsDNAPerChrom
ptPeaksTEsDNAPerChrom <- NULL
load(paste0(inDir1, "permTest_REC8_HA_Rep1_rangerPeaks_vs_TEsRNA.RData"))
pt1_RNA <- ptPeaksTEsRNAPerChrom
ptPeaksTEsRNAPerChrom <- NULL
pt1 <- c(pt1_DNA, pt1_RNA)

load(paste0(inDir2, "permTest_REC8_HA_Rep2_rangerPeaks_vs_TEsDNA.RData"))
pt2_DNA <- ptPeaksTEsDNAPerChrom
ptPeaksTEsDNAPerChrom <- NULL
load(paste0(inDir2, "permTest_REC8_HA_Rep2_rangerPeaks_vs_TEsRNA.RData"))
pt2_RNA <- ptPeaksTEsRNAPerChrom
ptPeaksTEsRNAPerChrom <- NULL
pt2 <- c(pt2_DNA, pt2_RNA)

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

pt1_log2ObsExp_sorted <- sort.int(unlist(pt1_log2ObsExp), decreasing = T)
pt1_log2alpha0.05_sorted <- unlist(pt1_log2alpha0.05[sort.int(unlist(pt1_log2ObsExp), decreasing = T, index.return = T)$ix])
pt1_famNames_sorted <- famNames[sort.int(unlist(pt1_log2ObsExp), decreasing = T, index.return = T)$ix]
pt1_famNamesPlot_sorted <- famNamesPlot[sort.int(unlist(pt1_log2ObsExp), decreasing = T, index.return = T)$ix]


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

pt2_log2ObsExp_sorted <- unlist(pt2_log2ObsExp[sort.int(unlist(pt1_log2ObsExp), decreasing = T, index.return = T)$ix])
pt2_log2alpha0.05_sorted <- unlist(pt2_log2alpha0.05[sort.int(unlist(pt1_log2ObsExp), decreasing = T, index.return = T)$ix])
pt2_famNames_sorted <- famNames[sort.int(unlist(pt1_log2ObsExp), decreasing = T, index.return = T)$ix]
pt2_famNamesPlot_sorted <- famNamesPlot[sort.int(unlist(pt1_log2ObsExp), decreasing = T, index.return = T)$ix]

df <- data.frame(Sample = rep(c(pt1Name, pt2Name),
                              each = length(pt1_log2ObsExp_sorted)),
                 Transposon_family = rep(pt1_famNamesPlot_sorted, 2),
                 log2ObsExp = c(pt1_log2ObsExp_sorted,
                                pt2_log2ObsExp_sorted),
                 log2alpha0.05 = c(pt1_log2alpha0.05_sorted, pt2_log2alpha0.05_sorted))

df$Transposon_family <- factor(df$Transposon_family,
                               levels = pt1_famNamesPlot_sorted)
df$Sample <- factor(df$Sample,
                    levels = c(pt1Name, pt2Name))

bp <- ggplot(data = df,
             mapping = aes(x = Transposon_family,
                           y = log2ObsExp,
                           fill = Sample)) +
      geom_bar(stat = "identity",
               position = position_dodge()) +
      scale_fill_manual(name = "Genotype",
                        values = c("black",
                                   "dodgerblue3"),
                        labels = c(pt1Name,
                                   pt2Name)) +
      geom_point(mapping = aes(Transposon_family, log2alpha0.05),
                 position = position_dodge(0.9),
                 shape = "-", colour  = "red", size = 4.25) +
      labs(x = "Transposon family",
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
ggsave(paste0(plotDir, "barplot_TE_families_permTestResults_",
              as.character(perms), "perms_",
              "log2_Observed_Expected_", pt1LibName, "_", pt2LibName, "_peaks.pdf"),
       plot = bp,
       height = 4, width = 5)

 
