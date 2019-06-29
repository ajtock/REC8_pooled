# Plot permutation test log2(observed:expected) peak overlap barplots of equal width

# Usage:
# Rscript TE_family_wt_vs_kss_REC8_HA_Rep1_vs_SPO11oligo_vs_H3K9me2_genome_wide_peaks_equalPlotWidths.R wt_REC8_HA_Rep1 kss_REC8_HA_Rep1 wt_SPO11oligo kss_SPO11oligo wt_H3K9me2 kss_H3K9me2 10000

library(ggplot2)
library(cowplot)
library(gridExtra)

args <- commandArgs(trailingOnly = T)
#pt1REC8LibName <- "wt_REC8_HA_Rep1"
#pt2REC8LibName <- "kss_REC8_HA_Rep1"
#pt1LibName <- "wt_SPO11oligo"
#pt2LibName <- "kss_SPO11oligo"
## Number of permutations (randomisations) performed
#perms <- 10000

pt1REC8LibName <- as.character(args[1])
pt2REC8LibName <- as.character(args[2])
pt1LibName <- as.character(args[3])
pt2LibName <- as.character(args[4])
pt3LibName <- as.character(args[5])
pt4LibName <- as.character(args[6])
# Number of permutations (randomisations) performed
perms <- as.numeric(args[7])

plotDir <- "/home/ajt200/analysis/180622_Chris_lambing_ChIP_REC8_HA_Col_kss/kss/peaks/PeakRanger1.18/ranger/rep_specific_input_p0.001_q0.01/kss_REC8_HA_Rep1/genome_wide/regioneR/noMinWidth_mergedOverlaps/plots/"

load(paste0(plotDir, "barplot_TE_families_permTestResults_",
            as.character(perms), "perms_",
            "log2_Observed_Expected_",
            pt1REC8LibName, "_",
            pt2REC8LibName, "_peaks_",
            pt1LibName, "_",
            pt2LibName, "_peaks.RData"))
REC8_SPO11oligo_bp <- bp
bp <- NULL

load(paste0(plotDir, "barplot_TE_families_permTestResults_",
            as.character(perms), "perms_",
            "log2_Observed_Expected_",
            pt1REC8LibName, "_",
            pt2REC8LibName, "_peaks_",
            pt3LibName, "_",
            pt4LibName, "_peaks.RData"))
REC8_H3K9me2_bp <- bp
bp <- NULL

all_bp <- align_plots(REC8_SPO11oligo_bp,
                      REC8_H3K9me2_bp,
                      align = "hv",
                      axis = "tblr")
REC8_SPO11oligo_bpx <- ggdraw(all_bp[[1]])
REC8_H3K9me2_bpx <- ggdraw(all_bp[[2]])

save_plot(paste0(plotDir,
                 "barplot_TE_families_permTestResults_",
                 as.character(perms), "perms_",
                 "log2_Observed_Expected_",
                 pt1REC8LibName, "_",
                 pt2REC8LibName, "_peaks_",
                 pt1LibName, "_",
                 pt2LibName, "_peaks_equalPlotWidth.pdf"),
                 REC8_SPO11oligo_bpx,
                 base_aspect_ratio = 1.3)
save_plot(paste0(plotDir,
                 "barplot_TE_families_permTestResults_",
                 as.character(perms), "perms_",
                 "log2_Observed_Expected_",
                 pt1REC8LibName, "_",
                 pt2REC8LibName, "_peaks_",
                 pt3LibName, "_",
                 pt4LibName, "_peaks_equalPlotWidth.pdf"),
                 REC8_H3K9me2_bpx,
                 base_aspect_ratio = 1.3)

all_bp_singleplot <- list(REC8_SPO11oligo_bpx,
                          REC8_H3K9me2_bpx)
pdf(paste0(plotDir,
           "barplot_TE_families_permTestResults_",
           as.character(perms), "perms_",
           "log2_Observed_Expected_",
           pt1REC8LibName, "_",
           pt2REC8LibName, "_peaks_",
           pt1LibName, "_",
           pt2LibName, "_peaks_",
           pt3LibName, "_",
           pt4LibName, "_peaks_equalPlotWidth.pdf"),
    onefile = TRUE,
    height = 8, width = 7)
grid.arrange(all_bp_singleplot[[1]],
             all_bp_singleplot[[2]],
             nrow = 2, ncol = 1)
dev.off()
