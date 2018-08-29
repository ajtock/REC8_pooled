# Plot permutation test log2(observed:expected) peak overlap barplots of equal width

library(ggplot2)
library(cowplot)
library(gridExtra)

perms <- 10000

REC8_MYC_Rep1_plotDir <- "/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_MYC_Rep1/genome_wide/regioneR/noMinWidth_mergedOverlaps/plots/"

load(paste0(REC8_MYC_Rep1_plotDir, "barplot_TE_families_permTestResults_",
            as.character(perms), "perms_",
            "log2_Observed_Expected_wt_REC8_MYC_Rep1_kss_REC8_HA_Rep1_peaks_wt_SPO11oligo_kss_SPO11oligo_peaks.RData"))
REC8_MYC_Rep1_SPO11oligo_bp <- bp
bp <- NULL

load(paste0(REC8_MYC_Rep1_plotDir, "barplot_TE_families_permTestResults_",
            as.character(perms), "perms_",
            "log2_Observed_Expected_wt_REC8_MYC_Rep1_kss_REC8_HA_Rep1_peaks_wt_H3K9me2_kss_H3K9me2_peaks.RData"))
REC8_MYC_Rep1_H3K9me2_bp <- bp
bp <- NULL

all_bp <- align_plots(REC8_MYC_Rep1_SPO11oligo_bp,
                      REC8_MYC_Rep1_H3K9me2_bp,
                      align = "hv",
                      axis = "tblr")
REC8_MYC_Rep1_SPO11oligo_bpx <- ggdraw(all_bp[[1]])
REC8_MYC_Rep1_H3K9me2_bpx <- ggdraw(all_bp[[2]])

save_plot(paste0(REC8_MYC_Rep1_plotDir,
                 "barplot_TE_families_permTestResults_",
                 as.character(perms), "perms_",
                 "log2_Observed_Expected_wt_REC8_MYC_Rep1_kss_REC8_HA_Rep1_peaks_wt_SPO11oligo_kss_SPO11oligo_peaks_equalPlotWidth.pdf"),
                 REC8_MYC_Rep1_SPO11oligo_bpx,
                 base_aspect_ratio = 1.3)
save_plot(paste0(REC8_MYC_Rep1_plotDir,
                 "barplot_TE_families_permTestResults_",
                 as.character(perms), "perms_",
                 "log2_Observed_Expected_wt_REC8_MYC_Rep1_kss_REC8_HA_Rep1_peaks_wt_H3K9me2_kss_H3K9me2_peaks_equalPlotWidth.pdf"),
                 REC8_MYC_Rep1_H3K9me2_bpx,
                 base_aspect_ratio = 1.3)

all_bp_singleplot <- list(REC8_MYC_Rep1_SPO11oligo_bpx,
                          REC8_MYC_Rep1_H3K9me2_bpx)
pdf(paste0(REC8_MYC_Rep1_plotDir,
           "barplot_TE_families_permTestResults_",
           as.character(perms), "perms_",
           "log2_Observed_Expected_wt_REC8_MYC_Rep1_kss_REC8_HA_Rep1_peaks_wt_SPO11oligo_kss_SPO11oligo_peaks_wt_H3K9me2_kss_H3K9me2_peaks_equalPlotWidth.pdf"),
    onefile = TRUE,
    height = 8, width = 7)
grid.arrange(all_bp_singleplot[[1]],
             all_bp_singleplot[[2]],
             nrow = 2, ncol = 1)
dev.off()
