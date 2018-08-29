# Plot permutation test log2(observed:expected) peak overlap barplots of equal width

library(ggplot2)
library(cowplot)
library(gridExtra)

perms <- 10000

REC8_HA_Rep1_plotDir <- "/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/genome_wide/regioneR/noMinWidth_mergedOverlaps/plots/"
load(paste0(REC8_HA_Rep1_plotDir, "barplot_TE_families_permTestResults_",
            as.character(perms), "perms_",
            "log2_Observed_Expected_wt_REC8_HA_Rep1_kss_REC8_HA_Rep1_peaks.RData"))
REC8_HA_Rep1_bp <- bp 
bp <- NULL

REC8_HA_Rep2_plotDir <- "/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep2/genome_wide/regioneR/noMinWidth_mergedOverlaps/plots/"
load(paste0(REC8_HA_Rep2_plotDir, "barplot_TE_families_permTestResults_",
            as.character(perms), "perms_",
            "log2_Observed_Expected_wt_REC8_HA_Rep2_kss_REC8_HA_Rep1_peaks.RData"))
REC8_HA_Rep2_bp <- bp 
bp <- NULL

REC8_MYC_Rep1_plotDir <- "/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_MYC_Rep1/genome_wide/regioneR/noMinWidth_mergedOverlaps/plots/"
load(paste0(REC8_MYC_Rep1_plotDir, "barplot_TE_families_permTestResults_",
            as.character(perms), "perms_",
            "log2_Observed_Expected_wt_REC8_MYC_Rep1_kss_REC8_HA_Rep1_peaks.RData"))
REC8_MYC_Rep1_bp <- bp 
bp <- NULL

SPO11oligos_plotDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/RPI1_RPI8_profiles_rangerPeaks_idr0.05/genome_wide/regioneR/plots/"
load(paste0(SPO11oligos_plotDir, "barplot_TE_families_permTestResults_",
            as.character(perms), "perms_",
            "log2_Observed_Expected_wt_SPO11oligo_kss_SPO11oligo_peaks.RData"))
SPO11oligo_bp <- bp
bp <- NULL

H3K9me2_plotDir <- "/projects/ajt200/BAM_masters/H3K9me2/WT/peaks/PeakRanger1.18/ranger/p0.05_q0.05/genome_wide/regioneR/noMinWidth_mergedOverlaps/plots/"
load(paste0(H3K9me2_plotDir, "barplot_TE_families_permTestResults_",
            as.character(perms), "perms_",
            "log2_Observed_Expected_wt_H3K9me2_kss_H3K9me2_peaks.RData"))
H3K9me2_bp <- bp
bp <- NULL

all_bp <- align_plots(REC8_HA_Rep1_bp,
                      REC8_HA_Rep2_bp,
                      REC8_MYC_Rep1_bp,
                      SPO11oligo_bp,
                      H3K9me2_bp,
                      align = "hv",
                      axis = "tblr")
REC8_HA_Rep1_bpx <- ggdraw(all_bp[[1]])
REC8_HA_Rep2_bpx <- ggdraw(all_bp[[2]])
REC8_MYC_Rep1_bpx <- ggdraw(all_bp[[3]])
SPO11oligo_bpx <- ggdraw(all_bp[[4]])
H3K9me2_bpx <- ggdraw(all_bp[[5]])

save_plot(paste0(REC8_HA_Rep1_plotDir,
                 "barplot_TE_families_permTestResults_",
                 as.character(perms), "perms_",
                 "log2_Observed_Expected_wt_REC8_HA_Rep1_kss_REC8_HA_Rep1_peaks_equalPlotWidth.pdf"),
                 REC8_HA_Rep1_bpx,
                 base_aspect_ratio = 1.3)
save_plot(paste0(REC8_HA_Rep2_plotDir,
                 "barplot_TE_families_permTestResults_",
                 as.character(perms), "perms_",
                 "log2_Observed_Expected_wt_REC8_HA_Rep2_kss_REC8_HA_Rep1_peaks_equalPlotWidth.pdf"),
                 REC8_HA_Rep2_bpx,
                 base_aspect_ratio = 1.3)
save_plot(paste0(REC8_MYC_Rep1_plotDir,
                 "barplot_TE_families_permTestResults_",
                 as.character(perms), "perms_",
                 "log2_Observed_Expected_wt_REC8_MYC_Rep1_kss_REC8_HA_Rep1_peaks_equalPlotWidth.pdf"),
                 REC8_MYC_Rep1_bpx,
                 base_aspect_ratio = 1.3)
save_plot(paste0(SPO11oligos_plotDir,
                 "barplot_TE_families_permTestResults_",
                 as.character(perms), "perms_",
                 "log2_Observed_Expected_wt_SPO11oligo_kss_SPO11oligo_peaks_equalPlotWidth.pdf"),
                 SPO11oligo_bpx,
                 base_aspect_ratio = 1.3)
save_plot(paste0(H3K9me2_plotDir,
                 "barplot_TE_families_permTestResults_",
                 as.character(perms), "perms_",
                 "log2_Observed_Expected_wt_H3K9me2_kss_H3K9me2_peaks_equalPlotWidth.pdf"),
                 H3K9me2_bpx,
                 base_aspect_ratio = 1.3)

all_bp_singleplot <- list(REC8_HA_Rep1_bpx,
                          REC8_HA_Rep2_bpx,
                          REC8_MYC_Rep1_bpx,
                          SPO11oligo_bpx,
                          H3K9me2_bpx)
pdf(paste0(REC8_HA_Rep1_plotDir,
           "barplot_TE_families_permTestResults_",
           as.character(perms), "perms_",
           "log2_Observed_Expected_wt_kss_peaks_equalPlotWidth.pdf"),
    onefile = TRUE,
    height = 20, width = 6)
grid.arrange(all_bp_singleplot[[1]],
             all_bp_singleplot[[2]],
             all_bp_singleplot[[3]],
             all_bp_singleplot[[4]],
             all_bp_singleplot[[5]],
             nrow = 5, ncol = 1)
dev.off()

## Doesn't work:
#all_bp_singleplot <- plot_grid(REC8_HA_Rep1_bp,
#                               REC8_HA_Rep2_bp,
#                               REC8_MYC_Rep1_bp,
#                               SPO11oligo_bp,
#                               ncol = 1,
#                               align = "v",
#                               rel_heights = c(10, 8, 6, 2))
#save_plot(paste0(REC8_HA_Rep1_plotDir,
#                 "barplot_TE_families_permTestResults_",
#                 as.character(perms), "perms_",
#                 "log2_Observed_Expected_wt_kss_peaks_equalPlotWidth.pdf"),
#                 all_bp_singleplot,
#                 base_aspect_ratio = 1.3)





