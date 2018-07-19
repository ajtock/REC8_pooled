
# Create Venn diagrams of the number of REC8-HA Rep1 peaks in each differential set

# Usage:
# Rscript peak_set_overlaps.R red blue green darkmagenta 

library(VennDiagram)
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

args <- commandArgs(trailingOnly = TRUE)
myColours <- c(args[1], args[2], args[3], args[4])

plotDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/genome_wide/peak_set_overlaps_and_width_dists/plots/"

REC8 <- rownames(read.table("/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/genome_wide/REC8_HA_Rep1_wt_kss_diff/REC8_HA_Rep1_peaks_mean_log2_REC8_HA_Rep1_wt_kss_diff_cov_top10percent.txt"))
RNAseq <- rownames(read.table("/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/genome_wide/RNAseq_Rep1_kss_wt_diff/REC8_HA_Rep1_peaks_mean_RNAseq_Rep1_kss_wt_diff_cov_top10percent.txt"))
SPO11oligo <- rownames(read.table("/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/genome_wide/SPO11oligo_Rep1_kss_wt_diff/REC8_HA_Rep1_peaks_mean_log2_SPO11oligo_Rep1_kss_wt_diff_cov_top10percent.txt"))
H3K9me2 <- rownames(read.table("/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/genome_wide/H3K9me2_wt_kss_diff/REC8_HA_Rep1_peaks_mean_log2_H3K9me2_wt_kss_diff_cov_top10percent.txt"))

peakRowNameList <- list(REC8,
                        RNAseq,
                        SPO11oligo,
                        H3K9me2)

names <- c("REC8-HA Rep1 wt-kss",
           "RNA-seq kss-wt",
           "SPO11-1-oligos kss-wt",
           "H3K9me2 wt-kss")

vennPlotFun <- function(peakRowNameList, myColours, mainTitle, catPos) {
  grid.draw(venn.diagram(x = peakRowNameList,
                         filename = NULL,
                         fill = c(myColours),
                         alpha = c(0.5),
                         cex = 0.75,
                         main.cex = 2,
                         cat.cex = 0.75,
                         fontfamily = "sans",
                         main.fontfamily = "sans",
                         cat.fontfamily = "sans",
                         main.fontface = "bold",
                         cat.fontface = "bold",
                         main = mainTitle,
                         category.names = names,
                         cat.pos = catPos))
}

pdf(paste0(plotDir, "VennDiagram_REC8_HA_Rep1_peaks_sets_differential_wt_kss.pdf"))
vennPlotFun(peakRowNameList = peakRowNameList,
            myColours = myColours,
            mainTitle = "REC8-HA Rep1 peaks\nwith greatest change in kss",
            catPos = c(-10, 10, 10, 10))
dev.off()



