#!/applications/R/R-3.3.2/bin/Rscript

# Convert txt file containing peaks to gff file

# Usage:
# ./txtTOgff.R REC8_HA_Rep1 REC8_HA_Rep1_wt_kss_diff

args <- commandArgs(trailingOnly = T)
libName <- args[1]
diffName <- args[2]

inDir <- paste0("/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/genome_wide/",
                diffName, "/")

peaks <- read.table(paste0(inDir,
                           libName, "_peaks_mean_log2_", diffName, "_cov_top10percent.txt"),
                    header = T)

peaks_sorted <- peaks[order(-peaks$diff_wt_kss_norm_peakCov),]
write.table(peaks_sorted,
            file = paste0(inDir,
                          libName, "_peaks_mean_log2_", diffName, "_cov_top10percent_sorted.txt"),
            sep = "\t", quote = F)

peaks_sorted_gff <- cbind(peaks_sorted$Chr,
                          rep("."),
                          rep(paste0(libName, "_peak")),
                          peaks_sorted$Start,
                          peaks_sorted$End,
                          peaks_sorted$diff_wt_kss_norm_peakCov,
                          rep("."),
                          rep("."),
                          rep("."))
colnames(peaks_sorted_gff) <- c("chr",
                                "source",
                                "feature",
                                "start",
                                "end",
                                "score",
                                "strand",
                                "frame",
                                "attribute")
write.table(peaks_sorted_gff,
            file = paste0(inDir,
                          libName, "_peaks_mean_log2_", diffName, "_cov_top10percent_sorted.gff"),
            sep = "\t", quote = F, row.names = F, col.names = F)

