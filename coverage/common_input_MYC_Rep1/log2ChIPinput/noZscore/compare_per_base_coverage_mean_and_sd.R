#!/applications/R/R-3.5.0/bin/Rscript

# Without Z-score standardisation
wt_REC8_MYC_Rep1 <- read.table("/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/noZscore/log2_REC8_MYC_Rep1_ChIP_REC8_MYC_Rep1_input_norm_allchrs_coverage_coord_tab_noZscore.bed",
                              colClasses = c(NA, NA, "NULL", NA))
wt_REC8_HA_Rep1 <- read.table("/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/noZscore/log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep1_input_norm_allchrs_coverage_coord_tab_noZscore.bed",
                              colClasses = c(NA, NA, "NULL", NA))
wt_REC8_HA_Rep2 <- read.table("/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/noZscore/log2_REC8_HA_Rep2_ChIP_REC8_MYC_Rep1_input_norm_allchrs_coverage_coord_tab_noZscore.bed",
                              colClasses = c(NA, NA, "NULL", NA))
kss_REC8_HA_Rep1 <- read.table("/home/ajt200/analysis/180622_Chris_lambing_ChIP_REC8_HA_Col_kss/kss/coverage/log2ChIPinput/replicate_specific_input/noZscore/log2_kss_REC8_HA_Rep1_ChIP_kss_REC8_HA_Rep1_input_norm_allchrs_coverage_coord_tab_noZscore.bed",
                               colClasses = c(NA, NA, "NULL", NA))
kss_REC8_HA_Rep2 <- read.table("/home/ajt200/analysis/180830_Chris_lambing_ChIP_rec8HA_kss_rep2/coverage/log2ChIPinput/kss_REC8_HA_Rep1_input/noZscore/log2_kss_REC8_HA_Rep2_ChIP_kss_REC8_HA_Rep1_input_norm_allchrs_coverage_coord_tab_noZscore.bed",
                               colClasses = c(NA, NA, "NULL", NA))

# With Z-score standardisation
wt_REC8_MYC_Rep1_Z <- read.table("/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/log2_REC8_MYC_Rep1_ChIP_REC8_MYC_Rep1_input_norm_allchrs_coverage_coord_tab.bed",
                              colClasses = c(NA, NA, "NULL", NA))
wt_REC8_HA_Rep1_Z <- read.table("/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep1_input_norm_allchrs_coverage_coord_tab.bed",
                              colClasses = c(NA, NA, "NULL", NA))
wt_REC8_HA_Rep2_Z <- read.table("/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/log2_REC8_HA_Rep2_ChIP_REC8_MYC_Rep1_input_norm_allchrs_coverage_coord_tab.bed",
                              colClasses = c(NA, NA, "NULL", NA))
kss_REC8_HA_Rep1_Z <- read.table("/home/ajt200/analysis/180622_Chris_lambing_ChIP_REC8_HA_Col_kss/kss/coverage/log2ChIPinput/replicate_specific_input/log2_kss_REC8_HA_Rep1_ChIP_kss_REC8_HA_Rep1_input_norm_allchrs_coverage_coord_tab.bed",
                               colClasses = c(NA, NA, "NULL", NA))
kss_REC8_HA_Rep2_Z <- read.table("/home/ajt200/analysis/180830_Chris_lambing_ChIP_rec8HA_kss_rep2/coverage/log2ChIPinput/kss_REC8_HA_Rep1_input/log2_kss_REC8_HA_Rep2_ChIP_kss_REC8_HA_Rep1_input_norm_allchrs_coverage_coord_tab.bed",
                               colClasses = c(NA, NA, "NULL", NA))

# Without Z-score standardisation
wt_SPO11oligos_Rep1 <- read.table("/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/noZscore/log2wtSPO11oligoRPI1NakedDNA_norm_allchrs_coverage_coord_tab_noZscore.bed",
                                  colClasses = c(NA, NA, "NULL", NA))
wt_SPO11oligos_Rep2 <- read.table("/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/noZscore/log2wtSPO11oligoRPI3NakedDNA_norm_allchrs_coverage_coord_tab_noZscore.bed",
                                  colClasses = c(NA, NA, "NULL", NA))
wt_SPO11oligos_Rep3 <- read.table("/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/noZscore/log2wtSPO11oligoRPI8NakedDNA_norm_allchrs_coverage_coord_tab_noZscore.bed",
                                  colClasses = c(NA, NA, "NULL", NA))
kss_SPO11oligos_Rep1 <- read.table("/projects/ajt200/BAM_masters/SPO11-oligo/suvh456/coverage/log2ChIPinput/noZscore/log2suvh456SPO11oligoRPI34NakedDNA_norm_allchrs_coverage_coord_tab_noZscore.bed",
                                   colClasses = c(NA, NA, "NULL", NA))
kss_SPO11oligos_Rep2 <- read.table("/projects/ajt200/BAM_masters/SPO11-oligo/suvh456/coverage/log2ChIPinput/noZscore/log2suvh456SPO11oligoRPI35NakedDNA_norm_allchrs_coverage_coord_tab_noZscore.bed",
                                   colClasses = c(NA, NA, "NULL", NA))

# With Z-score standardisation
wt_SPO11oligos_Rep1_Z <- read.table("/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/log2wtSPO11oligoRPI1NakedDNA_norm_allchrs_coverage_coord_tab.bed",
                                    colClasses = c(NA, NA, "NULL", NA))
wt_SPO11oligos_Rep2_Z <- read.table("/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/log2wtSPO11oligoRPI3NakedDNA_norm_allchrs_coverage_coord_tab.bed",
                                    colClasses = c(NA, NA, "NULL", NA))
wt_SPO11oligos_Rep3_Z <- read.table("/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/log2wtSPO11oligoRPI8NakedDNA_norm_allchrs_coverage_coord_tab.bed",
                                    colClasses = c(NA, NA, "NULL", NA))
kss_SPO11oligos_Rep1_Z <- read.table("/projects/ajt200/BAM_masters/SPO11-oligo/suvh456/coverage/log2ChIPinput/log2suvh456SPO11oligoRPI34NakedDNA_norm_allchrs_coverage_coord_tab.bed",
                                     colClasses = c(NA, NA, "NULL", NA))
kss_SPO11oligos_Rep2_Z <- read.table("/projects/ajt200/BAM_masters/SPO11-oligo/suvh456/coverage/log2ChIPinput/log2suvh456SPO11oligoRPI35NakedDNA_norm_allchrs_coverage_coord_tab.bed",
                                     colClasses = c(NA, NA, "NULL", NA))

# Without Z-score standardisation
wt_H3K9me2_Rep1 <- read.table("/projects/ajt200/BAM_masters/H3K9me2/WT/coverage/log2ChIPinput/noZscore/log2_WT_H3K9me2_ChIP_WT_H3K9me2_input_norm_allchrs_coverage_coord_tab_noZscore.bed",
                              colClasses = c(NA, NA, "NULL", NA))
kss_H3K9me2_Rep1 <- read.table("/projects/ajt200/BAM_masters/H3K9me2/kss/coverage/log2ChIPinput/noZscore/log2_kss_H3K9me2_ChIP_kss_H3K9me2_input_norm_allchrs_coverage_coord_tab_noZscore.bed",
                               colClasses = c(NA, NA, "NULL", NA))

# With Z-score standardisation
wt_H3K9me2_Rep1_Z <- read.table("/projects/ajt200/BAM_masters/H3K9me2/WT/coverage/log2ChIPinput/log2_WT_H3K9me2_ChIP_WT_H3K9me2_input_norm_allchrs_coverage_coord_tab.bed",
                                colClasses = c(NA, NA, "NULL", NA))
kss_H3K9me2_Rep1_Z <- read.table("/projects/ajt200/BAM_masters/H3K9me2/kss/coverage/log2ChIPinput/log2_kss_H3K9me2_ChIP_kss_H3K9me2_input_norm_allchrs_coverage_coord_tab.bed",
                                 colClasses = c(NA, NA, "NULL", NA))

# Histogram plotting function
plotHist <- function(dat,
                     name,
                     breaks,
                     Xlim,
                     Ylim) {
  hist(x = dat,
       breaks = breaks,
       col = "grey50",
       border = NA,
       lwd = 2,
       xlab = expression("Z-score log"[2]*"(ChIP/input)"),
       ylab = "Loci",
       xlim = Xlim,
       ylim = Ylim,
       main = name,
       cex.lab = 2, cex.axis = 2, cex.main = 2)
abline(v = mean(dat),
       col = "red", lty = 2, lwd = 2)
abline(v = c(mean(dat)-sd(dat),
             mean(dat)+sd(dat)),
       col = "blue", lty = 2, lwd = 2)
}

# Without Z-score standardisation
pdf("./wt_vs_kss_REC8_SPO11oligos_H3K9me2_coverage_histogram_noZscore.pdf",
    height = 20, width = 12)
par(mfrow = c(5, 3), mar = c(6, 6, 2, 2), mgp = c(4, 1.5, 0))
plotHist(dat = wt_REC8_MYC_Rep1[,3],
         name = "wt REC8-Myc Rep1",
         breaks = 250,
         Xlim = c(-6, 6),
#         Xlim = c(min(wt_REC8_MYC_Rep1[,3],
#                      wt_REC8_HA_Rep1[,3],
#                      wt_REC8_HA_Rep2[,3],
#                      kss_REC8_HA_Rep1[,3],
#                      kss_REC8_HA_Rep2[,3]),
#                  max(wt_REC8_MYC_Rep1[,3],
#                      wt_REC8_HA_Rep1[,3],
#                      wt_REC8_HA_Rep2[,3],
#                      kss_REC8_HA_Rep1[,3],
#                      kss_REC8_HA_Rep2[,3])),
         Ylim = c(0, 4e+06))
plotHist(dat = wt_REC8_HA_Rep1[,3],
         name = "wt REC8-HA Rep1",
         breaks = 250,
         Xlim = c(-6, 6),
         Ylim = c(0, 4e+06))
plotHist(dat = wt_REC8_HA_Rep2[,3],
         name = "wt REC8-HA Rep2",
         breaks = 250,
         Xlim = c(-6, 6),
         Ylim = c(0, 4e+06))
plotHist(dat = kss_REC8_HA_Rep1[,3],
         name = "kss REC8-HA Rep1",
         breaks = 250,
         Xlim = c(-6, 6),
         Ylim = c(0, 4e+06))
plotHist(dat = kss_REC8_HA_Rep2[,3],
         name = "kss REC8-HA Rep2",
         breaks = 250,
         Xlim = c(-6, 6),
         Ylim = c(0, 4e+06))
plot.new()
plotHist(dat = wt_SPO11oligos_Rep1[,3],
         name = "wt SPO11-1-oligos Rep1",
         breaks = 100,
         Xlim = c(-6, 6),
         Ylim = c(0, 6e+06))
plotHist(dat = wt_SPO11oligos_Rep2[,3],
         name = "wt SPO11-1-oligos Rep2",
         breaks = 100,
         Xlim = c(-6, 6),
         Ylim = c(0, 6e+06))
plotHist(dat = wt_SPO11oligos_Rep3[,3],
         name = "wt SPO11-1-oligos Rep3",
         breaks = 100,
         Xlim = c(-6, 6),
         Ylim = c(0, 6e+06))
plotHist(dat = kss_SPO11oligos_Rep1[,3],
         name = "kss SPO11-1-oligos Rep1",
         breaks = 100,
         Xlim = c(-6, 6),
         Ylim = c(0, 6e+06))
plotHist(dat = kss_SPO11oligos_Rep2[,3],
         name = "kss SPO11-1-oligos Rep2",
         breaks = 100,
         Xlim = c(-6, 6),
         Ylim = c(0, 6e+06))
plot.new()
plotHist(dat = wt_H3K9me2_Rep1[,3],
         name = "wt H3K9me2 Rep1",
         breaks = 250,
         Xlim = c(-6, 6),
         Ylim = c(0, 4e+06))
plotHist(dat = kss_H3K9me2_Rep1[,3],
         name = "kss H3K9me2 Rep1",
         breaks = 250,
         Xlim = c(-6, 6),
         Ylim = c(0, 4e+06))
dev.off()

# With Z-scrore standardisation
pdf("./wt_vs_kss_REC8_SPO11oligos_H3K9me2_coverage_histogram_Zscore.pdf",
    height = 20, width = 12)
par(mfrow = c(5, 3), mar = c(6, 6, 2, 2), mgp = c(4, 1.5, 0))
plotHist(dat = wt_REC8_MYC_Rep1_Z[,3],
         name = "wt REC8-Myc Rep1",
         breaks = 250,
         Xlim = c(-6, 6),
#         Xlim = c(min(wt_REC8_MYC_Rep1_Z[,3],
#                      wt_REC8_HA_Rep1_Z[,3],
#                      wt_REC8_HA_Rep2_Z[,3],
#                      kss_REC8_HA_Rep1_Z[,3],
#                      kss_REC8_HA_Rep2_Z[,3]),
#                  max(wt_REC8_MYC_Rep1_Z[,3],
#                      wt_REC8_HA_Rep1_Z[,3],
#                      wt_REC8_HA_Rep2_Z[,3],
#                      kss_REC8_HA_Rep1_Z[,3],
#                      kss_REC8_HA_Rep2_Z[,3])),
         Ylim = c(0, 4e+06))
plotHist(dat = wt_REC8_HA_Rep1_Z[,3],
         name = "wt REC8-HA Rep1",
         breaks = 250,
         Xlim = c(-6, 6),
         Ylim = c(0, 4e+06))
plotHist(dat = wt_REC8_HA_Rep2_Z[,3],
         name = "wt REC8-HA Rep2",
         breaks = 250,
         Xlim = c(-6, 6),
         Ylim = c(0, 4e+06))
plotHist(dat = kss_REC8_HA_Rep1_Z[,3],
         name = "kss REC8-HA Rep1",
         breaks = 250,
         Xlim = c(-6, 6),
         Ylim = c(0, 4e+06))
plotHist(dat = kss_REC8_HA_Rep2_Z[,3],
         name = "kss REC8-HA Rep2",
         breaks = 250,
         Xlim = c(-6, 6),
         Ylim = c(0, 4e+06))
plot.new()
plotHist(dat = wt_SPO11oligos_Rep1_Z[,3],
         name = "wt SPO11-1-oligos Rep1",
         breaks = 100,
         Xlim = c(-6, 6),
         Ylim = c(0, 6e+06))
plotHist(dat = wt_SPO11oligos_Rep2_Z[,3],
         name = "wt SPO11-1-oligos Rep2",
         breaks = 100,
         Xlim = c(-6, 6),
         Ylim = c(0, 6e+06))
plotHist(dat = wt_SPO11oligos_Rep3_Z[,3],
         name = "wt SPO11-1-oligos Rep3",
         breaks = 100,
         Xlim = c(-6, 6),
         Ylim = c(0, 6e+06))
plotHist(dat = kss_SPO11oligos_Rep1_Z[,3],
         name = "kss SPO11-1-oligos Rep1",
         breaks = 100,
         Xlim = c(-6, 6),
         Ylim = c(0, 6e+06))
plotHist(dat = kss_SPO11oligos_Rep2_Z[,3],
         name = "kss SPO11-1-oligos Rep2",
         breaks = 100,
         Xlim = c(-6, 6),
         Ylim = c(0, 6e+06))
plot.new()
plotHist(dat = wt_H3K9me2_Rep1_Z[,3],
         name = "wt H3K9me2 Rep1",
         breaks = 250,
         Xlim = c(-6, 6),
         Ylim = c(0, 4e+06))
plotHist(dat = kss_H3K9me2_Rep1_Z[,3],
         name = "kss H3K9me2 Rep1",
         breaks = 250,
         Xlim = c(-6, 6),
         Ylim = c(0, 4e+06))
dev.off()


# Create data.frame of unstandardised means and SDs
df <- data.frame(Sample = c("wt REC8-Myc Rep1",
                            "wt REC8-HA Rep1",
                            "wt REC8-HA Rep2",
                            "kss REC8-HA Rep1",
                            "kss REC8-HA Rep2",
                            "wt SPO11-1-oligos Rep1",
                            "wt SPO11-1-oligos Rep2",
                            "wt SPO11-1-oligos Rep3",
                            "kss SPO11-1-oligos Rep1",
                            "kss SPO11-1-oligos Rep2",
                            "wt H3K9me2 Rep1",
                            "kss H3K9me2 Rep1"),
                 Mean = c(mean(wt_REC8_MYC_Rep1[,3]),
                          mean(wt_REC8_HA_Rep1[,3]),
                          mean(wt_REC8_HA_Rep2[,3]),
                          mean(kss_REC8_HA_Rep1[,3]),
                          mean(kss_REC8_HA_Rep2[,3]),
                          mean(wt_SPO11oligos_Rep1[,3]),
                          mean(wt_SPO11oligos_Rep2[,3]),
                          mean(wt_SPO11oligos_Rep3[,3]),
                          mean(kss_SPO11oligos_Rep1[,3]),
                          mean(kss_SPO11oligos_Rep2[,3]),
                          mean(wt_H3K9me2_Rep1[,3]),
                          mean(kss_H3K9me2_Rep1[,3])),
                 SD = c(sd(wt_REC8_MYC_Rep1[,3]),
                        sd(wt_REC8_HA_Rep1[,3]),
                        sd(wt_REC8_HA_Rep2[,3]),
                        sd(kss_REC8_HA_Rep1[,3]),
                        sd(kss_REC8_HA_Rep2[,3]),
                        sd(wt_SPO11oligos_Rep1[,3]),
                        sd(wt_SPO11oligos_Rep2[,3]),
                        sd(wt_SPO11oligos_Rep3[,3]),
                        sd(kss_SPO11oligos_Rep1[,3]),
                        sd(kss_SPO11oligos_Rep2[,3]),
                        sd(wt_H3K9me2_Rep1[,3]),
                        sd(kss_H3K9me2_Rep1[,3])))
write.table(df,
            file = "./wt_vs_kss_REC8_SPO11oligos_H3K9me2_coverage_means_SDs.tsv",
            sep = "\t", quote = F, row.names = F, col.names = T)
