#!/applications/R/R-3.3.2/bin/Rscript

# Profile mean coverage of REC8_HA_Rep1 and other
# chromatin marks around target and random loci

# Source functions to be used in this script
source("/projects/ajt200/Rfunctions/plotAvgCov_plotAvgMeth_target_ranLoc_hexiles.R")

matDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/arms/peakHexiles/log2_REC8_HA_Rep1_ChIP_MYC_Rep2_input/analysis_01/matrices/"
plotDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/arms/peakHexiles/log2_REC8_HA_Rep1_ChIP_MYC_Rep2_input/analysis_01/plots/"

libNames <- c(
              "REC8_HA_Rep1",
              "REC8_MYC_Rep1",
              "REC8_MYC_Rep2",
              "WT_RNAseq_Chris_Rep1",
              "WT_RNAseq_Chris_Rep2",
              "WT_RNAseq_Kyuha_Rep1",
              "WT_RNAseq_Kyuha_Rep2",
              "WT_RNAseq_Kyuha_Rep3",
              "WT_RNAseq_Kyuha_Rep1_bowtie2",
              "WT_RNAseq_Kyuha_Rep2_bowtie2",
              "WT_RNAseq_Kyuha_Rep3_bowtie2",
              "WT_RNAseq_meiocyte_Rep1",
              "WT_RNAseq_meiocyte_Rep2",
              "WT_RNAseq_meiocyte_Rep3",
              "MNase",
              "H3K4me1",
              "H3K4me2",
              "H3K4me3_ChIP12",
              "H3K4me3_ChIP14",
              "H3K4me3_ChIP15",
              "H2AZ",
 
              "H2A",
              "H2AW",
              "H2AX",
              "H3K9me2",
              "H3K27me1",
              "H3K27me3",

              "SPO11_1_oligos_RPI1",
              "SPO11_1_oligos_RPI3",
              "SPO11_1_oligos_RPI8",
              "SPO11_1_ChIP4",
              "SPO11_1_ChIP13",
              "MSH4",
              "PolIV_Rep2",
              "PolV"
             )
DNAmethNames <- c(
                  "CGmeth",
                  "CHGmeth",
                  "CHHmeth"
                 )

# Define column mean coverage outfile (mean profiles)
outDFCM <- lapply(seq_along(libNames), function(x)
  list(paste0(matDir, libNames[[x]],
              "_norm_cov_peak_hexile", c(1:6), "_by_log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_mat1_smoothed_target_and_flank_dataframe_colMeans.txt"),
       paste0(matDir, libNames[[x]],
              "_norm_cov_ranLoc_hexile", c(1:6), "_mat2_smoothed_target_and_flank_dataframe_colMeans.txt")))

# Define column mean DNA methylation outfiles
DNAmethOutDFCM <- lapply(seq_along(DNAmethNames), function(x)
  list(paste0(matDir, DNAmethNames[[x]],
              "_norm_cov_peak_hexile", c(1:6), "_by_log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_mat1_smoothed_target_and_flank_dataframe_colMeans.txt"),
       paste0(matDir, DNAmethNames[[x]],
              "_norm_cov_ranLoc_hexile", c(1:6), "_mat2_smoothed_target_and_flank_dataframe_colMeans.txt")))

# Arrange target hexile coverage profiles into a list object
# in which each top-level element in the list corresponds to a different dataset (e.g., REC8)
targetHexileCovProfiles <- list()
for(x in 1:length(outDFCM)) {
  subTarget <- NULL
  for(y in 1:6) {
    subTarget <- cbind(subTarget, read.table(file = outDFCM[[x]][[1]][y])[,1])
  }
  targetHexileCovProfiles[[x]] <- subTarget
}
# Arrange random locus hexile coverage profiles into a list object
# in which each top-level element in the list corresponds to a different dataset (e.g., REC8)
ranLocHexileCovProfiles <- list()
for(x in 1:length(outDFCM)) {
  subRanLoc <- NULL
  for(y in 1:6) {
    subRanLoc <- cbind(subRanLoc, read.table(file = outDFCM[[x]][[2]][y])[,1])
  }
  ranLocHexileCovProfiles[[x]] <- subRanLoc
}

# Arrange target hexile DNA methylation profiles into a list object
# in which each top-level element in the list corresponds to a different dataset (e.g., REC8)
targetHexileDNAmethProfiles <- list()
for(x in 1:length(DNAmethOutDFCM)) {
  subTarget <- NULL
  for(y in 1:6) {
    subTarget <- cbind(subTarget, read.table(file = DNAmethOutDFCM[[x]][[1]][y])[,1])
  }
  targetHexileDNAmethProfiles[[x]] <- subTarget
}
# Arrange random locus hexile DNA methylation profiles into a list object
# in which each top-level element in the list corresponds to a different dataset (e.g., REC8)
ranLocHexileDNAmethProfiles <- list()
for(x in 1:length(DNAmethOutDFCM)) {
  subRanLoc <- NULL
  for(y in 1:6) {
    subRanLoc <- cbind(subRanLoc, read.table(file = DNAmethOutDFCM[[x]][[2]][y])[,1])
  }
  ranLocHexileDNAmethProfiles[[x]] <- subRanLoc
}

# Redefine names for use in plots
libNames <- c(
              "REC8-HA Rep1",
              "REC8-MYC Rep1",
              "REC8-MYC Rep2",
              "RNA-seq (buds; Chris) Rep1",
              "RNA-seq (buds; Chris) Rep2",
              "RNA-seq (buds; Kyuha) Rep1",
              "RNA-seq (buds; Kyuha) Rep2",
              "RNA-seq (buds; Kyuha) Rep3",
              "RNA-seq (buds; Kyuha) Rep1 bt2",
              "RNA-seq (buds; Kyuha) Rep2 bt2",
              "RNA-seq (buds; Kyuha) Rep3 bt2",
              "RNA-seq (meiocytes) Rep1",
              "RNA-seq (meiocytes) Rep2",
              "RNA-seq (meiocytes) Rep3",
              "MNase",
              "H3K4me1",
              "H3K4me2",
              "H3K4me3 ChIP12",
              "H3K4me3 ChIP14",
              "H3K4me3 ChIP15",
              "H2A.Z",
 
              "H2A",
              "H2A.W",
              "H2A.X",
              "H3K9me2",
              "H3K27me1",
              "H3K27me3",

              "SPO11-1-oligos RPI1",
              "SPO11-1-oligos RPI3",
              "SPO11-1-oligos RPI8",
              "SPO11-1 ChIP4",
              "SPO11-1 ChIP13",
              "MSH4",
              "Pol IV",
              "Pol V"
             )
DNAmethNames <- c(
                  "CG methylation",
                  "CHG methylation",
                  "CHH methylation"
                 )

# Plot mean REC8 vs other coverage profiles around target and random loci
pdf(paste0(plotDir, "REC8_HA_Rep1_arm_peakHexileProfiles_winSize20.pdf"), height = 95, width = 6)
par(mfrow = c(38, 2))
par(mar = c(2.1, 3.2, 2.1, 3.2))
par(mgp = c(2.25, 1, 0))
xplot <- seq(1, length(targetHexileCovProfiles[[1]][,1]), by = 1)
colfunc <- colorRampPalette(c("red", "blue"))
mycols <- colfunc(6)

for(x in 1:length(libNames)) {
  plotAvgCov_hexiles(xplot = xplot, dat1 = targetHexileCovProfiles[[x]], ranDat1 = ranLocHexileCovProfiles[[x]],
                     flankSize = 2000, winSize = 20, Ylabel1 = libNames[x],
                     flankLabL = "-2 kb", flankLabR = "+2 kb",
                     startLab1 = "Start", endLab1 = "End",
                     startLab2 = "Start", endLab2 = "End",
                     mycols = mycols)
}
for(x in 1:length(DNAmethNames)) {
  plotAvgCov_hexiles(xplot = xplot, dat1 = targetHexileDNAmethProfiles[[x]], ranDat1 = ranLocHexileDNAmethProfiles[[x]],
                     flankSize = 2000, winSize = 20, Ylabel1 = DNAmethNames[x],
                     flankLabL = "-2 kb", flankLabR = "+2 kb",
                     startLab1 = "Start", endLab1 = "End",
                     startLab2 = "Start", endLab2 = "End",
                     mycols = mycols)
}
dev.off()

# Plot on multiple pages (landscape)
pdf(paste0(plotDir, "REC8_HA_Rep1_arm_peakHexileProfiles_winSize20_landscape_4pages.pdf"), height = 10, width = 18)
par(mfrow = c(4, 6))
par(mar = c(2.1, 3.2, 2.1, 3.2))
par(mgp = c(2.25, 1, 0))
xplot <- seq(1, length(targetHexileCovProfiles[[1]][,1]), by = 1)
colfunc <- colorRampPalette(c("red", "blue"))
mycols <- colfunc(6)

for(x in 1:length(libNames)) {
  plotAvgCov_hexiles(xplot = xplot, dat1 = targetHexileCovProfiles[[x]], ranDat1 = ranLocHexileCovProfiles[[x]],
                     flankSize = 2000, winSize = 20, Ylabel1 = libNames[x],
                     flankLabL = "-2 kb", flankLabR = "+2 kb",
                     startLab1 = "Start", endLab1 = "End",
                     startLab2 = "Start", endLab2 = "End",
                     mycols = mycols)
}
for(x in 1:length(DNAmethNames)) {
  plotAvgCov_hexiles(xplot = xplot, dat1 = targetHexileDNAmethProfiles[[x]], ranDat1 = ranLocHexileDNAmethProfiles[[x]],
                     flankSize = 2000, winSize = 20, Ylabel1 = DNAmethNames[x],
                     flankLabL = "-2 kb", flankLabR = "+2 kb",
                     startLab1 = "Start", endLab1 = "End",
                     startLab2 = "Start", endLab2 = "End",
                     mycols = mycols)
}
dev.off()

