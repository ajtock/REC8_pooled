#!/applications/R/R-3.3.2/bin/Rscript

# Profile mean coverage of REC8_HA_Rep1 and other
# chromatin marks around target and random loci

# Source functions to be used in this script
source("/projects/ajt200/Rfunctions/plotAvgCov_plotAvgMeth_target_ranLoc_sameYaxis.R")

matDir_DNA <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/TEProfiles/analysis_01_DNAfams/matrices/"
matDir_RNA <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/TEProfiles/analysis_01_RNAfams/matrices/"
plotDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/TEProfiles/analysis_01_DNAfams/plots/"

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
outDFCM_DNA <- lapply(seq_along(libNames), function(x)
  list(paste0(matDir_DNA, libNames[[x]],
              "_norm_cov_dna_TEs_mat1_smoothed_target_and_flank_dataframe_colMeans.txt"),
       paste0(matDir_DNA, libNames[[x]],
              "_norm_cov_dna_ranLoc_mat2_smoothed_target_and_flank_dataframe_colMeans.txt")))

# Define column mean DNA methylation outfiles
DNAmethOutDFCM_DNA <- lapply(seq_along(DNAmethNames), function(x)
  list(paste0(matDir_DNA, DNAmethNames[[x]],
              "_norm_cov_dna_TEs_mat1_smoothed_target_and_flank_dataframe_colMeans.txt"),
       paste0(matDir_DNA, DNAmethNames[[x]],
              "_norm_cov_dna_ranLoc_mat2_smoothed_target_and_flank_dataframe_colMeans.txt")))

# Define column mean coverage outfile (mean profiles)
outDFCM_RNA <- lapply(seq_along(libNames), function(x)
  list(paste0(matDir_RNA, libNames[[x]],
              "_norm_cov_rna_TEs_mat1_smoothed_target_and_flank_dataframe_colMeans.txt"),
       paste0(matDir_RNA, libNames[[x]],
              "_norm_cov_rna_ranLoc_mat2_smoothed_target_and_flank_dataframe_colMeans.txt")))

# Define column mean DNA methylation outfiles
DNAmethOutDFCM_RNA <- lapply(seq_along(DNAmethNames), function(x)
  list(paste0(matDir_RNA, DNAmethNames[[x]],
              "_norm_cov_rna_TEs_mat1_smoothed_target_and_flank_dataframe_colMeans.txt"),
       paste0(matDir_RNA, DNAmethNames[[x]],
              "_norm_cov_rna_ranLoc_mat2_smoothed_target_and_flank_dataframe_colMeans.txt")))


# Read in target and ranLoc mean coverage profiles
target_covDat_DNA <- lapply(seq_along(libNames), function(x) {
  read.table(file = outDFCM_DNA[[x]][[1]])
})
ranLoc_covDat_DNA <- lapply(seq_along(libNames), function(x) {
  read.table(file = outDFCM_DNA[[x]][[2]])
})
# Read in target and ranLoc mean DNA methylation profiles
target_DNAmethDat_DNA <- lapply(seq_along(DNAmethNames), function(x) {
  read.table(file = DNAmethOutDFCM_DNA[[x]][[1]])
})
ranLoc_DNAmethDat_DNA <- lapply(seq_along(DNAmethNames), function(x) {
  read.table(file = DNAmethOutDFCM_DNA[[x]][[2]])
})


# Read in target and ranLoc mean coverage profiles
target_covDat_RNA <- lapply(seq_along(libNames), function(x) {
  read.table(file = outDFCM_RNA[[x]][[1]])
})
ranLoc_covDat_RNA <- lapply(seq_along(libNames), function(x) {
  read.table(file = outDFCM_RNA[[x]][[2]])
})
# Read in target and ranLoc mean DNA methylation profiles
target_DNAmethDat_RNA <- lapply(seq_along(DNAmethNames), function(x) {
  read.table(file = DNAmethOutDFCM_RNA[[x]][[1]])
})
ranLoc_DNAmethDat_RNA <- lapply(seq_along(DNAmethNames), function(x) {
  read.table(file = DNAmethOutDFCM_RNA[[x]][[2]])
})


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
pdf(paste0(plotDir, "REC8_HA_Rep1_dna_and_rna_TEProfiles_winSize20_sameYaxis.pdf"), height = 95, width = 12)
par(mfrow = c(38, 4))
par(mar = c(2.1, 3.2, 2.1, 3.2))
par(mgp = c(2.25, 1, 0))
xplot_DNA <- seq(1, length(target_covDat_DNA[[1]][,1]), by = 1)
xplot_RNA <- seq(1, length(target_covDat_RNA[[1]][,1]), by = 1)
mycols <- c("red", "blue")
mycolsMeth <- c("navy", "blue", "deepskyblue1")

for(x in 2:length(libNames)) {
  plotAvgCov_oneVanother_sameY(xplot = xplot_DNA,
                               dat1 = target_covDat_DNA[[1]][,1], dat2 = target_covDat_DNA[[x]][,1],
                               ranDat1 = ranLoc_covDat_DNA[[1]][,1], ranDat2 = ranLoc_covDat_DNA[[x]][,1],

                               dat3 = target_covDat_RNA[[1]][,1], dat4 = target_covDat_RNA[[x]][,1],
                               ranDat3 = ranLoc_covDat_RNA[[1]][,1], ranDat4 = ranLoc_covDat_RNA[[x]][,1],

                               flankSize = 2000, winSize = 20,
                               Ylabel1 = libNames[1], Ylabel2 = libNames[x],
                               flankLabL = "-2 kb", flankLabR = "+2 kb",
                               startLab1 = "Start", endLab1 = "End",
                               startLab2 = "Start", endLab2 = "End",
                               mycols = mycols)
  plotAvgCov_oneVanother_sameY(xplot = xplot_RNA,
                               dat1 = target_covDat_RNA[[1]][,1], dat2 = target_covDat_RNA[[x]][,1],
                               ranDat1 = ranLoc_covDat_RNA[[1]][,1], ranDat2 = ranLoc_covDat_RNA[[x]][,1],

                               dat3 = target_covDat_DNA[[1]][,1], dat4 = target_covDat_DNA[[x]][,1],
                               ranDat3 = ranLoc_covDat_DNA[[1]][,1], ranDat4 = ranLoc_covDat_DNA[[x]][,1],

                               flankSize = 2000, winSize = 20,
                               Ylabel1 = libNames[1], Ylabel2 = libNames[x],
                               flankLabL = "-2 kb", flankLabR = "+2 kb",
                               startLab1 = "Start", endLab1 = "End",
                               startLab2 = "Start", endLab2 = "End",
                               mycols = mycols)
}

for(x in 1:length(DNAmethNames)) {
  plotAvgCov_oneVanother_sameY(xplot = xplot_DNA,
                               dat1 = target_covDat_DNA[[1]][,1], dat2 = target_DNAmethDat_DNA[[x]][,1],
                               ranDat1 = ranLoc_covDat_DNA[[1]][,1], ranDat2 = ranLoc_DNAmethDat_DNA[[x]][,1],

                               dat3 = target_covDat_RNA[[1]][,1], dat4 = target_DNAmethDat_RNA[[x]][,1],
                               ranDat3 = ranLoc_covDat_RNA[[1]][,1], ranDat4 = ranLoc_DNAmethDat_RNA[[x]][,1], 

                               flankSize = 2000, winSize = 20,
                               Ylabel1 = libNames[1], Ylabel2 = DNAmethNames[x],
                               flankLabL = "-2 kb", flankLabR = "+2 kb",
                               startLab1 = "Start", endLab1 = "End",
                               startLab2 = "Start", endLab2 = "End",
                               mycols = mycols)
  plotAvgCov_oneVanother_sameY(xplot = xplot_RNA,
                               dat1 = target_covDat_RNA[[1]][,1], dat2 = target_DNAmethDat_RNA[[x]][,1],
                               ranDat1 = ranLoc_covDat_RNA[[1]][,1], ranDat2 = ranLoc_DNAmethDat_RNA[[x]][,1],

                               dat3 = target_covDat_DNA[[1]][,1], dat4 = target_DNAmethDat_DNA[[x]][,1],
                               ranDat3 = ranLoc_covDat_DNA[[1]][,1], ranDat4 = ranLoc_DNAmethDat_DNA[[x]][,1],

                               flankSize = 2000, winSize = 20,
                               Ylabel1 = libNames[1], Ylabel2 = DNAmethNames[x],
                               flankLabL = "-2 kb", flankLabR = "+2 kb",
                               startLab1 = "Start", endLab1 = "End",
                               startLab2 = "Start", endLab2 = "End",
                               mycols = mycols)
}

plotAvgCov_plotAvgMeth_sameY(xplot = xplot_DNA,
                             dat1 = target_covDat_DNA[[1]][,1],
                             CGmethDat1 = target_DNAmethDat_DNA[[1]][,1],
                             CHGmethDat1 = target_DNAmethDat_DNA[[2]][,1],
                             CHHmethDat1 = target_DNAmethDat_DNA[[3]][,1],
                             ranDat1 = ranLoc_covDat_DNA[[1]][,1],
                             CGmethRanDat1 = ranLoc_DNAmethDat_DNA[[1]][,1],
                             CHGmethRanDat1 = ranLoc_DNAmethDat_DNA[[2]][,1],
                             CHHmethRanDat1 = ranLoc_DNAmethDat_DNA[[3]][,1],

                             dat2 = target_covDat_RNA[[1]][,1],
                             CGmethDat2 = target_DNAmethDat_RNA[[1]][,1],
                             CHGmethDat2 = target_DNAmethDat_RNA[[2]][,1],
                             CHHmethDat2 = target_DNAmethDat_RNA[[3]][,1],
                             ranDat2 = ranLoc_covDat_RNA[[1]][,1],
                             CGmethRanDat2 = ranLoc_DNAmethDat_RNA[[1]][,1],
                             CHGmethRanDat2 = ranLoc_DNAmethDat_RNA[[2]][,1],
                             CHHmethRanDat2 = ranLoc_DNAmethDat_RNA[[3]][,1],

                             flankSize = 2000, winSize = 20,
                             Ylabel1 = libNames[1], Ylabel2 = "DNA methylation",
                             flankLabL = "-2 kb", flankLabR = "+2 kb",
                             startLab1 = "Start", endLab1 = "End",
                             startLab2 = "Start", endLab2 = "End",
                             legendLoc = "right",
                             mycols = mycols, mycolsMeth = mycolsMeth)

plotAvgCov_plotAvgMeth_sameY(xplot = xplot_RNA,
                             dat1 = target_covDat_RNA[[1]][,1],
                             CGmethDat1 = target_DNAmethDat_RNA[[1]][,1],
                             CHGmethDat1 = target_DNAmethDat_RNA[[2]][,1],
                             CHHmethDat1 = target_DNAmethDat_RNA[[3]][,1],
                             ranDat1 = ranLoc_covDat_RNA[[1]][,1],
                             CGmethRanDat1 = ranLoc_DNAmethDat_RNA[[1]][,1],
                             CHGmethRanDat1 = ranLoc_DNAmethDat_RNA[[2]][,1],
                             CHHmethRanDat1 = ranLoc_DNAmethDat_RNA[[3]][,1],

                             dat2 = target_covDat_DNA[[1]][,1],
                             CGmethDat2 = target_DNAmethDat_DNA[[1]][,1],
                             CHGmethDat2 = target_DNAmethDat_DNA[[2]][,1],
                             CHHmethDat2 = target_DNAmethDat_DNA[[3]][,1],
                             ranDat2 = ranLoc_covDat_DNA[[1]][,1],
                             CGmethRanDat2 = ranLoc_DNAmethDat_DNA[[1]][,1],
                             CHGmethRanDat2 = ranLoc_DNAmethDat_DNA[[2]][,1],
                             CHHmethRanDat2 = ranLoc_DNAmethDat_DNA[[3]][,1],

                             flankSize = 2000, winSize = 20,
                             Ylabel1 = libNames[1], Ylabel2 = "DNA methylation",
                             flankLabL = "-2 kb", flankLabR = "+2 kb",
                             startLab1 = "Start", endLab1 = "End",
                             startLab2 = "Start", endLab2 = "End",
                             legendLoc = "right",
                             mycols = mycols, mycolsMeth = mycolsMeth)

dev.off()
