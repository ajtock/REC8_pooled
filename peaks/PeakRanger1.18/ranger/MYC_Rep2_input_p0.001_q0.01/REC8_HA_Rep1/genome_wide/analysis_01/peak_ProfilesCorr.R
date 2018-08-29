#!/applications/R/R-3.3.2/bin/Rscript

# Profile mean coverage of REC8_HA_Rep1 and other
# chromatin marks around target and random loci

# Source functions to be used in this script
source("/projects/ajt200/Rfunctions/plotAvgCov_plotAvgMeth_target_ranLoc.R")
library(corrplot)

matDir <- "/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/genome_wide/analysis_01/matrices/"
plotDir <- "/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/genome_wide/analysis_01/plots/"

libNames <- c(
              "REC8_HA_Rep1",
              "REC8_MYC_Rep1",
              "REC8_MYC_Rep2",
              "WT_RNAseq_Chris_Rep1",
              "WT_RNAseq_Chris_Rep2",
              "WT_RNAseq_Kyuha_Rep1",
              "WT_RNAseq_Kyuha_Rep2",
              "WT_RNAseq_Kyuha_Rep3",
#              "WT_RNAseq_Kyuha_Rep1_bowtie2",
#              "WT_RNAseq_Kyuha_Rep2_bowtie2",
#              "WT_RNAseq_Kyuha_Rep3_bowtie2",
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
              "_norm_cov_peaks_mat1_smoothed_target_and_flank_dataframe_colMeans.txt"),
       paste0(matDir, libNames[[x]],
              "_norm_cov_ranLoc_mat2_smoothed_target_and_flank_dataframe_colMeans.txt")))

# Define column mean DNA methylation outfiles
DNAmethOutDFCM <- lapply(seq_along(DNAmethNames), function(x)
  list(paste0(matDir, DNAmethNames[[x]],
              "_norm_cov_peaks_mat1_smoothed_target_and_flank_dataframe_colMeans.txt"),
       paste0(matDir, DNAmethNames[[x]],
              "_norm_cov_ranLoc_mat2_smoothed_target_and_flank_dataframe_colMeans.txt")))

# Read in target and ranLoc mean coverage profiles
target_covDat <- lapply(seq_along(libNames), function(x) {
  read.table(file = outDFCM[[x]][[1]])[101:124,]
})
ranLoc_covDat <- lapply(seq_along(libNames), function(x) {
  read.table(file = outDFCM[[x]][[2]])[101:124,]
})
# Read in target and ranLoc mean DNA methylation profiles
target_DNAmethDat <- lapply(seq_along(DNAmethNames), function(x) {
  read.table(file = DNAmethOutDFCM[[x]][[1]])[101:124,]
})
ranLoc_DNAmethDat <- lapply(seq_along(DNAmethNames), function(x) {
  read.table(file = DNAmethOutDFCM[[x]][[2]])[101:124,]
})

target_Dat <- c(target_covDat, target_DNAmethDat)
ranLoc_Dat <- c(ranLoc_covDat, ranLoc_DNAmethDat)

# Redefine names for use in plots
libNames <- c(
              "REC8-HA Rep1",
              "REC8-MYC Rep1",
              "REC8-MYC Rep2",
              "RNA-seq (Chris) Rep1",
              "RNA-seq (Chris) Rep2",
              "RNA-seq (Kyuha) Rep1",
              "RNA-seq (Kyuha) Rep2",
              "RNA-seq (Kyuha) Rep3",
#              "RNA-seq (Kyuha) Rep1 bt2",
#              "RNA-seq (Kyuha) Rep2 bt2",
#              "RNA-seq (Kyuha) Rep3 bt2",
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
              "Pol V",
            
              "CG methylation",
              "CHG methylation",
              "CHH methylation"
             )

target_DatDF <- NULL
for(i in 1:length(target_Dat)) {
  target_DatDF <- cbind(target_DatDF, target_Dat[[i]])
}
colnames(target_DatDF) <- libNames
ranLoc_DatDF <- NULL
for(i in 1:length(ranLoc_Dat)) {
  ranLoc_DatDF <- cbind(ranLoc_DatDF, ranLoc_Dat[[i]])
}
colnames(ranLoc_DatDF) <- libNames

# Create and plot Spearman's rho correlation matrices for target and random loci mean profiles
pdf(paste0(plotDir, "REC8_HA_Rep1_peakProfileCorrMatrices_winSize20_Start_End.pdf"),
    height = 24, width = 12)
par(mfrow = c(2, 1))
par(mar = c(2.1, 3.2, 2.1, 3.2))
par(mgp = c(2.25, 1, 0))

# target loci
target_DatDF_corrMat <- cor(target_DatDF,
                            method = "spearman",
                            use = "pairwise.complete.obs")
col1 <- colorRampPalette(c("blue", "white", "red"))
corrplot(target_DatDF_corrMat,
         method = "color",
         type = "upper",
         col = col1(20),
         tl.col = "black",
         addgrid.col = "white",
         addCoef.col = "grey90",
         mar = c(0,0,1,0),
         tl.cex = 0.6,
         cl.cex = 0.6,
         number.cex = 0.5,
         title = "Peak profiles Spearman correlation matrix (Start to End)")
# random loci
ranLoc_DatDF_corrMat <- cor(ranLoc_DatDF,
                            method = "spearman",
                            use = "pairwise.complete.obs")
corrplot(ranLoc_DatDF_corrMat,
         method = "color",          
         type = "upper",
         col = col1(20),
         tl.col = "black",
         addgrid.col = "white",
         addCoef.col = "grey90",
         mar = c(0,0,1,0),
         tl.cex = 0.6,
         cl.cex = 0.6,
         number.cex = 0.5,
         title = "Random locus profiles Spearman correlation matrix (Start to End)")
dev.off()

