#!/applications/R/R-3.3.2/bin/Rscript

# Profile mean coverage of REC8_HA_Rep1 and other
# chromatin marks around target and random loci

# Source functions to be used in this script
source("/projects/ajt200/Rfunctions/plotAvgCov_plotAvgMeth_target_ranLoc.R")

matDir <- "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/DMRprofiles/kss_hypoCHG/analysis_01/matrices/"
plotDir <- "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/DMRprofiles/kss_hypoCHG/analysis_01/plots/"

libNames <- c(
              "REC8_HA_Rep1",
              "kss_REC8_HA_Rep1",

              "SPO11_1_oligos_RPI1",
              "kss_SPO11_1_oligos_RPI34",

              "CHGmeth",
              "kss_CHGmeth",

              "H3K9me2",
              "kss_H3K9me2",
              
              "SPO11_ChIP4",
              "MNase",

              "PolIV_Rep2",
              "PolV",

              "WT_RNAseq_Chris_Rep1",
              "kss_RNAseq_Chris_Rep1",
              "WT_RNAseq_Chris_Rep2",
              "kss_RNAseq_Chris_Rep2"
             )

# Define column mean coverage outfile (mean profiles)
outDFCM <- lapply(seq_along(libNames), function(x)
  list(paste0(matDir, libNames[[x]],
              "_norm_cov_kss_hypoCHG_DMRs_mat1_smoothed_target_and_flank_dataframe_colMeans.txt"),
       paste0(matDir, libNames[[x]],
              "_norm_cov_kss_hypoCHG_ranLoc_mat2_smoothed_target_and_flank_dataframe_colMeans.txt")))

# Read in target and ranLoc mean coverage profiles
target_covDat <- lapply(seq_along(libNames), function(x) {
  read.table(file = outDFCM[[x]][[1]])
})
ranLoc_covDat <- lapply(seq_along(libNames), function(x) {
  read.table(file = outDFCM[[x]][[2]])
})

# Redefine names for use in plots
libNames <- c(
              "Wild type",
              "kyp suvh5 suvh6",

              "Wild type",
              "kyp suvh5 suvh6",

              "Wild type",
              "kyp suvh5 suvh6",

              "Wild type",
              "kyp suvh5 suvh6",

              "Wild type",
              "Wild type",

              "Wild type",
              "Wild type",

              "Wild type",
              "kyp suvh5 suvh6",
              "Wild type",
              "kyp suvh5 suvh6"
             )

YlabelNames <- c(
                 "REC8",
                 "SPO11-1-oligos",
                 "CHG methylation",
                 "H3K9me2",
                 "SPO11-1 ChIP",
                 "MNase",
                 "Pol IV",
                 "Pol V",
                 "RNA-seq Rep1",
                 "RNA-seq Rep2"
                )

# Plot mean REC8 vs other coverage profiles around target and random loci
pdf(paste0(plotDir, "kss_hypoCHG_DMR_Profiles_wt_mutant_REC8_SPO11_1_oligos_CHGmeth_H3K9me2_RNAseq_SPO11_1_ChIP_MNase_winSize20_manuscript.pdf"), height = 17.5, width = 6)
par(mfrow = c(7, 2))
par(mar = c(2.1, 3.2, 2.1, 3.2))
par(mgp = c(2.25, 1, 0))
xplot <- seq(1, length(target_covDat[[1]][,1]), by = 1)

# REC8
plotAvgCov_WTvMutant(xplot = xplot,
    dat1 = target_covDat[[1]][,1],
    mutantDat1  = target_covDat[[2]][,1],
    ranDat1 = ranLoc_covDat[[1]][,1],
    mutantRanDat1 = ranLoc_covDat[[2]][,1],
    Ylabel1 = YlabelNames[1],
    flankSize = 2000, winSize = 20, 
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "Start", endLab1 = "End",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "topleft",
    legendLabs = c(libNames[1], libNames[2]),
    wtCol =  "blue", mutantCol = "red")
# SPO11-1-oligos
plotAvgCov_WTvMutant(xplot = xplot,
    dat1 = target_covDat[[3]][,1],
    mutantDat1  = target_covDat[[4]][,1],
    ranDat1 = ranLoc_covDat[[3]][,1],
    mutantRanDat1 = ranLoc_covDat[[4]][,1],
    Ylabel1 = YlabelNames[2],
    flankSize = 2000, winSize = 20, 
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "Start", endLab1 = "End",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "left",
    legendLabs = c(libNames[3], libNames[4]),
    wtCol =  "blue", mutantCol = "red")
# CHGmeth
plotAvgCov_WTvMutant(xplot = xplot,
    dat1 = target_covDat[[5]][,1],
    mutantDat1  = target_covDat[[6]][,1],
    ranDat1 = ranLoc_covDat[[5]][,1],
    mutantRanDat1 = ranLoc_covDat[[6]][,1],
    Ylabel1 = YlabelNames[3],
    flankSize = 2000, winSize = 20, 
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "Start", endLab1 = "End",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "left",
    legendLabs = c(libNames[5], libNames[6]),
    wtCol =  "blue", mutantCol = "red")
# H3K9me2
plotAvgCov_WTvMutant(xplot = xplot,
    dat1 = target_covDat[[7]][,1],
    mutantDat1  = target_covDat[[8]][,1],
    ranDat1 = ranLoc_covDat[[7]][,1],
    mutantRanDat1 = ranLoc_covDat[[8]][,1],
    Ylabel1 = YlabelNames[4],
    flankSize = 2000, winSize = 20, 
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "Start", endLab1 = "End",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "left",
    legendLabs = c(libNames[7], libNames[8]),
    wtCol =  "blue", mutantCol = "red")
# RNA-seq rep1
plotAvgCov_WTvMutant(xplot = xplot,
    dat1 = target_covDat[[13]][,1],
    mutantDat1  = target_covDat[[14]][,1],
    ranDat1 = ranLoc_covDat[[13]][,1],
    mutantRanDat1 = ranLoc_covDat[[14]][,1],
    Ylabel1 = YlabelNames[9],
    flankSize = 2000, winSize = 20, 
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "Start", endLab1 = "End",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "topright",
    legendLabs = c(libNames[13], libNames[14]),
    wtCol =  "blue", mutantCol = "red")
# RNA-seq rep2
plotAvgCov_WTvMutant(xplot = xplot,
    dat1 = target_covDat[[15]][,1],
    mutantDat1  = target_covDat[[16]][,1],
    ranDat1 = ranLoc_covDat[[15]][,1],
    mutantRanDat1 = ranLoc_covDat[[16]][,1],
    Ylabel1 = YlabelNames[10],
    flankSize = 2000, winSize = 20, 
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "Start", endLab1 = "End",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "topright",
    legendLabs = c(libNames[15], libNames[16]),
    wtCol =  "blue", mutantCol = "red")

# SPO11-1 ChIP, MNase, Pol IV and Pol V
plotAvgCov_oneVanother(xplot = xplot,
    dat1 = target_covDat[[9]][,1],
    dat2 = target_covDat[[10]][,1],
    ranDat1 = ranLoc_covDat[[9]][,1],
    ranDat2 = ranLoc_covDat[[10]][,1],
    flankSize = 2000, winSize = 20,
    Ylabel1 = "SPO11-1 ChIP",
    Ylabel2 = "MNase",
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "Start", endLab1 = "End",
    startLab2 = "Start", endLab2 = "End",
    mycols = c("green", "blue"))
dev.off()
