#!/applications/R/R-3.3.2/bin/Rscript

# Profile mean coverage of REC8_HA_Rep1 and other
# chromatin marks around target and random loci

# Source functions to be used in this script
source("/projects/ajt200/Rfunctions/plotAvgCov_plotAvgMeth_target_ranLoc.R")

matDir <- "/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/DiffBind/analysis_01_DB_sites_wt_exceeds_kss/matrices/"
plotDir <- "/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/DiffBind/analysis_01_DB_sites_wt_exceeds_kss/plots/"

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
              "PolV"
             )

# Define column mean coverage outfile (mean profiles)
outDFCM <- lapply(seq_along(libNames), function(x)
  list(paste0(matDir, libNames[[x]],
              "_norm_cov_peaks_mat1_smoothed_target_and_flank_dataframe_colMeans.txt"),
       paste0(matDir, libNames[[x]],
              "_norm_cov_ranLoc_mat2_smoothed_target_and_flank_dataframe_colMeans.txt")))

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
              "Wild type"
             )

YlabelNames <- c(
                 "REC8",
                 "SPO11-1-oligos",
                 "CHG methylation",
                 "H3K9me2",
                 "SPO11-1 ChIP",
                 "MNase",
                 "Pol IV",
                 "Pol V"
                )

# Plot mean REC8 vs other coverage profiles around target and random loci
pdf(paste0(plotDir, "DB_REC8_HA_peaks_wt_exceeds_kss_genome_wide_peakProfiles_wt_mutant_REC8_SPO11_1_oligos_CHGmeth_H3K9me2_SPO11_1_ChIP_MNase_PolIV_PolV_winSize20_manuscript.pdf"), height = 12.5, width = 6)
par(mfrow = c(5, 2))
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
    legendLoc = "topleft",
    legendLabs = c(libNames[7], libNames[8]),
    wtCol =  "blue", mutantCol = "red")

# SPO11-1 ChIP, MNase
plotAvgCov_WTvWTvsWTvsWT(xplot = xplot,
    dat1 = target_covDat[[9]][,1],
    dat2 = target_covDat[[10]][,1],
    dat3 = NULL,
    dat4 = NULL,
    ranDat1 = ranLoc_covDat[[9]][,1],
    ranDat2 = ranLoc_covDat[[10]][,1],
    ranDat3 = NULL,
    ranDat4 = NULL,
    flankSize = 2000, winSize = 20,
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "Start", endLab1 = "End",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "topleft",
    legendLabs = c("SPO11-1 ChIP", "MNase"),
    col1 = "green", col2 = "blue", col3 = NULL, col4 = NULL)
dev.off()

