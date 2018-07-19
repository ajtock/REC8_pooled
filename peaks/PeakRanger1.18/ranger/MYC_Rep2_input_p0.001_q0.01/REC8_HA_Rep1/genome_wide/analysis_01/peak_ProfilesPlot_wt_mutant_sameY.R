#!/applications/R/R-3.3.2/bin/Rscript

# Profile mean coverage of REC8_HA_Rep1 and other
# chromatin marks around target and random loci

# Source functions to be used in this script
source("/projects/ajt200/Rfunctions/plotAvgCov_plotAvgMeth_target_ranLoc.R")

matDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/genome_wide/analysis_01/matrices/"
plotDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/genome_wide/analysis_01/plots/"

libNames <- c(
              "REC8_HA_Rep1",

              "SPO11_1_oligos_RPI1",
              "SPO11_1_oligos_RPI3",
              "SPO11_1_oligos_RPI8",
              "kss_SPO11_1_oligos_RPI34",
              "kss_SPO11_1_oligos_RPI35",

              "H3K9me2",
              "kss_H3K9me2",
              "cmt3_H3K9me2"
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
              "REC8-HA Rep1",

              "wt RPI1",
              "wt RPI3",
              "wt RPI8",
              "kss RPI34",
              "kss RPI35",

              "wt",
              "kss",
              "cmt3"
             )

Ylabel2Names <- c(
                  "SPO11-1-oligos",
                  "H3K9me2"
                 )

# Plot mean REC8 vs other coverage profiles around target and random loci
pdf(paste0(plotDir, "REC8_HA_Rep1_genome_wide_peakProfiles_vs_wt_mutant_SPO11_1_oligos_H3K9me2_winSize20.pdf"), height = 20, width = 6)
par(mfrow = c(8, 2))
par(mar = c(2.1, 3.2, 2.1, 3.2))
par(mgp = c(2.25, 1, 0))
xplot <- seq(1, length(target_covDat[[1]][,1]), by = 1)
mycols <- c("red", "blue")
mycolsWTmutant <- c("blue", "deepskyblue1")

plotAvgCov_oneWTVanotherWTmutant(xplot = xplot,
                                 dat1 = target_covDat[[1]][,1],
                                 dat2 = target_covDat[[2]][,1],
                                 mutantdat2 = target_covDat[[5]][,1],
                                 ranDat1 = ranLoc_covDat[[1]][,1],
                                 ranDat2 = ranLoc_covDat[[2]][,1],
                                 mutantranDat2 = ranLoc_covDat[[5]][,1],
                                 Ylabel1 = libNames[1], Ylabel2 = Ylabel2Names[1],
                                 flankSize = 2000, winSize = 20,
                                 flankLabL = "-2 kb", flankLabR = "+2 kb",
                                 startLab1 = "Start", endLab1 = "End",
                                 startLab2 = "Start", endLab2 = "End",
                                 legendLoc = "topright",
                                 legendLabs = c(libNames[2], libNames[5]),
                                 mycolsDat1 = mycols,
                                 mycolsDat2 = mycolsWTmutant)
plotAvgCov_oneWTVanotherWTmutant(xplot = xplot,
                                 dat1 = target_covDat[[1]][,1],
                                 dat2 = target_covDat[[2]][,1],
                                 mutantdat2 = target_covDat[[6]][,1],
                                 ranDat1 = ranLoc_covDat[[1]][,1],
                                 ranDat2 = ranLoc_covDat[[2]][,1],
                                 mutantranDat2 = ranLoc_covDat[[6]][,1],
                                 Ylabel1 = libNames[1], Ylabel2 = Ylabel2Names[1],
                                 flankSize = 2000, winSize = 20,
                                 flankLabL = "-2 kb", flankLabR = "+2 kb",
                                 startLab1 = "Start", endLab1 = "End",
                                 startLab2 = "Start", endLab2 = "End",
                                 legendLoc = "topright",
                                 legendLabs = c(libNames[2], libNames[6]),
                                 mycolsDat1 = mycols,
                                 mycolsDat2 = mycolsWTmutant)
plotAvgCov_oneWTVanotherWTmutant(xplot = xplot,
                                 dat1 = target_covDat[[1]][,1],
                                 dat2 = target_covDat[[3]][,1],
                                 mutantdat2 = target_covDat[[5]][,1],
                                 ranDat1 = ranLoc_covDat[[1]][,1],
                                 ranDat2 = ranLoc_covDat[[3]][,1],
                                 mutantranDat2 = ranLoc_covDat[[5]][,1],
                                 Ylabel1 = libNames[1], Ylabel2 = Ylabel2Names[1],
                                 flankSize = 2000, winSize = 20,
                                 flankLabL = "-2 kb", flankLabR = "+2 kb",
                                 startLab1 = "Start", endLab1 = "End",
                                 startLab2 = "Start", endLab2 = "End",
                                 legendLoc = "topright",
                                 legendLabs = c(libNames[3], libNames[5]),
                                 mycolsDat1 = mycols,
                                 mycolsDat2 = mycolsWTmutant)
plotAvgCov_oneWTVanotherWTmutant(xplot = xplot,
                                 dat1 = target_covDat[[1]][,1],
                                 dat2 = target_covDat[[3]][,1],
                                 mutantdat2 = target_covDat[[6]][,1],
                                 ranDat1 = ranLoc_covDat[[1]][,1],
                                 ranDat2 = ranLoc_covDat[[3]][,1],
                                 mutantranDat2 = ranLoc_covDat[[6]][,1],
                                 Ylabel1 = libNames[1], Ylabel2 = Ylabel2Names[1],
                                 flankSize = 2000, winSize = 20,
                                 flankLabL = "-2 kb", flankLabR = "+2 kb",
                                 startLab1 = "Start", endLab1 = "End",
                                 startLab2 = "Start", endLab2 = "End",
                                 legendLoc = "topright",
                                 legendLabs = c(libNames[3], libNames[6]),
                                 mycolsDat1 = mycols,
                                 mycolsDat2 = mycolsWTmutant)
plotAvgCov_oneWTVanotherWTmutant(xplot = xplot,
                                 dat1 = target_covDat[[1]][,1],
                                 dat2 = target_covDat[[4]][,1],
                                 mutantdat2 = target_covDat[[5]][,1],
                                 ranDat1 = ranLoc_covDat[[1]][,1],
                                 ranDat2 = ranLoc_covDat[[4]][,1],
                                 mutantranDat2 = ranLoc_covDat[[5]][,1],
                                 Ylabel1 = libNames[1], Ylabel2 = Ylabel2Names[1],
                                 flankSize = 2000, winSize = 20,
                                 flankLabL = "-2 kb", flankLabR = "+2 kb",
                                 startLab1 = "Start", endLab1 = "End",
                                 startLab2 = "Start", endLab2 = "End",
                                 legendLoc = "topright",
                                 legendLabs = c(libNames[4], libNames[5]),
                                 mycolsDat1 = mycols,
                                 mycolsDat2 = mycolsWTmutant)
plotAvgCov_oneWTVanotherWTmutant(xplot = xplot,
                                 dat1 = target_covDat[[1]][,1],
                                 dat2 = target_covDat[[4]][,1],
                                 mutantdat2 = target_covDat[[6]][,1],
                                 ranDat1 = ranLoc_covDat[[1]][,1],
                                 ranDat2 = ranLoc_covDat[[4]][,1],
                                 mutantranDat2 = ranLoc_covDat[[6]][,1],
                                 Ylabel1 = libNames[1], Ylabel2 = Ylabel2Names[1],
                                 flankSize = 2000, winSize = 20,
                                 flankLabL = "-2 kb", flankLabR = "+2 kb",
                                 startLab1 = "Start", endLab1 = "End",
                                 startLab2 = "Start", endLab2 = "End",
                                 legendLoc = "topright",
                                 legendLabs = c(libNames[4], libNames[6]),
                                 mycolsDat1 = mycols,
                                 mycolsDat2 = mycolsWTmutant)


plotAvgCov_oneWTVanotherWTmutant(xplot = xplot,
                                 dat1 = target_covDat[[1]][,1],
                                 dat2 = target_covDat[[7]][,1],
                                 mutantdat2 = target_covDat[[8]][,1],
                                 ranDat1 = ranLoc_covDat[[1]][,1],
                                 ranDat2 = ranLoc_covDat[[7]][,1],
                                 mutantranDat2 = ranLoc_covDat[[8]][,1],
                                 Ylabel1 = libNames[1], Ylabel2 = Ylabel2Names[2],
                                 flankSize = 2000, winSize = 20,
                                 flankLabL = "-2 kb", flankLabR = "+2 kb",
                                 startLab1 = "Start", endLab1 = "End",
                                 startLab2 = "Start", endLab2 = "End",
                                 legendLoc = "topright",
                                 legendLabs = c(libNames[7], libNames[8]),
                                 mycolsDat1 = mycols,
                                 mycolsDat2 = mycolsWTmutant)
plotAvgCov_oneWTVanotherWTmutant(xplot = xplot,
                                 dat1 = target_covDat[[1]][,1],
                                 dat2 = target_covDat[[7]][,1],
                                 mutantdat2 = target_covDat[[9]][,1],
                                 ranDat1 = ranLoc_covDat[[1]][,1],
                                 ranDat2 = ranLoc_covDat[[7]][,1],
                                 mutantranDat2 = ranLoc_covDat[[9]][,1],
                                 Ylabel1 = libNames[1], Ylabel2 = Ylabel2Names[2],
                                 flankSize = 2000, winSize = 20,
                                 flankLabL = "-2 kb", flankLabR = "+2 kb",
                                 startLab1 = "Start", endLab1 = "End",
                                 startLab2 = "Start", endLab2 = "End",
                                 legendLoc = "topright",
                                 legendLabs = c(libNames[7], libNames[9]),
                                 mycolsDat1 = mycols,
                                 mycolsDat2 = mycolsWTmutant)
dev.off()


