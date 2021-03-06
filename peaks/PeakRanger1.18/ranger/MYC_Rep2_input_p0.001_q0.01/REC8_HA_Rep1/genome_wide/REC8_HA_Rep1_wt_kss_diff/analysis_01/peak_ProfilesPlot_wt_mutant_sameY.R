#!/applications/R/R-3.3.2/bin/Rscript

# Profile mean coverage of REC8_HA_Rep1 and other
# chromatin marks around target and random loci

# Source functions to be used in this script
source("/projects/ajt200/Rfunctions/plotAvgCov_plotAvgMeth_target_ranLoc.R")

matDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/genome_wide/REC8_HA_Rep1_wt_kss_diff/analysis_01/matrices/"
plotDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/genome_wide/REC8_HA_Rep1_wt_kss_diff/analysis_01/plots/"

libNames <- c(
              "REC8_HA_Rep1",
              "REC8_HA_Rep2",
              "REC8_MYC_Rep1",
              "kss_REC8_HA_Rep1",

              "SPO11_1_oligos_RPI1",
              "SPO11_1_oligos_RPI3",
              "SPO11_1_oligos_RPI8",
              "kss_SPO11_1_oligos_RPI34",
              "kss_SPO11_1_oligos_RPI35",

              "WT_RNAseq_Chris_Rep1",
              "WT_RNAseq_Chris_Rep2",
              "kss_RNAseq_Chris_Rep1",
              "kss_RNAseq_Chris_Rep2",

              "H3K9me2",
              "kss_H3K9me2",
              "cmt3_H3K9me2",

              "CGmeth",
              "kss_CGmeth",
              "cmt3_CGmeth",

              "CHGmeth",
              "kss_CHGmeth",
              "cmt3_CHGmeth",

              "CHHmeth",
              "kss_CHHmeth",
              "cmt3_CHHmeth",

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
              "wt HA Rep1",
              "wt HA Rep2",
              "wt MYC Rep1",
              "kss HA Rep1",

              "wt RPI1",
              "wt RPI3",
              "wt RPI8",
              "kss RPI34",
              "kss RPI35",

              "wt Rep1",
              "wt Rep2",
              "kss Rep1",
              "kss Rep2",

              "wt",
              "kss",
              "cmt3",

              "wt",
              "kss",
              "cmt3",

              "wt",
              "kss",
              "cmt3",

              "wt",
              "kss",
              "cmt3",

              "wt",
 
              "wt",
              "wt"
             )

YlabelNames <- c(
                 "REC8",
                 "SPO11-1-oligos",
                 "RNA-seq",
                 "H3K9me2",
                 "CG methylation proportion",
                 "CHG methylation proportion",
                 "CHH methylation proportion",
                 "MNase",
                 "Pol IV",
                 "Pol V"
                )

# Plot mean REC8 vs other coverage profiles around target and random loci
pdf(paste0(plotDir, "REC8_HA_Rep1_peaks_mean_log2_REC8_HA_Rep1_wt_kss_diff_cov_top10percent_genome_wide_peakProfiles_vs_wt_mutant_SPO11_1_oligos_RNAseq_H3K9me2_DNAmeth_MNase_PolIV_PolV_winSize20.pdf"), height = 40, width = 6)
par(mfrow = c(16, 2))
par(mar = c(2.1, 3.2, 2.1, 3.2))
par(mgp = c(2.25, 1, 0))
xplot <- seq(1, length(target_covDat[[1]][,1]), by = 1)
mycolsOneWTmutant <- c("red", "red4")
mycolsAnotherWTmutant <- c("blue", "navy")

plotAvgCov_oneWTmutantVanotherWTmutant(xplot = xplot,
    dat1 = target_covDat[[1]][,1],
    dat2 = target_covDat[[5]][,1],
    mutantDat1  = target_covDat[[4]][,1],
    mutantDat2 = target_covDat[[8]][,1],
    ranDat1 = ranLoc_covDat[[1]][,1],
    ranDat2 = ranLoc_covDat[[5]][,1],
    mutantRanDat1 = ranLoc_covDat[[4]][,1],
    mutantRanDat2 = ranLoc_covDat[[8]][,1],
    Ylabel1 = YlabelNames[1], Ylabel2 = YlabelNames[2],
    flankSize = 2000, winSize = 20,
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "Start", endLab1 = "End",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "topright",
    legendLabs = c(libNames[1], libNames[4], libNames[5], libNames[8]),
    mycolsDat1 = mycolsOneWTmutant,
    mycolsDat2 = mycolsAnotherWTmutant)
plotAvgCov_oneWTmutantVanotherWTmutant(xplot = xplot,
    dat1 = target_covDat[[2]][,1],
    dat2 = target_covDat[[6]][,1],
    mutantDat1  = target_covDat[[4]][,1],
    mutantDat2 = target_covDat[[8]][,1],
    ranDat1 = ranLoc_covDat[[2]][,1],
    ranDat2 = ranLoc_covDat[[6]][,1],
    mutantRanDat1 = ranLoc_covDat[[4]][,1],
    mutantRanDat2 = ranLoc_covDat[[8]][,1],
    Ylabel1 = YlabelNames[1], Ylabel2 = YlabelNames[2],
    flankSize = 2000, winSize = 20,
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "Start", endLab1 = "End",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "topright",
    legendLabs = c(libNames[2], libNames[4], libNames[6], libNames[8]),
    mycolsDat1 = mycolsOneWTmutant,
    mycolsDat2 = mycolsAnotherWTmutant)
plotAvgCov_oneWTmutantVanotherWTmutant(xplot = xplot,
    dat1 = target_covDat[[3]][,1],
    dat2 = target_covDat[[7]][,1],
    mutantDat1  = target_covDat[[4]][,1],
    mutantDat2 = target_covDat[[9]][,1],
    ranDat1 = ranLoc_covDat[[3]][,1],
    ranDat2 = ranLoc_covDat[[7]][,1],
    mutantRanDat1 = ranLoc_covDat[[4]][,1],
    mutantRanDat2 = ranLoc_covDat[[9]][,1],
    Ylabel1 = YlabelNames[1], Ylabel2 = YlabelNames[2],
    flankSize = 2000, winSize = 20, 
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "Start", endLab1 = "End",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "topright",
    legendLabs = c(libNames[3], libNames[4], libNames[7], libNames[9]),
    mycolsDat1 = mycolsOneWTmutant,
    mycolsDat2 = mycolsAnotherWTmutant)

plotAvgCov_oneWTmutantVanotherWT(xplot = xplot,
    dat1 = target_covDat[[1]][,1],
    dat2 = target_covDat[[26]][,1],
    mutantDat1 = target_covDat[[4]][,1],
    ranDat1 = ranLoc_covDat[[1]][,1],
    ranDat2 = ranLoc_covDat[[26]][,1],
    mutantRanDat1 = ranLoc_covDat[[4]][,1],
    Ylabel1 = YlabelNames[1], Ylabel2 = YlabelNames[8],
    flankSize = 2000, winSize = 20,
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "Start", endLab1 = "End",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "topright",
    legendLabs = c(libNames[1], libNames[4], libNames[26]),
    mycolsDat1 = mycolsOneWTmutant,
    mycolsDat2 = mycolsAnotherWTmutant[1])

plotAvgCov_oneWTmutantVanotherWT(xplot = xplot,
    dat1 = target_covDat[[1]][,1],
    dat2 = target_covDat[[27]][,1],
    mutantDat1 = target_covDat[[4]][,1],
    ranDat1 = ranLoc_covDat[[1]][,1],
    ranDat2 = ranLoc_covDat[[27]][,1],
    mutantRanDat1 = ranLoc_covDat[[4]][,1],
    Ylabel1 = YlabelNames[1], Ylabel2 = YlabelNames[9],
    flankSize = 2000, winSize = 20,
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "Start", endLab1 = "End",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "topright",
    legendLabs = c(libNames[1], libNames[4], libNames[27]),
    mycolsDat1 = mycolsOneWTmutant,
    mycolsDat2 = mycolsAnotherWTmutant[1])
plotAvgCov_oneWTmutantVanotherWT(xplot = xplot,
    dat1 = target_covDat[[1]][,1],
    dat2 = target_covDat[[28]][,1],
    mutantDat1 = target_covDat[[4]][,1],
    ranDat1 = ranLoc_covDat[[1]][,1],
    ranDat2 = ranLoc_covDat[[28]][,1],
    mutantRanDat1 = ranLoc_covDat[[4]][,1],
    Ylabel1 = YlabelNames[1], Ylabel2 = YlabelNames[10],
    flankSize = 2000, winSize = 20,
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "Start", endLab1 = "End",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "topright",
    legendLabs = c(libNames[1], libNames[4], libNames[28]),
    mycolsDat1 = mycolsOneWTmutant,
    mycolsDat2 = mycolsAnotherWTmutant[1])

plotAvgCov_oneWTmutantVanotherWTmutant(xplot = xplot,
    dat1 = target_covDat[[1]][,1],
    dat2 = target_covDat[[10]][,1],
    mutantDat1  = target_covDat[[4]][,1],
    mutantDat2 = target_covDat[[12]][,1],
    ranDat1 = ranLoc_covDat[[1]][,1],
    ranDat2 = ranLoc_covDat[[10]][,1],
    mutantRanDat1 = ranLoc_covDat[[4]][,1],
    mutantRanDat2 = ranLoc_covDat[[12]][,1],
    Ylabel1 = YlabelNames[1], Ylabel2 = YlabelNames[3],
    flankSize = 2000, winSize = 20,
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "Start", endLab1 = "End",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "topright",
    legendLabs = c(libNames[1], libNames[4], libNames[10], libNames[12]),
    mycolsDat1 = mycolsOneWTmutant,
    mycolsDat2 = mycolsAnotherWTmutant)
plotAvgCov_oneWTmutantVanotherWTmutant(xplot = xplot,
    dat1 = target_covDat[[1]][,1],
    dat2 = target_covDat[[11]][,1],
    mutantDat1  = target_covDat[[4]][,1],
    mutantDat2 = target_covDat[[13]][,1],
    ranDat1 = ranLoc_covDat[[1]][,1],
    ranDat2 = ranLoc_covDat[[11]][,1],
    mutantRanDat1 = ranLoc_covDat[[4]][,1],
    mutantRanDat2 = ranLoc_covDat[[13]][,1],
    Ylabel1 = YlabelNames[1], Ylabel2 = YlabelNames[3],
    flankSize = 2000, winSize = 20, 
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "Start", endLab1 = "End",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "topright",
    legendLabs = c(libNames[1], libNames[4], libNames[11], libNames[13]),
    mycolsDat1 = mycolsOneWTmutant,
    mycolsDat2 = mycolsAnotherWTmutant)
plotAvgCov_oneWTmutantVanotherWTmutant(xplot = xplot,
    dat1 = target_covDat[[1]][,1],
    dat2 = target_covDat[[14]][,1],
    mutantDat1  = target_covDat[[4]][,1],
    mutantDat2 = target_covDat[[15]][,1],
    ranDat1 = ranLoc_covDat[[1]][,1],
    ranDat2 = ranLoc_covDat[[14]][,1],
    mutantRanDat1 = ranLoc_covDat[[4]][,1],
    mutantRanDat2 = ranLoc_covDat[[15]][,1],
    Ylabel1 = YlabelNames[1], Ylabel2 = YlabelNames[4],
    flankSize = 2000, winSize = 20, 
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "Start", endLab1 = "End",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "topright",
    legendLabs = c(libNames[1], libNames[4], libNames[14], libNames[15]),
    mycolsDat1 = mycolsOneWTmutant,
    mycolsDat2 = mycolsAnotherWTmutant)
plotAvgCov_oneWTmutantVanotherWTmutant(xplot = xplot,
    dat1 = target_covDat[[1]][,1],
    dat2 = target_covDat[[14]][,1],
    mutantDat1  = target_covDat[[4]][,1],
    mutantDat2 = target_covDat[[16]][,1],
    ranDat1 = ranLoc_covDat[[1]][,1],
    ranDat2 = ranLoc_covDat[[14]][,1],
    mutantRanDat1 = ranLoc_covDat[[4]][,1],
    mutantRanDat2 = ranLoc_covDat[[16]][,1],
    Ylabel1 = YlabelNames[1], Ylabel2 = YlabelNames[4],
    flankSize = 2000, winSize = 20,
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "Start", endLab1 = "End",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "topright",
    legendLabs = c(libNames[1], libNames[4], libNames[14], libNames[16]),
    mycolsDat1 = mycolsOneWTmutant,
    mycolsDat2 = mycolsAnotherWTmutant)

plotAvgCov_oneWTmutantVanotherWTmutant(xplot = xplot,
    dat1 = target_covDat[[1]][,1],
    dat2 = target_covDat[[17]][,1],
    mutantDat1  = target_covDat[[4]][,1],
    mutantDat2 = target_covDat[[18]][,1],
    ranDat1 = ranLoc_covDat[[1]][,1],
    ranDat2 = ranLoc_covDat[[17]][,1],
    mutantRanDat1 = ranLoc_covDat[[4]][,1],
    mutantRanDat2 = ranLoc_covDat[[18]][,1],
    Ylabel1 = YlabelNames[1], Ylabel2 = YlabelNames[5],
    flankSize = 2000, winSize = 20,
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "Start", endLab1 = "End",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "topright",
    legendLabs = c(libNames[1], libNames[4], libNames[17], libNames[18]),
    mycolsDat1 = mycolsOneWTmutant,
    mycolsDat2 = mycolsAnotherWTmutant)
plotAvgCov_oneWTmutantVanotherWTmutant(xplot = xplot,
    dat1 = target_covDat[[1]][,1],
    dat2 = target_covDat[[17]][,1],
    mutantDat1  = target_covDat[[4]][,1],
    mutantDat2 = target_covDat[[19]][,1],
    ranDat1 = ranLoc_covDat[[1]][,1],
    ranDat2 = ranLoc_covDat[[17]][,1],
    mutantRanDat1 = ranLoc_covDat[[4]][,1],
    mutantRanDat2 = ranLoc_covDat[[19]][,1],
    Ylabel1 = YlabelNames[1], Ylabel2 = YlabelNames[5],
    flankSize = 2000, winSize = 20,
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "Start", endLab1 = "End",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "topright",
    legendLabs = c(libNames[1], libNames[4], libNames[17], libNames[19]),
    mycolsDat1 = mycolsOneWTmutant,
    mycolsDat2 = mycolsAnotherWTmutant)

plotAvgCov_oneWTmutantVanotherWTmutant(xplot = xplot,
    dat1 = target_covDat[[1]][,1],
    dat2 = target_covDat[[20]][,1],
    mutantDat1  = target_covDat[[4]][,1],
    mutantDat2 = target_covDat[[21]][,1],
    ranDat1 = ranLoc_covDat[[1]][,1],
    ranDat2 = ranLoc_covDat[[20]][,1],
    mutantRanDat1 = ranLoc_covDat[[4]][,1],
    mutantRanDat2 = ranLoc_covDat[[21]][,1],
    Ylabel1 = YlabelNames[1], Ylabel2 = YlabelNames[6],
    flankSize = 2000, winSize = 20,
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "Start", endLab1 = "End",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "topright",
    legendLabs = c(libNames[1], libNames[4], libNames[20], libNames[21]),
    mycolsDat1 = mycolsOneWTmutant,
    mycolsDat2 = mycolsAnotherWTmutant)
plotAvgCov_oneWTmutantVanotherWTmutant(xplot = xplot,
    dat1 = target_covDat[[1]][,1],
    dat2 = target_covDat[[20]][,1],
    mutantDat1  = target_covDat[[4]][,1],
    mutantDat2 = target_covDat[[22]][,1],
    ranDat1 = ranLoc_covDat[[1]][,1],
    ranDat2 = ranLoc_covDat[[20]][,1],
    mutantRanDat1 = ranLoc_covDat[[4]][,1],
    mutantRanDat2 = ranLoc_covDat[[22]][,1],
    Ylabel1 = YlabelNames[1], Ylabel2 = YlabelNames[6],
    flankSize = 2000, winSize = 20,
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "Start", endLab1 = "End",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "topright",
    legendLabs = c(libNames[1], libNames[4], libNames[20], libNames[22]),
    mycolsDat1 = mycolsOneWTmutant,
    mycolsDat2 = mycolsAnotherWTmutant)

plotAvgCov_oneWTmutantVanotherWTmutant(xplot = xplot,
    dat1 = target_covDat[[1]][,1],
    dat2 = target_covDat[[23]][,1],
    mutantDat1  = target_covDat[[4]][,1],
    mutantDat2 = target_covDat[[24]][,1],
    ranDat1 = ranLoc_covDat[[1]][,1],
    ranDat2 = ranLoc_covDat[[23]][,1],
    mutantRanDat1 = ranLoc_covDat[[4]][,1],
    mutantRanDat2 = ranLoc_covDat[[24]][,1],
    Ylabel1 = YlabelNames[1], Ylabel2 = YlabelNames[7],
    flankSize = 2000, winSize = 20,
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "Start", endLab1 = "End",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "topright",
    legendLabs = c(libNames[1], libNames[4], libNames[23], libNames[24]),
    mycolsDat1 = mycolsOneWTmutant,
    mycolsDat2 = mycolsAnotherWTmutant)
plotAvgCov_oneWTmutantVanotherWTmutant(xplot = xplot,
    dat1 = target_covDat[[1]][,1],
    dat2 = target_covDat[[23]][,1],
    mutantDat1  = target_covDat[[4]][,1],
    mutantDat2 = target_covDat[[25]][,1],
    ranDat1 = ranLoc_covDat[[1]][,1],
    ranDat2 = ranLoc_covDat[[23]][,1],
    mutantRanDat1 = ranLoc_covDat[[4]][,1],
    mutantRanDat2 = ranLoc_covDat[[25]][,1],
    Ylabel1 = YlabelNames[1], Ylabel2 = YlabelNames[7],
    flankSize = 2000, winSize = 20,
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "Start", endLab1 = "End",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "topright",
    legendLabs = c(libNames[1], libNames[4], libNames[23], libNames[25]),
    mycolsDat1 = mycolsOneWTmutant,
    mycolsDat2 = mycolsAnotherWTmutant)
dev.off()


