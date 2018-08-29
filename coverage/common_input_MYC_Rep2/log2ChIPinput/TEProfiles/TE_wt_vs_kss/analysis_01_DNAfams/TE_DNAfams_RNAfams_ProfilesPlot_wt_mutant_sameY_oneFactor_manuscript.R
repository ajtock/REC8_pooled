#!/applications/R/R-3.3.2/bin/Rscript

# Profile mean coverage of REC8_HA_Rep1 and other
# chromatin marks around target and random loci

# Source functions to be used in this script
source("/projects/ajt200/Rfunctions/plotAvgCov_plotAvgMeth_target_ranLoc.R")

matDir_DNA <- "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/TEProfiles/TE_wt_vs_kss/analysis_01_DNAfams/matrices/"
matDir_RNA <- "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/TEProfiles/TE_wt_vs_kss/analysis_01_RNAfams/matrices/"
plotDir <- "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/TEProfiles/TE_wt_vs_kss/analysis_01_DNAfams/plots/"

libNames <- c(
              "REC8_HA_Rep1",
              "kss_REC8_HA_Rep1",

              "SPO11_1_oligos_RPI1",
              "kss_SPO11_1_oligos_RPI34",

              "CGmeth",
              "kss_CGmeth",

              "CHGmeth",
              "kss_CHGmeth",

              "CHHmeth",
              "kss_CHHmeth",

              "H3K9me2",
              "kss_H3K9me2",
              
              "SPO11_ChIP4",
              "MNase",

              "PolIV_Rep2",
              "PolV"
             )

# Define column mean coverage outfile (mean profiles)
outDFCM_DNA <- lapply(seq_along(libNames), function(x)
  list(paste0(matDir_DNA, libNames[[x]],
              "_norm_cov_dna_TEs_mat1_smoothed_target_and_flank_dataframe_colMeans.txt"),
       paste0(matDir_DNA, libNames[[x]],
              "_norm_cov_dna_ranLoc_mat2_smoothed_target_and_flank_dataframe_colMeans.txt")))

# Read in target and ranLoc mean coverage profiles
target_covDat_DNA <- lapply(seq_along(libNames), function(x) {
  read.table(file = outDFCM_DNA[[x]][[1]])
})
ranLoc_covDat_DNA <- lapply(seq_along(libNames), function(x) {
  read.table(file = outDFCM_DNA[[x]][[2]])
})

# Define column mean coverage outfile (mean profiles)
outDFCM_RNA <- lapply(seq_along(libNames), function(x)
  list(paste0(matDir_RNA, libNames[[x]],
              "_norm_cov_rna_TEs_mat1_smoothed_target_and_flank_dataframe_colMeans.txt"),
       paste0(matDir_RNA, libNames[[x]],
              "_norm_cov_rna_ranLoc_mat2_smoothed_target_and_flank_dataframe_colMeans.txt")))

# Read in target and ranLoc mean coverage profiles
target_covDat_RNA <- lapply(seq_along(libNames), function(x) {
  read.table(file = outDFCM_RNA[[x]][[1]])
})
ranLoc_covDat_RNA <- lapply(seq_along(libNames), function(x) {
  read.table(file = outDFCM_RNA[[x]][[2]])
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
                 "CG methylation",
                 "CHG methylation",
                 "CHH methylation",
                 "H3K9me2",
                 "SPO11-1 ChIP",
                 "MNase",
                 "Pol IV",
                 "Pol V"
                )

# Plot mean REC8 vs other coverage profiles around target and random loci
pdf(paste0(plotDir, "dna_and_rna_TEProfiles_wt_mutant_REC8_SPO11_1_oligos_CGmeth_CHGmeth_CHHmeth_H3K9me2_SPO11_1_ChIP_MNase_PolIV_PolV_winSize20_manuscript.pdf"), height = 17.5, width = 12)
par(mfrow = c(7, 4))
par(mar = c(2.1, 3.2, 2.1, 3.2))
par(mgp = c(2.25, 1, 0))
xplot_DNA <- seq(1, length(target_covDat_DNA[[1]][,1]), by = 1)
xplot_RNA <- seq(1, length(target_covDat_RNA[[1]][,1]), by = 1)

# REC8
plotAvgCov_WTvMutant_DNA_RNA_TEs(xplot = xplot_DNA,
    DNAdat1 = target_covDat_DNA[[1]][,1],
    DNAmutantDat1  = target_covDat_DNA[[2]][,1],
    DNAranDat1 = ranLoc_covDat_DNA[[1]][,1],
    DNAmutantRanDat1 = ranLoc_covDat_DNA[[2]][,1],
    RNAdat1 = target_covDat_RNA[[1]][,1],
    RNAmutantDat1  = target_covDat_RNA[[2]][,1],
    RNAranDat1 = ranLoc_covDat_RNA[[1]][,1],
    RNAmutantRanDat1 = ranLoc_covDat_RNA[[2]][,1],
    Ylabel1 = YlabelNames[1],
    flankSize = 2000, winSize = 20, 
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "TSS", endLab1 = "TTS",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "bottomright",
    legendLabs = c(libNames[1], libNames[2]),
    wtCol =  "blue", mutantCol = "red")
plotAvgCov_WTvMutant_DNA_RNA_TEs(xplot = xplot_RNA,
    DNAdat1 = target_covDat_RNA[[1]][,1],
    DNAmutantDat1  = target_covDat_RNA[[2]][,1],
    DNAranDat1 = ranLoc_covDat_RNA[[1]][,1],
    DNAmutantRanDat1 = ranLoc_covDat_RNA[[2]][,1],
    RNAdat1 = target_covDat_DNA[[1]][,1],
    RNAmutantDat1  = target_covDat_DNA[[2]][,1],
    RNAranDat1 = ranLoc_covDat_DNA[[1]][,1],
    RNAmutantRanDat1 = ranLoc_covDat_DNA[[2]][,1],
    Ylabel1 = YlabelNames[1],
    flankSize = 2000, winSize = 20, 
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "TSS", endLab1 = "TTS",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "bottomright",
    legendLabs = c(libNames[1], libNames[2]),
    wtCol =  "blue", mutantCol = "red")
# SPO11-1-oligos
plotAvgCov_WTvMutant_DNA_RNA_TEs(xplot = xplot_DNA,
    DNAdat1 = target_covDat_DNA[[3]][,1],
    DNAmutantDat1  = target_covDat_DNA[[4]][,1],
    DNAranDat1 = ranLoc_covDat_DNA[[3]][,1],
    DNAmutantRanDat1 = ranLoc_covDat_DNA[[4]][,1],
    RNAdat1 = target_covDat_RNA[[3]][,1],
    RNAmutantDat1  = target_covDat_RNA[[4]][,1],
    RNAranDat1 = ranLoc_covDat_RNA[[3]][,1],
    RNAmutantRanDat1 = ranLoc_covDat_RNA[[4]][,1],
    Ylabel1 = YlabelNames[2],
    flankSize = 2000, winSize = 20, 
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "TSS", endLab1 = "TTS",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "bottomright",
    legendLabs = c(libNames[3], libNames[4]),
    wtCol =  "blue", mutantCol = "red")
plotAvgCov_WTvMutant_DNA_RNA_TEs(xplot = xplot_RNA,
    DNAdat1 = target_covDat_RNA[[3]][,1],
    DNAmutantDat1  = target_covDat_RNA[[4]][,1],
    DNAranDat1 = ranLoc_covDat_RNA[[3]][,1],
    DNAmutantRanDat1 = ranLoc_covDat_RNA[[4]][,1],
    RNAdat1 = target_covDat_DNA[[3]][,1],
    RNAmutantDat1  = target_covDat_DNA[[4]][,1],
    RNAranDat1 = ranLoc_covDat_DNA[[3]][,1],
    RNAmutantRanDat1 = ranLoc_covDat_DNA[[4]][,1],
    Ylabel1 = YlabelNames[2],
    flankSize = 2000, winSize = 20, 
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "TSS", endLab1 = "TTS",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "bottomright",
    legendLabs = c(libNames[3], libNames[4]),
    wtCol =  "blue", mutantCol = "red")
# CGmeth
plotAvgCov_WTvMutant_DNA_RNA_TEs(xplot = xplot_DNA,
    DNAdat1 = target_covDat_DNA[[5]][,1],
    DNAmutantDat1  = target_covDat_DNA[[6]][,1],
    DNAranDat1 = ranLoc_covDat_DNA[[5]][,1],
    DNAmutantRanDat1 = ranLoc_covDat_DNA[[6]][,1],
    RNAdat1 = target_covDat_RNA[[5]][,1],
    RNAmutantDat1  = target_covDat_RNA[[6]][,1],
    RNAranDat1 = ranLoc_covDat_RNA[[5]][,1],
    RNAmutantRanDat1 = ranLoc_covDat_RNA[[6]][,1],
    Ylabel1 = YlabelNames[3],
    flankSize = 2000, winSize = 20, 
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "TSS", endLab1 = "TTS",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "right",
    legendLabs = c(libNames[5], libNames[6]),
    wtCol =  "blue", mutantCol = "red")
plotAvgCov_WTvMutant_DNA_RNA_TEs(xplot = xplot_RNA,
    DNAdat1 = target_covDat_RNA[[5]][,1],
    DNAmutantDat1  = target_covDat_RNA[[6]][,1],
    DNAranDat1 = ranLoc_covDat_RNA[[5]][,1],
    DNAmutantRanDat1 = ranLoc_covDat_RNA[[6]][,1],
    RNAdat1 = target_covDat_DNA[[5]][,1],
    RNAmutantDat1  = target_covDat_DNA[[6]][,1],
    RNAranDat1 = ranLoc_covDat_DNA[[5]][,1],
    RNAmutantRanDat1 = ranLoc_covDat_DNA[[6]][,1],
    Ylabel1 = YlabelNames[3],
    flankSize = 2000, winSize = 20, 
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "TSS", endLab1 = "TTS",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "right",
    legendLabs = c(libNames[5], libNames[6]),
    wtCol =  "blue", mutantCol = "red")
# CHGmeth
plotAvgCov_WTvMutant_DNA_RNA_TEs(xplot = xplot_DNA,
    DNAdat1 = target_covDat_DNA[[7]][,1],
    DNAmutantDat1  = target_covDat_DNA[[8]][,1],
    DNAranDat1 = ranLoc_covDat_DNA[[7]][,1],
    DNAmutantRanDat1 = ranLoc_covDat_DNA[[8]][,1],
    RNAdat1 = target_covDat_RNA[[7]][,1],
    RNAmutantDat1  = target_covDat_RNA[[8]][,1],
    RNAranDat1 = ranLoc_covDat_RNA[[7]][,1],
    RNAmutantRanDat1 = ranLoc_covDat_RNA[[8]][,1],
    Ylabel1 = YlabelNames[4],
    flankSize = 2000, winSize = 20, 
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "TSS", endLab1 = "TTS",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "right",
    legendLabs = c(libNames[7], libNames[8]),
    wtCol =  "blue", mutantCol = "red")
plotAvgCov_WTvMutant_DNA_RNA_TEs(xplot = xplot_RNA,
    DNAdat1 = target_covDat_RNA[[7]][,1],
    DNAmutantDat1  = target_covDat_RNA[[8]][,1],
    DNAranDat1 = ranLoc_covDat_RNA[[7]][,1],
    DNAmutantRanDat1 = ranLoc_covDat_RNA[[8]][,1],
    RNAdat1 = target_covDat_DNA[[7]][,1],
    RNAmutantDat1  = target_covDat_DNA[[8]][,1],
    RNAranDat1 = ranLoc_covDat_DNA[[7]][,1],
    RNAmutantRanDat1 = ranLoc_covDat_DNA[[8]][,1],
    Ylabel1 = YlabelNames[4],
    flankSize = 2000, winSize = 20, 
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "TSS", endLab1 = "TTS",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "right",
    legendLabs = c(libNames[7], libNames[8]),
    wtCol =  "blue", mutantCol = "red")
# CHHmeth
plotAvgCov_WTvMutant_DNA_RNA_TEs(xplot = xplot_DNA,
    DNAdat1 = target_covDat_DNA[[9]][,1],
    DNAmutantDat1  = target_covDat_DNA[[10]][,1],
    DNAranDat1 = ranLoc_covDat_DNA[[9]][,1],
    DNAmutantRanDat1 = ranLoc_covDat_DNA[[10]][,1],
    RNAdat1 = target_covDat_RNA[[9]][,1],
    RNAmutantDat1  = target_covDat_RNA[[10]][,1],
    RNAranDat1 = ranLoc_covDat_RNA[[9]][,1],
    RNAmutantRanDat1 = ranLoc_covDat_RNA[[10]][,1],
    Ylabel1 = YlabelNames[5],
    flankSize = 2000, winSize = 20, 
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "TSS", endLab1 = "TTS",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "right",
    legendLabs = c(libNames[9], libNames[10]),
    wtCol =  "blue", mutantCol = "red")
plotAvgCov_WTvMutant_DNA_RNA_TEs(xplot = xplot_RNA,
    DNAdat1 = target_covDat_RNA[[9]][,1],
    DNAmutantDat1  = target_covDat_RNA[[10]][,1],
    DNAranDat1 = ranLoc_covDat_RNA[[9]][,1],
    DNAmutantRanDat1 = ranLoc_covDat_RNA[[10]][,1],
    RNAdat1 = target_covDat_DNA[[9]][,1],
    RNAmutantDat1  = target_covDat_DNA[[10]][,1],
    RNAranDat1 = ranLoc_covDat_DNA[[9]][,1],
    RNAmutantRanDat1 = ranLoc_covDat_DNA[[10]][,1],
    Ylabel1 = YlabelNames[5],
    flankSize = 2000, winSize = 20, 
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "TSS", endLab1 = "TTS",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "right",
    legendLabs = c(libNames[9], libNames[10]),
    wtCol =  "blue", mutantCol = "red")
# H3K9me2
plotAvgCov_WTvMutant_DNA_RNA_TEs(xplot = xplot_DNA,
    DNAdat1 = target_covDat_DNA[[11]][,1],
    DNAmutantDat1  = target_covDat_DNA[[12]][,1],
    DNAranDat1 = ranLoc_covDat_DNA[[11]][,1],
    DNAmutantRanDat1 = ranLoc_covDat_DNA[[12]][,1],
    RNAdat1 = target_covDat_RNA[[11]][,1],
    RNAmutantDat1  = target_covDat_RNA[[12]][,1],
    RNAranDat1 = ranLoc_covDat_RNA[[11]][,1],
    RNAmutantRanDat1 = ranLoc_covDat_RNA[[12]][,1],
    Ylabel1 = YlabelNames[6],
    flankSize = 2000, winSize = 20, 
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "TSS", endLab1 = "TTS",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "right",
    legendLabs = c(libNames[11], libNames[12]),
    wtCol =  "blue", mutantCol = "red")
plotAvgCov_WTvMutant_DNA_RNA_TEs(xplot = xplot_RNA,
    DNAdat1 = target_covDat_RNA[[11]][,1],
    DNAmutantDat1  = target_covDat_RNA[[12]][,1],
    DNAranDat1 = ranLoc_covDat_RNA[[11]][,1],
    DNAmutantRanDat1 = ranLoc_covDat_RNA[[12]][,1],
    RNAdat1 = target_covDat_DNA[[11]][,1],
    RNAmutantDat1  = target_covDat_DNA[[12]][,1],
    RNAranDat1 = ranLoc_covDat_DNA[[11]][,1],
    RNAmutantRanDat1 = ranLoc_covDat_DNA[[12]][,1],
    Ylabel1 = YlabelNames[6],
    flankSize = 2000, winSize = 20, 
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "TSS", endLab1 = "TTS",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "right",
    legendLabs = c(libNames[11], libNames[12]),
    wtCol =  "blue", mutantCol = "red")

# SPO11-1 ChIP, MNase, Pol IV and Pol V
plotAvgCov_WTvWTvsWTvsWT_DNA_RNA_TEs(xplot = xplot_DNA,
    DNAdat1 = target_covDat_DNA[[13]][,1],
    DNAdat2 = target_covDat_DNA[[14]][,1],
    DNAdat3 = target_covDat_DNA[[15]][,1],
    DNAdat4 = target_covDat_DNA[[16]][,1],
    DNAranDat1 = ranLoc_covDat_DNA[[13]][,1],
    DNAranDat2 = ranLoc_covDat_DNA[[14]][,1],
    DNAranDat3 = ranLoc_covDat_DNA[[15]][,1],
    DNAranDat4 = ranLoc_covDat_DNA[[16]][,1],
    RNAdat1 = target_covDat_RNA[[13]][,1],
    RNAdat2 = target_covDat_RNA[[14]][,1],
    RNAdat3 = target_covDat_RNA[[15]][,1],
    RNAdat4 = target_covDat_RNA[[16]][,1],
    RNAranDat1 = ranLoc_covDat_RNA[[13]][,1],
    RNAranDat2 = ranLoc_covDat_RNA[[14]][,1],
    RNAranDat3 = ranLoc_covDat_RNA[[15]][,1],
    RNAranDat4 = ranLoc_covDat_RNA[[16]][,1],
    flankSize = 2000, winSize = 20,
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "TSS", endLab1 = "TTS",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "bottomright",
    legendLabs = c("SPO11-1 ChIP", "MNase", "Pol IV", "Pol V"),
    col1 = "green", col2 = "blue", col3 = "red", col4 = "magenta")
plotAvgCov_WTvWTvsWTvsWT_DNA_RNA_TEs(xplot = xplot_RNA,
    DNAdat1 = target_covDat_RNA[[13]][,1],
    DNAdat2 = target_covDat_RNA[[14]][,1],
    DNAdat3 = target_covDat_RNA[[15]][,1],
    DNAdat4 = target_covDat_RNA[[16]][,1],
    DNAranDat1 = ranLoc_covDat_RNA[[13]][,1],
    DNAranDat2 = ranLoc_covDat_RNA[[14]][,1],
    DNAranDat3 = ranLoc_covDat_RNA[[15]][,1],
    DNAranDat4 = ranLoc_covDat_RNA[[16]][,1],
    RNAdat1 = target_covDat_DNA[[13]][,1],
    RNAdat2 = target_covDat_DNA[[14]][,1],
    RNAdat3 = target_covDat_DNA[[15]][,1],
    RNAdat4 = target_covDat_DNA[[16]][,1],
    RNAranDat1 = ranLoc_covDat_DNA[[13]][,1],
    RNAranDat2 = ranLoc_covDat_DNA[[14]][,1],
    RNAranDat3 = ranLoc_covDat_DNA[[15]][,1],
    RNAranDat4 = ranLoc_covDat_DNA[[16]][,1],
    flankSize = 2000, winSize = 20,
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "TSS", endLab1 = "TTS",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "bottomright",
    legendLabs = c("SPO11-1 ChIP", "MNase", "Pol IV", "Pol V"),
    col1 = "green", col2 = "blue", col3 = "red", col4 = "magenta")
dev.off()
