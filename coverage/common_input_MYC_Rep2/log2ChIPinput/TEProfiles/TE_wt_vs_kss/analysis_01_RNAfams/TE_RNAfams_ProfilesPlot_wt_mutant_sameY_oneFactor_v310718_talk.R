#!/applications/R/R-3.3.2/bin/Rscript

# Profile mean coverage of REC8_HA_Rep1 and other
# chromatin marks around target and random loci

# Source functions to be used in this script
source("/projects/ajt200/Rfunctions/plotAvgCov_plotAvgMeth_target_ranLoc.R")

matDir_RNA <- "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/TEProfiles/TE_wt_vs_kss/analysis_01_RNAfams/matrices/"
plotDir <- "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/TEProfiles/TE_wt_vs_kss/analysis_01_RNAfams/plots/"

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
pdf(paste0(plotDir, "rna_TEProfiles_vs_wt_mutant_REC8_SPO11_1_oligos_RNAseq_H3K9me2_DNAmeth_MNase_PolIV_PolV_winSize20_v310718_talk.pdf"), height = 10, width = 7.5)
par(mfcol = c(4, 3))
par(mar = c(2.1, 3.2, 2.1, 1))
par(mgp = c(2.25, 1, 0))
xplot <- seq(1, length(target_covDat_RNA[[1]][,1]), by = 1)

plotAvgCov_WTvMutant(xplot = xplot,
    dat1 = target_covDat_RNA[[1]][,1],
    mutantDat1  = target_covDat_RNA[[4]][,1],
    ranDat1 = ranLoc_covDat_RNA[[1]][,1],
    mutantRanDat1 = ranLoc_covDat_RNA[[4]][,1],
    Ylabel1 = YlabelNames[1],
    flankSize = 2000, winSize = 20, 
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "TSS", endLab1 = "TTS",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "topright",
    legendLabs = c(libNames[1], libNames[4]),
    wtCol =  "red", mutantCol = "red4")
plotAvgCov_WTvMutant(xplot = xplot,
    dat1 = target_covDat_RNA[[20]][,1],
    mutantDat1  = target_covDat_RNA[[21]][,1],
    ranDat1 = ranLoc_covDat_RNA[[20]][,1],
    mutantRanDat1 = ranLoc_covDat_RNA[[21]][,1],
    Ylabel1 = YlabelNames[6],
    flankSize = 2000, winSize = 20, 
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "TSS", endLab1 = "TTS",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "topright",
    legendLabs = c(libNames[20], libNames[21]),
    wtCol =  "green", mutantCol = "darkgreen")
plotAvgCov_WTvMutant(xplot = xplot,
    dat1 = target_covDat_RNA[[14]][,1],
    mutantDat1  = target_covDat_RNA[[15]][,1],
    ranDat1 = ranLoc_covDat_RNA[[14]][,1],
    mutantRanDat1 = ranLoc_covDat_RNA[[15]][,1],
    Ylabel1 = YlabelNames[4],
    flankSize = 2000, winSize = 20, 
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "TSS", endLab1 = "TTS",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "topright",
    legendLabs = c(libNames[14], libNames[15]),
    wtCol =  "deepskyblue", mutantCol = "deepskyblue4")
plotAvgCov_WTvMutant(xplot = xplot,
    dat1 = target_covDat_RNA[[23]][,1],
    mutantDat1  = target_covDat_RNA[[24]][,1],
    ranDat1 = ranLoc_covDat_RNA[[23]][,1],
    mutantRanDat1 = ranLoc_covDat_RNA[[24]][,1],
    Ylabel1 = YlabelNames[7],
    flankSize = 2000, winSize = 20, 
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "TSS", endLab1 = "TTS",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "topright",
    legendLabs = c(libNames[23], libNames[24]),
    wtCol =  "orange", mutantCol = "orange4")
plotAvgCov_WTvMutant(xplot = xplot,
    dat1 = target_covDat_RNA[[5]][,1],
    mutantDat1  = target_covDat_RNA[[8]][,1],
    ranDat1 = ranLoc_covDat_RNA[[5]][,1],
    mutantRanDat1 = ranLoc_covDat_RNA[[8]][,1],
    Ylabel1 = YlabelNames[2],
    flankSize = 2000, winSize = 20, 
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "TSS", endLab1 = "TTS",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "topright",
    legendLabs = c(libNames[5], libNames[8]),
    wtCol =  "blue", mutantCol = "navy")
plotAvgCov_WTvMutant(xplot = xplot,
    dat1 = target_covDat_RNA[[10]][,1],
    mutantDat1  = target_covDat_RNA[[12]][,1],
    ranDat1 = ranLoc_covDat_RNA[[10]][,1],
    mutantRanDat1 = ranLoc_covDat_RNA[[12]][,1],
    Ylabel1 = YlabelNames[3],
    flankSize = 2000, winSize = 20, 
    flankLabL = "-2 kb", flankLabR = "+2 kb",
    startLab1 = "TSS", endLab1 = "TTS",
    startLab2 = "Start", endLab2 = "End",
    legendLoc = "topright",
    legendLabs = c(libNames[10], libNames[12]),
    wtCol =  "magenta", mutantCol = "magenta4")
dev.off()
