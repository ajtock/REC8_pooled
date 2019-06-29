#!/applications/R/R-3.3.2/bin/Rscript

##################################################################################################
# Load log2 ratio of ChIP and input coverage values in genomic windows                           #
# and generate correlation matrices for chromosome arms, pericentromeric regions and genome-wide #
##################################################################################################

# REC8-HA Rep1
# REC8-Myc Rep1
# Control-Myc Rep1

# Usage:
# ./log2_ChIPinput_perwin_arm_peri_rho_corr_matrices_REC8_v250319.R 10kb

args <- commandArgs(trailingOnly = T)
winName <- args[1]

library(corrplot)
library(parallel)

# Genomic definitions
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
# pericentromeric regions are as defined in Supplemental Table S26
# of Ziolkowski et al. (2017) Genes Dev. 31
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)

# Make cumulative genomes
sumchr <- cumsum(c(0, chrLens))
print(sumchr)
sumchr_tot <- sumchr[length(sumchr)]
print(sumchr_tot)

centromeres <- sapply(seq_along(centromeres), function(x) {
  centromeres[x] + sumchr[x]
})
print(centromeres)
pericenStart <- sapply(seq_along(pericenStart), function(x) {
  pericenStart[x] + sumchr[x]
})
print(pericenStart)
pericenEnd <- sapply(seq_along(pericenEnd), function(x) {
  pericenEnd[x] + sumchr[x]
})
print(pericenEnd)

plotDir <- "./plots/corr_matrices/"
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# REC8
REC8HADir <- "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/genomeProfiles/genomewideZscore/"
REC8HAprefix <- "REC8_HA_Rep1_ChIP_REC8_MYC_Rep1_input"
REC8HAname <- "REC8-HA"
 
REC8MYCDir <- "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/genomeProfiles/genomewideZscore/"
REC8MYCprefix <- "REC8_MYC_Rep1_ChIP_REC8_MYC_Rep1_input"
REC8MYCname <- "REC8-Myc"

controlHADir <- "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/genomeProfiles/genomewideZscore/"
controlHAprefix <- "control_HA_Rep1_ChIP_REC8_MYC_Rep1_input"
controlHAname <- "control-HA"

controlMYCDir <- "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/genomeProfiles/genomewideZscore/"
controlMYCprefix <- "control_MYC_Rep1_ChIP_REC8_MYC_Rep1_input"
controlMYCname <- "control-Myc"

# REC8
inDirREC8 <- c(REC8HADir, REC8MYCDir, controlHADir, controlMYCDir)
prefixREC8 <- c(REC8HAprefix, REC8MYCprefix, controlHAprefix, controlMYCprefix)
nameREC8 <- c(REC8HAname, REC8MYCname, controlHAname, controlMYCname)

filt_REC8log2trans_list <- mclapply(seq_along(nameREC8), function(x) {
  read.table(paste0(inDirREC8[x], "filt_log2_", prefixREC8[x],
                    "_genome_norm_coverage_Zscore_", winName, ".txt"))
}, mc.cores = length(nameREC8))

# Create FILTERED, log2-transformed genome-wide Spearman's rho correlation matrix

REC8_HA <- filt_REC8log2trans_list[[1]]
REC8_Myc <- filt_REC8log2trans_list[[2]]
Control_HA <- filt_REC8log2trans_list[[3]]
Control_Myc <- filt_REC8log2trans_list[[4]]

allDF <- data.frame(REC8_HA[,3], REC8_Myc[,3], Control_HA[,3], Control_Myc[,3])
colnames(allDF) <- c("REC8-HA", "REC8-Myc", "Control-HA", "Control-Myc")
allDF_corMat <- cor(allDF, method = "spearman", use = "pairwise.complete.obs")
col1 <- colorRampPalette(c("blue", "white", "red"))
pdf(file = paste0(plotDir, "REC8_Control_filtered_genomewide_correlation_matrix_", winName, "_colouronly_genomewideZscore_v250319.pdf"),
                  height = 8.5, width = 8)
corrplot(allDF_corMat, method = "color", type = "upper", col = col1(20), tl.col = "black",
         addgrid.col = "white", addCoef.col = "grey90", mar = c(0,0,1,0), tl.cex = 1.5, cl.cex = 1.5, number.cex = 1.0,
         title = paste0("Filtered genome-wide Spearman correlation matrix (10-kb windows)"))
dev.off()


# Separate arms and pericentromeric regions to create separate Spearman's rho correlation matrices
REC8_HA_armL <- NULL
REC8_HA_armR <- NULL
REC8_HA_arm <- NULL
REC8_HA_peri <- NULL

REC8_Myc_armL <- NULL
REC8_Myc_armR <- NULL
REC8_Myc_arm <- NULL
REC8_Myc_peri <- NULL

Control_HA_armL <- NULL
Control_HA_armR <- NULL
Control_HA_arm <- NULL
Control_HA_peri <- NULL

Control_Myc_armL <- NULL
Control_Myc_armR <- NULL
Control_Myc_arm <- NULL
Control_Myc_peri <- NULL
for(i in 1:5) {
  print(i)
  REC8_HA_armL_tmp <- REC8_HA[REC8_HA[,2] > sumchr[i] & REC8_HA[,2] < pericenStart[i] & REC8_HA[,2] < sumchr[i+1],]
  REC8_HA_armR_tmp <- REC8_HA[REC8_HA[,2] > pericenEnd[i] & REC8_HA[,2] <= sumchr[i+1],]
  REC8_HA_arm_tmp <- REC8_HA[(REC8_HA[,2] > sumchr[i] & REC8_HA[,2] < pericenStart[i] & REC8_HA[,2] < sumchr[i+1]) |
                                           (REC8_HA[,2] > pericenEnd[i] & REC8_HA[,2] <= sumchr[i+1]),]
  REC8_HA_peri_tmp <- REC8_HA[REC8_HA[,2] >= pericenStart[i] & REC8_HA[,2] <= pericenEnd[i],]
  REC8_HA_peri <- rbind(REC8_HA_peri, REC8_HA_peri_tmp)
  REC8_HA_armL <- rbind(REC8_HA_armL, REC8_HA_armL_tmp)
  REC8_HA_armR <- rbind(REC8_HA_armR, REC8_HA_armR_tmp)
  REC8_HA_arm <- rbind(REC8_HA_arm, REC8_HA_arm_tmp)

  print(i)
  REC8_Myc_armL_tmp <- REC8_Myc[REC8_Myc[,2] > sumchr[i] & REC8_Myc[,2] < pericenStart[i] & REC8_Myc[,2] < sumchr[i+1],]
  REC8_Myc_armR_tmp <- REC8_Myc[REC8_Myc[,2] > pericenEnd[i] & REC8_Myc[,2] <= sumchr[i+1],]
  REC8_Myc_arm_tmp <- REC8_Myc[(REC8_Myc[,2] > sumchr[i] & REC8_Myc[,2] < pericenStart[i] & REC8_Myc[,2] < sumchr[i+1]) |
                                           (REC8_Myc[,2] > pericenEnd[i] & REC8_Myc[,2] <= sumchr[i+1]),]
  REC8_Myc_peri_tmp <- REC8_Myc[REC8_Myc[,2] >= pericenStart[i] & REC8_Myc[,2] <= pericenEnd[i],]
  REC8_Myc_peri <- rbind(REC8_Myc_peri, REC8_Myc_peri_tmp)
  REC8_Myc_armL <- rbind(REC8_Myc_armL, REC8_Myc_armL_tmp)
  REC8_Myc_armR <- rbind(REC8_Myc_armR, REC8_Myc_armR_tmp)
  REC8_Myc_arm <- rbind(REC8_Myc_arm, REC8_Myc_arm_tmp)

  print(i)
  Control_HA_armL_tmp <- Control_HA[Control_HA[,2] > sumchr[i] & Control_HA[,2] < pericenStart[i] & Control_HA[,2] < sumchr[i+1],]
  Control_HA_armR_tmp <- Control_HA[Control_HA[,2] > pericenEnd[i] & Control_HA[,2] <= sumchr[i+1],]
  Control_HA_arm_tmp <- Control_HA[(Control_HA[,2] > sumchr[i] & Control_HA[,2] < pericenStart[i] & Control_HA[,2] < sumchr[i+1]) |
                                           (Control_HA[,2] > pericenEnd[i] & Control_HA[,2] <= sumchr[i+1]),]
  Control_HA_peri_tmp <- Control_HA[Control_HA[,2] >= pericenStart[i] & Control_HA[,2] <= pericenEnd[i],]
  Control_HA_peri <- rbind(Control_HA_peri, Control_HA_peri_tmp)
  Control_HA_armL <- rbind(Control_HA_armL, Control_HA_armL_tmp)
  Control_HA_armR <- rbind(Control_HA_armR, Control_HA_armR_tmp)
  Control_HA_arm <- rbind(Control_HA_arm, Control_HA_arm_tmp)

  print(i)
  Control_Myc_armL_tmp <- Control_Myc[Control_Myc[,2] > sumchr[i] & Control_Myc[,2] < pericenStart[i] & Control_Myc[,2] < sumchr[i+1],]
  Control_Myc_armR_tmp <- Control_Myc[Control_Myc[,2] > pericenEnd[i] & Control_Myc[,2] <= sumchr[i+1],]
  Control_Myc_arm_tmp <- Control_Myc[(Control_Myc[,2] > sumchr[i] & Control_Myc[,2] < pericenStart[i] & Control_Myc[,2] < sumchr[i+1]) |
                                           (Control_Myc[,2] > pericenEnd[i] & Control_Myc[,2] <= sumchr[i+1]),]
  Control_Myc_peri_tmp <- Control_Myc[Control_Myc[,2] >= pericenStart[i] & Control_Myc[,2] <= pericenEnd[i],]
  Control_Myc_peri <- rbind(Control_Myc_peri, Control_Myc_peri_tmp)
  Control_Myc_armL <- rbind(Control_Myc_armL, Control_Myc_armL_tmp)
  Control_Myc_armR <- rbind(Control_Myc_armR, Control_Myc_armR_tmp)
  Control_Myc_arm <- rbind(Control_Myc_arm, Control_Myc_arm_tmp)
}

# Create arm and peri Spearman's rho correlation matrices
allDF_arm <- data.frame(REC8_HA_arm[,3], REC8_Myc_arm[,3], Control_HA_arm[,3], Control_Myc_arm[,3])
colnames(allDF_arm) <- c("REC8-HA", "REC8-Myc", "Control-HA", "Control-Myc")
allDF_arm_corMat <- cor(allDF_arm, method = "spearman", use = "pairwise.complete.obs")
col1 <- colorRampPalette(c("blue", "white", "red"))
pdf(file = paste0(plotDir, "REC8_control_filtered_arms_correlation_matrix_", winName, "_colouronly_genomewideZscore_v250319.pdf"),
                  height = 8.5, width = 8)
corrplot(allDF_arm_corMat, method = "color", type = "upper", col = col1(20), tl.col = "black",
         addgrid.col = "white", addCoef.col = "grey90", mar = c(0,0,1,0), tl.cex = 1.5, cl.cex = 1.5, number.cex = 1.0,
         title = paste0("Filtered arms Spearman correlation matrix (10-kb windows)"))
dev.off()

allDF_peri <- data.frame(REC8_HA_peri[,3], REC8_Myc_peri[,3], Control_HA_peri[,3], Control_Myc_peri[,3])
colnames(allDF_peri) <- c("REC8-HA", "REC8-Myc", "Control-HA", "Control-Myc")
allDF_peri_corMat <- cor(allDF_peri, method = "spearman", use = "pairwise.complete.obs")
col1 <- colorRampPalette(c("blue", "white", "red"))
pdf(file = paste0(plotDir, "REC8_control_filtered_pericentromeres_correlation_matrix_", winName, "_colouronly_genomewideZscore_v250319.pdf"),
                  height = 8.5, width = 8)
corrplot(allDF_peri_corMat, method = "color", type = "upper", col = col1(20), tl.col = "black",
         addgrid.col = "white", addCoef.col = "grey90", mar = c(0,0,1,0), tl.cex = 1.5, cl.cex = 1.5, number.cex = 1.0,
         title = paste0("Filtered pericentromeres Spearman correlation matrix (10-kb windows)"))      
dev.off()


pdf(file = paste0(plotDir, "REC8_control_filtered_genomewide_arms_pericentromeres_correlation_matrix_", winName, "_colouronly_genomewideZscore.pdf"),
                  height = 25.5, width = 8)
par(mfcol = c(3, 1))
corrplot(allDF_corMat, method = "color", type = "upper", col = col1(20), tl.col = "black",
         addgrid.col = "white", addCoef.col = "grey90", mar = c(0,0,1,0), tl.cex = 1.5, cl.cex = 1.5, number.cex = 1.0,
         title = expression(italic("r"[s])~~~"Whole genome"))
corrplot(allDF_arm_corMat, method = "color", type = "upper", col = col1(20), tl.col = "black",
         addgrid.col = "white", addCoef.col = "grey90", mar = c(0,0,1,0), tl.cex = 1.5, cl.cex = 1.5, number.cex = 1.0,
         title = expression(italic("r"[s])~~~"Chromosome arms"))
corrplot(allDF_peri_corMat, method = "color", type = "upper", col = col1(20), tl.col = "black",
         addgrid.col = "white", addCoef.col = "grey90", mar = c(0,0,1,0), tl.cex = 1.5, cl.cex = 1.5, number.cex = 1.0,
         title = expression(italic("r"[s])~~~"Pericentromeres & centromeres"))    
dev.off()
