#############################################################################################
# Load log2 ratio of ChIP and input coverage values in 10-kb windows                   #
# and generate correlation matrices for chromosome arms, pericentromeric regions and genome-wide generating                                                     #
#############################################################################################

# REC8-HA Rep1 vs:

## Euchromatin:
# Genes
# H3K4me1
# H3K4me2
# H3K4me3
# H2A.Z
# H3K27me3

## Heterochromatin:
# Transposons
# CG, CHG and CHH DNA methylation
# H2A.W
# H3K9me2
# H3K27me1
# Nucleosomes

library(corrplot)
library(parallel)

# genomic definitions
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
# pericentromeric regions are as defined in Supplemental Table S26 of Ziolkowski et al. (2017) Genes Dev. 31
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)

#################################
# make cumulative genomes       #
#################################

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

windows <- c(10000)
winNames <- c("10kb")

plotDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/genomeProfiles/plots/corr_matrices/"

outDir1 <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/genomeProfiles/"
ChIPnames1 <- "REC8_HA_Rep1_ChIP"
names1 <- "REC8-HA"

# Euchromatic marks
outDir2 <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K4me1/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames2 <- "WT_H3K4me1_Rep1_ChIP"
names2 <- "H3K4me1"

outDir3 <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K4me2/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames3 <- "WT_H3K4me2_Rep1_ChIP"
names3 <- "H3K4me2"

outDir4 <- "/projects/ajt200/BAM_masters/H3K4me3/replicates/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames4 <- "WT_H3K4me3_ChIP14"
names4 <- "H3K4me3"

outDir5 <- "/projects/ajt200/BAM_masters/H2A/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames5 <- "H2AZ_ChIP"
names5 <- "H2A.Z"

outDir6 <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K27me3/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames6 <- "WT_H3K27me3_Rep1_ChIP"
names6 <- "H3K27me3"

# Heterochromatic marks
outDir7 <- "/projects/ajt200/BAM_masters/H2A/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames7 <- "H2AW_ChIP"
names7 <- "H2A.W"

outDir8 <- "/projects/ajt200/BAM_masters/H3K9me2/WT/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames8 <- "WT_H3K9me2_ChIP"
names8 <- "H3K9me2"

outDir9 <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K27me1/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames9 <- "WT_H3K27me1_Rep1_ChIP"
names9 <- "H3K27me1"

outDir10 <- "/projects/ajt200/BAM_masters/nucleosomes/WT/coverage/nakedDNA_untrimmed_input/log2ChIPinput/genomeProfiles/"
ChIPnames10 <- "WT_nuc"
names10 <- "Nucleosomes"

outDirCombined <- c(outDir1, outDir2, outDir3, outDir4, outDir5,
                    outDir6, outDir7, outDir8, outDir9, outDir10)
ChIPnamesCombined <- c(ChIPnames1, ChIPnames2, ChIPnames3, ChIPnames4, ChIPnames5,
                       ChIPnames6, ChIPnames7, ChIPnames8, ChIPnames9, ChIPnames10)
namesCombined <- c(names1, names2, names3, names4, names5,
                   names6, names7, names8, names9, names10)

# Load concatenated filtered, log2-transformed coverage datasets
filt_log2trans_list <- mclapply(seq_along(ChIPnamesCombined), function(k) {
  read.table(file = paste0(outDirCombined[k], "filt_log2_", ChIPnamesCombined[k], "_genome_norm_coverage_", winNames[1], ".txt"))
}, mc.cores = length(ChIPnamesCombined))
filt_noNA_log2trans_list <- mclapply(seq_along(ChIPnamesCombined), function(k) {
  read.table(file = paste0(outDirCombined[k], "filt_noNA_log2_", ChIPnamesCombined[k], "_genome_norm_coverage_", winNames[1], ".txt"))
}, mc.cores = length(ChIPnamesCombined))

inDirFeatures <- "/projects/ajt200/REC8_MSH4/data_merged_fastq/coverage/log2ChIPinput/genomeProfiles/"

# Load concatenated filtered gene dataset
filt_GeneDat <- read.table(file = paste0(inDirFeatures, "filt_gene_density_genome_", winNames[1], ".txt"))
filt_GeneDat_noNA <- read.table(file = paste0(inDirFeatures, "filt_noNA_gene_density_genome_", winNames[1], ".txt"))

# Load concatenated filtered TE dataset
filt_TEDat <- read.table(file = paste0(inDirFeatures, "filt_TE_density_genome_", winNames[1], ".txt"))
filt_TEDat_noNA <- read.table(file = paste0(inDirFeatures, "filt_noNA_TE_density_genome_", winNames[1], ".txt"))

# Load concatenated filtered DNA meth datasets
filt_methDat <- read.table(file = paste0(inDirFeatures, "filt_wtMeth_genome_10kb.txt"))
filt_noNA_methDat <- read.table(file = paste0(inDirFeatures, "filt_noNA_wtMeth_genome_10kb.txt"))


# Create FILTERED, log2-transformed genome-wide Spearman's rho correlation matrix

REC8_HA <- filt_log2trans_list[[1]]
## Euchromatin:
Genes <- filt_GeneDat
H3K4me1 <- filt_log2trans_list[[2]]
H3K4me2 <- filt_log2trans_list[[3]]
H3K4me3 <- filt_log2trans_list[[4]]
H2AZ <- filt_log2trans_list[[5]]
H3K27me3 <- filt_log2trans_list[[6]]
## Heterochromatin:
Transposons <- filt_TEDat
DNAmeth <- filt_methDat
H2AW <- filt_log2trans_list[[7]]
H3K9me2 <- filt_log2trans_list[[8]]
H3K27me1 <- filt_log2trans_list[[9]]
Nucleosomes <- filt_log2trans_list[[10]]

allDF <- data.frame(REC8_HA[,2], Genes[,2], H3K4me1[,2], H3K4me2[,2], H3K4me3[,2], H2AZ[,2], H3K27me3[,2],
                    Transposons[,2], DNAmeth[,2], DNAmeth[,3], DNAmeth[,4], H2AW[,2], H3K9me2[,2], H3K27me1[,2], Nucleosomes[,2])
colnames(allDF) <- c("REC8-HA", "Genes", "H3K4me1", "H3K4me2", "H3K4me3", "H2A.Z", "H3K27me3",
                     "Transposons", "CG", "CHG", "CHH", "H2A.W", "H3K9me2", "H3K27me1", "Nucleosomes")
allDF_corMat <- cor(allDF, method = "spearman", use = "pairwise.complete.obs")
col1 <- colorRampPalette(c("blue", "white", "red"))
pdf(file = paste0(plotDir, "REC8_HA_Rep1_filtered_genome-wide_correlation_matrix_", winNames[1], "_colouronly.pdf"))
corrplot(allDF_corMat, method = "color", type = "upper", col = col1(20), tl.col = "black",
         addgrid.col = "white", addCoef.col = "grey90", mar = c(0,0,1,0), tl.cex = 0.7, cl.cex = 0.7, number.cex = 0.7,
         title = paste0("Filtered genome-wide Spearman correlation matrix (10-kb windows)"))
dev.off()

#### FUNCTION NEEDS WORK ####
# Funnction to separate arms and pericentromeric regions to create separate Spearman's rho correlation matrices
separateArmPeri <- function(dat, dat_armL, dat_armR, dat_arm, dat_peri, sumchr, pericenStart, pericenEnd) {
  for(i in 1:5) {
    print(i)
    dat_armL_tmp <- dat[dat[,1] > sumchr[i] & dat[,1] < pericenStart[i] & dat[,1] < sumchr[i+1],]
    dat_armR_tmp <- dat[dat[,1] > pericenEnd[i] & dat[,1] <= sumchr[i+1],]
    dat_arm_tmp <- dat[(dat[,1] > sumchr[i] & dat[,1] < pericenStart[i] & dat[,1] < sumchr[i+1]) | (dat[,1] > pericenEnd[i] & dat[,1] <= sumchr[i+1]),]
    dat_peri_tmp <- dat[dat[,1] >= pericenStart[i] & dat[,1] <= pericenEnd[i],]
    dat_peri <- rbind(dat_peri, dat_peri_tmp)
    dat_armL <- rbind(dat_armL, dat_armL_tmp)
    dat_armR <- rbind(dat_armR, dat_armR_tmp)
    dat_arm <- rbind(dat_arm, dat_arm_tmp)
  }
}

#filt_log2trans_list_separateArmPeri <- lapply(seq_along(filt_log2trans_list), function(x) {
#    separateArmPeri(dat = filt_log2trans_list[[x]],
#                    dat_armL = NULL, dat_armR = NULL,
#                    dat_arm = NULL, dat_peri = NULL,
#                    sumchr = sumchr,
#                    pericenStart = pericenStart, pericenEnd = pericenEnd)
#})


# Separate arms and pericentromeric regions to create separate Spearman's rho correlation matrices
REC8_HA_armL <- NULL
REC8_HA_armR <- NULL
REC8_HA_arm <- NULL
REC8_HA_peri <- NULL

H3K4me1_armL <- NULL
H3K4me1_armR <- NULL
H3K4me1_arm <- NULL
H3K4me1_peri <- NULL

H3K4me2_armL <- NULL
H3K4me2_armR <- NULL
H3K4me2_arm <- NULL
H3K4me2_peri <- NULL

H3K4me3_armL <- NULL
H3K4me3_armR <- NULL
H3K4me3_arm <- NULL
H3K4me3_peri <- NULL

H2AZ_armL <- NULL
H2AZ_armR <- NULL
H2AZ_arm <- NULL
H2AZ_peri <- NULL

H3K27me3_armL <- NULL
H3K27me3_armR <- NULL
H3K27me3_arm <- NULL
H3K27me3_peri <- NULL

H2AW_armL <- NULL
H2AW_armR <- NULL
H2AW_arm <- NULL
H2AW_peri <- NULL

H3K9me2_armL <- NULL
H3K9me2_armR <- NULL
H3K9me2_arm <- NULL
H3K9me2_peri <- NULL

H3K27me1_armL <- NULL
H3K27me1_armR <- NULL
H3K27me1_arm <- NULL
H3K27me1_peri <- NULL

Nucleosomes_armL <- NULL
Nucleosomes_armR <- NULL
Nucleosomes_arm <- NULL
Nucleosomes_peri <- NULL

Genes_armL <- NULL
Genes_armR <- NULL
Genes_arm <- NULL
Genes_peri <- NULL

Transposons_armL <- NULL
Transposons_armR <- NULL
Transposons_arm <- NULL
Transposons_peri <- NULL

DNAmeth_armL <- NULL
DNAmeth_armR <- NULL
DNAmeth_arm <- NULL
DNAmeth_peri <- NULL
for(i in 1:5) {
  print(i)
  REC8_HA_armL_tmp <- REC8_HA[REC8_HA[,1] > sumchr[i] & REC8_HA[,1] < pericenStart[i] & REC8_HA[,1] < sumchr[i+1],]
  REC8_HA_armR_tmp <- REC8_HA[REC8_HA[,1] > pericenEnd[i] & REC8_HA[,1] <= sumchr[i+1],]
  REC8_HA_arm_tmp <- REC8_HA[(REC8_HA[,1] > sumchr[i] & REC8_HA[,1] < pericenStart[i] & REC8_HA[,1] < sumchr[i+1]) |
                                           (REC8_HA[,1] > pericenEnd[i] & REC8_HA[,1] <= sumchr[i+1]),]
  REC8_HA_peri_tmp <- REC8_HA[REC8_HA[,1] >= pericenStart[i] & REC8_HA[,1] <= pericenEnd[i],]
  REC8_HA_peri <- rbind(REC8_HA_peri, REC8_HA_peri_tmp)
  REC8_HA_armL <- rbind(REC8_HA_armL, REC8_HA_armL_tmp)
  REC8_HA_armR <- rbind(REC8_HA_armR, REC8_HA_armR_tmp)
  REC8_HA_arm <- rbind(REC8_HA_arm, REC8_HA_arm_tmp)

  print(i)
  H3K4me1_armL_tmp <- H3K4me1[H3K4me1[,1] > sumchr[i] & H3K4me1[,1] < pericenStart[i] & H3K4me1[,1] < sumchr[i+1],]
  H3K4me1_armR_tmp <- H3K4me1[H3K4me1[,1] > pericenEnd[i] & H3K4me1[,1] <= sumchr[i+1],]
  H3K4me1_arm_tmp <- H3K4me1[(H3K4me1[,1] > sumchr[i] & H3K4me1[,1] < pericenStart[i] & H3K4me1[,1] < sumchr[i+1]) |
                                           (H3K4me1[,1] > pericenEnd[i] & H3K4me1[,1] <= sumchr[i+1]),]
  H3K4me1_peri_tmp <- H3K4me1[H3K4me1[,1] >= pericenStart[i] & H3K4me1[,1] <= pericenEnd[i],]
  H3K4me1_peri <- rbind(H3K4me1_peri, H3K4me1_peri_tmp)
  H3K4me1_armL <- rbind(H3K4me1_armL, H3K4me1_armL_tmp)
  H3K4me1_armR <- rbind(H3K4me1_armR, H3K4me1_armR_tmp)
  H3K4me1_arm <- rbind(H3K4me1_arm, H3K4me1_arm_tmp)

  print(i)
  H3K4me2_armL_tmp <- H3K4me2[H3K4me2[,1] > sumchr[i] & H3K4me2[,1] < pericenStart[i] & H3K4me2[,1] < sumchr[i+1],]
  H3K4me2_armR_tmp <- H3K4me2[H3K4me2[,1] > pericenEnd[i] & H3K4me2[,1] <= sumchr[i+1],]
  H3K4me2_arm_tmp <- H3K4me2[(H3K4me2[,1] > sumchr[i] & H3K4me2[,1] < pericenStart[i] & H3K4me2[,1] < sumchr[i+1]) |
                                           (H3K4me2[,1] > pericenEnd[i] & H3K4me2[,1] <= sumchr[i+1]),]
  H3K4me2_peri_tmp <- H3K4me2[H3K4me2[,1] >= pericenStart[i] & H3K4me2[,1] <= pericenEnd[i],]
  H3K4me2_peri <- rbind(H3K4me2_peri, H3K4me2_peri_tmp)
  H3K4me2_armL <- rbind(H3K4me2_armL, H3K4me2_armL_tmp)
  H3K4me2_armR <- rbind(H3K4me2_armR, H3K4me2_armR_tmp)
  H3K4me2_arm <- rbind(H3K4me2_arm, H3K4me2_arm_tmp)

  print(i)
  H3K4me3_armL_tmp <- H3K4me3[H3K4me3[,1] > sumchr[i] & H3K4me3[,1] < pericenStart[i] & H3K4me3[,1] < sumchr[i+1],]
  H3K4me3_armR_tmp <- H3K4me3[H3K4me3[,1] > pericenEnd[i] & H3K4me3[,1] <= sumchr[i+1],]
  H3K4me3_arm_tmp <- H3K4me3[(H3K4me3[,1] > sumchr[i] & H3K4me3[,1] < pericenStart[i] & H3K4me3[,1] < sumchr[i+1]) |
                                           (H3K4me3[,1] > pericenEnd[i] & H3K4me3[,1] <= sumchr[i+1]),]
  H3K4me3_peri_tmp <- H3K4me3[H3K4me3[,1] >= pericenStart[i] & H3K4me3[,1] <= pericenEnd[i],]
  H3K4me3_peri <- rbind(H3K4me3_peri, H3K4me3_peri_tmp)
  H3K4me3_armL <- rbind(H3K4me3_armL, H3K4me3_armL_tmp)
  H3K4me3_armR <- rbind(H3K4me3_armR, H3K4me3_armR_tmp)
  H3K4me3_arm <- rbind(H3K4me3_arm, H3K4me3_arm_tmp)

  print(i)
  H2AZ_armL_tmp <- H2AZ[H2AZ[,1] > sumchr[i] & H2AZ[,1] < pericenStart[i] & H2AZ[,1] < sumchr[i+1],]
  H2AZ_armR_tmp <- H2AZ[H2AZ[,1] > pericenEnd[i] & H2AZ[,1] <= sumchr[i+1],]
  H2AZ_arm_tmp <- H2AZ[(H2AZ[,1] > sumchr[i] & H2AZ[,1] < pericenStart[i] & H2AZ[,1] < sumchr[i+1]) |
                                           (H2AZ[,1] > pericenEnd[i] & H2AZ[,1] <= sumchr[i+1]),]
  H2AZ_peri_tmp <- H2AZ[H2AZ[,1] >= pericenStart[i] & H2AZ[,1] <= pericenEnd[i],]
  H2AZ_peri <- rbind(H2AZ_peri, H2AZ_peri_tmp)
  H2AZ_armL <- rbind(H2AZ_armL, H2AZ_armL_tmp)
  H2AZ_armR <- rbind(H2AZ_armR, H2AZ_armR_tmp)
  H2AZ_arm <- rbind(H2AZ_arm, H2AZ_arm_tmp)

  print(i)
  H3K27me3_armL_tmp <- H3K27me3[H3K27me3[,1] > sumchr[i] & H3K27me3[,1] < pericenStart[i] & H3K27me3[,1] < sumchr[i+1],]
  H3K27me3_armR_tmp <- H3K27me3[H3K27me3[,1] > pericenEnd[i] & H3K27me3[,1] <= sumchr[i+1],]
  H3K27me3_arm_tmp <- H3K27me3[(H3K27me3[,1] > sumchr[i] & H3K27me3[,1] < pericenStart[i] & H3K27me3[,1] < sumchr[i+1]) |
                                           (H3K27me3[,1] > pericenEnd[i] & H3K27me3[,1] <= sumchr[i+1]),]
  H3K27me3_peri_tmp <- H3K27me3[H3K27me3[,1] >= pericenStart[i] & H3K27me3[,1] <= pericenEnd[i],]
  H3K27me3_peri <- rbind(H3K27me3_peri, H3K27me3_peri_tmp)
  H3K27me3_armL <- rbind(H3K27me3_armL, H3K27me3_armL_tmp)
  H3K27me3_armR <- rbind(H3K27me3_armR, H3K27me3_armR_tmp)
  H3K27me3_arm <- rbind(H3K27me3_arm, H3K27me3_arm_tmp)

  print(i)
  H2AW_armL_tmp <- H2AW[H2AW[,1] > sumchr[i] & H2AW[,1] < pericenStart[i] & H2AW[,1] < sumchr[i+1],]
  H2AW_armR_tmp <- H2AW[H2AW[,1] > pericenEnd[i] & H2AW[,1] <= sumchr[i+1],]
  H2AW_arm_tmp <- H2AW[(H2AW[,1] > sumchr[i] & H2AW[,1] < pericenStart[i] & H2AW[,1] < sumchr[i+1]) |
                                           (H2AW[,1] > pericenEnd[i] & H2AW[,1] <= sumchr[i+1]),]
  H2AW_peri_tmp <- H2AW[H2AW[,1] >= pericenStart[i] & H2AW[,1] <= pericenEnd[i],]
  H2AW_peri <- rbind(H2AW_peri, H2AW_peri_tmp)
  H2AW_armL <- rbind(H2AW_armL, H2AW_armL_tmp)
  H2AW_armR <- rbind(H2AW_armR, H2AW_armR_tmp)
  H2AW_arm <- rbind(H2AW_arm, H2AW_arm_tmp)

  print(i)
  H3K9me2_armL_tmp <- H3K9me2[H3K9me2[,1] > sumchr[i] & H3K9me2[,1] < pericenStart[i] & H3K9me2[,1] < sumchr[i+1],]
  H3K9me2_armR_tmp <- H3K9me2[H3K9me2[,1] > pericenEnd[i] & H3K9me2[,1] <= sumchr[i+1],]
  H3K9me2_arm_tmp <- H3K9me2[(H3K9me2[,1] > sumchr[i] & H3K9me2[,1] < pericenStart[i] & H3K9me2[,1] < sumchr[i+1]) |
                                           (H3K9me2[,1] > pericenEnd[i] & H3K9me2[,1] <= sumchr[i+1]),]
  H3K9me2_peri_tmp <- H3K9me2[H3K9me2[,1] >= pericenStart[i] & H3K9me2[,1] <= pericenEnd[i],]
  H3K9me2_peri <- rbind(H3K9me2_peri, H3K9me2_peri_tmp)
  H3K9me2_armL <- rbind(H3K9me2_armL, H3K9me2_armL_tmp)
  H3K9me2_armR <- rbind(H3K9me2_armR, H3K9me2_armR_tmp)
  H3K9me2_arm <- rbind(H3K9me2_arm, H3K9me2_arm_tmp)

  print(i)
  H3K27me1_armL_tmp <- H3K27me1[H3K27me1[,1] > sumchr[i] & H3K27me1[,1] < pericenStart[i] & H3K27me1[,1] < sumchr[i+1],]
  H3K27me1_armR_tmp <- H3K27me1[H3K27me1[,1] > pericenEnd[i] & H3K27me1[,1] <= sumchr[i+1],]
  H3K27me1_arm_tmp <- H3K27me1[(H3K27me1[,1] > sumchr[i] & H3K27me1[,1] < pericenStart[i] & H3K27me1[,1] < sumchr[i+1]) |
                                           (H3K27me1[,1] > pericenEnd[i] & H3K27me1[,1] <= sumchr[i+1]),]
  H3K27me1_peri_tmp <- H3K27me1[H3K27me1[,1] >= pericenStart[i] & H3K27me1[,1] <= pericenEnd[i],]
  H3K27me1_peri <- rbind(H3K27me1_peri, H3K27me1_peri_tmp)
  H3K27me1_armL <- rbind(H3K27me1_armL, H3K27me1_armL_tmp)
  H3K27me1_armR <- rbind(H3K27me1_armR, H3K27me1_armR_tmp)
  H3K27me1_arm <- rbind(H3K27me1_arm, H3K27me1_arm_tmp)

  print(i)
  Nucleosomes_armL_tmp <- Nucleosomes[Nucleosomes[,1] > sumchr[i] & Nucleosomes[,1] < pericenStart[i] & Nucleosomes[,1] < sumchr[i+1],]
  Nucleosomes_armR_tmp <- Nucleosomes[Nucleosomes[,1] > pericenEnd[i] & Nucleosomes[,1] <= sumchr[i+1],]
  Nucleosomes_arm_tmp <- Nucleosomes[(Nucleosomes[,1] > sumchr[i] & Nucleosomes[,1] < pericenStart[i] & Nucleosomes[,1] < sumchr[i+1]) |
                                           (Nucleosomes[,1] > pericenEnd[i] & Nucleosomes[,1] <= sumchr[i+1]),]
  Nucleosomes_peri_tmp <- Nucleosomes[Nucleosomes[,1] >= pericenStart[i] & Nucleosomes[,1] <= pericenEnd[i],]
  Nucleosomes_peri <- rbind(Nucleosomes_peri, Nucleosomes_peri_tmp)
  Nucleosomes_armL <- rbind(Nucleosomes_armL, Nucleosomes_armL_tmp)
  Nucleosomes_armR <- rbind(Nucleosomes_armR, Nucleosomes_armR_tmp)
  Nucleosomes_arm <- rbind(Nucleosomes_arm, Nucleosomes_arm_tmp)

  print(i)
  Genes_armL_tmp <- Genes[Genes[,1] > sumchr[i] & Genes[,1] < pericenStart[i] & Genes[,1] < sumchr[i+1],]
  Genes_armR_tmp <- Genes[Genes[,1] > pericenEnd[i] & Genes[,1] <= sumchr[i+1],]
  Genes_arm_tmp <- Genes[(Genes[,1] > sumchr[i] & Genes[,1] < pericenStart[i] & Genes[,1] < sumchr[i+1]) |
                                           (Genes[,1] > pericenEnd[i] & Genes[,1] <= sumchr[i+1]),]
  Genes_peri_tmp <- Genes[Genes[,1] >= pericenStart[i] & Genes[,1] <= pericenEnd[i],]
  Genes_peri <- rbind(Genes_peri, Genes_peri_tmp)
  Genes_armL <- rbind(Genes_armL, Genes_armL_tmp)
  Genes_armR <- rbind(Genes_armR, Genes_armR_tmp)
  Genes_arm <- rbind(Genes_arm, Genes_arm_tmp)

  print(i)
  Transposons_armL_tmp <- Transposons[Transposons[,1] > sumchr[i] & Transposons[,1] < pericenStart[i] & Transposons[,1] < sumchr[i+1],]
  Transposons_armR_tmp <- Transposons[Transposons[,1] > pericenEnd[i] & Transposons[,1] <= sumchr[i+1],]
  Transposons_arm_tmp <- Transposons[(Transposons[,1] > sumchr[i] & Transposons[,1] < pericenStart[i] & Transposons[,1] < sumchr[i+1]) |
                                           (Transposons[,1] > pericenEnd[i] & Transposons[,1] <= sumchr[i+1]),]
  Transposons_peri_tmp <- Transposons[Transposons[,1] >= pericenStart[i] & Transposons[,1] <= pericenEnd[i],]
  Transposons_peri <- rbind(Transposons_peri, Transposons_peri_tmp)
  Transposons_armL <- rbind(Transposons_armL, Transposons_armL_tmp)
  Transposons_armR <- rbind(Transposons_armR, Transposons_armR_tmp)
  Transposons_arm <- rbind(Transposons_arm, Transposons_arm_tmp)

  print(i)
  DNAmeth_armL_tmp <- DNAmeth[DNAmeth[,1] > sumchr[i] & DNAmeth[,1] < pericenStart[i] & DNAmeth[,1] < sumchr[i+1],]
  DNAmeth_armR_tmp <- DNAmeth[DNAmeth[,1] > pericenEnd[i] & DNAmeth[,1] <= sumchr[i+1],]
  DNAmeth_arm_tmp <- DNAmeth[(DNAmeth[,1] > sumchr[i] & DNAmeth[,1] < pericenStart[i] & DNAmeth[,1] < sumchr[i+1]) |
                                           (DNAmeth[,1] > pericenEnd[i] & DNAmeth[,1] <= sumchr[i+1]),]
  DNAmeth_peri_tmp <- DNAmeth[DNAmeth[,1] >= pericenStart[i] & DNAmeth[,1] <= pericenEnd[i],]
  DNAmeth_peri <- rbind(DNAmeth_peri, DNAmeth_peri_tmp)
  DNAmeth_armL <- rbind(DNAmeth_armL, DNAmeth_armL_tmp)
  DNAmeth_armR <- rbind(DNAmeth_armR, DNAmeth_armR_tmp)
  DNAmeth_arm <- rbind(DNAmeth_arm, DNAmeth_arm_tmp)
}

# Create arm and peri Spearman's rho correlation matrices
allDF_arm <- data.frame(REC8_HA_arm[,2], Genes_arm[,2], H3K4me1_arm[,2], H3K4me2_arm[,2], H3K4me3_arm[,2], H2AZ_arm[,2], H3K27me3_arm[,2],
                        Transposons_arm[,2], DNAmeth_arm[,2], DNAmeth_arm[,3], DNAmeth_arm[,4], H2AW_arm[,2], H3K9me2_arm[,2], H3K27me1_arm[,2], Nucleosomes_arm[,2])
colnames(allDF_arm) <- c("REC8-HA", "Genes", "H3K4me1", "H3K4me2", "H3K4me3", "H2A.Z", "H3K27me3",
                         "Transposons", "CG", "CHG", "CHH", "H2A.W", "H3K9me2", "H3K27me1", "Nucleosomes")
allDF_arm_corMat <- cor(allDF_arm, method = "spearman", use = "pairwise.complete.obs")
col1 <- colorRampPalette(c("blue", "white", "red"))
pdf(file = paste0(plotDir, "REC8_HA_Rep1_filtered_arms_correlation_matrix_", winNames[1], "_colouronly.pdf"))
corrplot(allDF_arm_corMat, method = "color", type = "upper", col = col1(20), tl.col = "black",
         addgrid.col = "white", addCoef.col = "grey90", mar = c(0,0,1,0), tl.cex = 0.7, cl.cex = 0.7, number.cex = 0.7,
         title = paste0("Filtered arms Spearman correlation matrix (10-kb windows)"))
dev.off()

allDF_peri <- data.frame(REC8_HA_peri[,2], Genes_peri[,2], H3K4me1_peri[,2], H3K4me2_peri[,2], H3K4me3_peri[,2], H2AZ_peri[,2], H3K27me3_peri[,2],
                         Transposons_peri[,2], DNAmeth_peri[,2], DNAmeth_peri[,3], DNAmeth_peri[,4], H2AW_peri[,2], H3K9me2_peri[,2], H3K27me1_peri[,2], Nucleosomes_peri[,2])
colnames(allDF_peri) <- c("REC8-HA", "Genes", "H3K4me1", "H3K4me2", "H3K4me3", "H2A.Z", "H3K27me3",
                          "Transposons", "CG", "CHG", "CHH", "H2A.W", "H3K9me2", "H3K27me1", "Nucleosomes")
allDF_peri_corMat <- cor(allDF_peri, method = "spearman", use = "pairwise.complete.obs")
col1 <- colorRampPalette(c("blue", "white", "red"))
pdf(file = paste0(plotDir, "REC8_HA_Rep1_filtered_pericentromeres_correlation_matrix_", winNames[1], "_colouronly.pdf"))
corrplot(allDF_peri_corMat, method = "color", type = "upper", col = col1(20), tl.col = "black",
         addgrid.col = "white", addCoef.col = "grey90", mar = c(0,0,1,0), tl.cex = 0.7, cl.cex = 0.7, number.cex = 0.7,
         title = paste0("Filtered pericentromeres Spearman correlation matrix (10-kb windows)"))      
dev.off()


