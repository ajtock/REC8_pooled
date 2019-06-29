#!/applications/R/R-3.3.2/bin/Rscript

##################################################################################################
# Load log2 ratio of ChIP and input coverage values in genomic windows                           #
# and generate correlation matrices for chromosome arms, pericentromeric regions and genome-wide #
##################################################################################################

# REC8-HA Rep2 vs:

## Recombination:
# SPO11-1-oligos
# Crossovers

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

# Usage:
# ./log2_ChIPinput_perwin_arm_peri_rho_corr_matrices_v170619.R 10kb

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
REC8HAprefix <- "REC8_HA_Rep2_ChIP_REC8_MYC_Rep1_input"
REC8HAname <- "REC8-HA"

REC8MYCDir <- "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/genomeProfiles/genomewideZscore/"
REC8MYCprefix <- "REC8_MYC_Rep1_ChIP_REC8_MYC_Rep1_input"
REC8MYCname <- "REC8-Myc"

# Heterochromatin
inDir1 <- "/projects/ajt200/BAM_masters/H2A/coverage/log2ChIPinput/genomeProfiles/genomewideZscore/"
prefix1 <- "H2AW_input"
name1 <- "H2A.W"

inDir2 <- "/projects/ajt200/BAM_masters/H3K9me2/WT/coverage/log2ChIPinput/genomeProfiles/genomewideZscore/"
prefix2 <- "WT_H3K9me2_ChIP_WT_H3K9me2_input"
name2 <- "H3K9me2"

inDir3 <- "/home/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K27me1/coverage/log2ChIPinput/genomeProfiles/genomewideZscore/"
prefix3 <- "WT_H3K27me1_Rep1_ChIP_WT_H3K9me2_input"
name3 <- "H3K27me1"

inDir4 <- "/projects/ajt200/BAM_masters/nucleosomes/WT/coverage/nakedDNA_untrimmed_input/log2ChIPinput/genomeProfiles/genomewideZscore/"
prefix4 <- "WT_nuc_nakedDNAuntrimmed"
name4 <- "Nucleosomes"

# Euchromatin
inDir5 <- "/home/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K4me1/coverage/log2ChIPinput/genomeProfiles/genomewideZscore/"
prefix5 <- "WT_H3K4me1_Rep1_ChIP_WT_H3K9me2_input"
name5 <- "H3K4me1"

inDir6 <- "/home/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K4me2/coverage/log2ChIPinput/genomeProfiles/genomewideZscore/"
prefix6 <- "WT_H3K4me2_Rep1_ChIP_WT_H3K9me2_input"
name6 <- "H3K4me2"

inDir7 <- "/projects/ajt200/BAM_masters/H3K4me3/replicates/coverage/log2ChIPinput/genomeProfiles/genomewideZscore/"
prefix7 <- "WT_H3K4me3_ChIP14_WT_H3K9me2_input"
name7 <- "H3K4me3"

inDir8 <- "/projects/ajt200/BAM_masters/H2A/coverage/log2ChIPinput/genomeProfiles/genomewideZscore/"
prefix8 <- "H2AZ_input"
name8 <- "H2A.Z"

inDir9 <- "/home/ajt200/analysis/H3K27me3_bud_UWMadison_2015/coverage/log2ChIPinput/genomeProfiles/genomewideZscore/"
prefix9 <- "H3K27me3_ChIP_SRR1509478_WT_H3K9me2_input"
name9 <- "H3K27me3"

# Recombination
inDir10 <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/genomeProfiles/genomewideZscore/"
prefix10 <- "WT_SPO11oligo_RPI1_WT_nakedDNA_R1"
name10 <- "SPO11-1-oligos"

# REC8
inDirREC8 <- c(REC8HADir, REC8MYCDir)
prefixREC8 <- c(REC8HAprefix, REC8MYCprefix)
nameREC8 <- c(REC8HAname, REC8MYCname)

filt_REC8log2trans_list <- mclapply(seq_along(nameREC8), function(x) {
  read.table(paste0(inDirREC8[x], "filt_log2_", prefixREC8[x],
                    "_genome_norm_coverage_Zscore_", winName, ".txt"))
}, mc.cores = length(nameREC8))

# All
inDirAll <- c(inDir1, inDir2, inDir3, inDir4, inDir5, inDir6, inDir7, inDir8, inDir9, inDir10)
prefixAll <- c(prefix1, prefix2, prefix3, prefix4, prefix5, prefix6, prefix7, prefix8, prefix9, prefix10)
nameAll <- c(name1, name2, name3, name4, name5, name6, name7, name8, name9, name10)

filt_log2trans_list <- mclapply(seq_along(nameAll), function(x) {
  read.table(paste0(inDirAll[x], "filt_log2_", prefixAll[x],
                    "_genome_norm_coverage_Zscore_", winName, ".txt"))
}, mc.cores = length(nameAll))

# Feature frequency
inDirFeatures <- "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/genomeProfiles/genomewideZscore/"

# TEs
filt_TEs <- read.table(paste0(inDirFeatures, "filt_TE_frequency_genome_",
                              winName, ".txt"))
# genes
filt_genes <- read.table(paste0(inDirFeatures, "filt_gene_frequency_genome_",
                                winName, ".txt"))
# COs
filt_COs <- read.table(paste0(inDirFeatures, "filt_CO_frequency_genome_",
                              winName, ".txt"))
# DNA methylation in 200-kb windows
filt_meth <- read.table(paste0(inDirFeatures, "filt_DNAmeth_GSM980986_WT_rep2_genome_",
                               winName, ".txt"))


# Create FILTERED, log2-transformed genome-wide Spearman's rho correlation matrix

REC8_HA <- filt_REC8log2trans_list[[1]]
REC8_MYC <- filt_REC8log2trans_list[[2]]
## Recombination:
SPO11_1_oligos <- filt_log2trans_list[[10]]
Crossovers <- filt_COs
## Euchromatin:
Genes <- filt_genes
H3K4me1 <- filt_log2trans_list[[5]]
H3K4me2 <- filt_log2trans_list[[6]]
H3K4me3 <- filt_log2trans_list[[7]]
H2AZ <- filt_log2trans_list[[8]]
H3K27me3 <- filt_log2trans_list[[9]]
## Heterochromatin:
Transposons <- filt_TEs
DNAmeth <- filt_meth
H2AW <- filt_log2trans_list[[1]]
H3K9me2 <- filt_log2trans_list[[2]]
H3K27me1 <- filt_log2trans_list[[3]]
Nucleosomes <- filt_log2trans_list[[4]]

allDF <- data.frame(REC8_HA[,3], REC8_MYC[,3], SPO11_1_oligos[,3], Crossovers[,3],
                    Genes[,3], H3K4me1[,3], H3K4me2[,3], H3K4me3[,3], H2AZ[,3], H3K27me3[,3],
                    Transposons[,3], DNAmeth[,3], DNAmeth[,4], DNAmeth[,5], H2AW[,3], H3K9me2[,3], H3K27me1[,3], Nucleosomes[,3])
colnames(allDF) <- c("REC8-HA", "REC8-Myc", "SPO11-1-oligos", "Crossovers",
                     "Genes", "H3K4me1", "H3K4me2", "H3K4me3", "H2A.Z", "H3K27me3",
                     "Transposons", "mCG", "mCHG", "mCHH", "H2A.W", "H3K9me2", "H3K27me1", "Nucleosomes")
allDF_corMat <- cor(allDF, method = "spearman", use = "pairwise.complete.obs")
col1 <- colorRampPalette(c("blue", "white", "red"))
pdf(file = paste0(plotDir, "REC8_HA_Rep2_filtered_genomewide_correlation_matrix_", winName, "_colouronly_genomewideZscore_v170619.pdf"),
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

REC8_MYC_armL <- NULL
REC8_MYC_armR <- NULL
REC8_MYC_arm <- NULL
REC8_MYC_peri <- NULL

SPO11_1_oligos_armL <- NULL
SPO11_1_oligos_armR <- NULL
SPO11_1_oligos_arm <- NULL
SPO11_1_oligos_peri <- NULL

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

Crossovers_armL <- NULL
Crossovers_armR <- NULL
Crossovers_arm <- NULL
Crossovers_peri <- NULL

DNAmeth_armL <- NULL
DNAmeth_armR <- NULL
DNAmeth_arm <- NULL
DNAmeth_peri <- NULL
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
  REC8_MYC_armL_tmp <- REC8_MYC[REC8_MYC[,2] > sumchr[i] & REC8_MYC[,2] < pericenStart[i] & REC8_MYC[,2] < sumchr[i+1],]
  REC8_MYC_armR_tmp <- REC8_MYC[REC8_MYC[,2] > pericenEnd[i] & REC8_MYC[,2] <= sumchr[i+1],]
  REC8_MYC_arm_tmp <- REC8_MYC[(REC8_MYC[,2] > sumchr[i] & REC8_MYC[,2] < pericenStart[i] & REC8_MYC[,2] < sumchr[i+1]) |
                                           (REC8_MYC[,2] > pericenEnd[i] & REC8_MYC[,2] <= sumchr[i+1]),]
  REC8_MYC_peri_tmp <- REC8_MYC[REC8_MYC[,2] >= pericenStart[i] & REC8_MYC[,2] <= pericenEnd[i],]
  REC8_MYC_peri <- rbind(REC8_MYC_peri, REC8_MYC_peri_tmp)
  REC8_MYC_armL <- rbind(REC8_MYC_armL, REC8_MYC_armL_tmp)
  REC8_MYC_armR <- rbind(REC8_MYC_armR, REC8_MYC_armR_tmp)
  REC8_MYC_arm <- rbind(REC8_MYC_arm, REC8_MYC_arm_tmp)

  print(i)
  SPO11_1_oligos_armL_tmp <- SPO11_1_oligos[SPO11_1_oligos[,2] > sumchr[i] & SPO11_1_oligos[,2] < pericenStart[i] & SPO11_1_oligos[,2] < sumchr[i+1],]
  SPO11_1_oligos_armR_tmp <- SPO11_1_oligos[SPO11_1_oligos[,2] > pericenEnd[i] & SPO11_1_oligos[,2] <= sumchr[i+1],]
  SPO11_1_oligos_arm_tmp <- SPO11_1_oligos[(SPO11_1_oligos[,2] > sumchr[i] & SPO11_1_oligos[,2] < pericenStart[i] & SPO11_1_oligos[,2] < sumchr[i+1]) |
                                           (SPO11_1_oligos[,2] > pericenEnd[i] & SPO11_1_oligos[,2] <= sumchr[i+1]),]
  SPO11_1_oligos_peri_tmp <- SPO11_1_oligos[SPO11_1_oligos[,2] >= pericenStart[i] & SPO11_1_oligos[,2] <= pericenEnd[i],]
  SPO11_1_oligos_peri <- rbind(SPO11_1_oligos_peri, SPO11_1_oligos_peri_tmp)
  SPO11_1_oligos_armL <- rbind(SPO11_1_oligos_armL, SPO11_1_oligos_armL_tmp)
  SPO11_1_oligos_armR <- rbind(SPO11_1_oligos_armR, SPO11_1_oligos_armR_tmp)
  SPO11_1_oligos_arm <- rbind(SPO11_1_oligos_arm, SPO11_1_oligos_arm_tmp)

  print(i)
  H3K4me1_armL_tmp <- H3K4me1[H3K4me1[,2] > sumchr[i] & H3K4me1[,2] < pericenStart[i] & H3K4me1[,2] < sumchr[i+1],]
  H3K4me1_armR_tmp <- H3K4me1[H3K4me1[,2] > pericenEnd[i] & H3K4me1[,2] <= sumchr[i+1],]
  H3K4me1_arm_tmp <- H3K4me1[(H3K4me1[,2] > sumchr[i] & H3K4me1[,2] < pericenStart[i] & H3K4me1[,2] < sumchr[i+1]) |
                                           (H3K4me1[,2] > pericenEnd[i] & H3K4me1[,2] <= sumchr[i+1]),]
  H3K4me1_peri_tmp <- H3K4me1[H3K4me1[,2] >= pericenStart[i] & H3K4me1[,2] <= pericenEnd[i],]
  H3K4me1_peri <- rbind(H3K4me1_peri, H3K4me1_peri_tmp)
  H3K4me1_armL <- rbind(H3K4me1_armL, H3K4me1_armL_tmp)
  H3K4me1_armR <- rbind(H3K4me1_armR, H3K4me1_armR_tmp)
  H3K4me1_arm <- rbind(H3K4me1_arm, H3K4me1_arm_tmp)

  print(i)
  H3K4me2_armL_tmp <- H3K4me2[H3K4me2[,2] > sumchr[i] & H3K4me2[,2] < pericenStart[i] & H3K4me2[,2] < sumchr[i+1],]
  H3K4me2_armR_tmp <- H3K4me2[H3K4me2[,2] > pericenEnd[i] & H3K4me2[,2] <= sumchr[i+1],]
  H3K4me2_arm_tmp <- H3K4me2[(H3K4me2[,2] > sumchr[i] & H3K4me2[,2] < pericenStart[i] & H3K4me2[,2] < sumchr[i+1]) |
                                           (H3K4me2[,2] > pericenEnd[i] & H3K4me2[,2] <= sumchr[i+1]),]
  H3K4me2_peri_tmp <- H3K4me2[H3K4me2[,2] >= pericenStart[i] & H3K4me2[,2] <= pericenEnd[i],]
  H3K4me2_peri <- rbind(H3K4me2_peri, H3K4me2_peri_tmp)
  H3K4me2_armL <- rbind(H3K4me2_armL, H3K4me2_armL_tmp)
  H3K4me2_armR <- rbind(H3K4me2_armR, H3K4me2_armR_tmp)
  H3K4me2_arm <- rbind(H3K4me2_arm, H3K4me2_arm_tmp)

  print(i)
  H3K4me3_armL_tmp <- H3K4me3[H3K4me3[,2] > sumchr[i] & H3K4me3[,2] < pericenStart[i] & H3K4me3[,2] < sumchr[i+1],]
  H3K4me3_armR_tmp <- H3K4me3[H3K4me3[,2] > pericenEnd[i] & H3K4me3[,2] <= sumchr[i+1],]
  H3K4me3_arm_tmp <- H3K4me3[(H3K4me3[,2] > sumchr[i] & H3K4me3[,2] < pericenStart[i] & H3K4me3[,2] < sumchr[i+1]) |
                                           (H3K4me3[,2] > pericenEnd[i] & H3K4me3[,2] <= sumchr[i+1]),]
  H3K4me3_peri_tmp <- H3K4me3[H3K4me3[,2] >= pericenStart[i] & H3K4me3[,2] <= pericenEnd[i],]
  H3K4me3_peri <- rbind(H3K4me3_peri, H3K4me3_peri_tmp)
  H3K4me3_armL <- rbind(H3K4me3_armL, H3K4me3_armL_tmp)
  H3K4me3_armR <- rbind(H3K4me3_armR, H3K4me3_armR_tmp)
  H3K4me3_arm <- rbind(H3K4me3_arm, H3K4me3_arm_tmp)

  print(i)
  H2AZ_armL_tmp <- H2AZ[H2AZ[,2] > sumchr[i] & H2AZ[,2] < pericenStart[i] & H2AZ[,2] < sumchr[i+1],]
  H2AZ_armR_tmp <- H2AZ[H2AZ[,2] > pericenEnd[i] & H2AZ[,2] <= sumchr[i+1],]
  H2AZ_arm_tmp <- H2AZ[(H2AZ[,2] > sumchr[i] & H2AZ[,2] < pericenStart[i] & H2AZ[,2] < sumchr[i+1]) |
                                           (H2AZ[,2] > pericenEnd[i] & H2AZ[,2] <= sumchr[i+1]),]
  H2AZ_peri_tmp <- H2AZ[H2AZ[,2] >= pericenStart[i] & H2AZ[,2] <= pericenEnd[i],]
  H2AZ_peri <- rbind(H2AZ_peri, H2AZ_peri_tmp)
  H2AZ_armL <- rbind(H2AZ_armL, H2AZ_armL_tmp)
  H2AZ_armR <- rbind(H2AZ_armR, H2AZ_armR_tmp)
  H2AZ_arm <- rbind(H2AZ_arm, H2AZ_arm_tmp)

  print(i)
  H3K27me3_armL_tmp <- H3K27me3[H3K27me3[,2] > sumchr[i] & H3K27me3[,2] < pericenStart[i] & H3K27me3[,2] < sumchr[i+1],]
  H3K27me3_armR_tmp <- H3K27me3[H3K27me3[,2] > pericenEnd[i] & H3K27me3[,2] <= sumchr[i+1],]
  H3K27me3_arm_tmp <- H3K27me3[(H3K27me3[,2] > sumchr[i] & H3K27me3[,2] < pericenStart[i] & H3K27me3[,2] < sumchr[i+1]) |
                                           (H3K27me3[,2] > pericenEnd[i] & H3K27me3[,2] <= sumchr[i+1]),]
  H3K27me3_peri_tmp <- H3K27me3[H3K27me3[,2] >= pericenStart[i] & H3K27me3[,2] <= pericenEnd[i],]
  H3K27me3_peri <- rbind(H3K27me3_peri, H3K27me3_peri_tmp)
  H3K27me3_armL <- rbind(H3K27me3_armL, H3K27me3_armL_tmp)
  H3K27me3_armR <- rbind(H3K27me3_armR, H3K27me3_armR_tmp)
  H3K27me3_arm <- rbind(H3K27me3_arm, H3K27me3_arm_tmp)

  print(i)
  H2AW_armL_tmp <- H2AW[H2AW[,2] > sumchr[i] & H2AW[,2] < pericenStart[i] & H2AW[,2] < sumchr[i+1],]
  H2AW_armR_tmp <- H2AW[H2AW[,2] > pericenEnd[i] & H2AW[,2] <= sumchr[i+1],]
  H2AW_arm_tmp <- H2AW[(H2AW[,2] > sumchr[i] & H2AW[,2] < pericenStart[i] & H2AW[,2] < sumchr[i+1]) |
                                           (H2AW[,2] > pericenEnd[i] & H2AW[,2] <= sumchr[i+1]),]
  H2AW_peri_tmp <- H2AW[H2AW[,2] >= pericenStart[i] & H2AW[,2] <= pericenEnd[i],]
  H2AW_peri <- rbind(H2AW_peri, H2AW_peri_tmp)
  H2AW_armL <- rbind(H2AW_armL, H2AW_armL_tmp)
  H2AW_armR <- rbind(H2AW_armR, H2AW_armR_tmp)
  H2AW_arm <- rbind(H2AW_arm, H2AW_arm_tmp)

  print(i)
  H3K9me2_armL_tmp <- H3K9me2[H3K9me2[,2] > sumchr[i] & H3K9me2[,2] < pericenStart[i] & H3K9me2[,2] < sumchr[i+1],]
  H3K9me2_armR_tmp <- H3K9me2[H3K9me2[,2] > pericenEnd[i] & H3K9me2[,2] <= sumchr[i+1],]
  H3K9me2_arm_tmp <- H3K9me2[(H3K9me2[,2] > sumchr[i] & H3K9me2[,2] < pericenStart[i] & H3K9me2[,2] < sumchr[i+1]) |
                                           (H3K9me2[,2] > pericenEnd[i] & H3K9me2[,2] <= sumchr[i+1]),]
  H3K9me2_peri_tmp <- H3K9me2[H3K9me2[,2] >= pericenStart[i] & H3K9me2[,2] <= pericenEnd[i],]
  H3K9me2_peri <- rbind(H3K9me2_peri, H3K9me2_peri_tmp)
  H3K9me2_armL <- rbind(H3K9me2_armL, H3K9me2_armL_tmp)
  H3K9me2_armR <- rbind(H3K9me2_armR, H3K9me2_armR_tmp)
  H3K9me2_arm <- rbind(H3K9me2_arm, H3K9me2_arm_tmp)

  print(i)
  H3K27me1_armL_tmp <- H3K27me1[H3K27me1[,2] > sumchr[i] & H3K27me1[,2] < pericenStart[i] & H3K27me1[,2] < sumchr[i+1],]
  H3K27me1_armR_tmp <- H3K27me1[H3K27me1[,2] > pericenEnd[i] & H3K27me1[,2] <= sumchr[i+1],]
  H3K27me1_arm_tmp <- H3K27me1[(H3K27me1[,2] > sumchr[i] & H3K27me1[,2] < pericenStart[i] & H3K27me1[,2] < sumchr[i+1]) |
                                           (H3K27me1[,2] > pericenEnd[i] & H3K27me1[,2] <= sumchr[i+1]),]
  H3K27me1_peri_tmp <- H3K27me1[H3K27me1[,2] >= pericenStart[i] & H3K27me1[,2] <= pericenEnd[i],]
  H3K27me1_peri <- rbind(H3K27me1_peri, H3K27me1_peri_tmp)
  H3K27me1_armL <- rbind(H3K27me1_armL, H3K27me1_armL_tmp)
  H3K27me1_armR <- rbind(H3K27me1_armR, H3K27me1_armR_tmp)
  H3K27me1_arm <- rbind(H3K27me1_arm, H3K27me1_arm_tmp)

  print(i)
  Nucleosomes_armL_tmp <- Nucleosomes[Nucleosomes[,2] > sumchr[i] & Nucleosomes[,2] < pericenStart[i] & Nucleosomes[,2] < sumchr[i+1],]
  Nucleosomes_armR_tmp <- Nucleosomes[Nucleosomes[,2] > pericenEnd[i] & Nucleosomes[,2] <= sumchr[i+1],]
  Nucleosomes_arm_tmp <- Nucleosomes[(Nucleosomes[,2] > sumchr[i] & Nucleosomes[,2] < pericenStart[i] & Nucleosomes[,2] < sumchr[i+1]) |
                                           (Nucleosomes[,2] > pericenEnd[i] & Nucleosomes[,2] <= sumchr[i+1]),]
  Nucleosomes_peri_tmp <- Nucleosomes[Nucleosomes[,2] >= pericenStart[i] & Nucleosomes[,2] <= pericenEnd[i],]
  Nucleosomes_peri <- rbind(Nucleosomes_peri, Nucleosomes_peri_tmp)
  Nucleosomes_armL <- rbind(Nucleosomes_armL, Nucleosomes_armL_tmp)
  Nucleosomes_armR <- rbind(Nucleosomes_armR, Nucleosomes_armR_tmp)
  Nucleosomes_arm <- rbind(Nucleosomes_arm, Nucleosomes_arm_tmp)

  print(i)
  Genes_armL_tmp <- Genes[Genes[,2] > sumchr[i] & Genes[,2] < pericenStart[i] & Genes[,2] < sumchr[i+1],]
  Genes_armR_tmp <- Genes[Genes[,2] > pericenEnd[i] & Genes[,2] <= sumchr[i+1],]
  Genes_arm_tmp <- Genes[(Genes[,2] > sumchr[i] & Genes[,2] < pericenStart[i] & Genes[,2] < sumchr[i+1]) |
                                           (Genes[,2] > pericenEnd[i] & Genes[,2] <= sumchr[i+1]),]
  Genes_peri_tmp <- Genes[Genes[,2] >= pericenStart[i] & Genes[,2] <= pericenEnd[i],]
  Genes_peri <- rbind(Genes_peri, Genes_peri_tmp)
  Genes_armL <- rbind(Genes_armL, Genes_armL_tmp)
  Genes_armR <- rbind(Genes_armR, Genes_armR_tmp)
  Genes_arm <- rbind(Genes_arm, Genes_arm_tmp)

  print(i)
  Transposons_armL_tmp <- Transposons[Transposons[,2] > sumchr[i] & Transposons[,2] < pericenStart[i] & Transposons[,2] < sumchr[i+1],]
  Transposons_armR_tmp <- Transposons[Transposons[,2] > pericenEnd[i] & Transposons[,2] <= sumchr[i+1],]
  Transposons_arm_tmp <- Transposons[(Transposons[,2] > sumchr[i] & Transposons[,2] < pericenStart[i] & Transposons[,2] < sumchr[i+1]) |
                                           (Transposons[,2] > pericenEnd[i] & Transposons[,2] <= sumchr[i+1]),]
  Transposons_peri_tmp <- Transposons[Transposons[,2] >= pericenStart[i] & Transposons[,2] <= pericenEnd[i],]
  Transposons_peri <- rbind(Transposons_peri, Transposons_peri_tmp)
  Transposons_armL <- rbind(Transposons_armL, Transposons_armL_tmp)
  Transposons_armR <- rbind(Transposons_armR, Transposons_armR_tmp)
  Transposons_arm <- rbind(Transposons_arm, Transposons_arm_tmp)

  print(i)
  Crossovers_armL_tmp <- Crossovers[Crossovers[,2] > sumchr[i] & Crossovers[,2] < pericenStart[i] & Crossovers[,2] < sumchr[i+1],]
  Crossovers_armR_tmp <- Crossovers[Crossovers[,2] > pericenEnd[i] & Crossovers[,2] <= sumchr[i+1],]
  Crossovers_arm_tmp <- Crossovers[(Crossovers[,2] > sumchr[i] & Crossovers[,2] < pericenStart[i] & Crossovers[,2] < sumchr[i+1]) |
                                           (Crossovers[,2] > pericenEnd[i] & Crossovers[,2] <= sumchr[i+1]),]
  Crossovers_peri_tmp <- Crossovers[Crossovers[,2] >= pericenStart[i] & Crossovers[,2] <= pericenEnd[i],]
  Crossovers_peri <- rbind(Crossovers_peri, Crossovers_peri_tmp)
  Crossovers_armL <- rbind(Crossovers_armL, Crossovers_armL_tmp)
  Crossovers_armR <- rbind(Crossovers_armR, Crossovers_armR_tmp)
  Crossovers_arm <- rbind(Crossovers_arm, Crossovers_arm_tmp)

  print(i)
  DNAmeth_armL_tmp <- DNAmeth[DNAmeth[,2] > sumchr[i] & DNAmeth[,2] < pericenStart[i] & DNAmeth[,2] < sumchr[i+1],]
  DNAmeth_armR_tmp <- DNAmeth[DNAmeth[,2] > pericenEnd[i] & DNAmeth[,2] <= sumchr[i+1],]
  DNAmeth_arm_tmp <- DNAmeth[(DNAmeth[,2] > sumchr[i] & DNAmeth[,2] < pericenStart[i] & DNAmeth[,2] < sumchr[i+1]) |
                                           (DNAmeth[,2] > pericenEnd[i] & DNAmeth[,2] <= sumchr[i+1]),]
  DNAmeth_peri_tmp <- DNAmeth[DNAmeth[,2] >= pericenStart[i] & DNAmeth[,2] <= pericenEnd[i],]
  DNAmeth_peri <- rbind(DNAmeth_peri, DNAmeth_peri_tmp)
  DNAmeth_armL <- rbind(DNAmeth_armL, DNAmeth_armL_tmp)
  DNAmeth_armR <- rbind(DNAmeth_armR, DNAmeth_armR_tmp)
  DNAmeth_arm <- rbind(DNAmeth_arm, DNAmeth_arm_tmp)
}

# Create arm and peri Spearman's rho correlation matrices
allDF_arm <- data.frame(REC8_HA_arm[,3], REC8_MYC_arm[,3], SPO11_1_oligos_arm[,3], Crossovers_arm[,3],
                        Genes_arm[,3], H3K4me1_arm[,3], H3K4me2_arm[,3], H3K4me3_arm[,3], H2AZ_arm[,3], H3K27me3_arm[,3],
                        Transposons_arm[,3], DNAmeth_arm[,3], DNAmeth_arm[,4], DNAmeth_arm[,5], H2AW_arm[,3], H3K9me2_arm[,3], H3K27me1_arm[,3], Nucleosomes_arm[,3])
colnames(allDF_arm) <- c("REC8-HA", "REC8-Myc", "SPO11-1-oligos", "Crossovers",
                         "Genes", "H3K4me1", "H3K4me2", "H3K4me3", "H2A.Z", "H3K27me3",
                         "Transposons", "mCG", "mCHG", "mCHH", "H2A.W", "H3K9me2", "H3K27me1", "Nucleosomes")
allDF_arm_corMat <- cor(allDF_arm, method = "spearman", use = "pairwise.complete.obs")
col1 <- colorRampPalette(c("blue", "white", "red"))
pdf(file = paste0(plotDir, "REC8_HA_Rep2_filtered_arms_correlation_matrix_", winName, "_colouronly_genomewideZscore_v170619.pdf"),
                  height = 8.5, width = 8)
corrplot(allDF_arm_corMat, method = "color", type = "upper", col = col1(20), tl.col = "black",
         addgrid.col = "white", addCoef.col = "grey90", mar = c(0,0,1,0), tl.cex = 1.5, cl.cex = 1.5, number.cex = 1.0,
         title = paste0("Filtered arms Spearman correlation matrix (10-kb windows)"))
dev.off()

allDF_peri <- data.frame(REC8_HA_peri[,3], REC8_MYC_peri[,3], SPO11_1_oligos_peri[,3], Crossovers_peri[,3],
                         Genes_peri[,3], H3K4me1_peri[,3], H3K4me2_peri[,3], H3K4me3_peri[,3], H2AZ_peri[,3], H3K27me3_peri[,3],
                         Transposons_peri[,3], DNAmeth_peri[,3], DNAmeth_peri[,4], DNAmeth_peri[,5], H2AW_peri[,3], H3K9me2_peri[,3], H3K27me1_peri[,3], Nucleosomes_peri[,3])
colnames(allDF_peri) <- c("REC8-HA", "REC8-Myc", "SPO11-1-oligos", "Crossovers",
                          "Genes", "H3K4me1", "H3K4me2", "H3K4me3", "H2A.Z", "H3K27me3",
                          "Transposons", "mCG", "mCHG", "mCHH", "H2A.W", "H3K9me2", "H3K27me1", "Nucleosomes")
allDF_peri_corMat <- cor(allDF_peri, method = "spearman", use = "pairwise.complete.obs")
col1 <- colorRampPalette(c("blue", "white", "red"))
pdf(file = paste0(plotDir, "REC8_HA_Rep2_filtered_pericentromeres_correlation_matrix_", winName, "_colouronly_genomewideZscore_v170619.pdf"),
                  height = 8.5, width = 8)
corrplot(allDF_peri_corMat, method = "color", type = "upper", col = col1(20), tl.col = "black",
         addgrid.col = "white", addCoef.col = "grey90", mar = c(0,0,1,0), tl.cex = 1.5, cl.cex = 1.5, number.cex = 1.0,
         title = paste0("Filtered pericentromeres Spearman correlation matrix (10-kb windows)"))      
dev.off()


pdf(file = paste0(plotDir, "REC8_HA_Rep2_filtered_genomewide_arms_pericentromeres_correlation_matrix_", winName, "_colouronly_genomewideZscore_v170619.pdf"),
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
