#!/applications/R/R-3.3.2/bin/Rscript

# Plot bar chart of log2(observed:expected) peaks overlapping TEs within each family

# Usage:
# Rscript TE_family_wt_vs_kss_REC8_HA_Rep2_SPO11oligo_H3K9me2_genome_wide_peaks.R "REC8 peaks, SPO11-1-oligo hotspots and H3K9me2 peaks" "wt REC8-HA Rep2" "kss REC8-HA Rep1" "wt SPO11-1-oligos" "kss SPO11-1-oligos" "wt H3K9me2" "kss H3K9me2" wt_REC8_HA_Rep2 kss_REC8_HA_Rep1 wt_SPO11oligo kss_SPO11oligo wt_H3K9me2 kss_H3K9me2 10000

library(ggplot2)
library(ggthemes)

dataName <- "REC8 peaks, SPO11-1-oligo hotspots and H3K9me2 peaks"
pt1REC8Name <- "wt REC8-HA Rep2"
pt2REC8Name <- "kss REC8-HA Rep1"
pt1Name <- "wt SPO11-1-oligos" 
pt2Name <-  "kss SPO11-1-oligos" 
pt3Name <- "wt H3K9me2" 
pt4Name <-  "kss H3K9me2" 
pt1REC8LibName <- "wt_REC8_HA_Rep2"
pt2REC8LibName <- "kss_REC8_HA_Rep1"
pt1LibName <- "wt_SPO11oligo"
pt2LibName <- "kss_SPO11oligo"
pt3LibName <- "wt_H3K9me2"
pt4LibName <- "kss_H3K9me2"
# Number of permutations (randomisations) performed
perms <- 10000

args <- commandArgs(trailingOnly = T)
dataName <- as.character(args[1])
pt1REC8Name <- as.character(args[2])
pt2REC8Name <- as.character(args[3])
pt1Name <- as.character(args[4])
pt2Name <- as.character(args[5])
pt3Name <- as.character(args[6])
pt4Name <- as.character(args[7])
pt1REC8LibName <- as.character(args[8])
pt2REC8LibName <- as.character(args[9])
pt1LibName <- as.character(args[10])
pt2LibName <- as.character(args[11])
pt3LibName <- as.character(args[12])
pt4LibName <- as.character(args[13])
# Number of permutations (randomisations) performed
perms <- as.numeric(args[14])

REC8Dir <- "/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/genome_wide/regioneR/noMinWidth_mergedOverlaps/"
REC8Dir1 <- "/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep2/genome_wide/regioneR/noMinWidth_mergedOverlaps/"
REC8Dir2 <- "/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/kss_REC8_HA_Rep1/genome_wide/regioneR/noMinWidth_mergedOverlaps/"
inDir1 <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/RPI1_RPI8_profiles_rangerPeaks_idr0.05/genome_wide/regioneR/"
inDir2 <- "/projects/ajt200/BAM_masters/SPO11-oligo/suvh456/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/RPI34_RPI35_profiles_rangerPeaks_idr0.05/genome_wide/regioneR/"
inDir3 <- "/projects/ajt200/BAM_masters/H3K9me2/WT/peaks/PeakRanger1.18/ranger/p0.05_q0.05/genome_wide/regioneR/noMinWidth_mergedOverlaps/"
inDir4 <- "/projects/ajt200/BAM_masters/H3K9me2/kss/peaks/PeakRanger1.18/ranger/p0.05_q0.05/genome_wide/regioneR/noMinWidth_mergedOverlaps/"

plotDir <- "/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep2/genome_wide/regioneR/noMinWidth_mergedOverlaps/plots/"

famNames <- c("dna", "heli", "ptmari", "mudr", "enspm", "hat", "harbinger",
              "rna", "gypsy", "copia", "linel1", "sine")
famNamesPlot <- c("DNA", "Helitrons", "Pogo/Tc1/Mariner", "MuDR", "EnSpm", "hAT", "Harbinger",
                  "RNA", "Gypsy LTR", "Copia LTR", "LINE-1", "SINE")

load(paste0(REC8Dir, "permTest_REC8_HA_Rep1_rangerPeaks_vs_TEsDNA.RData"))
pt_REC8_DNA <- ptPeaksTEsDNAPerChrom
ptPeaksTEsDNAPerChrom <- NULL
load(paste0(REC8Dir, "permTest_REC8_HA_Rep1_rangerPeaks_vs_TEsRNA.RData"))
pt_REC8_RNA <- ptPeaksTEsRNAPerChrom
ptPeaksTEsRNAPerChrom <- NULL
pt_REC8 <- c(pt_REC8_DNA, pt_REC8_RNA)

load(paste0(REC8Dir1, "permTest_REC8_HA_Rep2_rangerPeaks_vs_TEsDNA.RData"))
pt1_REC8_DNA <- ptPeaksTEsDNAPerChrom
ptPeaksTEsDNAPerChrom <- NULL
load(paste0(REC8Dir1, "permTest_REC8_HA_Rep2_rangerPeaks_vs_TEsRNA.RData"))
pt1_REC8_RNA <- ptPeaksTEsRNAPerChrom
ptPeaksTEsRNAPerChrom <- NULL
pt1_REC8 <- c(pt1_REC8_DNA, pt1_REC8_RNA)

load(paste0(REC8Dir2, "permTest_kss_REC8_HA_Rep1_rangerPeaks_vs_TEsDNA.RData"))
pt2_REC8_DNA <- ptPeaksTEsDNAPerChrom
ptPeaksTEsDNAPerChrom <- NULL
load(paste0(REC8Dir2, "permTest_kss_REC8_HA_Rep1_rangerPeaks_vs_TEsRNA.RData"))
pt2_REC8_RNA <- ptPeaksTEsRNAPerChrom
ptPeaksTEsRNAPerChrom <- NULL
pt2_REC8 <- c(pt2_REC8_DNA, pt2_REC8_RNA)

load(paste0(inDir1, "permTest_SPO11_RPI1_RPI8_rangerPeaks_vs_TEsDNA.RData"))
pt1_DNA <- ptPeaksTEsDNAPerChrom
ptPeaksTEsDNAPerChrom <- NULL
load(paste0(inDir1, "permTest_SPO11_RPI1_RPI8_rangerPeaks_vs_TEsRNA.RData"))
pt1_RNA <- ptPeaksTEsRNAPerChrom
ptPeaksTEsRNAPerChrom <- NULL
pt1 <- c(pt1_DNA, pt1_RNA)

load(paste0(inDir2, "permTest_kss_SPO11_RPI34_RPI35_rangerPeaks_vs_TEsDNA.RData"))
pt2_DNA <- ptPeaksTEsDNAPerChrom
ptPeaksTEsDNAPerChrom <- NULL
load(paste0(inDir2, "permTest_kss_SPO11_RPI34_RPI35_rangerPeaks_vs_TEsRNA.RData"))
pt2_RNA <- ptPeaksTEsRNAPerChrom
ptPeaksTEsRNAPerChrom <- NULL
pt2 <- c(pt2_DNA, pt2_RNA)

load(paste0(inDir3, "permTest_H3K9me2_rangerPeaks_vs_TEsDNA.RData"))
pt3_DNA <- ptPeaksTEsDNAPerChrom
ptPeaksTEsDNAPerChrom <- NULL
load(paste0(inDir3, "permTest_H3K9me2_rangerPeaks_vs_TEsRNA.RData"))
pt3_RNA <- ptPeaksTEsRNAPerChrom
ptPeaksTEsRNAPerChrom <- NULL
pt3 <- c(pt3_DNA, pt3_RNA)

load(paste0(inDir4, "permTest_kss_H3K9me2_rangerPeaks_vs_TEsDNA.RData"))
pt4_DNA <- ptPeaksTEsDNAPerChrom
ptPeaksTEsDNAPerChrom <- NULL
load(paste0(inDir4, "permTest_kss_H3K9me2_rangerPeaks_vs_TEsRNA.RData"))
pt4_RNA <- ptPeaksTEsRNAPerChrom
ptPeaksTEsRNAPerChrom <- NULL
pt4 <- c(pt4_DNA, pt4_RNA)

# pt_REC8
pt_REC8_Pval <- lapply(seq_along(pt_REC8), function(x) {
  pt_REC8[[x]]$numOverlaps$pval
})
pt_REC8_Obs <- lapply(seq_along(pt_REC8), function(x) {
  pt_REC8[[x]]$numOverlaps$observed
})
pt_REC8_Perm <- lapply(seq_along(pt_REC8), function(x) {
  pt_REC8[[x]]$numOverlaps$permuted
})
pt_REC8_Exp <- lapply(seq_along(pt_REC8), function(x) {
  mean(pt_REC8[[x]]$numOverlaps$permuted)
})
pt_REC8_log2ObsExp <- lapply(seq_along(pt_REC8_Obs), function(x) {
  log2(pt_REC8_Obs[[x]]/pt_REC8_Exp[[x]])
})
pt_REC8_Zscore <- lapply(seq_along(pt_REC8), function(x) {
  pt_REC8[[x]]$numOverlaps$zscore
})
pt_REC8_AltHyp <- lapply(seq_along(pt_REC8), function(x) {
  pt_REC8[[x]]$numOverlaps$alternative
})
pt_REC8_alpha0.05 <- lapply(seq_along(pt_REC8_Perm), function(x) {
  if(pt_REC8_AltHyp[[x]] == "greater") {
    quantile(pt_REC8_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(pt_REC8_Perm[[x]], 0.05)[[1]]
  }
})
pt_REC8_log2alpha0.05 <- lapply(seq_along(pt_REC8_alpha0.05), function(x) {
  log2(pt_REC8_alpha0.05[[x]]/pt_REC8_Exp[[x]])
})

pt_REC8_log2ObsExp_sorted <- sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T)
pt_REC8_log2alpha0.05_sorted <- unlist(pt_REC8_log2alpha0.05[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix])
pt_REC8_famNames_sorted <- famNames[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix]
pt_REC8_famNamesPlot_sorted <- famNamesPlot[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix]

# pt1_REC8
pt1_REC8_Pval <- lapply(seq_along(pt1_REC8), function(x) {
  pt1_REC8[[x]]$numOverlaps$pval
})
pt1_REC8_Obs <- lapply(seq_along(pt1_REC8), function(x) {
  pt1_REC8[[x]]$numOverlaps$observed
})
pt1_REC8_Perm <- lapply(seq_along(pt1_REC8), function(x) {
  pt1_REC8[[x]]$numOverlaps$permuted
})
pt1_REC8_Exp <- lapply(seq_along(pt1_REC8), function(x) {
  mean(pt1_REC8[[x]]$numOverlaps$permuted)
})
pt1_REC8_log2ObsExp <- lapply(seq_along(pt1_REC8_Obs), function(x) {
  log2(pt1_REC8_Obs[[x]]/pt1_REC8_Exp[[x]])
})
pt1_REC8_Zscore <- lapply(seq_along(pt1_REC8), function(x) {
  pt1_REC8[[x]]$numOverlaps$zscore
})
pt1_REC8_AltHyp <- lapply(seq_along(pt1_REC8), function(x) {
  pt1_REC8[[x]]$numOverlaps$alternative
})
pt1_REC8_alpha0.05 <- lapply(seq_along(pt1_REC8_Perm), function(x) {
  if(pt1_REC8_AltHyp[[x]] == "greater") {
    quantile(pt1_REC8_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(pt1_REC8_Perm[[x]], 0.05)[[1]]
  }
})
pt1_REC8_log2alpha0.05 <- lapply(seq_along(pt1_REC8_alpha0.05), function(x) {
  log2(pt1_REC8_alpha0.05[[x]]/pt1_REC8_Exp[[x]])
})

pt1_REC8_log2ObsExp_sorted <- unlist(pt1_REC8_log2ObsExp[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix])
pt1_REC8_log2alpha0.05_sorted <- unlist(pt1_REC8_log2alpha0.05[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix])
pt1_REC8_famNames_sorted <- famNames[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix]
pt1_REC8_famNamesPlot_sorted <- famNamesPlot[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix]

# pt2_REC8
pt2_REC8_Pval <- lapply(seq_along(pt2_REC8), function(x) {
  pt2_REC8[[x]]$numOverlaps$pval
})
pt2_REC8_Obs <- lapply(seq_along(pt2_REC8), function(x) {
  pt2_REC8[[x]]$numOverlaps$observed
})
pt2_REC8_Perm <- lapply(seq_along(pt2_REC8), function(x) {
  pt2_REC8[[x]]$numOverlaps$permuted
})
pt2_REC8_Exp <- lapply(seq_along(pt2_REC8), function(x) {
  mean(pt2_REC8[[x]]$numOverlaps$permuted)
})
pt2_REC8_log2ObsExp <- lapply(seq_along(pt2_REC8_Obs), function(x) {
  log2(pt2_REC8_Obs[[x]]/pt2_REC8_Exp[[x]])
})
pt2_REC8_Zscore <- lapply(seq_along(pt2_REC8), function(x) {
  pt2_REC8[[x]]$numOverlaps$zscore
})
pt2_REC8_AltHyp <- lapply(seq_along(pt2_REC8), function(x) {
  pt2_REC8[[x]]$numOverlaps$alternative
})
pt2_REC8_alpha0.05 <- lapply(seq_along(pt2_REC8_Perm), function(x) {
  if(pt2_REC8_AltHyp[[x]] == "greater") {
    quantile(pt2_REC8_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(pt2_REC8_Perm[[x]], 0.05)[[1]]
  }
})
pt2_REC8_log2alpha0.05 <- lapply(seq_along(pt2_REC8_alpha0.05), function(x) {
  log2(pt2_REC8_alpha0.05[[x]]/pt2_REC8_Exp[[x]])
})

pt2_REC8_log2ObsExp_sorted <- unlist(pt2_REC8_log2ObsExp[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix])
pt2_REC8_log2alpha0.05_sorted <- unlist(pt2_REC8_log2alpha0.05[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix])
pt2_REC8_famNames_sorted <- famNames[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix]
pt2_REC8_famNamesPlot_sorted <- famNamesPlot[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix]

# pt1
pt1_Pval <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$pval
})
pt1_Obs <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$observed
})
pt1_Perm <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$permuted
})
pt1_Exp <- lapply(seq_along(pt1), function(x) {
  mean(pt1[[x]]$numOverlaps$permuted)
})
pt1_log2ObsExp <- lapply(seq_along(pt1_Obs), function(x) {
  log2(pt1_Obs[[x]]/pt1_Exp[[x]])
})
pt1_Zscore <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$zscore
})
pt1_AltHyp <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$alternative
})
pt1_alpha0.05 <- lapply(seq_along(pt1_Perm), function(x) {
  if(pt1_AltHyp[[x]] == "greater") {
    quantile(pt1_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(pt1_Perm[[x]], 0.05)[[1]]
  }
})
pt1_log2alpha0.05 <- lapply(seq_along(pt1_alpha0.05), function(x) {
  log2(pt1_alpha0.05[[x]]/pt1_Exp[[x]])
})

pt1_log2ObsExp_sorted <- unlist(pt1_log2ObsExp[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix])
pt1_log2alpha0.05_sorted <- unlist(pt1_log2alpha0.05[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix])
pt1_famNames_sorted <- famNames[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix]
pt1_famNamesPlot_sorted <- famNamesPlot[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix]

# pt2
pt2_Pval <- lapply(seq_along(pt2), function(x) {
  pt2[[x]]$numOverlaps$pval
})
pt2_Obs <- lapply(seq_along(pt2), function(x) {
  pt2[[x]]$numOverlaps$observed
})
pt2_Perm <- lapply(seq_along(pt2), function(x) {
  pt2[[x]]$numOverlaps$permuted
})
pt2_Exp <- lapply(seq_along(pt2), function(x) {
  mean(pt2[[x]]$numOverlaps$permuted)
})
pt2_log2ObsExp <- lapply(seq_along(pt2_Obs), function(x) {
  log2(pt2_Obs[[x]]/pt2_Exp[[x]])
})
pt2_Zscore <- lapply(seq_along(pt2), function(x) {
  pt2[[x]]$numOverlaps$zscore
})
pt2_AltHyp <- lapply(seq_along(pt2), function(x) {
  pt2[[x]]$numOverlaps$alternative
})
pt2_alpha0.05 <- lapply(seq_along(pt2_Perm), function(x) {
  if(pt2_AltHyp[[x]] == "greater") {
    quantile(pt2_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(pt2_Perm[[x]], 0.05)[[1]]
  }
})
pt2_log2alpha0.05 <- lapply(seq_along(pt2_alpha0.05), function(x) {
  log2(pt2_alpha0.05[[x]]/pt2_Exp[[x]])
})

pt2_log2ObsExp_sorted <- unlist(pt2_log2ObsExp[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix])
pt2_log2alpha0.05_sorted <- unlist(pt2_log2alpha0.05[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix])
pt2_famNames_sorted <- famNames[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix]
pt2_famNamesPlot_sorted <- famNamesPlot[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix]

# pt3
pt3_Pval <- lapply(seq_along(pt3), function(x) {
  pt3[[x]]$numOverlaps$pval
})
pt3_Obs <- lapply(seq_along(pt3), function(x) {
  pt3[[x]]$numOverlaps$observed
})
pt3_Perm <- lapply(seq_along(pt3), function(x) {
  pt3[[x]]$numOverlaps$permuted
})
pt3_Exp <- lapply(seq_along(pt3), function(x) {
  mean(pt3[[x]]$numOverlaps$permuted)
})
pt3_log2ObsExp <- lapply(seq_along(pt3_Obs), function(x) {
  log2(pt3_Obs[[x]]/pt3_Exp[[x]])
})
pt3_Zscore <- lapply(seq_along(pt3), function(x) {
  pt3[[x]]$numOverlaps$zscore
})
pt3_AltHyp <- lapply(seq_along(pt3), function(x) {
  pt3[[x]]$numOverlaps$alternative
})
pt3_alpha0.05 <- lapply(seq_along(pt3_Perm), function(x) {
  if(pt3_AltHyp[[x]] == "greater") {
    quantile(pt3_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(pt3_Perm[[x]], 0.05)[[1]]
  }
})
pt3_log2alpha0.05 <- lapply(seq_along(pt3_alpha0.05), function(x) {
  log2(pt3_alpha0.05[[x]]/pt3_Exp[[x]])
})

pt3_log2ObsExp_sorted <- unlist(pt3_log2ObsExp[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix])
pt3_log2alpha0.05_sorted <- unlist(pt3_log2alpha0.05[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix])
pt3_famNames_sorted <- famNames[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix]
pt3_famNamesPlot_sorted <- famNamesPlot[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix]

# pt4
pt4_Pval <- lapply(seq_along(pt4), function(x) {
  pt4[[x]]$numOverlaps$pval
})
pt4_Obs <- lapply(seq_along(pt4), function(x) {
  pt4[[x]]$numOverlaps$observed
})
pt4_Perm <- lapply(seq_along(pt4), function(x) {
  pt4[[x]]$numOverlaps$permuted
})
pt4_Exp <- lapply(seq_along(pt4), function(x) {
  mean(pt4[[x]]$numOverlaps$permuted)
})
pt4_log2ObsExp <- lapply(seq_along(pt4_Obs), function(x) {
  log2(pt4_Obs[[x]]/pt4_Exp[[x]])
})
pt4_Zscore <- lapply(seq_along(pt4), function(x) {
  pt4[[x]]$numOverlaps$zscore
})
pt4_AltHyp <- lapply(seq_along(pt4), function(x) {
  pt4[[x]]$numOverlaps$alternative
})
pt4_alpha0.05 <- lapply(seq_along(pt4_Perm), function(x) {
  if(pt4_AltHyp[[x]] == "greater") {
    quantile(pt4_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(pt4_Perm[[x]], 0.05)[[1]]
  }
})
pt4_log2alpha0.05 <- lapply(seq_along(pt4_alpha0.05), function(x) {
  log2(pt4_alpha0.05[[x]]/pt4_Exp[[x]])
})

pt4_log2ObsExp_sorted <- unlist(pt4_log2ObsExp[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix])
pt4_log2alpha0.05_sorted <- unlist(pt4_log2alpha0.05[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix])
pt4_famNames_sorted <- famNames[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix]
pt4_famNamesPlot_sorted <- famNamesPlot[sort.int(unlist(pt_REC8_log2ObsExp), decreasing = T, index.return = T)$ix]

df <- data.frame(Sample = rep(c(pt1REC8Name, pt2REC8Name,
                                pt1Name, pt2Name,
                                pt3Name, pt4Name),
                              each = length(pt1_log2ObsExp_sorted)),
                 Transposon_family = rep(pt1_famNamesPlot_sorted, 2),
                 log2ObsExp = c(pt1_REC8_log2ObsExp_sorted,
                                pt2_REC8_log2ObsExp_sorted,
                                pt1_log2ObsExp_sorted,
                                pt2_log2ObsExp_sorted,
                                pt3_log2ObsExp_sorted,
                                pt4_log2ObsExp_sorted),
                 log2alpha0.05 = c(pt1_REC8_log2alpha0.05_sorted,
                                   pt2_REC8_log2alpha0.05_sorted,
                                   pt1_log2alpha0.05_sorted,
                                   pt2_log2alpha0.05_sorted,
                                   pt3_log2alpha0.05_sorted,
                                   pt4_log2alpha0.05_sorted))

df$Transposon_family <- factor(df$Transposon_family,
                               levels = pt1_famNamesPlot_sorted)
df$Sample <- factor(df$Sample,
                    levels = c(pt1REC8Name, pt2REC8Name,
                               pt1Name, pt2Name,
                               pt3Name, pt4Name))

bp <- ggplot(data = df,
             mapping = aes(x = Transposon_family,
                           y = log2ObsExp,
                           fill = Sample)) +
      geom_bar(stat = "identity",
               position = position_dodge()) +
      scale_fill_manual(name = "Genotype",
                        values = c("red",
                                   "red4",
                                   "blue",
                                   "navy",
                                   "green",
                                   "darkgreen"),
                        labels = c(pt1REC8Name,
                                   pt2REC8Name,
                                   pt1Name,
                                   pt2Name,
                                   pt3Name,
                                   pt4Name)) +
      geom_point(mapping = aes(Transposon_family, log2alpha0.05),
                 position = position_dodge(0.9),
                 shape = "-", colour  = "grey70", size = 3.5) +
      labs(x = "Transposon family",
           y = expression("Log"[2]*"(observed:expected) peak overlap")) +
      theme_bw() +
      theme(axis.line.y = element_line(size = 0.5, colour = "black"),
            axis.ticks.y = element_line(size = 0.25, colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.ticks.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 7),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      ggtitle(paste0(dataName, " (", as.character(perms), " permutations)"))
ggsave(paste0(plotDir, "barplot_TE_families_permTestResults_",
              as.character(perms), "perms_",
              "log2_Observed_Expected_",
               pt1REC8LibName, "_", pt2REC8LibName, "_peaks_",
               pt1LibName, "_", pt2LibName, "_peaks_",
               pt3LibName, "_", pt4LibName, "_peaks.pdf"),
       plot = bp,
       height = 4, width = 9)
save(bp,
     file = paste0(plotDir, "barplot_TE_families_permTestResults_",
                   as.character(perms), "perms_",
                  "log2_Observed_Expected_",
                   pt1REC8LibName, "_", pt2REC8LibName, "_peaks_",
                   pt1LibName, "_", pt2LibName, "_peaks_",
                   pt3LibName, "_", pt4LibName, "_peaks.RData"))

