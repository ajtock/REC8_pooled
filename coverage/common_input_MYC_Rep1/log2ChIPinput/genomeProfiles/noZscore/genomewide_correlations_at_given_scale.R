#!/applications/R/R-3.5.0/bin/Rscript

# Calculate genome-wide Spearman's rank-order correlations between libraries using
# unsmoothed windowed log2(ChIP/input) coverage values 

winNames <- c("2kb", "5kb", "10kb", "20kb", "50kb", "100kb")

setwd("/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/genomeProfiles/noZscore/")

ASY1_Rep1_list <- lapply(seq_along(winNames), function(x) {
  read.table(paste0("ASY1_Rep1_ChIP_REC8_MYC_Rep1_input_genome_norm_coverage_",
                    winNames[x], "_noZscore.tsv"),
             header = T)[,6]
})
ASY1_Rep2_list <- lapply(seq_along(winNames), function(x) {
  read.table(paste0("ASY1_Rep2_ChIP_REC8_MYC_Rep1_input_genome_norm_coverage_",
                    winNames[x], "_noZscore.tsv"),
             header = T)[,6]
})
REC8_HA_Rep1_list <- lapply(seq_along(winNames), function(x) {
  read.table(paste0("REC8_HA_Rep1_ChIP_REC8_MYC_Rep1_input_genome_norm_coverage_",
                    winNames[x], "_noZscore.tsv"),
             header = T)[,6]
})
REC8_HA_Rep2_list <- lapply(seq_along(winNames), function(x) {
  read.table(paste0("REC8_HA_Rep2_ChIP_REC8_MYC_Rep1_input_genome_norm_coverage_",
                    winNames[x], "_noZscore.tsv"),
             header = T)[,6]
})
REC8_MYC_Rep1_list <- lapply(seq_along(winNames), function(x) {
  read.table(paste0("REC8_MYC_Rep1_ChIP_REC8_MYC_Rep1_input_genome_norm_coverage_",
                    winNames[x], "_noZscore.tsv"),
             header = T)[,6]
})

ASY1_Rep1_v_ASY1_Rep2 <- sapply(seq_along(winNames), function(x) {
  cor(x = ASY1_Rep1_list[[x]],
      y = ASY1_Rep2_list[[x]],
      method = "spearman")
})
# ASY1_Rep1
ASY1_Rep1_v_REC8_HA_Rep1 <- sapply(seq_along(winNames), function(x) {
  cor(x = ASY1_Rep1_list[[x]],
      y = REC8_HA_Rep1_list[[x]],
      method = "spearman")
})
ASY1_Rep1_v_REC8_HA_Rep2 <- sapply(seq_along(winNames), function(x) {
  cor(x = ASY1_Rep1_list[[x]],
      y = REC8_HA_Rep2_list[[x]],
      method = "spearman")
})
ASY1_Rep1_v_REC8_MYC_Rep1 <- sapply(seq_along(winNames), function(x) {
  cor(x = ASY1_Rep1_list[[x]],
      y = REC8_MYC_Rep1_list[[x]],
      method = "spearman")
})
# ASY1_Rep2
ASY1_Rep2_v_REC8_HA_Rep1 <- sapply(seq_along(winNames), function(x) {
  cor(x = ASY1_Rep2_list[[x]],
      y = REC8_HA_Rep1_list[[x]],
      method = "spearman")
})
ASY1_Rep2_v_REC8_HA_Rep2 <- sapply(seq_along(winNames), function(x) {
  cor(x = ASY1_Rep2_list[[x]],
      y = REC8_HA_Rep2_list[[x]],
      method = "spearman")
})
ASY1_Rep2_v_REC8_MYC_Rep1 <- sapply(seq_along(winNames), function(x) {
  cor(x = ASY1_Rep2_list[[x]],
      y = REC8_MYC_Rep1_list[[x]],
      method = "spearman")
})
# REC8_HA_Rep1
REC8_HA_Rep1_v_REC8_HA_Rep2 <- sapply(seq_along(winNames), function(x) {
  cor(x = REC8_HA_Rep1_list[[x]],
      y = REC8_HA_Rep2_list[[x]],
      method = "spearman")
})
REC8_HA_Rep1_v_REC8_MYC_Rep1 <- sapply(seq_along(winNames), function(x) {
  cor(x = REC8_HA_Rep1_list[[x]],
      y = REC8_MYC_Rep1_list[[x]],
      method = "spearman")
})
# REC8_HA_Rep2
REC8_HA_Rep2_v_REC8_MYC_Rep1 <- sapply(seq_along(winNames), function(x) {
  cor(x = REC8_HA_Rep2_list[[x]],
      y = REC8_MYC_Rep1_list[[x]],
      method = "spearman")
})

corrDF <- data.frame(ASY1_Rep1_v_ASY1_Rep2 = round(ASY1_Rep1_v_ASY1_Rep2, digits = 2),
                     ASY1_Rep1_v_REC8_HA_Rep1 = round(ASY1_Rep1_v_REC8_HA_Rep1, digits = 2),
                     ASY1_Rep1_v_REC8_HA_Rep2 = round(ASY1_Rep1_v_REC8_HA_Rep2, digits = 2),
                     ASY1_Rep1_v_REC8_MYC_Rep1 = round(ASY1_Rep1_v_REC8_MYC_Rep1, digits = 2),
                     ASY1_Rep2_v_REC8_HA_Rep1 = round(ASY1_Rep2_v_REC8_HA_Rep1, digits = 2),
                     ASY1_Rep2_v_REC8_HA_Rep2 = round(ASY1_Rep2_v_REC8_HA_Rep2, digits = 2),
                     ASY1_Rep2_v_REC8_MYC_Rep1 = round(ASY1_Rep2_v_REC8_MYC_Rep1, digits = 2),
                     REC8_HA_Rep1_v_REC8_HA_Rep2 = round(REC8_HA_Rep1_v_REC8_HA_Rep2, digits = 2),
                     REC8_HA_Rep1_v_REC8_MYC_Rep1 = round(REC8_HA_Rep1_v_REC8_MYC_Rep1, digits = 2),
                     REC8_HA_Rep2_v_REC8_MYC_Rep1 = round(REC8_HA_Rep2_v_REC8_MYC_Rep1, digits = 2))

