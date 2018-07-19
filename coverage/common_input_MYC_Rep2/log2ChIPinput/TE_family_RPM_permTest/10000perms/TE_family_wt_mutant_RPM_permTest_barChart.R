#!/applications/R/R-3.3.2/bin/Rscript

# Plot bar chart of log2(observed:expected) reads in wt and mutant
# mapping to transposons within each transposon family

# Usage via Condor submission system on node7:
# csmit -m 1G -c 1 "Rscript TE_family_wt_mutant_RPM_permTest_barChart.R H3K9me2 "wild type" cmt3 wt_H3K9me2_ChIP cmt3_H3K9me2_ChIP 10000"

library(ggplot2)
library(ggthemes)

args <- commandArgs(trailingOnly = T)
dataName <- as.character(args[1])
wtName <- as.character(args[2])
mutantName <- as.character(args[3])
wtLibName <- as.character(args[4])
mutantLibName <- as.character(args[5])
# Number of permutations (randomisations) performed
perms <- as.numeric(args[6])

outDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/TE_family_RPM_permTest/10000perms/"
plotDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/TE_family_RPM_permTest/10000perms/plots/"

famNames <- c("heli", "ptmari", "mudr", "enspm", "hat", "harbinger",
              "gypsy", "copia", "linel1", "sine")
famNamesPlot <- c("Helitrons", "Pogo/Tc1/Mariner", "MuDR", "EnSpm", "hAT", "Harbinger",
                  "Gypsy LTR", "Copia LTR", "LINE-1", "SINE")

# wt
wtpermTestResultsList <- lapply(seq_along(famNames), function(x) {
  load(file = paste0(outDir,
                     famNames[x], "_TEs_permTestResults_",
                     as.character(perms), "perms_",
                     wtLibName, "_RPM.RData"))
  permTestResults <- permTestResults
})

wtlog2ObsExp <- sapply(wtpermTestResultsList, function(x) {
  as.vector(x@log2_Observed_Expected)
})
wtlog2alpha0.05 <- sapply(wtpermTestResultsList, function(x) {
  as.vector(log2((x@alpha0.05)/(x@Expected_RPM)))
})

wtlog2ObsExp_sorted <- sort.int(wtlog2ObsExp, decreasing = T)
wtlog2alpha0.05_sorted <- wtlog2alpha0.05[sort.int(wtlog2ObsExp, decreasing = T, index.return = T)$ix]
wtfamNames_sorted <- famNames[sort.int(wtlog2ObsExp, decreasing = T, index.return = T)$ix]
wtfamNamesPlot_sorted <- famNamesPlot[sort.int(wtlog2ObsExp, decreasing = T, index.return = T)$ix]

# mutant
mutantpermTestResultsList <- lapply(seq_along(famNames), function(x) {
  load(file = paste0(outDir,
                     famNames[x], "_TEs_permTestResults_",
                     as.character(perms), "perms_",
                     mutantLibName, "_RPM.RData"))
  permTestResults <- permTestResults
})

mutantlog2ObsExp <- sapply(mutantpermTestResultsList, function(x) {
  as.vector(x@log2_Observed_Expected)
})
mutantlog2alpha0.05 <- sapply(mutantpermTestResultsList, function(x) {
  as.vector(log2((x@alpha0.05)/(x@Expected_RPM)))
})

mutantlog2ObsExp_sorted <- mutantlog2ObsExp[sort.int(wtlog2ObsExp, decreasing = T, index.return = T)$ix]
mutantlog2alpha0.05_sorted <- mutantlog2alpha0.05[sort.int(wtlog2ObsExp, decreasing = T, index.return = T)$ix]
mutantfamNames_sorted <- famNames[sort.int(wtlog2ObsExp, decreasing = T, index.return = T)$ix]
mutantfamNamesPlot_sorted <- famNamesPlot[sort.int(wtlog2ObsExp, decreasing = T, index.return = T)$ix]

df <- data.frame(Genotype = rep(c(wtName, mutantName),
                                each = length(wtlog2ObsExp_sorted)),
                 Transposon_family = rep(wtfamNamesPlot_sorted, 2),
                 log2ObsExp = c(wtlog2ObsExp_sorted,
                                mutantlog2ObsExp_sorted),
                 log2alpha0.05 = c(wtlog2alpha0.05_sorted, mutantlog2alpha0.05_sorted))

df$Transposon_family <- factor(df$Transposon_family,
                               levels = wtfamNamesPlot_sorted)
df$Genotype <- factor(df$Genotype,
                      levels = c(wtName, mutantName))

bp <- ggplot(data = df,
             mapping = aes(x = Transposon_family,
                           y = log2ObsExp,
                           fill = Genotype)) +
      geom_bar(stat = "identity",
               position = position_dodge()) +
      scale_fill_manual(name = "Genotype",
                        values = c("black",
                                   "dodgerblue3"),
                        labels = c(wtName,
                                   bquote(italic(.(as.character(mutantName)))))) +
      geom_point(mapping = aes(Transposon_family, log2alpha0.05),
                 position = position_dodge(0.9),
                 shape = "-", colour = "red", size = 5) +
      labs(x = "Transposon family",
           y = expression("Log"[2]*"(observed:expected) RPM")) +
      theme_bw() +
      theme(axis.line.y = element_line(size = 0.5, colour = "black"),
            axis.ticks.y = element_line(size = 0.25, colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.ticks.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
            #legend.text = element_text(face = c("plain", "italic"))) +
      ggtitle(paste0(dataName, " (", as.character(perms), " permutations)")) 
ggsave(paste0(plotDir, "barplot_TE_families_permTestResults_",
              as.character(perms), "perms_",
              "log2_Observed_Expected_", wtLibName, "_", mutantLibName, "_RPM.pdf"),
       plot = bp,
       height = 4, width = 5)

