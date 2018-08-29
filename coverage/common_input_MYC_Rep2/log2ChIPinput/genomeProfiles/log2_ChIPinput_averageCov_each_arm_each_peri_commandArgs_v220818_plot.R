#!/applications/R/R-3.4.0/bin/Rscript

#############################################################################
# Sum all per base coverage values in each of 10 chromosome arm regions and #
# in each of 5 pericentromeric regions in wt and mutant                     #
############################################################################# 

# Usage on hydrogen node7
# csmit -m 20G -c 1 "./log2_ChIPinput_averageCov_each_arm_each_peri_commandArgs_v220818_plot.R REC8_HA_Rep1 kss_REC8_HA_Rep1"

library(ggbeeswarm)
library(reshape)
library(gridExtra)

args <- commandArgs(trailingOnly = TRUE)
wt_libName <- args[1]
mutant_libName <- args[2]

outDir <- "./" 

# Genomic definitions
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
# Pericentromeric regions are as defined in Supplemental Table S26
# of Ziolkowski et al. (2017) Genes Dev. 31
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)

wt_chrCov_armL_all <- read.table(paste0(outDir, wt_libName, "_chrCov_armL_all.txt"))
wt_chrCov_armR_all <- read.table(paste0(outDir, wt_libName, "_chrCov_armR_all.txt"))
wt_chrCov_peri_all <- read.table(paste0(outDir, wt_libName, "_chrCov_peri_all.txt"))
mutant_chrCov_armL_all <- read.table(paste0(outDir, mutant_libName, "_chrCov_armL_all.txt"))
mutant_chrCov_armR_all <- read.table(paste0(outDir, mutant_libName, "_chrCov_armR_all.txt"))
mutant_chrCov_peri_all <- read.table(paste0(outDir, mutant_libName, "_chrCov_peri_all.txt"))

wt_chrCov_all <- t(cbind(wt_chrCov_armL_all, wt_chrCov_peri_all, wt_chrCov_armR_all))
rownames(wt_chrCov_all) <- NULL
wt_chrCov_all <- data.frame(cbind(wt_chrCov_all, c("Left arm", "Pericentromere", "Right arm")))
colnames(wt_chrCov_all) <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Region")
wt_chrCov_all <- melt(wt_chrCov_all, id = "Region")
wt_chrCov_all$value <- as.numeric(as.character(wt_chrCov_all$value))
colnames(wt_chrCov_all) <- c("Region", "Chromosome", "Coverage")

mutant_chrCov_all <- t(cbind(mutant_chrCov_armL_all, mutant_chrCov_peri_all, mutant_chrCov_armR_all))
rownames(mutant_chrCov_all) <- NULL
mutant_chrCov_all <- data.frame(cbind(mutant_chrCov_all, c("Left arm", "Pericentromere", "Right arm")))
colnames(mutant_chrCov_all) <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Region")
mutant_chrCov_all <- melt(mutant_chrCov_all, id = "Region")
mutant_chrCov_all$value <- as.numeric(as.character(mutant_chrCov_all$value))
colnames(mutant_chrCov_all) <- c("Region", "Chromosome", "Coverage")

wt_chrCovPlot <- ggplot(wt_chrCov_all,
                     aes(x = Region, y = Coverage, colour = factor(Chromosome), group = factor(Chromosome))) +
                 geom_point(cex = 3) +
                 geom_line() +
                 xlab("Region") +
                 ylab("Normalized coverage") +
                 ylim(min(wt_chrCov_all$Coverage, mutant_chrCov_all$Coverage),
                      max(wt_chrCov_all$Coverage, mutant_chrCov_all$Coverage)) +
                 labs(colour = NULL) +
                 theme_bw() +
                 theme(axis.line.y = element_line(size = 0.5, colour = "black"),
                       axis.ticks.y = element_line(size = 0.25, colour = "black"),
                       axis.text.y = element_text(colour = "black"),
                       axis.ticks.x = element_blank(),
                       axis.text.x = element_text(colour = "black"),
                       panel.grid = element_blank(),
                       panel.border = element_blank(),
                       panel.background = element_blank(),
#                       plot.margin = grid::unit(c(0,0,0,0), "mm"),
                       plot.title = element_text(hjust = 0.5)) +
                 ggtitle(paste0(wt_libName, " average chromosome coverage"))

mutant_chrCovPlot <- ggplot(mutant_chrCov_all,
                            aes(x = Region, y = Coverage, colour = factor(Chromosome), group = factor(Chromosome))) +
                     geom_point(cex = 3) +
                     geom_line() +
                     xlab("Region") +
                     ylab("Normalized coverage") +
                     ylim(min(wt_chrCov_all$Coverage, mutant_chrCov_all$Coverage),
                          max(wt_chrCov_all$Coverage, mutant_chrCov_all$Coverage)) +
                     labs(colour = NULL) +
                     theme_bw() +
                     theme(axis.line.y = element_line(size = 0.5, colour = "black"),
                           axis.ticks.y = element_line(size = 0.25, colour = "black"),
                           axis.text.y = element_text(colour = "black"),
                           axis.ticks.x = element_blank(),
                           axis.text.x = element_text(colour = "black"),
                           panel.grid = element_blank(),
                           panel.border = element_blank(),
                           panel.background = element_blank(),
#                           plot.margin = grid::unit(c(0,0,0,0), "mm"),
                           plot.title = element_text(hjust = 0.5)) +
                     ggtitle(paste0(mutant_libName, " average chromosome coverage"))
singleplot <- list(wt_chrCovPlot,
                   mutant_chrCovPlot)
pdf(paste0(outDir, wt_libName, "_", mutant_libName, "_average_chromosome_coverage_revised.pdf"),
    onefile = T,
    width = 12, height = 5)
grid.arrange(singleplot[[1]],
             singleplot[[2]],
             nrow = 1, ncol = 2)
dev.off()

