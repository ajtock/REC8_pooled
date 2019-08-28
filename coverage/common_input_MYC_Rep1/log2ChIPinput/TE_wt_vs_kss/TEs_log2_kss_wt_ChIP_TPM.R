#!/applications/R/R-3.5.0/bin/Rscript

# Rank TEs by differential wt vs kss ChIP TPM

# Usage:
# Rscript ./TEs_log2_kss_wt_ChIP_TPM.R '/home/ajt200/analysis/180622_Chris_lambing_ChIP_REC8_HA_Col_kss/WT/' '/home/ajt200/analysis/180622_Chris_lambing_ChIP_REC8_HA_Col_kss/kss/' REC8_HA_Rep2 kss_REC8_HA_Rep1 REC8 PE 2

library(GenomicAlignments)
library(dplyr)
library(ggplot2)
library(hexbin)

#geno1Dir <- "/home/ajt200/analysis/180622_Chris_lambing_ChIP_REC8_HA_Col_kss/WT/"
#geno2Dir <- "/home/ajt200/analysis/180622_Chris_lambing_ChIP_REC8_HA_Col_kss/kss/"
#geno1Name <- "REC8_HA_Rep2"
#geno2Name <- "kss_REC8_HA_Rep1"
#ChIP <- "REC8"
#libType <- "PE"
#L2FC <- 2

args <- commandArgs(trailingOnly = T)
geno1Dir <- args[1]
geno2Dir <- args[2]
geno1Name <- args[3]
geno2Name <- args[4]
ChIP <- args[5]
libType <- args[6]
L2FC <- as.numeric(args[7])

outDir <- paste0(ChIP, "/")
profileDir <- paste0(outDir, "heatmaps_analysis_01_ChIPonly/")
gainDir <- paste0(profileDir, "gain/")
lossDir <- paste0(profileDir, "loss/")
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
system(paste0("[ -d ", profileDir, " ] || mkdir ", profileDir))
system(paste0("[ -d ", gainDir, " ] || mkdir ", gainDir))
system(paste0("[ -d ", lossDir, " ] || mkdir ", lossDir))

# Load TEs and convert into GRanges
TEs <- read.table("/projects/ajt200/TAIR10/TAIR10_Buisine_TEs_strand_tab_ann.txt",
                  header = T)
colnames(TEs) <- c("chr", "start", "end", "strand", "TE_ID", "TE_family", "TE_superfamily")
featuresGR <- GRanges(seqnames = TEs$chr,
                      ranges = IRanges(start = TEs$start,
                                       end = TEs$end),
                      strand = TEs$strand,
                      TE_ID = TEs$TE_ID,
                      TE_family = TEs$TE_family,
                      TE_superfamily = TEs$TE_superfamily)

# Load BAMs and create GAlignmentPairs objects
if(libType == "SE") {
  geno1_readGAlignments <- readGAlignments(paste0(geno1Dir,
                                                  geno1Name,
                                                  "_k10_bt2_mapped_lowmiss_unique_both_sort_rmdup_new_sort.bam"))
  geno2_readGAlignments <- readGAlignments(paste0(geno2Dir,
                                                  geno2Name,
                                                  "_k10_bt2_mapped_lowmiss_unique_both_sort_rmdup_new_sort.bam"))
} else if(libType == "PE") {
  geno1_readGAlignments <- readGAlignmentPairs(paste0(geno1Dir,
                                                      geno1Name,
                                                      "_ChIP_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam"))
  geno2_readGAlignments <- readGAlignmentPairs(paste0(geno2Dir,
                                                      geno2Name,
                                                      "_ChIP_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam"))
} else {
  stop("libType must be SE or PE")
}

# Convert into GRanges
geno1GR <- as(geno1_readGAlignments[seqnames(geno1_readGAlignments) != "mitochondria" &
                                    seqnames(geno1_readGAlignments) != "chloroplast"],
              "GRanges")
seqlevels(geno1GR) <- sub("", "Chr", seqlevels(geno1GR))
geno2GR <- as(geno2_readGAlignments[seqnames(geno2_readGAlignments) != "mitochondria" &
                                    seqnames(geno2_readGAlignments) != "chloroplast"],
              "GRanges")
seqlevels(geno2GR) <- sub("", "Chr", seqlevels(geno2GR))

# Function to calculate feature TPM for each library
featureCovCalc <- function(features, reads, libName) {
  feature_reads <- countOverlaps(query = features,
                                 subject = reads,
                                 type = "any",
                                 ignore.strand = TRUE)
  feature_RPK <- feature_reads/(width(features)/1e+03)
  RPKM_scaling_factor <- sum(feature_RPK)/1e+06
  feature_TPM <- feature_RPK/RPKM_scaling_factor
  data.frame(reads = as.integer(feature_reads),
             RPK = as.numeric(feature_RPK),
             TPM = as.numeric(feature_TPM))
}

# Apply featureTPMcalc() to each library and rename colnames
featureCov_geno1 <- featureCovCalc(features = featuresGR,
                                   reads = geno1GR,
                                   libName = geno1Name)
featureCov_geno2 <- featureCovCalc(features = featuresGR,
                                   reads = geno2GR,
                                   libName = geno2Name)

# Plot feature TPM means against SDs
feature_mean_TPM <- sapply(seq_along(featuresGR), function(x) {
  mean(c(featureCov_geno1$TPM[x], featureCov_geno2$TPM[x]))
})
feature_SD_TPM <- sapply(seq_along(featuresGR), function(x) {

  sd(c(featureCov_geno1$TPM[x], featureCov_geno2$TPM[x]))
})

feature_mean_TPM_SDs_df <- as_data_frame(cbind(TPM_mean = feature_mean_TPM,
                                               TPM_SD = feature_SD_TPM))
featureTPM_df <- bind_rows(
  as_data_frame(cbind(featureCov_geno1$TPM,
                      featureCov_geno2$TPM)) %>%
    mutate(transformation = "none"),
  as_data_frame(cbind(log2(featureCov_geno1$TPM+1),
                      log2(featureCov_geno2$TPM+1))) %>%
    mutate(transformation = "log2(TPM+1)"),
)
colnames(featureTPM_df)[1:2] <- c(geno1Name, geno2Name)

TPM_means_SDs_plot <- ggplot(feature_mean_TPM_SDs_df,
                             aes(x = TPM_mean,
                                 y = TPM_SD)) +
                      geom_hex(bins = 80) +
                      xlim(0, 500) + 
                      ylim(0, 500) + 
                      coord_fixed() +
                      labs(fill = "TEs") +
                      theme_classic()
ggsave(TPM_means_SDs_plot,
       file = paste0(outDir,
                     geno1Name, "_vs_", geno2Name, "_TE_TPM_mean_vs_SD.pdf"))

#TPM_plot <- ggplot(featureTPM_df,
#                   aes(x = geno1Name,
#                       y = geno2Name)) +
#            geom_hex(bins = 80) +
#            coord_fixed() +
#            facet_grid(. ~ transformation) +
#            labs(fill = "Occurrences") +
#            theme_classic()
#ggsave(TPM_plot, height = 5, width = 10,
#       file = paste0(outDir,
#                     geno1Name, "_vs_", geno2Name, "_TE_TPM.pdf"))

# Select features with greatest loss and gain of ChIP coverage in kss
featuresGR <- as(data.frame(featuresGR,
                            wt_TPM = featureCov_geno1$TPM,
                            kss_TPM = featureCov_geno2$TPM,
                            mean_TPM = feature_mean_TPM,
                            log2_kss_wt_TPM = log2((featureCov_geno2$TPM+1) /
                                                        (featureCov_geno1$TPM+1))),
                 "GRanges")

featuresGR_kssLoss <- featuresGR[featuresGR$log2_kss_wt_TPM <= -L2FC]
featuresGR_kssGain <- featuresGR[featuresGR$log2_kss_wt_TPM >= L2FC]
save(featuresGR_kssLoss,
     file = paste0(outDir,
                   "TAIR10_Buisine_TEs_log2_", geno2Name, "_", geno1Name,
                   "_TPM_lessThanOrEqualToMinus", as.character(L2FC),
                   ".RData"))
save(featuresGR_kssGain,
     file = paste0(outDir,
                   "TAIR10_Buisine_TEs_log2_", geno2Name, "_", geno1Name,
                   "_TPM_moreThanOrEqualToPlus", as.character(L2FC),
                   ".RData"))

# Make MA plot           
feature_meanTPM_log2TPM_df <- as_data_frame(cbind(mean_TPM = featuresGR$mean_TPM,
                                                  log2_TPM = featuresGR$log2_kss_wt_TPM))
# Plotting function
MAplotFun <- function(L2FC, SCMvalues) {
 ggplot(data = feature_meanTPM_log2TPM_df,
        mapping = aes(x = mean_TPM,
                      y = log2_TPM)) +
 geom_point(aes(colour = cut(log2_TPM,
                             c(-Inf, -L2FC, L2FC-0.001, Inf))),
            shape = 20) +
 scale_color_manual(name = expression("Log"[2] * "(fold change)"),
                    values = SCMvalues,
                    labels = c(bquote("<= -" * .(L2FC)),
                               bquote("-" * .(L2FC) * " < change < " * .(L2FC)),
                               bquote(">= " * .(L2FC)))) +
 geom_hline(yintercept = c(-L2FC, L2FC),
            linetype = "dashed",
            colour = "black",
            size = 1.0) +
 xlab("Mean TPM") +
 ylab(bquote("Log"[2] * "(" *
             .(geno2Name) * "/" * .(geno1Name) *
             ") TPM")) +
 xlim(0, max(feature_meanTPM_log2TPM_df$mean_TPM)) +
 ylim(-(max(abs(min(feature_meanTPM_log2TPM_df$log2_TPM)),
            max(feature_meanTPM_log2TPM_df$log2_TPM))),
      max(abs(min(feature_meanTPM_log2TPM_df$log2_TPM)),
          max(feature_meanTPM_log2TPM_df$log2_TPM))) +
 theme_bw() +
 theme(axis.ticks.length = unit(5, "pt"),
       axis.ticks.x = element_line(size = 1, colour = "black"),
       axis.ticks.y = element_line(size = 1, colour = "black"),
       axis.text.x = element_text(size = 14, colour = "black"),
       axis.text.y = element_text(size = 14, colour = "black"),
       axis.title = element_text(size = 14, colour = "black"),
       legend.text = element_text(size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black"),
       panel.border = element_rect(size = 1.5, colour = "black"),
       panel.grid = element_blank(),
       panel.background = element_blank(),
       plot.margin = unit(c(0.3, 0.9, 0.9, 0.3), "cm"),
       plot.title = element_text(hjust = 0.5, size = 14)) +
 ggtitle(bquote("MA plot of log"[2] * "(" *
                .(geno2Name) * "/" * .(geno1Name) *
                ") TPM at TEs"))
}

MAplot <- MAplotFun(L2FC = L2FC,
                    SCMvalues = c("dodgerblue2", "grey40", "orangered"))
ggsave(MAplot,
       file = paste0(outDir,
                     "TAIR10_Buisine_TEs_log2_",
                     geno2Name, "_", geno1Name,
                     "_TPM_L2FC", as.character(L2FC), "_MAplot.pdf"),
       height = 6, width = 10)
