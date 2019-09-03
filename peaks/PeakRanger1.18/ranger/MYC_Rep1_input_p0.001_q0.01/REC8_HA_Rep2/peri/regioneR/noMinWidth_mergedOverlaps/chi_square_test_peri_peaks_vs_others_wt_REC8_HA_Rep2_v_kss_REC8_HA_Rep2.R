#!/applications/R/R-3.5.0/bin/Rscript
#
# Compare proportions of wt and kss peaks overlapping features using chi-square tests of independence

# Usage: 
# /applications/R/R-3.5.0/bin/Rscript chi_square_test_peri_peaks_vs_others_wt_REC8_HA_Rep2_v_kss_REC8_HA_Rep2.R REC8_HA_Rep2_ChIP_rangerPeaksGR_peri_mergedOverlaps_noMinWidth.RData kss_REC8_HA_Rep2_ChIP_rangerPeaksGR_peri_mergedOverlaps_noMinWidth.RData "wt REC8-HA_Rep2" "kss REC8-HA_Rep2" "wt SPO11-1" "kss SPO11-1" REC8_HA_Rep2 kss_REC8_HA_Rep2 wt_SPO11oligo kss_SPO11oligo "Pericentromeric REC8-HA peaks and SPO11-1-oligo hotspots" peri

#peakFile1 <- "REC8_HA_Rep2_ChIP_rangerPeaksGR_peri_mergedOverlaps_noMinWidth.RData"
#peakFile2 <- "kss_REC8_HA_Rep2_ChIP_rangerPeaksGR_peri_mergedOverlaps_noMinWidth.RData"
#peakName1 <- "wt REC8-HA_Rep2"
#peakName2 <- "kss REC8-HA_Rep2"
#peakName3 <- "wt SPO11-1"
#peakName4 <- "kss SPO11-1"
#peakLibName1 <- "REC8_HA_Rep2"
#peakLibName2 <- "kss_REC8_HA_Rep2"
#peakLibName3 <- "wt_SPO11oligo"
#peakLibName4 <- "kss_SPO11oligo"
#plotTitle <- "Pericentromeric REC8-HA peaks and SPO11-1-oligo hotspots"
#region <- "peri"

args <- commandArgs(trailingOnly = T)
peakFile1 <- args[1]
peakFile2 <- args[2]
peakName1 <- args[3]
peakName2 <- args[4]
peakName3 <- args[5]
peakName4 <- args[6]
peakLibName1 <- args[7]
peakLibName2 <- args[8]
peakLibName3 <- args[9]
peakLibName4 <- args[10]
plotTitle <- args[11]
region <- args[12]

library(GenomicRanges)
library(regioneR)
library(ggplot2)
library(ggthemes)
library(ggsignif)

# Genomic definitions
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chrStart <- c(1, 1, 1, 1, 1)
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))
mask <- toGRanges(data.frame(rep(chrs, 2),
                             c(chrStart, pericenEnd),
                             c(pericenStart, chrLens)))

inDir <- "/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep1_input_p0.001_q0.01/"
outDir <- "./"

plotDir <- "./bar_graphs/"
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# wt REC8 peaks
load(paste0(inDir,
            peakFile1))
wtREC8GR <- rangerPeaksGR_peri_mergedOverlaps
rangerPeaksGR_peri_mergedOverlaps <- NULL
strand(wtREC8GR) <- "*"
print("***********wt_REC8_HA_Rep2_peaks***********")
print(wtREC8GR)
print(length(wtREC8GR))

# kss REC8 peaks
load(paste0("/home/ajt200/analysis/180830_Chris_lambing_ChIP_rec8HA_kss_rep2/peaks/PeakRanger1.18/ranger/kss_REC8_HA_Rep1_input_p0.001_q0.01/",
            peakFile2))
kssREC8GR <- rangerPeaksGR_peri_mergedOverlaps
rangerPeaksGR_peri_mergedOverlaps <- NULL
strand(kssREC8GR) <- "*"
print("***********kss_REC8_HA_Rep2_peaks***********")
print(kssREC8GR)
print(length(kssREC8GR))

# wt SPO11-1-oligo hotspots
load("/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/perirangerPeaksGRmerge_RPI1_RPI8_idr0.05_noMinWidth.RData")
wtSPO11GR <- perirangerPeaksGRmerge
perirangerPeaksGRmerge <- NULL
strand(wtSPO11GR) <- "*"
print(length(wtSPO11GR))

# kss SPO11-1-oligo hotspots
load("/projects/ajt200/BAM_masters/SPO11-oligo/suvh456/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/perirangerPeaksGRmerge_RPI34_RPI35_idr0.05_noMinWidth.RData")
kssSPO11GR <- perirangerPeaksGRmerge
perirangerPeaksGRmerge <- NULL
strand(kssSPO11GR) <- "*"
print(length(kssSPO11GR))

# Others
load("/projects/ajt200/REC8_MSH4/nuc_peaks/log2ChIPinput/nucleR/trim/analysis_01/periPeaksSH99GRmerge.RData")
nucleRnucsGR <- periPeaksSH99GRmerge
periPeaksSH99GRmerge <- NULL
strand(nucleRnucsGR) <- "*"
print(length(nucleRnucsGR))

load("/projects/ajt200/BAM_masters/H3K9me2/WT/peaks/PeakRanger1.18/ranger/REC8_MYC_Rep1_input_p0.05_q0.05/WT_H3K9me2_ChIP_rangerPeaksGR_peri_mergedOverlaps_noMinWidth.RData")
H3K9me2GR <- rangerPeaksGR_peri_mergedOverlaps
rangerPeaksGR_peri_mergedOverlaps <- NULL
strand(H3K9me2GR) <- "*"
print(length(H3K9me2GR))

kss_hypoCHG_DMRs <- read.table("/home/ajt200/BS_Seq/Stroud_2013/DMRs/suvh456_hypoCHG_DMR_vs3reps_min4filter_mg200.bed", header = F)
kss_hypoCHG_DMRsGR <- sort(GRanges(seqnames = kss_hypoCHG_DMRs[,1],
                                   ranges = IRanges(start = kss_hypoCHG_DMRs[,2],
                                                    end = kss_hypoCHG_DMRs[,3]),
                                   strand = "*"))
seqlevels(kss_hypoCHG_DMRsGR) <- sub("chr", "Chr", seqlevels(kss_hypoCHG_DMRsGR))
mask_kss_hypoCHG_DMRsOverlaps <- findOverlaps(query = mask,
                                              subject = kss_hypoCHG_DMRsGR,
                                              ignore.strand = TRUE,
                                              select = "all")
kss_hypoCHG_DMRsGR <- kss_hypoCHG_DMRsGR[-subjectHits(mask_kss_hypoCHG_DMRsOverlaps)]
print(length(kss_hypoCHG_DMRsGR))

kss_hypoCHH_DMRs <- read.table("/home/ajt200/BS_Seq/Stroud_2013/DMRs/suvh456_hypoCHH_DMR_vs3reps_min4filter_mg200.bed", header = F)
kss_hypoCHH_DMRsGR <- sort(GRanges(seqnames = kss_hypoCHH_DMRs[,1],
                                   ranges = IRanges(start = kss_hypoCHH_DMRs[,2],
                                                    end = kss_hypoCHH_DMRs[,3]),
                                   strand = "*"))
seqlevels(kss_hypoCHH_DMRsGR) <- sub("chr", "Chr", seqlevels(kss_hypoCHH_DMRsGR))
mask_kss_hypoCHH_DMRsOverlaps <- findOverlaps(query = mask,
                                              subject = kss_hypoCHH_DMRsGR,
                                              ignore.strand = TRUE,
                                              select = "all")
kss_hypoCHH_DMRsGR <- kss_hypoCHH_DMRsGR[-subjectHits(mask_kss_hypoCHH_DMRsOverlaps)]
print(length(kss_hypoCHH_DMRsGR))

# TEs
DNAfamNames <- c("dna", "heli", "ptmari", "mudr", "enspm", "hat", "harbinger")
DNAfamNames <- DNAfamNames[5]
RNAfamNames <- c("rna", "gypsy", "copia", "linel1", "sine")
RNAfamNames <- RNAfamNames[2]
DNAfamNamesPlot <- c("DNA", "Helitron", "Pogo/Tc1/Mariner", "MuDR", "EnSpm/CACTA", "hAT", "Harbinger")
DNAfamNamesPlot <- DNAfamNamesPlot[5]
RNAfamNamesPlot <- c("RNA", "Gypsy LTR", "Copia LTR", "LINE-1", "SINE")
RNAfamNamesPlot <- RNAfamNamesPlot[2]

DNAdir <- "/projects/ajt200/TAIR10/TE_classes/DNA/"
RNAdir <- "/projects/ajt200/TAIR10/TE_classes/RNA/"

### DNA TEs
TEsDNAGR <- lapply(seq_along(DNAfamNames), function(x) {
  TEsDNA <- read.table(paste0(DNAdir,
                              "TAIR10_Buisine_TEs_strand_tab_ann_",
                              DNAfamNames[x], ".txt"),
                       header = T)
  TEsDNAGR_all <- GRanges(seqnames = TEsDNA$chr,
                          ranges = IRanges(start = TEsDNA$start,
                                           end = TEsDNA$end),
                          strand = "*")
  maskTEsDNAoverlaps <- findOverlaps(query = mask,
                                     subject = TEsDNAGR_all,
                                     ignore.strand = TRUE,
                                     select = "all")
  TEsDNAGR_all[-subjectHits(maskTEsDNAoverlaps)]
})
### RNA TEs
TEsRNAGR <- lapply(seq_along(RNAfamNames), function(x) {
  TEsRNA <- read.table(paste0(RNAdir,
                              "TAIR10_Buisine_TEs_strand_tab_ann_",
                              RNAfamNames[x], ".txt"),
                       header = T)
  TEsRNAGR_all <- GRanges(seqnames = TEsRNA$chr,
                          ranges = IRanges(start = TEsRNA$start,
                                           end = TEsRNA$end),
                          strand = "*")
  maskTEsRNAoverlaps <- findOverlaps(query = mask,
                                     subject = TEsRNAGR_all,
                                     ignore.strand = TRUE,
                                     select = "all")
  TEsRNAGR_all[-subjectHits(maskTEsRNAoverlaps)]
})

otherNames <- c(
                "nucleRnucsGR",
                "kss_hypoCHH_DMRsGR",
                "kss_hypoCHG_DMRsGR",
                "H3K9me2GR",
                "enspm",
                "gypsy"
               )
otherNamesPlot <- c(
                    "Nucleosomes",
                    "Hypo-mCHH",
                    "Hypo-mCHG",
                    "H3K9me2",
                    "EnSpm/CACTA",
                    "Gypsy LTR"
                   )
grl <- c(
         "nucleRnucsGR" = nucleRnucsGR,
         "kss_hypoCHH_DMRsGR" = kss_hypoCHH_DMRsGR,
         "kss_hypoCHG_DMRsGR" = kss_hypoCHG_DMRsGR,
         "H3K9me2GR" = H3K9me2GR,
         "enspm" = TEsDNAGR[[1]],
         "gypsy" = TEsRNAGR[[1]]
        )

# Function to create 2x2 contingency table where wtPeaksIn and wtPeaksOut are the
# number of wild type peaks overlapping a given set of features
contingencyTable <- function(wtPeaksIn, wtPeaksOut, kssPeaksIn, kssPeaksOut) {
  conTab <- matrix(c(wtPeaksIn, wtPeaksOut, kssPeaksIn, kssPeaksOut), ncol = 2)
  rownames(conTab) <- c("inFeature", "outFeature")
  colnames(conTab) <- c("wt", "kss")
  conTab
}

# REC8
# Count number of peaks overlapping given sets of features
wtREC8overlaps <- sapply(seq_along(grl), function(x) {
  sum(countOverlaps(query = wtREC8GR,
                    subject = grl[[x]],
                    ignore.strand = TRUE) > 0)
})
kssREC8overlaps <- sapply(seq_along(grl), function(x) {
  sum(countOverlaps(query = kssREC8GR,
                    subject = grl[[x]],
                    ignore.strand = TRUE) > 0)
})

# Calculate proportion of peaks overlapping given sets of features
wtREC8prop <- sapply(seq_along(wtREC8overlaps), function(x) {
  (wtREC8overlaps[x]/(length(wtREC8GR)-wtREC8overlaps[x]))*100
})
kssREC8prop <- sapply(seq_along(kssREC8overlaps), function(x) {
  (kssREC8overlaps[x]/(length(kssREC8GR)-kssREC8overlaps[x]))*100
})

# Perform chi-square test on each contingency table and extract P-values
chisqPvalsREC8 <- sapply(1:length(grl), function(x) {
  chisq.test(x = contingencyTable(wtPeaksIn = wtREC8overlaps[x],
                                  wtPeaksOut = length(wtREC8GR)-wtREC8overlaps[x],
                                  kssPeaksIn = kssREC8overlaps[[x]], 
                                  kssPeaksOut = length(kssREC8GR)-kssREC8overlaps[x]))$p.value
})
  
# SPO11
# Count number of peaks overlapping given sets of features
wtSPO11overlaps <- sapply(seq_along(grl), function(x) {
  sum(countOverlaps(query = wtSPO11GR,
                    subject = grl[[x]],
                    ignore.strand = TRUE) > 0)
})
kssSPO11overlaps <- sapply(seq_along(grl), function(x) {
  sum(countOverlaps(query = kssSPO11GR,
                    subject = grl[[x]],
                    ignore.strand = TRUE) > 0)
})

# Calculate proportion of peaks overlapping given sets of features
wtSPO11prop <- sapply(seq_along(wtSPO11overlaps), function(x) {
  (wtSPO11overlaps[x]/(length(wtSPO11GR)-wtSPO11overlaps[x]))*100
})
kssSPO11prop <- sapply(seq_along(kssSPO11overlaps), function(x) {
  (kssSPO11overlaps[x]/(length(kssSPO11GR)-kssSPO11overlaps[x]))*100
})

# Perform chi-square test on each contingency table and extract P-values
chisqPvalsSPO11 <- sapply(1:length(grl), function(x) {
  chisq.test(x = contingencyTable(wtPeaksIn = wtSPO11overlaps[x],
                                  wtPeaksOut = length(wtSPO11GR)-wtSPO11overlaps[x],
                                  kssPeaksIn = kssSPO11overlaps[[x]],
                                  kssPeaksOut = length(kssSPO11GR)-kssSPO11overlaps[x]))$p.value
})

# Disable scientific notation
options(scipen = 999)

# Simplify P-values
chisqPvalsREC8char <- sapply(seq_along(chisqPvalsREC8),
  function(x) {
  if(chisqPvalsREC8[x] < 0.0001) {
    "< 0.0001"
  } else {
    paste0("= ", as.character(round(chisqPvalsREC8[x], digits = 4)))
  }
})
chisqPvalsSPO11char <- sapply(seq_along(chisqPvalsSPO11),
  function(x) {
  if(chisqPvalsSPO11[x] < 0.0001) {
    "< 0.0001"
  } else {
    paste0("= ", as.character(round(chisqPvalsSPO11[x], digits = 4)))
  }
})
  
# Combine in data.frame
df <- data.frame(Sample = rep(c(peakName1, peakName2, peakName3, peakName4),
                              each = length(grl)),
                 Feature = rep(otherNamesPlot, 4),
                 Overlap_proportion = c(wtREC8prop,
                                        kssREC8prop,
                                        wtSPO11prop,
                                        kssSPO11prop))

df$Feature <- factor(df$Feature,
                     levels = otherNamesPlot)
df$Sample <- factor(df$Sample,
                    levels = c(peakName1, peakName2, peakName3, peakName4))

bp <- ggplot(data = df,
             mapping = aes(x = Feature,
                           y = Overlap_proportion,
                           fill = Sample)) +
  geom_bar(stat = "identity",
           position = position_dodge()) +
  scale_fill_manual(name = "",
                    values = c("red",
                               "red4",
                               "dodgerblue2",
                               "navy"),
                    labels = c(peakName1,
                               bquote(italic(.(unlist(strsplit(peakName2,
                                                               split = " "))[1])) ~
                                     .(unlist(strsplit(peakName2,
                                                       split = " "))[2])),
                               peakName3,
                               bquote(italic(.(unlist(strsplit(peakName4,
                                                               split = " "))[1])) ~
                                      .(unlist(strsplit(peakName4,
                                                        split = " "))[2])))) +
  labs(y =  expression("% peaks overlapping features"[ ])) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_x_discrete(position = "top") +
  guides(fill = guide_legend(direction = "horizontal",
                             label.position = "top",
                             label.theme = element_text(size = 20, hjust = 0, vjust = 0.5, angle = 90),
                             nrow = 2,
                             byrow = TRUE)) +
  theme_bw() +
  theme(axis.line.y = element_line(size = 1, colour = "black"),
        axis.ticks.y = element_line(size = 1, colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black", hjust = 0.5, vjust = 0.5, angle = 90),
        axis.title.y = element_text(size = 20, colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 20, colour = "black", hjust = 0, vjust = 0.5, angle = 90),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(colour = "transparent",
                                  fill = "transparent"),
        legend.key.size = unit(1, "cm"),
        plot.margin = unit(c(5.5, 5.5, 10.5, 5.5), "pt"),
        plot.title = element_text(size = 17, colour = "black", hjust = 0.5)) +
  ggtitle(plotTitle) +
  annotate(geom = "text",
           size = 6,
           x = c(0.8, 1.2,
                 1.8, 2.2,
                 2.8, 3.2,
                 3.8, 4.2,
                 4.8, 5.2,
                 5.8, 6.2),
           y = c(90), angle = 90,
           label = c(paste0("P ", chisqPvalsREC8char[1]),
                     paste0("P ", chisqPvalsSPO11char[1]),
                     paste0("P ", chisqPvalsREC8char[2]),
                     paste0("P ", chisqPvalsSPO11char[2]), 
                     paste0("P ", chisqPvalsREC8char[3]),
                     paste0("P ", chisqPvalsSPO11char[3]),
                     paste0("P ", chisqPvalsREC8char[4]),
                     paste0("P ", chisqPvalsSPO11char[4]),
                     paste0("P ", chisqPvalsREC8char[5]),
                     paste0("P ", chisqPvalsSPO11char[5]),
                     paste0("P ", chisqPvalsREC8char[6]),
                     paste0("P ", chisqPvalsSPO11char[6])))
#  geom_signif(y_position = c(70), xmin = c(0.4), xmax = c(0.8),
#              annotation = c(paste0("P = ", chisqPvalsREC8[1])),
#              tip_length = 0) 
ggsave(paste0(plotDir, "barplot_selected_other_features_TE_families_proportion_of_overlapping_",
              peakLibName1, "_", peakLibName2, "_",
              peakLibName3, "_", peakLibName4, "_",
              region, "_peaks.pdf"),
       plot = bp,
       height = 8, width = 10)
