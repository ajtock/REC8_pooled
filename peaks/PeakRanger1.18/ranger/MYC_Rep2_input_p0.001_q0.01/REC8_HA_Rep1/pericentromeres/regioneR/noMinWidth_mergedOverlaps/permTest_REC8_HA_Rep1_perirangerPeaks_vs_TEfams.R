#!/applications/R/R-3.3.2/bin/Rscript

# Use permutation test function in regioneR to determine if peaks overlap features
# of interest (e.g., other peaks, genes, TEs) more or less than expected by chance

# Usage on hydrogen node7:
# csmit -m 50G -c 48 "Rscript permTest_REC8_HA_Rep1_perirangerPeaks_vs_TEfams.R"

library(regioneR)

# Genomic definitions
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chrStart <- c(1, 1, 1, 1, 1)
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))
mask <- toGRanges(data.frame(rep(chrs, 2), c(chrStart, pericenEnd), c(pericenStart, chrLens)))

inDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/"
outDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/pericentromeres/regioneR/noMinWidth_mergedOverlaps/"
DNAplotDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/pericentromeres/regioneR/noMinWidth_mergedOverlaps/plots/TEsDNA/"
RNAplotDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/pericentromeres/regioneR/noMinWidth_mergedOverlaps/plots/TEsRNA/"

# REC8 peaks
load(paste0(inDir,
            "REC8_HA_Rep1_perirangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData"))
peaksGR <- perirangerPeaksGRmergedOverlaps
perirangerPeaksGRmergedOverlaps <- NULL
maskPeaksOverlaps <- findOverlaps(mask, peaksGR,
                                  ignore.strand = T,
                                  select = "all")
print("Peaks located within masked regions:")
print(maskPeaksOverlaps)
#peaksGR <- peaksGR[width(peaksGR) <= 2000]
print(length(peaksGR))
#[1] 24044 

DNAfamNames <- c("dna", "heli", "ptmari", "mudr", "enspm", "hat", "harbinger")
RNAfamNames <- c("rna", "gypsy", "copia", "linel1", "sine")
DNAdir <- "/projects/ajt200/TAIR10/TE_classes/DNA/"
RNAdir <- "/projects/ajt200/TAIR10/TE_classes/RNA/"

# ### DNA TEs
# 
# TEsDNAGR <- lapply(seq_along(DNAfamNames), function(x) {
#   TEsDNA <- read.table(file = paste0(DNAdir, "TAIR10_Buisine_TEs_strand_tab_ann_", DNAfamNames[x], ".txt"), header = T)
#   GRanges(seqnames = TEsDNA$chr, ranges = IRanges(start = TEsDNA$start, end = TEsDNA$end), strand = "*")
# })
# 
# # Remove TEs located within masked regions
# TEsDNAGRmasked <- lapply(seq_along(TEsDNAGR), function(x) {
#   maskTEsOverlaps <- findOverlaps(mask, TEsDNAGR[[x]], ignore.strand = TRUE, select = "all")
#   TEsDNAGR[[x]][-subjectHits(maskTEsOverlaps)]
# })
# 
# # Perform permutation tests with randomized regions generated on a per chromosome basis;
# # same per-chromosome number and size of regions in B as in A
# set.seed(123)
# ptPeaksTEsDNAPerChrom <- lapply(seq_along(TEsDNAGRmasked), function(x) {
#   permTest(A = peaksGR, B = TEsDNAGRmasked[[x]], genome = genome, mask = mask,
#            randomize.function = randomizeRegions,
#            allow.overlaps = TRUE, per.chromosome = TRUE,
#            evaluate.function = numOverlaps, count.once = TRUE,
#            ntimes = 10000, mc.set.seed = FALSE, mc.cores = 20)
# })
# 
# for(i in 1:length(ptPeaksTEsDNAPerChrom)) {
#   assign(paste0(DNAfamNames[i]), ptPeaksTEsDNAPerChrom[[i]])
# }
# save(ptPeaksTEsDNAPerChrom,
#      file = paste0(outDir,
#                    "permTest_REC8_HA_Rep1_perirangerPeaks_vs_TEsDNA.RData"))
# 
# # Summarise results in a table
# noOfFeatures <- NULL
# expected <- NULL
# observed <- NULL
# pval <- NULL
# zscore <- NULL
# for(i in 1:length(ptPeaksTEsDNAPerChrom)) {
#   noOfFeaturesi <- print(length(TEsDNAGRmasked[[i]]))
#   noOfFeatures <- c(noOfFeatures, noOfFeaturesi)
#   expectedi <- print(round(mean(ptPeaksTEsDNAPerChrom[[i]]$numOverlaps$permuted)))
#   expected <- c(expected, expectedi)
#   observedi <- print(ptPeaksTEsDNAPerChrom[[i]]$numOverlaps$observed)
#   observed <- c(observed, observedi)
#   pvali <- print(round(ptPeaksTEsDNAPerChrom[[i]]$numOverlaps$pval, 4))
#   pval <- c(pval, pvali)
#   zscorei <- print(round(ptPeaksTEsDNAPerChrom[[i]]$numOverlaps$zscore, 4))
#   zscore <- c(zscore, zscorei)
# }
# ptPeaksTEsDNAPerChromDataFrame <- cbind(noOfFeatures, expected, observed, pval, zscore)
# write.table(ptPeaksTEsDNAPerChromDataFrame,
#             file = paste0(outDir,
#                           "permTest_REC8_HA_Rep1_perirangerPeaks_vs_TEsDNA_DataFrame.txt"),
#             sep = "\t", row.names = F)
# 
# # plot graphical summaries of results
# for(i in 1:length(ptPeaksTEsDNAPerChrom)) {
#   pdf(file = paste0(DNAplotDir, DNAfamNames[i], "_permTest_nperm10000_REC8_HA_Rep1_perirangerPeaks_perChrom.pdf"), width = 10, height = 7)
#   plot(ptPeaksTEsDNAPerChrom[[i]], main = paste0("REC8_HA_Rep1_perirangerPeaks vs ", DNAfamNames[i], " in pericentromeric regions"), xlab = "Number of overlaps", ylab = "Density")
#   dev.off()
# 
#   # Using the localZScore() function, evaluate whether the association between peaks and other is highly dependent on their exact position
#   lz_1kb <- localZScore(pt = ptPeaksTEsDNAPerChrom[[i]], A = peaksGR, B = TEsDNAGRmasked[[i]],
#                         window = 1000, step = 50, count.once = TRUE)
#   lz_10kb <- localZScore(pt = ptPeaksTEsDNAPerChrom[[i]], A = peaksGR, B = TEsDNAGRmasked[[i]],
#                          window = 10000, step = 500, count.once = TRUE)
#   lz_custom <- localZScore(pt = ptPeaksTEsDNAPerChrom[[i]], A = peaksGR, B = TEsDNAGRmasked[[i]],
#                            window = 10*mean(width(peaksGR)), step = mean(width(peaksGR))/2, count.once = TRUE)
#   win <- as.character(round((10*mean(width(peaksGR)))/1000))
#   step <- as.character(round(mean(width(peaksGR))/2))
#   pdf(file = paste0(DNAplotDir, DNAfamNames[i], "_localZscore_permTest_nperm10000_REC8_HA_Rep1_perirangerPeaks_w1kb_s50bp_w10kb_s500bp_w", win ,"kb_s", step, "bp_perChrom.pdf"))
#   par(mar=c(5.1, 4.1, 4.1, 2.1))
#   plot(lz_1kb, main = paste0("REC8_HA_Rep1_perirangerPeaks vs ", DNAfamNames[i], " in pericentromeric regions (1-kb shift)"))
#   mtext(side = 3, at = 2, text = paste0("REC8_HA_Rep1_perirangerPeaks vs ", DNAfamNames[i], " in pericentromeric regions (1-kb shift)"))
#   plot(lz_10kb, main = paste0("REC8_HA_Rep1_perirangerPeaks vs ", DNAfamNames[i], " in pericentromeric regions (10-kb shift)"))
#   mtext(side = 3, at = 2, text = paste0("REC8_HA_Rep1_perirangerPeaks vs ", DNAfamNames[i], " in pericentromeric regions (10-kb shift)"))
#   plot(lz_custom, main = paste0("REC8_HA_Rep1_perirangerPeaks vs ", DNAfamNames[i], " in pericentromeric regions (~", win, "-kb shift)"))
#   mtext(side = 3, at = 2, text = paste0("REC8_HA_Rep1_perirangerPeaks vs ", DNAfamNames[i], " in pericentromeric regions (~", win, "-kb shift)"))
#   dev.off()
# }
 

### RNA TEs

TEsRNAGR <- lapply(seq_along(RNAfamNames), function(x) {
  TEsRNA <- read.table(file = paste0(RNAdir, "TAIR10_Buisine_TEs_strand_tab_ann_", RNAfamNames[x], ".txt"), header = T)
  GRanges(seqnames = TEsRNA$chr, ranges = IRanges(start = TEsRNA$start, end = TEsRNA$end), strand = "*")
})

# Remove TEs located within masked regions
TEsRNAGRmasked <- lapply(seq_along(TEsRNAGR), function(x) {
  maskTEsOverlaps <- findOverlaps(mask, TEsRNAGR[[x]], ignore.strand = TRUE, select = "all")
  TEsRNAGR[[x]][-subjectHits(maskTEsOverlaps)]
})

# Perform permutation tests with randomized regions generated on a per chromosome basis;
# same per-chromosome number and size of regions in B as in A
set.seed(123)
ptPeaksTEsRNAPerChrom <- lapply(seq_along(TEsRNAGRmasked), function(x) {
  permTest(A = peaksGR, B = TEsRNAGRmasked[[x]], genome = genome, mask = mask,
           randomize.function = randomizeRegions,
           allow.overlaps = TRUE, per.chromosome = TRUE,
           evaluate.function = numOverlaps, count.once = TRUE,
           ntimes = 10000, mc.set.seed = FALSE, mc.cores = 20)
})

for(i in 1:length(ptPeaksTEsRNAPerChrom)) {
  assign(paste0(RNAfamNames[i]), ptPeaksTEsRNAPerChrom[[i]])
}
save(ptPeaksTEsRNAPerChrom,
     file = paste0(outDir,
                   "permTest_REC8_HA_Rep1_perirangerPeaks_vs_TEsRNA.RData"))

# Summarise results in a table
noOfFeatures <- NULL
expected <- NULL
observed <- NULL
pval <- NULL
zscore <- NULL
for(i in 1:length(ptPeaksTEsRNAPerChrom)) {
  noOfFeaturesi <- print(length(TEsRNAGRmasked[[i]]))
  noOfFeatures <- c(noOfFeatures, noOfFeaturesi)
  expectedi <- print(round(mean(ptPeaksTEsRNAPerChrom[[i]]$numOverlaps$permuted)))
  expected <- c(expected, expectedi)
  observedi <- print(ptPeaksTEsRNAPerChrom[[i]]$numOverlaps$observed)
  observed <- c(observed, observedi)
  pvali <- print(round(ptPeaksTEsRNAPerChrom[[i]]$numOverlaps$pval, 4))
  pval <- c(pval, pvali)
  zscorei <- print(round(ptPeaksTEsRNAPerChrom[[i]]$numOverlaps$zscore, 4))
  zscore <- c(zscore, zscorei)
}
ptPeaksTEsRNAPerChromDataFrame <- cbind(noOfFeatures, expected, observed, pval, zscore)
write.table(ptPeaksTEsRNAPerChromDataFrame,
            file = paste0(outDir,
                          "permTest_REC8_HA_Rep1_perirangerPeaks_vs_TEsRNA_DataFrame.txt"),
            sep = "\t", row.names = F)

# plot graphical summaries of results
for(i in 1:length(ptPeaksTEsRNAPerChrom)) {
  pdf(file = paste0(RNAplotDir, RNAfamNames[i], "_permTest_nperm10000_REC8_HA_Rep1_perirangerPeaks_perChrom.pdf"), width = 10, height = 7)
  plot(ptPeaksTEsRNAPerChrom[[i]], main = paste0("REC8_HA_Rep1_perirangerPeaks vs ", RNAfamNames[i], " in pericentromeric regions"), xlab = "Number of overlaps", ylab = "Density")
  dev.off()

  # Using the localZScore() function, evaluate whether the association between peaks and other is highly dependent on their exact position
  lz_1kb <- localZScore(pt = ptPeaksTEsRNAPerChrom[[i]], A = peaksGR, B = TEsRNAGRmasked[[i]],
                        window = 1000, step = 50, count.once = TRUE)
  lz_10kb <- localZScore(pt = ptPeaksTEsRNAPerChrom[[i]], A = peaksGR, B = TEsRNAGRmasked[[i]],
                         window = 10000, step = 500, count.once = TRUE)
  lz_custom <- localZScore(pt = ptPeaksTEsRNAPerChrom[[i]], A = peaksGR, B = TEsRNAGRmasked[[i]],
                           window = 10*mean(width(peaksGR)), step = mean(width(peaksGR))/2, count.once = TRUE)
  win <- as.character(round((10*mean(width(peaksGR)))/1000))
  step <- as.character(round(mean(width(peaksGR))/2))
  pdf(file = paste0(RNAplotDir, RNAfamNames[i], "_localZscore_permTest_nperm10000_REC8_HA_Rep1_perirangerPeaks_w1kb_s50bp_w10kb_s500bp_w", win ,"kb_s", step, "bp_perChrom.pdf"))
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  plot(lz_1kb, main = paste0("REC8_HA_Rep1_perirangerPeaks vs ", RNAfamNames[i], " in pericentromeric regions (1-kb shift)"))
  mtext(side = 3, at = 2, text = paste0("REC8_HA_Rep1_perirangerPeaks vs ", RNAfamNames[i], " in pericentromeric regions (1-kb shift)"))
  plot(lz_10kb, main = paste0("REC8_HA_Rep1_perirangerPeaks vs ", RNAfamNames[i], " in pericentromeric regions (10-kb shift)"))
  mtext(side = 3, at = 2, text = paste0("REC8_HA_Rep1_perirangerPeaks vs ", RNAfamNames[i], " in pericentromeric regions (10-kb shift)"))
  plot(lz_custom, main = paste0("REC8_HA_Rep1_perirangerPeaks vs ", RNAfamNames[i], " in pericentromeric regions (~", win, "-kb shift)"))
  mtext(side = 3, at = 2, text = paste0("REC8_HA_Rep1_perirangerPeaks vs ", RNAfamNames[i], " in pericentromeric regions (~", win, "-kb shift)"))
  dev.off()
}

