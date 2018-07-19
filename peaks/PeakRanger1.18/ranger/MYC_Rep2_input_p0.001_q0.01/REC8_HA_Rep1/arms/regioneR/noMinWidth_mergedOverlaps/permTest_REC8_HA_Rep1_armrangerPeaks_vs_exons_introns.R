#!/applications/R/R-3.3.2/bin/Rscript

# Use permutation test function in regioneR to determine if peaks overlap features
# of interest (e.g., other peaks, genes, TEs) more or less than expected by chance

# Usage on hydrogen node7:
# csmit -m 50G -c 48 "Rscript permTest_REC8_HA_Rep1_armrangerPeaks_vs_exons_introns.R"

library(regioneR)

# Genomic definitions
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chrStart <- c(1, 1, 1, 1, 1)
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))
mask <- toGRanges(data.frame(chrs, pericenStart, pericenEnd))

inDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/"
outDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/arms/regioneR/noMinWidth_mergedOverlaps/"
plotDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/arms/regioneR/noMinWidth_mergedOverlaps/plots/"

# REC8 peaks
load(paste0(inDir,
            "REC8_HA_Rep1_armrangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData"))
peaksGR <- armrangerPeaksGRmergedOverlaps
armrangerPeaksGRmergedOverlaps <- NULL
maskPeaksOverlaps <- findOverlaps(mask, peaksGR,
                                  ignore.strand = T,
                                  select = "all")
print("Peaks located within masked regions:")
print(maskPeaksOverlaps)
#peaksGR <- peaksGR[width(peaksGR) <= 2000]
print(length(peaksGR))
#[1] 63694

# Import exons as GRanges object
genes_exons <- read.table("/projects/ajt200/TAIR10/all_exons.txt",
                          header = T)
exons <- genes_exons[genes_exons$ge.ex == "exon",]
levels(exons$strand) <- c("-", "*", "+")
exonsGR <- GRanges(seqnames = exons$chr,
                   ranges = IRanges(start = exons$start, end = exons$end),
                   strand = "*",
                   gene_model = exons$name)
# Remove exons located within pericentromeric regions
maskexonsOverlaps <- findOverlaps(mask, exonsGR, ignore.strand = TRUE, select = "all")
exonsGR <- exonsGR[-subjectHits(maskexonsOverlaps)]
print(length(exonsGR))
print("***********exons***********")
print(exonsGR)

# Import introns tables and convert to GRanges object
intronsPlus <- read.table("/projects/ajt200/TAIR10/all_plus_introns.txt", header = T)
intronsMinus <- read.table("/projects/ajt200/TAIR10/all_minus_introns.txt", header = T)
intronsPlusGR <- GRanges(seqnames = paste0("Chr", intronsPlus$chr),
                         ranges = (IRanges(start = intronsPlus$all.intron.starts+1,
                                           end = intronsPlus$all.intron.stops-1)),
                         strand = "*", gene_model = intronsPlus$all.atgs)
intronsMinusGR <- GRanges(seqnames = paste0("Chr", intronsMinus$chr),
                          ranges = (IRanges(start = intronsMinus$all.intron.stops+1,
                                            end = intronsMinus$all.intron.starts-1)),
                          strand = "*", gene_model = intronsMinus$all.atgs)
intronsGR <- sort(append(intronsPlusGR, intronsMinusGR), by = ~ seqnames + start + end)
# Remove introns located within pericentromeric regions
maskintronsOverlaps <- findOverlaps(mask, intronsGR, ignore.strand = TRUE, select = "all")
intronsGR <- intronsGR[-subjectHits(maskintronsOverlaps)]
print(length(intronsGR))
print("***********introns***********")
print(intronsGR)

otherNames <- c("exonsGR", "intronsGR")

grl <- GRangesList("exonsGR" = exonsGR, "intronsGR" = intronsGR)

# Perform permutation tests with randomized regions generated on a per chromosome basis;
# same per-chromosome number and size of regions in B as in A
set.seed(123)
ptPeaksOtherPerChrom <- lapply(seq_along(grl), function(x) {
  permTest(A = peaksGR, B = grl[[x]], genome = genome, mask = mask,
           randomize.function = randomizeRegions,
           allow.overlaps = TRUE, per.chromosome = TRUE,
           evaluate.function = numOverlaps, count.once = TRUE,
           ntimes = 10000, mc.set.seed = FALSE, mc.cores = 48)
})

for(i in 1:length(ptPeaksOtherPerChrom)) {
  assign(paste0(otherNames[i]), ptPeaksOtherPerChrom[[i]])
}
save(ptPeaksOtherPerChrom,
     file = paste0(outDir,
                   "permTest_REC8_HA_Rep1_armrangerPeaks_vs_exons_introns.RData"))

# Summarise results in a table
noOfFeatures <- NULL
expected <- NULL
observed <- NULL
pval <- NULL
zscore <- NULL
for(i in 1:length(ptPeaksOtherPerChrom)) {
  noOfFeaturesi <- print(length(grl[[i]]))
  noOfFeatures <- c(noOfFeatures, noOfFeaturesi)
  expectedi <- print(round(mean(ptPeaksOtherPerChrom[[i]]$numOverlaps$permuted)))
  expected <- c(expected, expectedi)
  observedi <- print(ptPeaksOtherPerChrom[[i]]$numOverlaps$observed)
  observed <- c(observed, observedi)
  pvali <- print(round(ptPeaksOtherPerChrom[[i]]$numOverlaps$pval, 4))
  pval <- c(pval, pvali)
  zscorei <- print(round(ptPeaksOtherPerChrom[[i]]$numOverlaps$zscore, 4))
  zscore <- c(zscore, zscorei)
}
ptPeaksOtherPerChromDataFrame <- cbind(noOfFeatures, expected, observed, pval, zscore)
write.table(ptPeaksOtherPerChromDataFrame,
            file = paste0(outDir,
                          "permTest_REC8_HA_Rep1_armrangerPeaks_vs_exons_introns_DataFrame.txt"),
            sep = "\t", row.names = F)

# plot graphical summaries of results
for(i in 1:length(ptPeaksOtherPerChrom)) {
  pdf(file = paste0(plotDir, otherNames[i], "_permTest_nperm10000_REC8_HA_Rep1_armrangerPeaks_perChrom.pdf"), width = 10, height = 7)
  plot(ptPeaksOtherPerChrom[[i]], main = paste0("REC8_HA_Rep1_armrangerPeaks vs ", otherNames[i], " on chromosome arms"), xlab = "Number of overlaps", ylab = "Density")
  dev.off()

  # Using the localZScore() function, evaluate whether the association between peaks and other is highly dependent on their exact position
  lz_1kb <- localZScore(pt = ptPeaksOtherPerChrom[[i]], A = peaksGR, B = grl[[i]],
                        window = 1000, step = 50, count.once = TRUE)
  lz_10kb <- localZScore(pt = ptPeaksOtherPerChrom[[i]], A = peaksGR, B = grl[[i]],
                         window = 10000, step = 500, count.once = TRUE)
  lz_custom <- localZScore(pt = ptPeaksOtherPerChrom[[i]], A = peaksGR, B = grl[[i]],
                           window = 10*mean(width(peaksGR)), step = mean(width(peaksGR))/2, count.once = TRUE)
  win <- as.character(round((10*mean(width(peaksGR)))/1000))
  step <- as.character(round(mean(width(peaksGR))/2))
  pdf(file = paste0(plotDir, otherNames[i], "_localZscore_permTest_nperm10000_REC8_HA_Rep1_armrangerPeaks_w1kb_s50bp_w10kb_s500bp_w", win ,"kb_s", step, "bp_perChrom.pdf"))
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  plot(lz_1kb, main = paste0("REC8_HA_Rep1_armrangerPeaks vs ", otherNames[i], " on chromosome arms (1-kb shift)"))
  mtext(side = 3, at = 2, text = paste0("REC8_HA_Rep1_armrangerPeaks vs ", otherNames[i], " on chromosome arms (1-kb shift)"))
  plot(lz_10kb, main = paste0("REC8_HA_Rep1_armrangerPeaks vs ", otherNames[i], " on chromosome arms (10-kb shift)"))
  mtext(side = 3, at = 2, text = paste0("REC8_HA_Rep1_armrangerPeaks vs ", otherNames[i], " on chromosome arms (10-kb shift)"))
  plot(lz_custom, main = paste0("REC8_HA_Rep1_armrangerPeaks vs ", otherNames[i], " on chromosome arms (~", win, "-kb shift)"))
  mtext(side = 3, at = 2, text = paste0("REC8_HA_Rep1_armrangerPeaks vs ", otherNames[i], " on chromosome arms (~", win, "-kb shift)"))
  dev.off()
}


