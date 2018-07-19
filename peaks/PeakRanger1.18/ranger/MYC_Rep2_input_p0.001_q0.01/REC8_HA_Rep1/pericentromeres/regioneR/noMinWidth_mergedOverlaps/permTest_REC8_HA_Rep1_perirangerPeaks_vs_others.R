#!/applications/R/R-3.3.2/bin/Rscript

# Use permutation test function in regioneR to determine if peaks overlap features
# of interest (e.g., other peaks, genes, TEs) more or less than expected by chance

# Usage on hydrogen node7:
# csmit -m 50G -c 48 "Rscript permTest_REC8_HA_Rep1_perirangerPeaks_vs_others.R"

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
plotDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/pericentromeres/regioneR/noMinWidth_mergedOverlaps/plots/"

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
#[1] 

load(paste0(inDir,
            "REC8_MYC_Rep1_perirangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData"))
REC8_MYC_Rep1GR <- perirangerPeaksGRmergedOverlaps
perirangerPeaksGRmergedOverlaps <- NULL
#REC8_MYC_Rep1GR <- REC8_MYC_Rep1GR[width(REC8_MYC_Rep1GR) <= 2000]
print(length(REC8_MYC_Rep1GR))
#[1] 

load(paste0(inDir,
            "REC8_MYC_Rep2_perirangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData"))
REC8_MYC_Rep2GR <- perirangerPeaksGRmergedOverlaps
perirangerPeaksGRmergedOverlaps <- NULL
#REC8_MYC_Rep2GR <- REC8_MYC_Rep2GR[width(REC8_MYC_Rep2GR) <= 2000]
print(length(REC8_MYC_Rep2GR))
#[1] 

# Others
load("/projects/ajt200/REC8_MSH4/nuc_peaks/log2ChIPinput/nucleR/trim/analysis_01/periPeaksSH99GRmerge.RData")
nucleRnucsGR <- periPeaksSH99GRmerge
periPeaksSH99GRmerge <- NULL
print(length(nucleRnucsGR))
#[1] 20205 

load("/projects/ajt200/BAM_masters/nucleosomes/WT/peaks/PeakRanger1.18/ranger/nakedDNA_untrimmed_input_p0.05_q0.05_l147/perirangerPeaksGRmerge_WT_nuc_p0.05_q0.05_noMinWidth.RData")
rangernucsGR <- perirangerPeaksGRmerge
perirangerPeaksGRmerge <- NULL
rangernucsGR <- rangernucsGR[width(rangernucsGR) <= 500]
print(length(rangernucsGR))
#[1] 12773

load("/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/perirangerPeaksGRmerge_RPI1_RPI8_idr0.05_noMinWidth.RData")
SPO11GR <- perirangerPeaksGRmerge
perirangerPeaksGRmerge <- NULL
print(length(SPO11GR))
#[1] 1574

load("/projects/ajt200/BAM_masters/H3K4me3/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/perirangerPeaksGRmerge_WT_H3K4me3_ChIP14_WT_H3K4me3_ChIP15_idr0.05_noMinWidth.RData")
H3K4me3GR <- perirangerPeaksGRmerge
perirangerPeaksGRmerge <- NULL
print(length(H3K4me3GR))
#[1] 1235

load("/projects/ajt200/BAM_masters/H3K9me2/WT/peaks/PeakRanger1.18/ranger/p0.05_q0.05/perirangerPeaksGRmerge_WT_H3K9me2_ChIP_p0.05_q0.05_noMinWidth.RData")
H3K9me2GR <- perirangerPeaksGRmerge
perirangerPeaksGRmerge <- NULL
print(length(H3K9me2GR))
#[1] 8054

load("/projects/ajt200/BAM_masters/H3K9me2/WT/peaks/PeakRanger1.18/BCP/p0.05/peribcpPeaksGRmerge_WT_H3K9me2_ChIP_p0.05_noMinWidth.RData")
H3K9me2GRbcp <- peribcpPeaksGRmerge
peribcpPeaksGRmerge <- NULL
print(length(H3K9me2GRbcp))
#[1] 344

load("/projects/ajt200/GBS_CO/HS_CU_080617/wt/COsGRcoords.RData")
COsGR <- COsGRcoords
print(length(COsGR))
#[1] 3320
# Remove COs located within pericentromeric regions
maskCOsOverlaps <- findOverlaps(mask, COsGR, ignore.strand = TRUE, select = "all")
COsGR <- COsGR[-subjectHits(maskCOsOverlaps)]
print(length(COsGR))
#[1] 868

genes <- read.table(file = "/projects/ajt200/TAIR10/representative_genes/representative_genes_uniq_fmt_strand.txt", header = T)
genesGR <- GRanges(seqnames = genes$chr, ranges = IRanges(start = genes$start, end = genes$end), strand = "*")
seqlevels(genesGR) <- sub("", "Chr", seqlevels(genesGR))
print(length(genesGR))
#[1] 27204
# Remove genes located within pericentromeric regions
maskgenesOverlaps <- findOverlaps(mask, genesGR, ignore.strand = TRUE, select = "all")
genesGR <- genesGR[-subjectHits(maskgenesOverlaps)]
print(length(genesGR))
#[1] 3723

genesGRprom <- GRanges(seqnames = genes$chr, ranges = IRanges(start = genes$start, end = genes$end), strand = genes$strand)
print(length(genesGRprom))
#[1] 27204
promotersGR <- promoters(genesGRprom, upstream = 500, downstream = 0)
seqlevels(promotersGR) <- sub("", "Chr", seqlevels(promotersGR))
strand(promotersGR) <- "*"
# Remove promoters located within pericentromeric regions
maskpromotersOverlaps <- findOverlaps(mask, promotersGR, ignore.strand = TRUE, select = "all")
promotersGR <- promotersGR[-subjectHits(maskpromotersOverlaps)]
print(length(promotersGR))
#[1] 3729

TSSdownstream500GR <- promoters(genesGRprom, upstream = 0, downstream = 500)
seqlevels(TSSdownstream500GR) <- sub("", "Chr", seqlevels(TSSdownstream500GR))
strand(TSSdownstream500GR) <- "*"
# Remove promoters located within pericentromeric regions
maskTSSdownstreamOverlaps <- findOverlaps(mask, TSSdownstream500GR, ignore.strand = TRUE, select = "all")
TSSdownstream500GR <- TSSdownstream500GR[-subjectHits(maskTSSdownstreamOverlaps)]
print(length(TSSdownstream500GR))
#[1] 3727

source("/projects/ajt200/Rfunctions/TTSplus.R")
genesGRterm <- GRanges(seqnames = genes$chr, ranges = IRanges(start = genes$start, end = genes$end), strand = genes$strand)
print(length(genesGRterm))
#[1] 27204
terminatorsGR <- TTSplus(genesGRterm, upstream = -1, downstream = 500)
seqlevels(terminatorsGR) <- sub("", "Chr", seqlevels(terminatorsGR))
strand(terminatorsGR) <- "*"
# Remove terminators located within pericentromeric regions
maskterminatorsOverlaps <- findOverlaps(mask, terminatorsGR, ignore.strand = TRUE, select = "all")
terminatorsGR <- terminatorsGR[-subjectHits(maskterminatorsOverlaps)]
print(length(terminatorsGR))
#[1] 3725

TTSupstream500GR <- TTSplus(genesGRterm, upstream = 499, downstream = 0)
seqlevels(TTSupstream500GR) <- sub("", "Chr", seqlevels(TTSupstream500GR))
strand(TTSupstream500GR) <- "*"
# Remove promoters located within pericentromeric regions
maskTTSupstreamOverlaps <- findOverlaps(mask, TTSupstream500GR, ignore.strand = TRUE, select = "all")
TTSupstream500GR <- TTSupstream500GR[-subjectHits(maskTTSupstreamOverlaps)]
print(length(TTSupstream500GR))
#[1] 3725

TEs <- read.table(file = "/projects/ajt200/TAIR10/TAIR10_Buisine_TEs_strand_tab_ann.txt", header=T)
TEsGR <- GRanges(seqnames = TEs$Chr, ranges = IRanges(start = TEs$start, end = TEs$end), strand = "*")
print(length(TEsGR))
#[1] 31189
# Remove TEs located within pericentromeric regions
maskTEsOverlaps <- findOverlaps(mask, TEsGR, ignore.strand = TRUE, select = "all")
TEsGR <- TEsGR[-subjectHits(maskTEsOverlaps)]
print(length(TEsGR))
#[1] 18171

otherNames <- c("REC8_MYC_Rep1GR", "REC8_MYC_Rep2GR",
                "nucleRnucsGR", "rangernucsGR", "SPO11GR",
                "H3K4me3GR", "H3K9me2GR", "H3K9me2GRbcp", "COsGR",
                "genesGR", "promotersGR", "terminatorsGR",
                "TSSdownstream500GR", "TTSupstream500GR", "TEsGR")

grl <- GRangesList("REC8_MYC_Rep1GR" = REC8_MYC_Rep1GR, "REC8_MYC_Rep2GR" = REC8_MYC_Rep2GR,
                   "nucleRnucsGR" = nucleRnucsGR, "rangernucsGR" = rangernucsGR, "SPO11GR" = SPO11GR,
                   "H3K4me3GR" = H3K4me3GR, "H3K9me2GR" = H3K9me2GR, "H3K9me2GRbcp" = H3K9me2GRbcp, "COsGR" = COsGR,
                   "genesGR" = genesGR, "promotersGR" = promotersGR, "terminatorsGR" = terminatorsGR,
                   "TSSdownstream500GR" = TSSdownstream500GR, "TTSupstream500GR" = TTSupstream500GR, "TEsGR" = TEsGR)


# Perform permutation tests with randomized regions generated on a per chromosome basis;
# same per-chromosome number and size of regions in B as in A
set.seed(123)
ptPeaksOtherPerChrom <- lapply(seq_along(grl), function(x) {
  permTest(A = peaksGR, B = grl[[x]], genome = genome, mask = mask,
           randomize.function = randomizeRegions,
           allow.overlaps = TRUE, per.chromosome = TRUE,
           evaluate.function = numOverlaps, count.once = TRUE,
           ntimes = 10000, mc.set.seed = FALSE, mc.cores = 20)
})

for(i in 1:length(ptPeaksOtherPerChrom)) {
  assign(paste0(otherNames[i]), ptPeaksOtherPerChrom[[i]])
}
save(ptPeaksOtherPerChrom,
     file = paste0(outDir,
                   "permTest_REC8_HA_Rep1_perirangerPeaks_vs_others.RData"))

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
                          "permTest_REC8_HA_Rep1_perirangerPeaks_vs_others_DataFrame.txt"),
            sep = "\t", row.names = F)

# plot graphical summaries of results
for(i in 1:length(ptPeaksOtherPerChrom)) {
  pdf(file = paste0(plotDir, otherNames[i], "_permTest_nperm10000_REC8_HA_Rep1_perirangerPeaks_perChrom.pdf"), width = 10, height = 7)
  plot(ptPeaksOtherPerChrom[[i]], main = paste0("REC8_HA_Rep1_perirangerPeaks vs ", otherNames[i], " in pericentromeric regions"), xlab = "Number of overlaps", ylab = "Density")
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
  pdf(file = paste0(plotDir, otherNames[i], "_localZscore_permTest_nperm10000_REC8_HA_Rep1_perirangerPeaks_w1kb_s50bp_w10kb_s500bp_w", win ,"kb_s", step, "bp_perChrom.pdf"))
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  plot(lz_1kb, main = paste0("REC8_HA_Rep1_perirangerPeaks vs ", otherNames[i], " in pericentromeric regions (1-kb shift)"))
  mtext(side = 3, at = 2, text = paste0("REC8_HA_Rep1_perirangerPeaks vs ", otherNames[i], " in pericentromeric regions (1-kb shift)"))
  plot(lz_10kb, main = paste0("REC8_HA_Rep1_perirangerPeaks vs ", otherNames[i], " in pericentromeric regions (10-kb shift)"))
  mtext(side = 3, at = 2, text = paste0("REC8_HA_Rep1_perirangerPeaks vs ", otherNames[i], " in pericentromeric regions (10-kb shift)"))
  plot(lz_custom, main = paste0("REC8_HA_Rep1_perirangerPeaks vs ", otherNames[i], " in pericentromeric regions (~", win, "-kb shift)"))
  mtext(side = 3, at = 2, text = paste0("REC8_HA_Rep1_perirangerPeaks vs ", otherNames[i], " in pericentromeric regions (~", win, "-kb shift)"))
  dev.off()
}


# Plot histograms of peri peak widths

histDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/pericentromeres/regioneR/noMinWidth_mergedOverlaps/hist/"

load(paste0(inDir,
            "REC8_HA_Rep1_perirangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData"))
REC8_HA_Rep1GR <- perirangerPeaksGRmergedOverlaps
perirangerPeaksGRmergedOverlaps <- NULL
load(paste0(inDir,
            "REC8_MYC_Rep1_perirangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData"))
REC8_MYC_Rep1GR <- perirangerPeaksGRmergedOverlaps
perirangerPeaksGRmergedOverlaps <- NULL
load(paste0(inDir,
            "REC8_MYC_Rep2_perirangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData"))
REC8_MYC_Rep2GR <- perirangerPeaksGRmergedOverlaps
perirangerPeaksGRmergedOverlaps <- NULL
load("/projects/ajt200/REC8_MSH4/nuc_peaks/log2ChIPinput/nucleR/trim/analysis_01/periPeaksSH99GRmerge.RData")
nucleRnucsGR <- periPeaksSH99GRmerge
periPeaksSH99GRmerge <- NULL
load("/projects/ajt200/BAM_masters/nucleosomes/WT/peaks/PeakRanger1.18/ranger/nakedDNA_untrimmed_input_p0.05_q0.05_l147/perirangerPeaksGRmerge_WT_nuc_p0.05_q0.05_noMinWidth.RData")
rangernucsGR <- perirangerPeaksGRmerge
perirangerPeaksGRmerge <- NULL
load("/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/perirangerPeaksGRmerge_RPI1_RPI8_idr0.05_noMinWidth.RData")
SPO11GR <- perirangerPeaksGRmerge
perirangerPeaksGRmerge <- NULL
load("/projects/ajt200/BAM_masters/H3K4me3/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/perirangerPeaksGRmerge_WT_H3K4me3_ChIP14_WT_H3K4me3_ChIP15_idr0.05_noMinWidth.RData")
H3K4me3GR <- perirangerPeaksGRmerge
perirangerPeaksGRmerge <- NULL
load("/projects/ajt200/BAM_masters/H3K9me2/WT/peaks/PeakRanger1.18/ranger/p0.05_q0.05/perirangerPeaksGRmerge_WT_H3K9me2_ChIP_p0.05_q0.05_noMinWidth.RData")
H3K9me2GR <- perirangerPeaksGRmerge
perirangerPeaksGRmerge <- NULL
load("/projects/ajt200/BAM_masters/H3K9me2/WT/peaks/PeakRanger1.18/BCP/p0.05/peribcpPeaksGRmerge_WT_H3K9me2_ChIP_p0.05_noMinWidth.RData")
H3K9me2GRbcp <- peribcpPeaksGRmerge
peribcpPeaksGRmerge <- NULL

grlHist <- GRangesList(REC8_HA_Rep1GR, REC8_MYC_Rep1GR, REC8_MYC_Rep2GR, nucleRnucsGR, rangernucsGR, SPO11GR, H3K4me3GR, H3K9me2GR, H3K9me2GRbcp)
namesHist <- c("REC8_HA_Rep1 (ranger)", "REC8_MYC_Rep1 (ranger)", "REC8_MYC_Rep2 (ranger)", "Nucleosome (nucleR)", "Nucleosome (ranger)", "SPO11-1-oligo (ranger)", "H3K4me3 (ranger)", "H3K9me2 (ranger)", "H3K9me2 (bcp)")

pdf(file = paste0(histDir, "hist_peak_widths_noMinWidth_noMaxWidth_brks250.pdf"), height = 12, width = 12)
par(mfrow = c(3, 3), mar = c(3, 3, 2, 2), mgp = c(2, 0.75, 0))
for(i in 1:length(grlHist)) {
  xname <- paste0(namesHist[i], " peri peak width (bp)")
  hist(width(grlHist[[i]]), breaks = 250, col = "grey60", border = NA, lwd = 2,
       xlab = paste0(namesHist[i], " peri peak width (bp)"), ylab = "Peaks", main = "", cex.lab = 1.25, cex.axis = 1.25)
  abline(v = mean(width(grlHist[[i]])), col = "black", lty = 2, lwd = 1)
}
dev.off()

