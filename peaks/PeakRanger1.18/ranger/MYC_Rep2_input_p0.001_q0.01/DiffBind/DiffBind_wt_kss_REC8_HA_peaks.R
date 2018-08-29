#!/applications/R/R-3.5.0/bin/Rscript

# R version 3.5.0 (invoke this version by specifying path at command line:
# /applications/R/R-3.5.0/bin/R 
# DiffBind version 2.8.0
# DESeq2 version 1.16.1
# Note that R version 3.4.0 or later is required for DESeq2 version 1.16 or later,
# in which "the log2 fold change shrinkage is no longer default for the DESeq and
# nbinomWaldTest functions". DESeq2 version 1.16 introduces "a separate function
# lfcShrink, which performs log2 fold change shrinkage for visualization and ranking
# of genes." (see https://support.bioconductor.org/p/95695/ and
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#changes)

# Use DiffBind to identify REC8-HA peaks that are differentially bound by REC8
# in kss vs wild type

# Usage:
# DiffBind_wt_kss_REC8_HA_peaks.R sample_sheet_wt_kss_REC8_HA_peaks.csv 0.05

library(DiffBind)
print(packageVersion("DiffBind"))
#[1] ‘2.8.0’
print(packageVersion("DESeq2"))
#[1] ‘1.16.1’

#sampleSheet <- "/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/DiffBind/sample_sheet_wt_kss_REC8_HA_peaks.csv"
#FDRnum <- 0.05
#FDRchar <- "0.05"

args <- commandArgs(trailingOnly = T)
sampleSheet <- args[1]
FDRnum <- as.numeric(args[2])
FDRchar <- as.character(args[2])

sampleSheetBasename <- basename(sampleSheet)
len <- nchar(sampleSheetBasename)
featureName <- substr(sampleSheetBasename, 14, len-4)

setwd(dirname(sampleSheet))
outDir <- paste0(dirname(sampleSheet), "/")
print(outDir)
#[1] "/projects/ajt200/BAM_masters/SPO11-oligo/arp6/gene_promoters/DiffBind/"

# Read in sampleSheet
samples <- read.csv(sampleSheet)
print(names(samples))
#[1] "SampleID"   "Condition"  "Replicate"  "bamReads"   "Peaks"     
#[6] "PeakCaller"
print(samples)
#          SampleID Condition Replicate                        bamReads
#1  wt_REC8_HA_Rep1        wt         1  reads/WT_REC8_HA_Rep1_ChIP.bam
#2  wt_REC8_HA_Rep2        wt         2  reads/WT_REC8_HA_Rep2_ChIP.bam
#3 kss_REC8_HA_Rep1       kss         1 reads/kss_REC8_HA_Rep1_ChIP.bam
#                        Peaks PeakCaller
#1  WT_REC8_HA_Rep1.narrowPeak     narrow
#2  WT_REC8_HA_Rep2.narrowPeak     narrow
#3 kss_REC8_HA_Rep1.narrowPeak     narrow

dbaObj <- dba(sampleSheet = sampleSheet)
print(dbaObj)
#3 Samples, 61280 sites in matrix (117997 total):
#                ID Condition Replicate Caller Intervals
#1  wt_REC8_HA_Rep1        wt         1 narrow     88007
#2  wt_REC8_HA_Rep2        wt         2 narrow     82952
#3 kss_REC8_HA_Rep1       kss         1 narrow    103945


# 3.2 Counting reads
# Calculate a binding matrix with scores based on read counts for every sample
# (affinity scores), rather than confidence scores for only those peaks called in
# a specific sample (occupancy scores)
# If working with "punctate", narrow peaks, use of the "summits" option is
# recommended to control for the widening effect of merging overlapping peaks. 
# This re-centres each peak around the point of greatest enrichment and keeps
# the peaks at a consistent width (in this case, with summits = 250, the peaks will
# be 500 bp, extending 250 bp upstream and downstream of the summit):
dbaObjCounts <- dba.count(dbaObj, summits = 250)
print(dbaObjCounts)
#3 Samples, 60236 sites in matrix:
#                ID Condition Replicate Caller Intervals FRiP
#1  wt_REC8_HA_Rep1        wt         1 counts     60236 0.39
#2  wt_REC8_HA_Rep2        wt         2 counts     60236 0.37
#3 kss_REC8_HA_Rep1       kss         1 counts     60236 0.40

# FRiP (Fraction of Reads in Peaks) - the proportion of reads for that sample that
# overlap a peak in the consensus peakset.
# FRiP can be used to indicate which samples show more enrichment overall

system(paste0("[ -d ./plots ] || mkdir ./plots"))
plotDir <- paste0(outDir, "plots/")

# Plot correlation heatmap to give a clustering of the samples using the
# cross-correlation of each row of the dbaObj$binding matrix (affinity scores)
pdf(paste0(plotDir, "sample_correlation_heatmap_binding_affinity_", featureName, "_" , FDRchar, ".pdf"))
dba.plotHeatmap(DBA = dbaObjCounts, margin = 20)
dev.off()

# 3.3 Establishing a contrast
# Specify which treatments fall in which groups
#dbaObjContrast <- dba.contrast(DBA = dbaObjCounts,
#                               categories = DBA_CONDITION,
#                               minMembers = 2)
dbaObjContrast <- dba.contrast(DBA = dbaObjCounts,
                               group1 = dbaObjCounts$masks$wt,
                               group2 = dbaObjCounts$masks$kss,
                               name1 = "wild type",
                               name2 = "kyp suvh5 suvh6")
print(dbaObjContrast)
#3 Samples, 60236 sites in matrix:
#                ID Condition Replicate Caller Intervals FRiP
#1  wt_REC8_HA_Rep1        wt         1 counts     60236 0.39
#2  wt_REC8_HA_Rep2        wt         2 counts     60236 0.37
#3 kss_REC8_HA_Rep1       kss         1 counts     60236 0.40
#
#1 Contrast:
#     Group1 Members1          Group2 Members2
#1 wild type        2 kyp suvh5 suvh6        1


# 3.4 Performing the differential analysis
# Run DESeq2 analysis using the default binding matrix of binding affinity scores
# Set bFullLibrarySize parameter to FALSE in dba.analyze() so that
# estimateSizeFactors is invoked, whereby the total number of reads mapping to
# target loci (the sum of each column)  is used to normalise by "library" size
# Can specify FDR threshold with, e.g.:
#dbaObjContrast$config$th <- 0.05
dbaObjAnalyze <- dba.analyze(dbaObjContrast,
                             method = DBA_DESEQ2,
                             bFullLibrarySize = FALSE,
# bTagwise - logical indicating if dispersion should be calculated on a tagwise
# (or per-condition) basis. If there are only a very few members of each group
# in a contrast (e.g. no replicates), this should be set to FALSE.
                             bTagwise = FALSE)
#converting counts to integer mode
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#Warning message:
#Some groups have no replicates. Results may be unreliable. 

print(dbaObjAnalyze)
#                ID Condition Replicate Caller Intervals FRiP
#1  wt_REC8_HA_Rep1        wt         1 counts     60236 0.39
#2  wt_REC8_HA_Rep2        wt         2 counts     60236 0.37
#3 kss_REC8_HA_Rep1       kss         1 counts     60236 0.40
#
#1 Contrast:
#     Group1 Members1          Group2 Members2 DB.DESeq2
#1 wild type        2 kyp suvh5 suvh6        1      7829

# dbaObjAnalyze is a DBA object
# 7829 of the 13492 sites are identified as being significantly differentially bound (DB)
# using the default threshold of FDR <= 0.05

# Plot correlation heatmap to cluster the samples using only significantly DB sites
# based on cross-correlation of each row of the dbaObjAnalyze$binding matrix
# (affinity scores) that corresponds to a DB site
pdf(paste0(plotDir, "sample_correlation_heatmap_binding_affinity_DB_", featureName, "_" , FDRchar, ".pdf"))
dba.plotHeatmap(DBA = dbaObjAnalyze, contrast = 1, margin = 20)
dev.off()

# 3.5 Retrieving the DB sites
dbaObjAnalyzeDB <- dba.report(dbaObjAnalyze, th = 0.05)
save(dbaObjAnalyzeDB,
     file = paste0(outDir, "DB_sites_", featureName, "_", FDRchar, ".RData"))
print(dbaObjAnalyzeDB)
#GRanges object with 7829 ranges and 6 metadata columns:
#        seqnames               ranges strand |      Conc Conc_wild type
#           <Rle>            <IRanges>  <Rle> | <numeric>      <numeric>
#  40126     chr4 [ 5541579,  5542079]      * |      9.48           5.97
#  12726     chr1 [25170753, 25171253]      * |      9.95           7.08
#  51664     chr5 [ 9644192,  9644692]      * |     10.07              7
#  59895     chr5 [26314276, 26314776]      * |      7.72          -0.27
#  59948     chr5 [26400895, 26401859]      * |      7.69          -0.27
#    ...      ...                  ...    ... .       ...            ...
#  15105     chr1 [29822302, 29822802]      * |      7.46           7.82
#  23376     chr2 [15482725, 15483225]      * |      7.24           7.61
#  38692     chr4 [ 2561577,  2562077]      * |      7.36           7.71
#  23971     chr2 [16580440, 16580940]      * |      7.89           8.19
#  52844     chr5 [12109089, 12109589]      * |      7.93           7.44
#        Conc_kyp suvh5 suvh6      Fold   p-value       FDR
#                   <numeric> <numeric> <numeric> <numeric>
#  40126                10.97        -5  1.14e-22  6.06e-18
#  12726                11.39     -4.32  2.01e-22  6.06e-18
#  51664                11.54     -4.54  2.17e-20  4.35e-16
#  59895                  9.3     -9.57  1.47e-17  2.22e-13
#  59948                 9.27     -9.53  2.26e-17  2.72e-13
#    ...                  ...       ...       ...       ...
#  15105                 6.26      1.56   0.00648    0.0499
#  23376                 6.01       1.6   0.00649    0.0499
#  38692                 6.18      1.53   0.00649    0.0499
#  23971                 6.98      1.21   0.00649    0.0499
#  52844                 8.58     -1.14    0.0065      0.05
#  -------
#  seqinfo: 5 sequences from an unspecified genome; no seqlengths
dbaObjAnalyzeDB_wt_exceeds_kss <- dbaObjAnalyzeDB[dbaObjAnalyzeDB$Fold > 0]
dbaObjAnalyzeDB_kss_exceeds_wt <- dbaObjAnalyzeDB[dbaObjAnalyzeDB$Fold < 0]
print(dbaObjAnalyzeDB_wt_exceeds_kss)
#GRanges object with 5048 ranges and 6 metadata columns:
#        seqnames               ranges strand |      Conc Conc_wild type
#           <Rle>            <IRanges>  <Rle> | <numeric>      <numeric>
#  39583     chr4 [ 4448917,  4449417]      * |      8.15           8.69
#  53662     chr5 [13803164, 13803664]      * |         8           8.55
#  56526     chr5 [19576720, 19577220]      * |      7.89           8.44
#  53275     chr5 [13020977, 13021477]      * |      8.02           8.56
#  44534     chr4 [14152415, 14152915]      * |      8.03           8.56
#    ...      ...                  ...    ... .       ...            ...
#   8542     chr1 [17163033, 17163533]      * |      7.62           7.94
#  15105     chr1 [29822302, 29822802]      * |      7.46           7.82
#  23376     chr2 [15482725, 15483225]      * |      7.24           7.61
#  38692     chr4 [ 2561577,  2562077]      * |      7.36           7.71
#  23971     chr2 [16580440, 16580940]      * |      7.89           8.19
#        Conc_kyp suvh5 suvh6      Fold   p-value       FDR
#                   <numeric> <numeric> <numeric> <numeric>
#  39583                  4.7      3.99  7.59e-15  2.55e-11
#  53662                 4.36      4.19  3.16e-14  6.35e-11
#  56526                 4.15      4.28  1.01e-13  1.48e-10
#  53275                 4.54      4.02  1.83e-13  2.35e-10
#  44534                 4.84      3.71  2.46e-13   2.9e-10
#    ...                  ...       ...       ...       ...
#   8542                  6.6      1.34   0.00648    0.0499
#  15105                 6.26      1.56   0.00648    0.0499
#  23376                 6.01       1.6   0.00649    0.0499
#  38692                 6.18      1.53   0.00649    0.0499
#  23971                 6.98      1.21   0.00649    0.0499
#  -------
#  seqinfo: 5 sequences from an unspecified genome; no seqlengths
print(dbaObjAnalyzeDB_kss_exceeds_wt)
#GRanges object with 2781 ranges and 6 metadata columns:
#        seqnames               ranges strand |      Conc Conc_wild type
#           <Rle>            <IRanges>  <Rle> | <numeric>      <numeric>
#  40126     chr4 [ 5541579,  5542079]      * |      9.48           5.97
#  12726     chr1 [25170753, 25171253]      * |      9.95           7.08
#  51664     chr5 [ 9644192,  9644692]      * |     10.07              7
#  59895     chr5 [26314276, 26314776]      * |      7.72          -0.27
#  59948     chr5 [26400895, 26401859]      * |      7.69          -0.27
#    ...      ...                  ...    ... .       ...            ...
#  14211     chr1 [28058973, 28059473]      * |      7.67           7.12
#  25342     chr2 [19222083, 19222583]      * |      7.88           7.36
#  37070     chr3 [22870446, 22870946]      * |      7.76           7.24
#  20516     chr2 [10189438, 10190415]      * |      8.12           7.63
#  52844     chr5 [12109089, 12109589]      * |      7.93           7.44
#        Conc_kyp suvh5 suvh6      Fold   p-value       FDR
#                   <numeric> <numeric> <numeric> <numeric>
#  40126                10.97        -5  1.14e-22  6.06e-18
#  12726                11.39     -4.32  2.01e-22  6.06e-18
#  51664                11.54     -4.54  2.17e-20  4.35e-16
#  59895                  9.3     -9.57  1.47e-17  2.22e-13
#  59948                 9.27     -9.53  2.26e-17  2.72e-13
#    ...                  ...       ...       ...       ...
#  14211                 8.38     -1.26   0.00647    0.0498
#  25342                 8.57      -1.2   0.00647    0.0499
#  37070                 8.45     -1.22   0.00648    0.0499
#  20516                 8.77     -1.14   0.00648    0.0499
#  52844                 8.58     -1.14    0.0065      0.05
#  -------
#  seqinfo: 5 sequences from an unspecified genome; no seqlengths
save(dbaObjAnalyzeDB_wt_exceeds_kss,
     file = paste0(outDir, "DB_sites_wt_exceeds_kss_", featureName, "_", FDRchar, ".RData"))
save(dbaObjAnalyzeDB_kss_exceeds_wt,
     file = paste0(outDir, "DB_sites_kss_exceeds_wt_", featureName, "_", FDRchar, ".RData"))


## 4. Example: Plotting

# 4.2 PCA plots
# Generate PCA plot using binding affinity data for all sites in the consensus peakset
pdf(paste0(plotDir, "PCA_", featureName, "_" , FDRchar, ".pdf"))
dba.plotPCA(dbaObjAnalyze, attributes = DBA_CONDITION, label = DBA_REPLICATE)
dev.off()

# Generate PCA plot using binding affinity data for only the DB sites
# Differential analysis identifies sites that can be used to separate the sample
# groups along the second component, indicating why the hierarchical clustering
# shown in the corresponding heatmap is imperfect
pdf(paste0(plotDir, "PCA_DB_", featureName, "_" , FDRchar, ".pdf"))
dba.plotPCA(dbaObjAnalyze, contrast = 1,
            attributes = DBA_CONDITION, label = DBA_REPLICATE)
dev.off()

# 4.3 MA plots of contrasts
# Visualise the effect of normalization on the data, and datapoints identified as DB
pdf(paste0(plotDir, "MAplot_", featureName, "_" , FDRchar, ".pdf"))
dba.plotMA(dbaObjAnalyze)
dev.off()

# Plot log concentrations for each condition in the contrast against one another
pdf(paste0(plotDir, "MAplot_bXY_", featureName, "_" , FDRchar, ".pdf"))
dba.plotMA(dbaObjAnalyze, bXY = TRUE)
dev.off()

# 4.4 Volcano plots of contrasts
pdf(paste0(plotDir, "volcano_", featureName, "_" , FDRchar, ".pdf"))
dba.plotVolcano(dbaObjAnalyze)
dev.off()

# 4.5 Boxplots
# Visualise how read distributions differ between classes of binding sites and sample groups
print(sum(dbaObjAnalyzeDB$Fold < 0))
#[1] 5048
print(sum(dbaObjAnalyzeDB$Fold > 0))
#[1] 2781
pdf(paste0(plotDir, "boxplot_read_distributions_DB_", featureName, "_" , FDRchar, ".pdf"))
pvals <- dba.plotBox(dbaObjAnalyze)
#Error in apply(report[, idx], 1, mean) : 
#  dim(X) must have a positive length
dev.off()
# dba.plotBox returns a matrix of p-values (computed using a two-sided Wilcoxon
# 'Mann-Whitney' test, paired where appropriate) indicating which of these distributions
# are significantly different from another distribution
print(pvals)

# 4.6 Heatmaps
# Visualise patterns of binding affinity directly in the DB sites
# via a binding affinity heatmap
pdf(paste0(plotDir, "heatmap_binding_affinity_DB_", featureName, "_" , FDRchar, ".pdf"))
corvals <- dba.plotHeatmap(dbaObjAnalyze, contrast = 1, correlations = FALSE,
                           margin = 20)
dev.off()
pdf(paste0(plotDir, "heatmap_binding_affinity_DB_", featureName, "_row_scaling_" , FDRchar, ".pdf"))
corvals <- dba.plotHeatmap(dbaObjAnalyze, contrast = 1, correlations = FALSE,
                           margin = 20, scale = "row")
dev.off()


 
