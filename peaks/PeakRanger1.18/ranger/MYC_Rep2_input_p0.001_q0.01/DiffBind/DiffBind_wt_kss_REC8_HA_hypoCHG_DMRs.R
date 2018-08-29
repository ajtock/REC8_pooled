#!/applications/R/R-3.4.0/bin/Rscript

# R version 3.4.0 (invoke this version by specifying path at command line:
# /applications/R/R-3.4.0/bin/R 
# DiffBind version 2.4.8
# DESeq2 version 1.16.1
# Note that R version 3.4.0 or later is required for DESeq2 version 1.16 or later,
# in which "the log2 fold change shrinkage is no longer default for the DESeq and
# nbinomWaldTest functions". DESeq2 version 1.16 introduces "a separate function
# lfcShrink, which performs log2 fold change shrinkage for visualization and ranking
# of genes." (see https://support.bioconductor.org/p/95695/ and
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#changes)

# Use DiffBind to identify suvh456 hypoCHG DMRs that are differentially bound by
# REC8 in kss vs wild type

# Usage:
# DiffBind_wt_kss_REC8_HA_hypoCHG_DMRs.R sample_sheet_wt_kss_REC8_HA_hypoCHG_DMRs.csv

library(DiffBind)
print(packageVersion("DiffBind"))
#[1] ‘2.4.8’
print(packageVersion("DESeq2"))
#[1] ‘1.16.1’

#sampleSheet <- "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/DMRprofiles/kss_hypoCHG/DiffBind/sample_sheet_wt_kss_REC8_HA_hypoCHG_DMRs.csv"

args <- commandArgs(trailingOnly = T)
sampleSheet <- args[1]

sampleSheetBasename <- basename(sampleSheet)
len <- nchar(sampleSheetBasename)
featureName <- substr(sampleSheetBasename, 14, len-4)

setwd(dirname(sampleSheet))
outDir <- paste0(dirname(sampleSheet), "/")
print(outDir)
#[1] "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/DMRprofiles/kss_hypoCHG/DiffBind/"

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
#                    Peaks PeakCaller
#1 suvh456_hypoCHG_DMR.bed        bed
#2 suvh456_hypoCHG_DMR.bed        bed
#3 suvh456_hypoCHG_DMR.bed        bed

dbaObj <- dba(sampleSheet = sampleSheet)
print(dbaObj)
#3 Samples, 15562 sites in matrix:
#                ID Condition Replicate Caller Intervals
#1  wt_REC8_HA_Rep1        wt         1    bed     15562
#2  wt_REC8_HA_Rep2        wt         2    bed     15562
#3 kss_REC8_HA_Rep1       kss         1    bed     15562


# 3.2 Counting reads
# Calculate a binding matrix with scores based on read counts for every sample
# (affinity scores), rather than confidence scores for only those peaks called in
# a specific sample (occupancy scores)
# If working with "punctate", narrow peaks, use of the "summits" option in
# dba.count() is recommended to control for the widening effect of merging
# overlapping peaks. 
# This re-centres each peak around the point of greatest enrichment and keeps
# the peaks at a consistent width (in this case, with summits = 250, the peaks will
# be 500 bp, extending 250 bp upstream and downstream of the summit):
dbaObjCounts <- dba.count(dbaObj)
print(dbaObjCounts)
#3 Samples, 15562 sites in matrix:
#                ID Condition Replicate Caller Intervals FRiP
#1  wt_REC8_HA_Rep1        wt         1 counts     15562 0.13
#2  wt_REC8_HA_Rep2        wt         2 counts     15562 0.11
#3 kss_REC8_HA_Rep1       kss         1 counts     15562 0.11

# FRiP (Fraction of Reads in Peaks) - the proportion of reads for that sample that
# overlap a peak in the consensus peakset.
# FRiP can be used to indicate which samples show more enrichment overall

system(paste0("[ -d ./plots ] || mkdir ./plots"))
plotDir <- paste0(outDir, "plots/")

# Plot correlation heatmap to give a clustering of the samples using the
# cross-correlation of each row of the dbaObj$binding matrix (affinity scores)
pdf(paste0(plotDir, "sample_correlation_heatmap_binding_affinity_", featureName, ".pdf"))
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
#3 Samples, 15562 sites in matrix:
#                ID Condition Replicate Caller Intervals FRiP
#1  wt_REC8_HA_Rep1        wt         1 counts     15562 0.13
#2  wt_REC8_HA_Rep2        wt         2 counts     15562 0.11
#3 kss_REC8_HA_Rep1       kss         1 counts     15562 0.11
#
#1 Contrast:
#     Group1 Members1          Group2 Members2
#1 wild type        2 kyp suvh5 suvh6        1


# 3.4 Performing the differential analysis
# Run DESeq2 analysis using the default binding matrix of binding affinity scores
# Set bFullLibrarySize parameter to FALSE in dba.analyze() so that
# estimateSizeFactors is invoked, whereby the total number of reads mapping to
# target loci (the sum of each column)  is used to normalise by "library" size
# Can specify FDR threshold with "dbaObjContrast$config$th <- 0.05"
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
#3 Samples, 15562 sites in matrix:
#                ID Condition Replicate Caller Intervals FRiP
#1  wt_REC8_HA_Rep1        wt         1 counts     15562 0.13
#2  wt_REC8_HA_Rep2        wt         2 counts     15562 0.11
#3 kss_REC8_HA_Rep1       kss         1 counts     15562 0.11
#
#1 Contrast:
#     Group1 Members1          Group2 Members2 DB.DESeq2
#1 wild type        2 kyp suvh5 suvh6        1      2165

# dbaObjAnalyze is a DBA object
# 2165 of the 15562 sites are identified as being significantly differentially bound (DB)
# using the default threshold of FDR <= 0.05. 3124 of 15562 when FDR <= 0.1

# Plot correlation heatmap to cluster the samples using only significantly DB sites
# based on cross-correlation of each row of the dbaObjAnalyze$binding matrix
# (affinity scores) that corresponds to a DB site
pdf(paste0(plotDir, "sample_correlation_heatmap_binding_affinity_DB_", featureName, ".pdf"))
dba.plotHeatmap(DBA = dbaObjAnalyze, contrast = 1, margin = 20)
dev.off()

# 3.5 Retrieving the DB sites
dbaObjAnalyzeDB <- dba.report(dbaObjAnalyze, th = 0.05)
save(dbaObjAnalyzeDB,
     file = paste0(outDir, "DB_sites_", featureName, ".RData"))
print(dbaObjAnalyzeDB)
#GRanges object with 2165 ranges and 6 metadata columns:
#        seqnames               ranges strand |      Conc Conc_wild type
#           <Rle>            <IRanges>  <Rle> | <numeric>      <numeric>
#  12789     chr5 [ 9644501,  9644700]      * |      9.91           5.43
#  15339     chr5 [18207401, 18208500]      * |      8.33           8.87
#  15553     chr5 [26258601, 26259900]      * |      7.08          -0.35
#  15550     chr5 [25967801, 25968600]      * |      6.83          -0.35
#  15551     chr5 [26070001, 26070100]      * |      6.51          -0.35
#    ...      ...                  ...    ... .       ...            ...
#   5059     chr2 [ 5209401,  5210800]      * |      7.45           7.78
#   3368     chr2 [ 2074201,  2074900]      * |      8.46           8.08
#  11766     chr4 [ 5474201,  5474300]      * |      5.16            5.7
#  12352     chr4 [16054601, 16054800]      * |      5.12           5.68
#   9819     chr4 [ 1755801,  1755900]      * |      4.97           5.53
#        Conc_kyp suvh5 suvh6      Fold   p-value       FDR
#                   <numeric> <numeric> <numeric> <numeric>
#  12789                11.45     -6.03 1.74e-111 2.71e-107
#  15339                 5.03      3.84  4.55e-16  3.54e-12
#  15553                 8.66     -9.02  9.37e-15  4.86e-11
#  15550                 8.41     -8.76  1.63e-13  6.35e-10
#  15551                 8.09     -8.44  5.48e-12  1.71e-08
#    ...                  ...       ...       ...       ...
#   5059                 6.42      1.37   0.00685    0.0493
#   3368                 9.01     -0.93   0.00685    0.0493
#  11766                 1.64      4.06   0.00687    0.0494
#  12352                 0.64      5.05   0.00689    0.0496
#   9819                 0.64      4.89   0.00695    0.0499
#  -------
#  seqinfo: 5 sequences from an unspecified genome; no seqlengths


## 4. Example: Plotting

# 4.2 PCA plots
# Generate PCA plot using binding affinity data for all sites in the consensus peakset
pdf(paste0(plotDir, "PCA_", featureName, ".pdf"))
dba.plotPCA(dbaObjAnalyze, attributes = DBA_CONDITION, label = DBA_REPLICATE)
dev.off()

# Generate PCA plot using binding affinity data for only the DB sites
# Differential analysis identifies sites that can be used to separate the sample
# groups along the second component, indicating why the hierarchical clustering
# shown in the corresponding heatmap is imperfect
pdf(paste0(plotDir, "PCA_DB_", featureName, ".pdf"))
dba.plotPCA(dbaObjAnalyze, contrast = 1,
            attributes = DBA_CONDITION, label = DBA_REPLICATE)
dev.off()

# 4.3 MA plots of contrasts
# Visualise the effect of normalization on the data, and datapoints identified as DB
pdf(paste0(plotDir, "MAplot_", featureName, ".pdf"))
dba.plotMA(dbaObjAnalyze)
dev.off()

# Plot log concentrations for each condition in the contrast against one another
pdf(paste0(plotDir, "MAplot_bXY_", featureName, ".pdf"))
dba.plotMA(dbaObjAnalyze, bXY = TRUE)
dev.off()

# 4.4 Volcano plots of contrasts
pdf(paste0(plotDir, "volcano_", featureName, ".pdf"))
dba.plotVolcano(dbaObjAnalyze)
dev.off()

# 4.5 Boxplots
# Visualise how read distributions differ between classes of binding sites and sample groups
print(sum(dbaObjAnalyzeDB$Fold > 0))
#[1] 1353
print(sum(dbaObjAnalyzeDB$Fold < 0))
#[1] 812
pdf(paste0(plotDir, "boxplot_read_distributions_DB_", featureName, ".pdf"))
pvals <- dba.plotBox(dbaObjAnalyze)
dev.off()
# dba.plotBox returns a matrix of p-values (computed using a two-sided Wilcoxon
# 'Mann-Whitney' test, paired where appropriate) indicating which of these distributions
# are significantly different from another distribution
print(pvals)

# 4.6 Heatmaps
# Visualise patterns of binding affinity directly in the DB sites
# via a binding affinity heatmap
pdf(paste0(plotDir, "heatmap_binding_affinity_DB_", featureName, ".pdf"))
corvals <- dba.plotHeatmap(dbaObjAnalyze, contrast = 1, correlations = FALSE,
                           margin = 20)
dev.off()
pdf(paste0(plotDir, "heatmap_binding_affinity_DB_", featureName, "_row_scaling.pdf"))
corvals <- dba.plotHeatmap(dbaObjAnalyze, contrast = 1, correlations = FALSE,
                           margin = 20, scale = "row")
dev.off()


 
