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
# H3K9me2 in kss vs wild type

# Usage:
# DiffBind_wt_kss_H3K9me2_hypoCHG_DMRs.R sample_sheet_wt_kss_H3K9me2_hypoCHG_DMRs.csv

library(DiffBind)
print(packageVersion("DiffBind"))
#[1] ‘2.4.8’
print(packageVersion("DESeq2"))
#[1] ‘1.16.1’

#sampleSheet <- "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/DMRprofiles/kss_hypoCHG/DiffBind/sample_sheet_wt_kss_H3K9me2_hypoCHG_DMRs.csv"

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
#[1] "SampleID"   "Condition"  "Replicate"  "bamReads"   "ControlID" 
#[6] "bamControl" "Peaks"      "PeakCaller"
print(samples)
#          SampleID Condition Replicate                   bamReads
#1  wt_H3K9me2_ChIP        wt         1  reads/WT_H3K9me2_ChIP.bam
#2 kss_H3K9me2_ChIP       kss         1 reads/kss_H3K9me2_ChIP.bam
#          ControlID                  bamControl                   Peaks
#1  wt_H3K9me2_input  reads/WT_H3K9me2_input.bam suvh456_hypoCHG_DMR.bed
#2 kss_H3K9me2_input reads/kss_H3K9me2_input.bam suvh456_hypoCHG_DMR.bed
#  PeakCaller
#1        bed
#2        bed

dbaObj <- dba(sampleSheet = sampleSheet)
print(dbaObj)
#2 Samples, 15562 sites in matrix:
#                ID Condition Replicate Caller Intervals
#1  wt_H3K9me2_ChIP        wt         1    bed     15562
#2 kss_H3K9me2_ChIP       kss         1    bed     15562


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
#2 Samples, 15562 sites in matrix:
#                ID Condition Replicate Caller Intervals FRiP
#1  wt_H3K9me2_ChIP        wt         1 counts     15562 0.29
#2 kss_H3K9me2_ChIP       kss         1 counts     15562 0.10

# FRiP (Fraction of Reads in Peaks) - the proportion of reads for that sample that
# overlap a peak in the consensus peakset.
# FRiP can be used to indicate which samples show more enrichment overall

system(paste0("[ -d ./plots ] || mkdir ./plots"))
plotDir <- paste0(outDir, "plots/")

# Plot correlation heatmap to give a clustering of the samples using the
# cross-correlation of each row of the dbaObj$binding matrix (affinity scores)
pdf(paste0(plotDir, "sample_correlation_heatmap_binding_affinity_", featureName, ".pdf"))
dba.plotHeatmap(DBA = dbaObjCounts, margin = 26)
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
#2 Samples, 15562 sites in matrix:
#                ID Condition Replicate Caller Intervals FRiP
#1  wt_H3K9me2_ChIP        wt         1 counts     15562 0.29
#2 kss_H3K9me2_ChIP       kss         1 counts     15562 0.10
#
#1 Contrast:
#     Group1 Members1          Group2 Members2
#1 wild type        1 kyp suvh5 suvh6        1


# 3.4 Performing the differential analysis
# Run DESeq2 analysis using the default binding matrix of binding affinity scores
# Set bFullLibrarySize parameter to FALSE in dba.analyze() so that
# estimateSizeFactors is invoked, whereby the total number of reads mapping to
# target loci (the sum of each column)  is used to normalise by "library" size
# Can specify FDR threshold with "dbaObjContrast$config$th <- 0.05"
dbaObjContrast$config$th <- 0.1
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
#1: Some groups have no replicates. Results may be unreliable. 
#2: In checkForExperimentalReplicates(object, modelMatrix) :
#  same number of samples and coefficients to fit,
#  estimating dispersion by treating samples as replicates.
#  read the ?DESeq section on 'Experiments without replicates'
print(dbaObjAnalyze)
#2 Samples, 15562 sites in matrix:
#                ID Condition Replicate Caller Intervals FRiP
#1  wt_H3K9me2_ChIP        wt         1 counts     15562 0.29
#2 kss_H3K9me2_ChIP       kss         1 counts     15562 0.10
#
#1 Contrast:
#     Group1 Members1          Group2 Members2 DB.DESeq2
#1 wild type        1 kyp suvh5 suvh6        1       554

# dbaObjAnalyze is a DBA object
# 554 of the 15562 sites are identified as being significantly differentially bound (DB)
# using a threshold of FDR <= 0.1. 0 of 15562 when FDR <= 0.05

# Plot correlation heatmap to cluster the samples using only significantly DB sites
# based on cross-correlation of each row of the dbaObjAnalyze$binding matrix
# (affinity scores) that corresponds to a DB site
pdf(paste0(plotDir, "sample_correlation_heatmap_binding_affinity_DB_", featureName, "_FDR0.1.pdf"))
dba.plotHeatmap(DBA = dbaObjAnalyze, contrast = 1, margin = 20)
dev.off()

# 3.5 Retrieving the DB sites
dbaObjAnalyzeDB <- dba.report(dbaObjAnalyze, th = 0.1)
save(dbaObjAnalyzeDB,
     file = paste0(outDir, "DB_sites_", featureName, "_FDR0.1.RData"))
print(dbaObjAnalyzeDB)
#GRanges object with 554 ranges and 6 metadata columns:
#        seqnames               ranges strand |      Conc Conc_wild type
#           <Rle>            <IRanges>  <Rle> | <numeric>      <numeric>
#  15171     chr5 [15653301, 15653700]      * |      6.85          -1.24
#   2418     chr1 [19722601, 19723400]      * |      6.76          -2.24
#   2534     chr1 [21691501, 21691800]      * |      6.76          -2.24
#   6360     chr3 [ 6431701,  6432200]      * |      6.76          -2.24
#  12156     chr4 [ 8393201,  8393600]      * |      6.76          -2.24
#    ...      ...                  ...    ... .       ...            ...
#  14548     chr5 [13163701, 13165400]      * |      7.81           8.79
#   4152     chr2 [ 3734201,  3735100]      * |      7.81           8.79
#   8275     chr3 [14167601, 14168700]      * |      7.81           8.79
#   5776     chr2 [ 7232201,  7236100]      * |      9.15           10.1
#   1949     chr1 [16374901, 16377100]      * |      7.79           8.78
#        Conc_kyp suvh5 suvh6      Fold   p-value       FDR
#                   <numeric> <numeric> <numeric> <numeric>
#  15171                 7.85     -9.09  0.000286    0.0869
#   2418                 7.76     -9.99   0.00031    0.0869
#   2534                 7.76     -9.99   0.00031    0.0869
#   6360                 7.76     -9.99   0.00031    0.0869
#  12156                 7.76     -9.99   0.00031    0.0869
#    ...                  ...       ...       ...       ...
#  14548                 2.24      6.56   0.00536    0.0983
#   4152                 2.24      6.56   0.00536    0.0983
#   8275                 2.24      6.56   0.00537    0.0983
#   5776                 5.41      4.69   0.00537    0.0983
#   1949                 2.24      6.54   0.00546    0.0999
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
pdf(paste0(plotDir, "PCA_DB_", featureName, "_FDR0.1.pdf"))
dba.plotPCA(dbaObjAnalyze, contrast = 1,
            attributes = DBA_CONDITION, label = DBA_REPLICATE)
dev.off()

# 4.3 MA plots of contrasts
# Visualise the effect of normalization on the data, and datapoints identified as DB
pdf(paste0(plotDir, "MAplot_", featureName, "_FDR0.1.pdf"))
dba.plotMA(dbaObjAnalyze)
dev.off()

# Plot log concentrations for each condition in the contrast against one another
pdf(paste0(plotDir, "MAplot_bXY_", featureName, "_FDR0.1.pdf"))
dba.plotMA(dbaObjAnalyze, bXY = TRUE)
dev.off()

# 4.4 Volcano plots of contrasts
pdf(paste0(plotDir, "volcano_", featureName, "_FDR0.1.pdf"))
dba.plotVolcano(dbaObjAnalyze)
dev.off()

# 4.5 Boxplots
# Visualise how read distributions differ between classes of binding sites and sample groups
print(sum(dbaObjAnalyzeDB$Fold > 0))
#[1] 117
print(sum(dbaObjAnalyzeDB$Fold < 0))
#[1] 437
pdf(paste0(plotDir, "boxplot_read_distributions_DB_", featureName, "_FDR0.1.pdf"))
pvals <- dba.plotBox(dbaObjAnalyze)
dev.off()
# dba.plotBox returns a matrix of p-values (computed using a two-sided Wilcoxon
# 'Mann-Whitney' test, paired where appropriate) indicating which of these distributions
# are significantly different from another distribution
print(pvals)

# 4.6 Heatmaps
# Visualise patterns of binding affinity directly in the DB sites
# via a binding affinity heatmap
pdf(paste0(plotDir, "heatmap_binding_affinity_DB_", featureName, "_FDR0.1.pdf"))
corvals <- dba.plotHeatmap(dbaObjAnalyze, contrast = 1, correlations = FALSE,
                           margin = 26)
dev.off()
pdf(paste0(plotDir, "heatmap_binding_affinity_DB_", featureName, "_row_scaling_FDR0.1.pdf"))
corvals <- dba.plotHeatmap(dbaObjAnalyze, contrast = 1, correlations = FALSE,
                           margin = 26, scale = "row")
dev.off()


 
