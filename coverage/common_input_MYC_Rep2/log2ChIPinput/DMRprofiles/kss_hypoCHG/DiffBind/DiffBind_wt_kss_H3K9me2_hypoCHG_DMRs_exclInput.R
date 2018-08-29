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
# DiffBind_wt_kss_H3K9me2_hypoCHG_DMRs_exclInput.R sample_sheet_wt_kss_H3K9me2_hypoCHG_DMRs_exclInput.csv

library(DiffBind)
print(packageVersion("DiffBind"))
#[1] ‘2.4.8’
print(packageVersion("DESeq2"))
#[1] ‘1.16.1’

#sampleSheet <- "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/DMRprofiles/kss_hypoCHG/DiffBind/sample_sheet_wt_kss_H3K9me2_hypoCHG_DMRs_exclInput.csv"

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
#          SampleID Condition Replicate                   bamReads
#1  wt_H3K9me2_ChIP        wt         1  reads/WT_H3K9me2_ChIP.bam
#2 kss_H3K9me2_ChIP       kss         1 reads/kss_H3K9me2_ChIP.bam
#                    Peaks PeakCaller
#1 suvh456_hypoCHG_DMR.bed        bed
#2 suvh456_hypoCHG_DMR.bed        bed

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
pdf(paste0(plotDir, "sample_correlation_heatmap_binding_affinity_", featureName, "_exclInput.pdf"))
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
#1 wild type        1 kyp suvh5 suvh6        1       0

# dbaObjAnalyze is a DBA object
# 0 of the 15562 sites are identified as being significantly differentially bound (DB)
# using a threshold of FDR <= 0.1.

