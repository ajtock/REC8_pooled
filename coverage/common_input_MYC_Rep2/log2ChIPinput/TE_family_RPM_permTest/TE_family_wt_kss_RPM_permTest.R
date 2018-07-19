#!/applications/R/R-3.3.2/bin/Rscript

# Calculate mean z-score standardised log2-transformed library-size-normalised
# wt and kss coverage levels at TEs (mean over TE width) within each family
# To do with separate script:
# Calculate correlation coefficients for levels in wt vs kss and for
# dataset vs dataset

# Usage via Condor submission system on node7:
# csmit -m 500G -c 48 "Rscript TE_family_wt_kss_RPM_permTest_TEST.R /projects/ajt200/BAM_masters/SPO11-oligo/WT/WT_SPO11-oligo_RPI1_k10_bt2_mapped_lowmiss_unique_both_sort_rmdup_new_sort.bam wt_SPO11_1_oligos_RPI1 10000 0.0001"

library(regioneR)
library(GenomicRanges)
library(GenomicAlignments)
library(parallel)
library(plotrix)
#library(doParallel)
#registerDoParallel(cores=7)
#print("Currently registered parallel backend name, version and cores")
#print(getDoParName())
#print(getDoParVersion())
#print(getDoParWorkers())

#libPath <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/WT_SPO11-oligo_RPI1_k10_bt2_mapped_lowmiss_unique_both_sort_rmdup_new_sort.bam/"
#libName <- "SPO11_1_oligos_RPI1"
#perms <- as.numeric(1000)
#minPval <- as.numeric(0.001)

args <- commandArgs(trailingOnly = T)
libPath <- as.character(args[1])
libName <- as.character(args[2])
# Number of permutations (randomisations) to perform
perms <- as.numeric(args[3])
# Corresponding minimum P-value
minPval <- as.numeric(args[4])

outDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/TE_family_wt_kss_libNorm_read_counts_permTest/"
plotDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/TE_family_wt_kss_libNorm_read_counts_permTest/plots/"

chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chrStart <- c(1, 1, 1, 1, 1)
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))

DNAfamNames <- c("dna", "heli", "ptmari", "mudr", "enspm", "hat", "harbinger")
RNAfamNames <- c("rna", "gypsy", "copia", "linel1", "sine")
DNAdir <- "/projects/ajt200/TAIR10/TE_classes/DNA/"
RNAdir <- "/projects/ajt200/TAIR10/TE_classes/RNA/"

### DNA transposons
TEsDNAGR <- lapply(seq_along(DNAfamNames), function(x) {
  TEsDNA <- read.table(file = paste0(DNAdir, "TAIR10_Buisine_TEs_strand_tab_ann_", DNAfamNames[x], ".txt"), header = T)
  GRanges(seqnames = TEsDNA$chr, ranges = IRanges(start = TEsDNA$start, end = TEsDNA$end), strand = "*")
})

### RNA transposons
TEsRNAGR <- lapply(seq_along(RNAfamNames), function(x) {
  TEsRNA <- read.table(file = paste0(RNAdir, "TAIR10_Buisine_TEs_strand_tab_ann_", RNAfamNames[x], ".txt"), header = T)
  GRanges(seqnames = TEsRNA$chr, ranges = IRanges(start = TEsRNA$start, end = TEsRNA$end), strand = "*")
})

# Load BAM and create RangedData object
lib_readGAlignments <- readGAlignments(libPath)
lib_ranges <- ranges(lib_readGAlignments)
lib_chrs <- as.character(as.data.frame(seqnames(lib_readGAlignments))[,1])
lib_strand <- as.character(strand(lib_readGAlignments))
lib_ranged <- RangedData(space = lib_chrs,
                         ranges = lib_ranges,
                         strand = "*")
lib_ranged <- lib_ranged[space(lib_ranged) != "chloroplast" & space(lib_ranged) != "mitochondria",]
names(lib_ranged) <- sub("", "Chr", names(lib_ranged))
save(lib_ranged, file = paste0(outDir, libName, "_RangedData.RData"))

# Calculate library size
lib_size <- length(space(lib_ranged))

# Calculate "per million" scaling factor
RPM_scaling_factor <- lib_size/1e+06

# Convert library ranged data into GRanges object
lib_rangedGR <- as(lib_ranged, "GRanges")

# Set class for permutation test results object
setClass("permTest",
         representation(Pval = "numeric",
                        feature_RPM_sum = "numeric",
                        ranLoc_RPM_sum = "numeric",
                        Expected_RPM = "numeric",
                        log2_Observed_Expected = "numeric"))

# Function to define perms sets of random loci based on target features and
# calculate RPM at each target and random loci 
locRPMcalc <- function(targets, targetsName, reads) {
  # Calculate RPM for each target feature
  feature_RPM <- (countOverlaps(targets, reads)/RPM_scaling_factor)
  save(feature_RPM,  
       file = paste0(outDir, targetsName, "_TEs_", libName, "_RPM.RData"))

  # Generate GRangesList object containing perms sets of random loci of the same
  # number and size distribution as for targets
  set.seed(9346)
  ranLocGRL <- GRangesList( mclapply(1:perms, function(i) {
    randomizeRegions(targets,
                     genome = genome,
                     per.chromosome = TRUE,
                     allow.overlaps = TRUE)
  }, mc.cores = 48) )
  ranLoc_RPM <- mclapply(seq_along(ranLocGRL), function(i) {
     (countOverlaps(ranLocGRL[[i]], reads)/RPM_scaling_factor)
  }, mc.cores = 48)
  save(ranLoc_RPM,
       file = paste0(outDir, targetsName, "_TEs_ranLoc_", libName, "_RPM.RData"))
  
  # Test whether more reads map to target features than to random loci
  ranLoc_RPM_sum <- mclapply(seq_along(ranLoc_RPM), function(i) {
    sum(ranLoc_RPM[[i]])
  })
  ranLocReads_lessThan_featureReads_Bool <- mclapply(seq_along(ranLoc_RPM_sum),
    function(i) {
      ranLoc_RPM_sum[[i]] < sum(feature_RPM)
    }, mc.cores = 48)
   
  # Calculate P-value
  Pval <- 1-(sum(unlist(ranLocReads_lessThan_featureReads_Bool))/perms)

  if(Pval == 0) {
    Pval <- minPval
  }

  # Create permutation test results object
  permTestResults <- new("permTest",
                         Pval = Pval,
                         feature_RPM_sum = sum(feature_RPM),
                         ranLoc_RPM_sum = unlist(ranLoc_RPM_sum),
                         Expected_RPM = mean(unlist(ranLoc_RPM_sum)),
                         log2_Observed_Expected = log2((sum(feature_RPM))/(mean(unlist(ranLoc_RPM_sum)))))
  save(permTestResults,
       file = paste0(outDir,
                     targetsName, "_TEs_permTestResults_",
                     as.character(perms), "perms_",
                     libName, "_RPM.RData"))

  options(scipen = 100)
  # Generate histogram
  pdf(paste0(plotDir, "hist_",
             targetsName, "_TEs_permTestResults_",
             as.character(perms), "perms_",
             libName, "_RPM.pdf"),
      height = 4, width = 5)
  # Disable scientific notation (e.g., 0.0001 rather than 1e-04)
  # Calculate max density
  maxDensityPlus <- max(density(permTestResults@ranLoc_RPM_sum)$y)*1.2
  alpha0.05 <- quantile(permTestResults@ranLoc_RPM_sum, 0.95)[[1]]
  hist(permTestResults@ranLoc_RPM_sum,
       freq = FALSE,
       col = "grey70",
       border = NA,
       lwd = 2,
       xlim = c(pmax(0, min(permTestResults@ranLoc_RPM_sum)-.1),
                pmax(permTestResults@feature_RPM_sum+.1, alpha0.05+.1)),
       ylim = c(0,
                maxDensityPlus),
       xlab = paste0(libName, " RPM"),
       ylab = "Density",
       main = bquote(atop(italic("P")~" = "~.(as.character(round(permTestResults@Pval,
                                                                 digits = 6))),
                     "Permutations = "~.(as.character(perms)))),
       cex.main = 1, cex.lab = 1, cex.axis = 1)
  lines(density(permTestResults@ranLoc_RPM_sum), lwd = 1.5)
  ablineclip(v = permTestResults@Expected_RPM,
             y1 = 0, y2 = maxDensityPlus*.92, lwd = 2)
  ablineclip(v = permTestResults@feature_RPM_sum,
             y1 = 0, y2 = maxDensityPlus*.92, lwd = 2, col = "forestgreen")
  ablineclip(v = alpha0.05,
             y1 = 0, y2 = maxDensityPlus*.92, lwd = 2, lty = 5, col = "red")
  text(x = c(pmax(0.05, min(permTestResults@ranLoc_RPM_sum)-.05),
             permTestResults@Expected_RPM,
             permTestResults@feature_RPM_sum,
             alpha0.05),
       y = c(maxDensityPlus*.95,
             maxDensityPlus,
             maxDensityPlus,
             maxDensityPlus*.95),
       labels = c("Permuted loci",
                  "Expected",
                  "Observed",
                  expression(alpha~" = 0.05")),
       col = c("grey70",
               "black",
               "forestgreen",
               "red"),
       cex = 0.7)
  box(lwd = 2)
  dev.off()

}

lapply(seq_along(TEsDNAGR), function(x) {
  locRPMcalc(targets = TEsDNAGR[[x]],
             targetsName = DNAfamNames[[x]],
             reads = lib_rangedGR)
})

lapply(seq_along(TEsRNAGR), function(x) {
  locRPMcalc(targets = TEsRNAGR[[x]],
             targetsName = RNAfamNames[[x]],
             reads = lib_rangedGR)
})


