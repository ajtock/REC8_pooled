#!/applications/R/R-3.5.0/bin/Rscript

# Plot feature average coverage profiles with 95% CIs

# Usage:
# /applications/R/R-3.5.0/bin/Rscript features_avgProfile_noCorr_2yaxes.R Genes 'Genes' 2000 2kb '2 kb' 20 20bp 'ASY1_Rep2,REC8_HA_Rep2,MNase,SPO11_1_oligos_RPI1' 'ASY1,REC8-HA,Nucleosomes,SPO11-1-oligos' 'red,green2,purple4,dodgerblue2'

#featureName <- "Genes"
#featureNamePlot <- "Genes"
#upstream <- 2000
#downstream <- 2000
#flankName <- "2kb"
#flankNamePlot <- "2 kb"
#binSize <- 20
#binName <- "20bp"
#libNames <- unlist(strsplit("ASY1_Rep2,REC8_HA_Rep2,MNase,SPO11_1_oligos_RPI1",
#                            split = ","))
#libNamesPlot <- unlist(strsplit("ASY1,REC8-HA,Nucleosomes,SPO11-1-oligos",
#                                split = ","))
#colours <- unlist(strsplit("red,green2,purple4,dodgerblue2",
#                           split = ","))

args <- commandArgs(trailingOnly = T)
featureName <- args[1]
featureNamePlot <- args[2]
upstream <- as.numeric(args[3])
downstream <- as.numeric(args[3])
flankName <- args[4]
flankNamePlot <- args[5]
binSize <- as.numeric(args[6])
binName <- args[7]
libNames <- unlist(strsplit(args[8],
                            split = ","))
libNamesPlot <- unlist(strsplit(args[9],
                                split = ","))
colours <- unlist(strsplit(args[10],
                           split = ","))

library(parallel)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(grid)
library(gridExtra)
library(extrafont)

matDir <- "matrices/"
plotDir <- "plots/"
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Load feature coverage matrix for each ChIP-seq dataset
featureMats <- mclapply(seq_along(libNames), function(x) {
  as.matrix(read.table(paste0(matDir,
                              libNames[x],                            
                              "_norm_cov_feature_smoothed_target_and_",
                              flankName, "_flank_dataframe.txt"),
                       header = T))
}, mc.cores = length(libNames))

# Load ranLoc coverage matrix for each ChIP-seq dataset
ranLocMats <- mclapply(seq_along(libNames), function(x) {
  as.matrix(read.table(paste0(matDir,
                              libNames[x],                            
                              "_norm_cov_ranLoc_smoothed_target_and_",
                              flankName, "_flank_dataframe.txt"),
                       header = T))
}, mc.cores = length(libNames))

## feature
# Transpose matrix and convert to dataframe
# in which first column is window name
wideDFfeature_list <- mclapply(seq_along(featureMats), function(x) {
  data.frame(window = colnames(featureMats[[x]]),
             t(featureMats[[x]]))
}, mc.cores = length(featureMats))

# Convert into tidy data.frame (long format)
tidyDFfeature_list  <- mclapply(seq_along(wideDFfeature_list), function(x) {
  gather(data  = wideDFfeature_list[[x]],
         key   = feature,
         value = coverage,
         -window)
}, mc.cores = length(wideDFfeature_list))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list)) {
  tidyDFfeature_list[[x]]$window <- factor(tidyDFfeature_list[[x]]$window,
                                           levels = as.character(wideDFfeature_list[[x]]$window))
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list  <- mclapply(seq_along(tidyDFfeature_list), function(x) {
  data.frame(window = as.character(wideDFfeature_list[[x]]$window),
             n      = tapply(X     = tidyDFfeature_list[[x]]$coverage,
                             INDEX = tidyDFfeature_list[[x]]$window,
                             FUN   = length),
             mean   = tapply(X     = tidyDFfeature_list[[x]]$coverage,
                             INDEX = tidyDFfeature_list[[x]]$window,
                             FUN   = mean),
             sd     = tapply(X     = tidyDFfeature_list[[x]]$coverage,
                             INDEX = tidyDFfeature_list[[x]]$window,
                             FUN   = sd))
}, mc.cores = length(tidyDFfeature_list))

for(x in seq_along(summaryDFfeature_list)) {
  summaryDFfeature_list[[x]]$window <- factor(summaryDFfeature_list[[x]]$window,
                                              levels = as.character(wideDFfeature_list[[x]]$window))
  summaryDFfeature_list[[x]]$winNo <- factor(1:dim(summaryDFfeature_list[[x]])[1])
  summaryDFfeature_list[[x]]$sem <- summaryDFfeature_list[[x]]$sd/sqrt(summaryDFfeature_list[[x]]$n-1)
  summaryDFfeature_list[[x]]$CI_lower <- summaryDFfeature_list[[x]]$mean -
    qt(0.975, df = summaryDFfeature_list[[x]]$n-1)*summaryDFfeature_list[[x]]$sem
  summaryDFfeature_list[[x]]$CI_upper <- summaryDFfeature_list[[x]]$mean +
    qt(0.975, df = summaryDFfeature_list[[x]]$n-1)*summaryDFfeature_list[[x]]$sem
}

names(summaryDFfeature_list) <- libNamesPlot

# Convert list summaryDFfeature_list into a single data.frame for plotting
summaryDFfeature <- bind_rows(summaryDFfeature_list, .id = "libName")
summaryDFfeature$libName <- factor(summaryDFfeature$libName,
                                   levels = names(summaryDFfeature_list))


## ranLoc
# Transpose matrix and convert to dataframe
# in which first column is window name
wideDFranLoc_list <- mclapply(seq_along(ranLocMats), function(x) {
  data.frame(window = colnames(ranLocMats[[x]]),
             t(ranLocMats[[x]]))
}, mc.cores = length(ranLocMats))

# Convert into tidy data.frame (long format)
tidyDFranLoc_list  <- mclapply(seq_along(wideDFranLoc_list), function(x) {
  gather(data  = wideDFranLoc_list[[x]],
         key   = ranLoc,
         value = coverage,
         -window)
}, mc.cores = length(wideDFranLoc_list))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFranLoc_list)) {
  tidyDFranLoc_list[[x]]$window <- factor(tidyDFranLoc_list[[x]]$window,
                                          levels = as.character(wideDFranLoc_list[[x]]$window))
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (ranLocs) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFranLoc_list  <- mclapply(seq_along(tidyDFranLoc_list), function(x) {
  data.frame(window = as.character(wideDFranLoc_list[[x]]$window),
             n      = tapply(X     = tidyDFranLoc_list[[x]]$coverage,
                             INDEX = tidyDFranLoc_list[[x]]$window,
                             FUN   = length),
             mean   = tapply(X     = tidyDFranLoc_list[[x]]$coverage,
                             INDEX = tidyDFranLoc_list[[x]]$window,
                             FUN   = mean),
             sd     = tapply(X     = tidyDFranLoc_list[[x]]$coverage,
                             INDEX = tidyDFranLoc_list[[x]]$window,
                             FUN   = sd))
}, mc.cores = length(tidyDFranLoc_list))

for(x in seq_along(summaryDFranLoc_list)) {
  summaryDFranLoc_list[[x]]$window <- factor(summaryDFranLoc_list[[x]]$window,
                                             levels = as.character(wideDFranLoc_list[[x]]$window))
  summaryDFranLoc_list[[x]]$winNo <- factor(1:dim(summaryDFranLoc_list[[x]])[1])
  summaryDFranLoc_list[[x]]$sem <- summaryDFranLoc_list[[x]]$sd/sqrt(summaryDFranLoc_list[[x]]$n-1)
  summaryDFranLoc_list[[x]]$CI_lower <- summaryDFranLoc_list[[x]]$mean -
    qt(0.975, df = summaryDFranLoc_list[[x]]$n-1)*summaryDFranLoc_list[[x]]$sem
  summaryDFranLoc_list[[x]]$CI_upper <- summaryDFranLoc_list[[x]]$mean +
    qt(0.975, df = summaryDFranLoc_list[[x]]$n-1)*summaryDFranLoc_list[[x]]$sem
}

names(summaryDFranLoc_list) <- libNamesPlot

# Convert list summaryDFranLoc_list into a single data.frame for plotting
summaryDFranLoc <- bind_rows(summaryDFranLoc_list, .id = "libName")
summaryDFranLoc$libName <- factor(summaryDFranLoc$libName,
                                  levels = names(summaryDFranLoc_list))

if(featureName == "Genes") {
  featureStartLab <- "TSS"
  featureEndLab <- "TTS"
} else {
  featureStartLab <- "Start"
  featureEndLab <- "End"
}

# Source plotting functions
source("/projects/ajt200/Rfunctions/plotAvgCov_plotAvgMeth_target_ranLoc.R")

# Plot one vs another average coverage profile around target and random loci
pdf(paste0(plotDir, paste(libNames, collapse = "_"),
           "_around_", featureName, "_basePlot.pdf"),
    height = (length(libNames)/2)*2.5, width = 6)
par(mfrow = c(length(libNames)/2, 2),
    mar = c(2.1, 3.2, 2.1, 3.2),
    mgp = c(3, 0.5, 0))
for(x in seq(from = 1, to = length(libNames), by = 2)) {
  plotAvgCov_oneVanother(xplot = seq(from = 1,
                                     to = dim(summaryDFfeature_list[[1]])[1], 
                                     by = 1),
                         dat1 = summaryDFfeature_list[[x]]$mean,
                         dat2 = summaryDFfeature_list[[x+1]]$mean,
                         ranDat1 = summaryDFranLoc_list[[x]]$mean,
                         ranDat2 = summaryDFranLoc_list[[x+1]]$mean,
                         flankSize = upstream, winSize = binSize,
                         Ylabel1 = libNamesPlot[x], Ylabel2 = libNamesPlot[x+1],
                         flankLabL = paste0("-", flankNamePlot), paste0("+", flankNamePlot),
                         startLab1 = featureStartLab, endLab1 = featureEndLab,
                         startLab2 = "Start", endLab2 = "End",
                         mycols = colours[x:(x+1)])
}
dev.off()
