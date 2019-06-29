#!/applications/R/R-3.4.0/bin/Rscript

# Plot heatmaps of features sorted by coverage levels between start and end sites

# Usage:
# /applications/R/R-3.4.0/bin/Rscript features_heatmap_sorted.R REC8_MYC_Rep1_ChIP genes 2000 2kb '2 kb' 20 20bp bodies

#libName <- "REC8_MYC_Rep1_ChIP"
#featureName <- "genes"
#upstream <- 2000
#downstream <- 2000
#flankName <- "2kb"
#flankNamePlot <- "2 kb"
#binSize <- 20
#binName <- "20bp"
#region <- "bodies"

args <- commandArgs(trailingOnly = T)
libName <- args[1]
featureName <- args[2]
upstream <- as.numeric(args[3])
downstream <- as.numeric(args[3])
flankName <- args[4]
flankNamePlot <- args[5]
binSize <- as.numeric(args[6])
binName <- args[7]
region <- args[8]

library(EnrichedHeatmap)
library(parallel)
library(circlize)
library(RColorBrewer)
library(tidyr)
library(ggplot2)
library(doParallel)

matDir <- "./matrices/"
plotDir <- "./plots/"
regionPlotDir <- paste0(plotDir, region,
                        "_by_", libName, "/")
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))
system(paste0("[ -d ", regionPlotDir, " ] || mkdir ", regionPlotDir))
plotDir <- regionPlotDir

ChIPseqNames <- c( 
                  "REC8_MYC_Rep1_ChIP",
                  "REC8_HA_Rep1_ChIP",
                  "control_MYC_Rep1_ChIP",
                  "control_HA_Rep1_ChIP",
                  "REC8_MYC_Rep1_input"
             )
ChIPseqNamesPlot <- c( 
                      "REC8-Myc Rep1 ChIP",
                      "REC8-HA Rep1 ChIP",
                      "Control-Myc Rep1 ChIP",
                      "Control-HA Rep1 ChIP",
                      "REC8-Myc Rep1 input"
                     )

ChIPseqColours <- c(
                    "red",
                    "red",
                    "black",
                    "black",
                    "grey60"
                   )

# Load matrix and extract region for sorting of features
mat1 <- as.matrix(read.table(paste0(matDir,
                                    libName,
                                    "_norm_cov_feature_smoothed_target_and_",
                                    flankName, "_flank_dataframe.txt"),
                             header = T))
bodyLength <- (dim(mat1)[2]-((upstream+downstream)/binSize))*binSize
if( region == "promoters" ) {
  mat1Region <- mat1[,(((upstream-500)/binSize)+1):(upstream/binSize)]
} else if ( region == "terminators" ) {
  mat1Region <- mat1[,(((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(500/binSize))]
} else if ( region == "bodies" ) {
  mat1Region <- mat1[,((upstream/binSize)+1):((upstream+bodyLength)/binSize)]
} else {
  print("The region name provided does not match 'promoter', 'terminator', or 'bodies'")
}
mat1RegionRowMeans <- rowMeans(mat1Region, na.rm = T)
mat1RegionRowMeansSorted <- sort.int(mat1RegionRowMeans,
                                     decreasing = T,
                                     index.return = T,
                                     na.last = T)
mat1RegionRowSums <- rowSums(mat1Region, na.rm = T)
mat1RegionRowSumsSorted <- sort.int(mat1RegionRowSums,
                                    decreasing = T,
                                    index.return = T,
                                    na.last = T)
mat1RegionRowMedians <- apply(X = mat1Region,
                              MARGIN = 1,
                              FUN = median)
mat1RegionRowMediansSorted <- sort.int(mat1RegionRowMedians,
                                       decreasing = T,
                                       index.return = T,
                                       na.last = T)

# Load feature coverage matrix for each ChIP-seq dataset
mats <- mclapply(seq_along(ChIPseqNames), function(x) {
  as.matrix(read.table(paste0(matDir,
                              ChIPseqNames[x],                            
                              "_norm_cov_feature_smoothed_target_and_",
                              flankName, "_flank_dataframe.txt"),
                       header = T))
}, mc.cores = length(ChIPseqNames))

# Sort by decreasing mat1RegionRowMeans
matsSorted <- mclapply(seq_along(mats), function(x) {
  mats[[x]][sort.int(mat1RegionRowMeans,
                     decreasing = T,
                     index.return = T,
                     na.last = T)$ix,]
}, mc.cores = length(ChIPseqNames))

# Convert matrices to class "normalizedMatrix" with associated attributes
# for use by EnrichedHeatmap function
for(x in seq_along(matsSorted)) {
  attr(matsSorted[[x]], "upstream_index")         = 1:(upstream/binSize)
  attr(matsSorted[[x]], "target_index")           = ((upstream/binSize)+1):((upstream+bodyLength)/binSize)
  attr(matsSorted[[x]], "downstream_index")       = (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))
  attr(matsSorted[[x]], "extend")                 = c(upstream, downstream)
  attr(matsSorted[[x]], "smooth")                 = TRUE
  attr(matsSorted[[x]], "signal_name")            = ChIPseqNamesPlot[x]
  attr(matsSorted[[x]], "target_name")            = featureName
  attr(matsSorted[[x]], "target_is_single_point") = FALSE
  attr(matsSorted[[x]], "background")             = 0
  attr(matsSorted[[x]], "signal_is_categorical")  = FALSE
  class(matsSorted[[x]]) = c("normalizedMatrix", "matrix")
}

if(featureName == "genes") {
  featureStartLab <- "TSS"
  featureEndLab <- "TTS"
} else {
  featureStartLab <- "Start"
  featureEndLab <- "End"
}

# Heatmap plotting function
featureHeatmap <- function(matSorted,
                           col_fun,
                           colour,
                           datName,
                           rowOrder) {
  print(EnrichedHeatmap(mat = matSorted,
                        col = col_fun,
                        top_annotation = HeatmapAnnotation(enriched = anno_enriched(gp = gpar(col = colour,
                                                                                              lwd = 3),
                                                                                    yaxis_side = "right",
                                                                                    yaxis_facing = "right",
                                                                                    yaxis_gp = gpar(fontsize = 10),
                                                                                    pos_line_gp = gpar(col = "black",
                                                                                                       lty = 2,
                                                                                                       lwd = 2))),
                        top_annotation_height = unit(2, "cm"),
                        name = "Coverage",
                        row_order = rowOrder,
                        column_title = datName,
                        axis_name = c(paste0("-", flankNamePlot),
                                      featureStartLab, featureEndLab,
                                      paste0("+", flankNamePlot)),
                        axis_name_gp = gpar(fontsize = 12),
                        border = FALSE,
                        pos_line_gp = gpar(col = "black", lty = 2, lwd = 2),
                        use_raster = FALSE))
}

# Define heatmap colours
rich8to6equal <- c("#000041", "#0000CB", "#0081FF", "#FDEE02", "#FFAB00", "#FF3300")

library(doParallel)
registerDoParallel(cores = length(ChIPseqNames))
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

# Plot
# ChIP-seq
foreach(x = 1:length(ChIPseqNames)) %dopar% {
  ChIPseq_col_fun <- colorRamp2(quantile(matsSorted[[x]],
                                         c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
                                         na.rm = T),
                                rich8to6equal)
  pdf(paste0(plotDir, ChIPseqNames[x], "_around_", featureName,
             "_heatmap_ordered_by_", libName, "_in_", region, ".pdf"),
      width = 5, height = 10)
  featureHeatmap(matSorted = matsSorted[[x]],
                 col_fun = ChIPseq_col_fun,
                 colour = ChIPseqColours[x],
                 datName = ChIPseqNamesPlot[x],
                 rowOrder = c(1:dim(matsSorted[[x]])[1]))
  dev.off()
}

