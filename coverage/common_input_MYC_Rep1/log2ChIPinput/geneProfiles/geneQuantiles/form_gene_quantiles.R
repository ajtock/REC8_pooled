#!/applications/R/R-3.5.0/bin/Rscript

##################################################################
# Select genes in TAIR10 annotation for quantile formation based #
# on mean normalised coverage levels within given gene "region"  #
##################################################################

# Usage:
# ./form_gene_quantiles.R '/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/' 'log2_REC8_HA_Rep2_ChIP_REC8_MYC_Rep1_input' 'bodies' 6 

#covDir <- "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/"
#covName <- "log2_REC8_HA_Rep2_ChIP_REC8_MYC_Rep1_input"
#region <- "bodies"
#groups <- 6

args <- commandArgs(trailingOnly = T)
covDir <- args[1]
covName <- args[2]
region <- args[3]
groups <- as.integer(args[4])

quantDir <- paste0(covName, "_in_gene_", region, "/")
system(paste0("[ -d ", quantDir, " ] || mkdir ", quantDir))

library(segmentSeq)

genes <- read.table("/projects/ajt200/TAIR10/representative_genes/representative_genes_uniq_fmt_strand.txt",
                    header = T)
genes <- data.frame(chr = as.character(paste0("Chr", genes$chr)),
                    genes[,-1])
print(dim(genes))
#[1] 27204     5
# Convert into GRanges
genesGR <- GRanges(seqnames = genes$chr,
                   ranges = IRanges(start = genes$start,
                                    end = genes$end),
                   strand = genes$strand,
                   gene_model = genes$gene_model)
# Get promoters (TSS-1 to TSS-500) and terminators (TTS+1 to TTS+500)
promotersGR <- promoters(genesGR, upstream = 500, downstream = 0)
source("/projects/ajt200/Rfunctions/TTSplus.R")
terminatorsGR <- TTSplus(genesGR, upstream = -1, downstream = 500)
# Convert into data.frames
promoters <- data.frame(promotersGR)[,-4]
colnames(promoters) <- c("chr", colnames(promoters[,2:5]))
terminators <- data.frame(terminatorsGR)[,-4]
colnames(terminators) <- c("chr", colnames(terminators[,2:5]))

if(region == "bodies") {
  targets <- genes 
  targetsGR <- genesGR
} else if(region == "promoters") {
  targets <- promoters
  targetsGR <- promotersGR
} else if(region == "terminators") {
  targets = terminators
  targetsGR = terminatorsGR
} else {
  stop("region is not bodies, promoters or terminators")
} 

# Calculate mean coverage levels within given region
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")

normCov <- read.table(paste0(covDir, covName,
                             "_norm_allchrs_coverage_coord_tab.bed"))
allCovTargets <- NULL
for(i in 1:5) {
  print(i)
  chrCov <- normCov[normCov[,1] == chrs[i],]
  chrCov <- chrCov[,4]
  covCoords <- seq(from = 1, to = length(chrCov), by = 1)
  covCoordsGR <- GRanges(seqnames = chrs[i],
                         ranges = IRanges(start = covCoords,
                                          width = 1),
                         strand = "*")
  chrTargets <- targets[targets$chr == chrs[i],]
  chrTargetsGR <- targetsGR[seqnames(targetsGR) == chrs[i]]
  strand(chrTargetsGR) <- "*"
  targetsOverlaps <- getOverlaps(coordinates = chrTargetsGR,
                                 segments = covCoordsGR,
                                 whichOverlaps = TRUE,
                                 ignoreStrand = TRUE)
  normCov_targetCov <- sapply(targetsOverlaps,
                              function(x) mean(chrCov[x]))
  chrTargets <- cbind(chrTargets, normCov_targetCov)
  allCovTargets <- rbind(allCovTargets, chrTargets)
}
write.table(allCovTargets,
            file = paste0(quantDir,
                          covName, "_in_gene_", region, "_meanCov.txt"),
            row.names = F, sep = "\t", quote = F)

# Group genes into quantiles based on normalised coverage
print(head(allCovTargets))
allCovTargets <- allCovTargets[order(allCovTargets$normCov_targetCov,
                                     decreasing = T),]
print(head(allCovTargets))
print(tail(allCovTargets))
quant <- round(length(allCovTargets[,1])/groups)
quant <- cumsum(c(1, rep(quant, times = groups)))
quant <- c(quant, length(allCovTargets[,1]))
quantDat <- NULL
for(j in 1:length(quant)) {
  print(j)
  if(j == 1) {
    quantjTargets <- allCovTargets[quant[j]:quant[j+1]-1,]
    print(paste0("condition 1: quantile ", j))
    print(dim(quantjTargets))
  }
  if(j > 1 & j < groups) {
    quantjTargets <- allCovTargets[(quant[j]):(quant[j+1]-1),]
    print(paste0("condition 2: quantile ", j))
    print(dim(quantjTargets))
  }
  if(j == groups) {
    quantjTargets <- allCovTargets[(quant[j]):(quant[length(quant)]),]
    print(paste0("condition 3: quantile ", j))
    print(dim(quantjTargets))
  } 
  if(j <= groups) {    
    dat <- c(j, length(quantjTargets[,1]),
             round(mean(quantjTargets[,3]-quantjTargets[,2])),
             sum(quantjTargets[,3]-quantjTargets[,2]),
             mean(quantjTargets[,6]))
    quantDat <- rbind(quantDat, dat)
    write.table(quantjTargets,
                file = paste0(quantDir,
                              "quantile_", j, "_",
                              covName, "_in_gene_", region, "_meanCov.txt"),
                row.names = F, sep = "\t", quote = F)
  }
}
colnames(quantDat) <- c("quantile",
                        "n",
                        "mean_bp",
                        "total_bp",
                        "normCov_meanCov")
write.table(quantDat,
            file = paste0(quantDir,
                          covName, "_in_gene_", region, "_meanCov_quantile_stats.txt"),
            row.names = F, sep = "\t", quote = F)
