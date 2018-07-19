#!/applications/R/R-3.4.0/bin/Rscript

###########################################################################
# Plot heatmaps of log2-transformed coverage levels
# around targets                         
###########################################################################

# Usage via Condor submission system on node7:
# csmit -m 100G -c 1 "geneLoops_EnrichedHeatmap_commandArgs.R /home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab.bed REC8_HA_Rep1"

## Source functions to be used in this script
#source("/projects/ajt200/Rfunctions/covMatrix_DNAmethMatrix_target_ranLoc.R")

#covDatPath <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab_head100000.bed"
#libName <- "REC8_HA_Rep1"

args <- commandArgs(trailingOnly = T)
covDatPath <- as.character(args[1])
libName <- as.character(args[2])

library(parallel)
library(doParallel)
library(EnrichedHeatmap)
library(genomation)
library(circlize)
library(RColorBrewer)

inDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/geneLoopHeatmaps/"
matDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/geneLoopHeatmaps/analysis_01/matrices/"
plotDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/geneLoopHeatmaps/analysis_01/plots/"
row_orderDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/geneLoopHeatmaps/analysis_01/row_order/"

# Import gene and transposon annotation in GTF format
# to obtain strand for correct orientation in heatmap
genes <- read.table(paste0(inDir, "Araport11_GFF3_genes_transposons.201606.gtf"),
                    header = F)
genes <- as.data.frame(cbind(as.vector(genes[,1]),
                             as.vector(genes[,7]),
                             as.vector(genes[,13])))
genes <- unique(genes)
colnames(genes) <- c("Chromosome", "Strand", "Gene_ID")

# Import gene loops as two GRanges objects:
# One with focal point start and end coordinates
# and another with interacting partner start and end coordinates
geneLoops <- read.table(paste0(inDir, "Liu_2016_Supplemental_Table_S5.txt"),
                        header = T)
geneLoopsStrand <- merge(x = geneLoops,
                         y = genes,
                         by.x = "Gene_ID",
                         by.y = "Gene_ID")
geneLoopsPlus <- geneLoopsStrand[geneLoopsStrand$Strand == "+",]
geneLoopsMinus <- geneLoopsStrand[geneLoopsStrand$Strand == "-",]

focusPlusGR <- GRanges(seqnames = geneLoopsPlus$Chromosome.y,
                       ranges = IRanges(start = geneLoopsPlus$Region_a_from,
                                        end = geneLoopsPlus$Region_a_to),
                       strand = geneLoopsPlus$Strand,
                       gene_ID  = geneLoopsPlus$Gene_ID,
                       q_value = geneLoopsPlus$q_value,
                       center = geneLoopsPlus$Center_a,
                       distance = geneLoopsPlus$Distance)
focusPlusGR <- sort(focusPlusGR,
                    by = ~ q_value, 
                    decreasing = F)
partnerPlusGR <- GRanges(seqnames = geneLoopsPlus$Chromosome.y,
                         ranges = IRanges(start = geneLoopsPlus$Region_b_from,
                                          end = geneLoopsPlus$Region_b_to),
                         strand = geneLoopsPlus$Strand,
                         gene_ID  = geneLoopsPlus$Gene_ID,
                         q_value = geneLoopsPlus$q_value,
                         center = geneLoopsPlus$Center_b,
                         distance = geneLoopsPlus$Distance)
partnerPlusGR <- sort(partnerPlusGR,
                      by = ~ q_value,
                      decreasing = F)

focusMinusGR <- GRanges(seqnames = geneLoopsMinus$Chromosome.y,
                        ranges = IRanges(start = geneLoopsMinus$Region_b_from,
                                         end = geneLoopsMinus$Region_b_to),
                        strand = geneLoopsMinus$Strand,
                        gene_ID  = geneLoopsMinus$Gene_ID,
                        q_value = geneLoopsMinus$q_value,
                        center = geneLoopsMinus$Center_b,
                        distance = geneLoopsMinus$Distance)
focusMinusGR <- sort(focusMinusGR,
                     by = ~ q_value,
                     decreasing = F)
partnerMinusGR <- GRanges(seqnames = geneLoopsMinus$Chromosome.y,
                          ranges = IRanges(start = geneLoopsMinus$Region_a_from,
                                           end = geneLoopsMinus$Region_a_to),
                          strand = geneLoopsMinus$Strand,
                          gene_ID  = geneLoopsMinus$Gene_ID,
                          q_value = geneLoopsMinus$q_value,
                          center = geneLoopsMinus$Center_a,
                          distance = geneLoopsMinus$Distance)
partnerMinusGR <- sort(partnerMinusGR,
                       by = ~ q_value,
                       decreasing = F)

focusGR <- sort(c(focusPlusGR, focusMinusGR),
                by = ~ q_value,
                decreasing = F)
partnerGR <- sort(c(partnerPlusGR, partnerMinusGR),
                  by = ~ q_value,
                  decreasing = F)
#focusGR <- focusGR[start(focusGR) <= 100000 & end(focusGR) <= 100000]
#partnerGR <- partnerGR[start(partnerGR) <= 100000 & end(partnerGR) <= 100000]
#partnerGR <- partnerGR[partnerGR$gene_ID != "AT3G01310"]

# Specify locations of normalised per base coverage files
libPath <- system(paste0("ls ", covDatPath), intern = T)

# Import coverage files as GRanges objects and assign to library names
gr <- readGeneric(libPath, meta.col = list(coverage = 4))

## Define window size for calculation of average coverage
## values within windows
#w <- 20
## Define lower and upper covvalue quantiles to be trimmed
#keep <- c(0, 0.01)

# function to create coverage matrices and heatmaps for genes
## and flanking regions for each dataset, and to extract row order
geneLoopHeatmap <- function(signal,
                            focus,
                            focusSize,
                            partner,
                            partnerSize,
                            flankSize,
                            flankLabL,
                            flankLabR,
                            winSize) {
  set.seed(4823)
  focus_smoothed <- normalizeToMatrix(signal = signal,
                                      target = focus,
                                      value_column = "coverage",
                                      extend = flankSize,
                                      mean_mode = "w0",
                                      w = winSize,
                                      background = 0,
                                      smooth = TRUE,
                                      include_target = TRUE,
                                      target_ratio = focusSize/(focusSize+(flankSize*2)))
  print("focus_smoothed")
  print(focus_smoothed)
  print(length(focus_smoothed))
  save(focus_smoothed,
       file = paste0(matDir, libName,
                     "_norm_cov_focus_and_flank_smoothed.RData"))

  partner_smoothed <- normalizeToMatrix(signal = signal,
                                        target = partner,
                                        value_column = "coverage",
                                        extend = flankSize,
                                        mean_mode = "w0",
                                        w = winSize,
                                        background = 0,
                                        smooth = TRUE,
                                        include_target = TRUE,
                                        target_ratio = partnerSize/(partnerSize+(flankSize*2)))
  print("partner_smoothed")
  print(partner_smoothed)
  print(length(partner_smoothed))
  save(partner_smoothed,
       file = paste0(matDir, libName,
                     "_norm_cov_partner_and_flank_smoothed.RData"))

  col_fun1 <- colorRamp2(quantile(c(focus_smoothed, partner_smoothed), c(0.5, 0.75, 0.95)), c("blue", "white", "red"))
  ht_list  <- print(EnrichedHeatmap(focus_smoothed,
                                    top_annotation = NULL,
                                    col = col_fun1,
                                    name = paste0(libName, " F"),
                                    row_order = c(1:length(focus)),
                                    column_title = "Foci ordered by gene loop FDR",
                                    axis_name = c(flankLabL, "Start", "End", flankLabR),
                                    border = FALSE)) +
              print(EnrichedHeatmap(partner_smoothed,
                                    top_annotation = NULL,
                                    col = col_fun1,
                                    name = paste0(libName, " P"),
                                    row_order = c(1:length(partner)),
                                    column_title = "Partners ordered by gene loop FDR",
                                    axis_name = c(flankLabL, "Start", "End", flankLabR),
                                    border = FALSE))
  pdf(paste0(plotDir, libName, "_heatmap_norm_cov_focus_partner_and_flank_smoothed.pdf"))
  draw(ht_list, gap = unit(10, "mm"))
  dev.off()
}

geneLoopHeatmap(signal = gr,
                focus = focusGR,
                focusSize = mean(width(focusGR)),
                partner = partnerGR,
                partnerSize = mean(width(partnerGR)),
                flankSize = 1000,
                flankLabL = "-1 kb",
                flankLabR = "+1 kb",
                winSize = 10)
 

