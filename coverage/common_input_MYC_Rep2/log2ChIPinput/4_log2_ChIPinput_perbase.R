#!/applications/R/R-3.3.2/bin/Rscript

####################################################################
# Calculate log2 ratio of ChIP and input coverage values           #
# and write to 1-based BED-like file suitable for EnrichedHeatmap  #
####################################################################

# Usage:
# csmit -m 20G -c 1 "./4_log2_ChIPinput_perbase.R REC8_HA_Rep1_ChIP REC8_MYC_Rep2_input"

args <- commandArgs(trailingOnly = TRUE)
ChIPname <- args[1]
inputname <- args[2]
outname <- paste0("log2_", ChIPname, "_", inputname)

inDir <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/"
outDir <- paste0(inDir,
                 "common_input_MYC_Rep2/log2ChIPinput/")

ChIP <- read.table(paste0(inDir,
                          ChIPname,
                          "_norm_allchrs_coverage_coord_tab.bed"))
print(head(ChIP))
input <- read.table(paste0(inDir,
                           inputname,
                           "_norm_allchrs_coverage_coord_tab.bed"))
print(head(input))
ChIP_input_offset <- cbind(ChIP, ChIP$V4+1, input$V4+1)  
ChIP_input_offset <- ChIP_input_offset[,-4]
colnames(ChIP_input_offset) <- c("V1", "V2", "V3", "V4", "V5")
print(head(ChIP_input_offset))
norm <- log2(ChIP_input_offset$V4/ChIP_input_offset$V5)
norm <- (norm-mean(norm, na.rm = T))/sd(norm, na.rm = T)
ChIPnormInput <- cbind(ChIP_input_offset, norm)
log2ChIPinput <- ChIPnormInput[,-4:-5] 
write.table(log2ChIPinput,
            file = paste0(outDir,
                          outname,
                          "_norm_allchrs_coverage_coord_tab.bed"),
            sep = "\t", quote = F,
            row.names = F, col.names = F)
rm(ChIPnormInput)
gc()

