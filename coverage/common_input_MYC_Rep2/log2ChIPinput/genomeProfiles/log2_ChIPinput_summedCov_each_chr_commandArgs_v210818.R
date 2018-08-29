#!/applications/R/R-3.4.0/bin/Rscript

# Sum all Z-score-standardised log2-transformed library-size-normalised
# per base coverage values in each of the 5 chromosomes

# Usage on cluster node7:
# csmit -m 20G -c 1 "./log2_ChIPinput_summedCov_each_chr_commandArgs_v210818.R /home/ajt200/bedlike_per_base_cov_files/log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab.bed REC8_HA_Rep1"

args <- commandArgs(trailingOnly = TRUE)
covFile <- args[1]
libName <- args[2]

outDir <- "./" 

covDat <- read.table(covFile,
                     colClasses = c(NA, NA, "NULL", NA))

chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")

chrCov_all <- NULL
for(i in 1:5) {
  chrCov <- sum(covDat[covDat$V1 == chrs[i],]$V4)
  chrCov_all <- c(chrCov_all, chrCov)
}
write.table(chrCov_all,
            file = paste0(outDir, libName, "_summed_Zscore_log2_ChIPinput_chrCov_all.txt"))

