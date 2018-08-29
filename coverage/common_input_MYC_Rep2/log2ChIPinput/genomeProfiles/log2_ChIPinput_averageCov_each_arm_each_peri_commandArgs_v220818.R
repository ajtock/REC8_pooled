#!/applications/R/R-3.4.0/bin/Rscript

#######################################################################################
# Calculate average of all per base coverage values in each of                        #
# 10 chromosome arm regions and in each of 5 pericentromeric regions in wt and mutant #
####################################################################################### 

# Usage on hydrogen node7
# csmit -m 20G -c 1 "./log2_ChIPinput_averageCov_each_arm_each_peri_commandArgs_v220818.R REC8_HA_Rep1 kss_REC8_HA_Rep1"

args <- commandArgs(trailingOnly = TRUE)
wt_libName <- args[1]
mutant_libName <- args[2]

wt_covFile <- paste0("../log2_", wt_libName,
                     "_ChIP_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab.bed")
mutant_covFile <- paste0("../log2_", mutant_libName,
                         "_ChIP_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab.bed")
outDir <- "./" 

wt_covDat <- read.table(wt_covFile,
                        colClasses = c(NA, NA, "NULL", NA))
mutant_covDat <- read.table(mutant_covFile,
                            colClasses = c(NA, NA, "NULL", NA))

# Genomic definitions
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
# Pericentromeric regions are as defined in Supplemental Table S26
# of Ziolkowski et al. (2017) Genes Dev. 31
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)

wt_chrCov_armL_all <- NULL
wt_chrCov_armR_all <- NULL
wt_chrCov_peri_all <- NULL
mutant_chrCov_armL_all <- NULL
mutant_chrCov_armR_all <- NULL
mutant_chrCov_peri_all <- NULL
for(i in 1:5) {
  wt_chrCov <- wt_covDat[wt_covDat$V1 == chrs[i],]
  wt_chrCov_armL <- sum(wt_chrCov[wt_chrCov$V2 < pericenStart[i],]$V4)/(pericenStart[i]-1)
  wt_chrCov_armL_all <- c(wt_chrCov_armL_all, wt_chrCov_armL)
  wt_chrCov_armR <- sum(wt_chrCov[wt_chrCov$V2 > pericenEnd[i],]$V4)/(chrLens[i]-pericenEnd[i])
  wt_chrCov_armR_all <- c(wt_chrCov_armR_all, wt_chrCov_armR)
  wt_chrCov_peri <- sum(wt_chrCov[wt_chrCov$V2 >= pericenStart[i] & wt_chrCov$V2 <= pericenEnd[i],]$V4)/(pericenEnd[i]-pericenStart[i]+1)
  wt_chrCov_peri_all <- c(wt_chrCov_peri_all, wt_chrCov_peri)

  mutant_chrCov <- mutant_covDat[mutant_covDat$V1 == chrs[i],]
  mutant_chrCov_armL <- sum(mutant_chrCov[mutant_chrCov$V2 < pericenStart[i],]$V4)/(pericenStart[i]-1)
  mutant_chrCov_armL_all <- c(mutant_chrCov_armL_all, mutant_chrCov_armL)
  mutant_chrCov_armR <- sum(mutant_chrCov[mutant_chrCov$V2 > pericenEnd[i],]$V4)/(chrLens[i]-pericenEnd[i])
  mutant_chrCov_armR_all <- c(mutant_chrCov_armR_all, mutant_chrCov_armR)
  mutant_chrCov_peri <- sum(mutant_chrCov[mutant_chrCov$V2 >= pericenStart[i] & mutant_chrCov$V2 <= pericenEnd[i],]$V4)/(pericenEnd[i]-pericenStart[i]+1)
  mutant_chrCov_peri_all <- c(mutant_chrCov_peri_all, mutant_chrCov_peri)
}
write.table(wt_chrCov_armL_all,
            file = paste0(outDir, wt_libName, "_chrCov_armL_all.txt"))
write.table(wt_chrCov_armR_all,
            file = paste0(outDir, wt_libName, "_chrCov_armR_all.txt"))
write.table(wt_chrCov_peri_all,
            file = paste0(outDir, wt_libName, "_chrCov_peri_all.txt"))
write.table(mutant_chrCov_armL_all,
            file = paste0(outDir, mutant_libName, "_chrCov_armL_all.txt"))
write.table(mutant_chrCov_armR_all,
            file = paste0(outDir, mutant_libName, "_chrCov_armR_all.txt"))
write.table(mutant_chrCov_peri_all,
            file = paste0(outDir, mutant_libName, "_chrCov_peri_all.txt"))

