#!/applications/R/R-3.4.0/bin/Rscript

# Perform U tests using average of all per base coverage values in each of
# 10 chromosome arm regions and in each of 5 pericentromeric regions in wt and mutant

#############################################################################
# Sum all per base coverage values in each of 10 chromosome arm regions and #
# in each of 5 pericentromeric regions in wt and mutant                     #
############################################################################# 

# Usage on hydrogen node7
# csmit -m 20G -c 1 "./log2_ChIPinput_averageCov_each_arm_each_peri_commandArgs_v220818_u_test.R REC8_HA_Rep1 kss_REC8_HA_Rep1"

args <- commandArgs(trailingOnly = TRUE)
wt_libName <- args[1]
mutant_libName <- args[2]

outDir <- "./" 

wt_chrCov_armL_all <- read.table(paste0(outDir, wt_libName, "_chrCov_armL_all.txt"))
wt_chrCov_armR_all <- read.table(paste0(outDir, wt_libName, "_chrCov_armR_all.txt"))
wt_chrCov_peri_all <- read.table(paste0(outDir, wt_libName, "_chrCov_peri_all.txt"))
mutant_chrCov_armL_all <- read.table(paste0(outDir, mutant_libName, "_chrCov_armL_all.txt"))
mutant_chrCov_armR_all <- read.table(paste0(outDir, mutant_libName, "_chrCov_armR_all.txt"))
mutant_chrCov_peri_all <- read.table(paste0(outDir, mutant_libName, "_chrCov_peri_all.txt"))

arm_u_test <- wilcox.test(x = c(wt_chrCov_armL_all$x, wt_chrCov_armR_all$x),
                          y = c(mutant_chrCov_armL_all$x, mutant_chrCov_armR_all$x),
                          alternative = "less")
armL_u_test <- wilcox.test(x = wt_chrCov_armL_all$x,
                           y = mutant_chrCov_armL_all$x,
                           alternative = "less")
armR_u_test <- wilcox.test(x = wt_chrCov_armR_all$x,
                           y = mutant_chrCov_armR_all$x,
                           alternative = "less")
peri_u_test <- wilcox.test(x = wt_chrCov_peri_all$x,
                           y = mutant_chrCov_peri_all$x,
                           alternative = "greater")
save(arm_u_test,
     file = paste0(outDir,
                   wt_libName, "_",
                   mutant_libName,
                   "_chrCov_arm_u_test.RData"))
save(armL_u_test,
     file = paste0(outDir,
                   wt_libName, "_",
                   mutant_libName,
                   "_chrCov_armL_u_test.RData"))
save(armR_u_test,
     file = paste0(outDir,
                   wt_libName, "_",
                   mutant_libName,
                   "_chrCov_armR_u_test.RData"))
save(peri_u_test,
     file = paste0(outDir,
                   wt_libName, "_",
                   mutant_libName,
                   "_chrCov_peri_u_test.RData"))

load(
     file = paste0(outDir,
                   wt_libName, "_",
                   mutant_libName,
                   "_chrCov_arm_u_test.RData"))
load(
     file = paste0(outDir,
                   wt_libName, "_",
                   mutant_libName,
                   "_chrCov_armL_u_test.RData"))
load(
     file = paste0(outDir,
                   wt_libName, "_",
                   mutant_libName,
                   "_chrCov_armR_u_test.RData"))
load(
     file = paste0(outDir,
                   wt_libName, "_",
                   mutant_libName,
                   "_chrCov_peri_u_test.RData"))
