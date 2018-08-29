#!/applications/R/R-3.4.0/bin/Rscript

#############################################################################
# Sum all per base coverage values in each of 10 chromosome arm regions and #
# in each of 5 pericentromeric regions in wt and mutant                     #
############################################################################# 

# Usage on hydrogen node7
# csmit -m 20G -c 1 "./log2_ChIPinput_averageCov_each_arm_each_peri_commandArgs_v220818_t_test.R REC8_HA_Rep1 kss_REC8_HA_Rep1"

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

arm_f_test <- var.test(x = c(wt_chrCov_armL_all$x, wt_chrCov_armR_all$x),
                       y = c(mutant_chrCov_armL_all$x, mutant_chrCov_armR_all$x),
                       ratio = 1,
                       alternative = "two.sided")
armL_f_test <- var.test(x = wt_chrCov_armL_all$x,
                        y = mutant_chrCov_armL_all$x,
                        ratio = 1,
                        alternative = "two.sided")
armR_f_test <- var.test(x = wt_chrCov_armR_all$x,
                        y = mutant_chrCov_armR_all$x,
                        ratio = 1,
                        alternative = "two.sided")
peri_f_test <- var.test(x = wt_chrCov_peri_all$x,
                        y = mutant_chrCov_peri_all$x,
                        ratio = 1,
                        alternative = "two.sided")
save(arm_f_test,
     file = paste0(outDir,
                   wt_libName, "_",
                   mutant_libName,
                   "_chrCov_arm_f_test.RData"))
save(armL_f_test,
     file = paste0(outDir,
                   wt_libName, "_",
                   mutant_libName,
                   "_chrCov_armL_f_test.RData"))
save(armR_f_test,
     file = paste0(outDir,
                   wt_libName, "_",
                   mutant_libName,
                   "_chrCov_armR_f_test.RData"))
save(peri_f_test,
     file = paste0(outDir,
                   wt_libName, "_",
                   mutant_libName,
                   "_chrCov_peri_f_test.RData"))

arm_t_test <- t.test(x = c(wt_chrCov_armL_all$x, wt_chrCov_armR_all$x),
                     y = c(mutant_chrCov_armL_all$x, mutant_chrCov_armR_all$x),
                     alternative = "less",
                     paired = FALSE,
                     var.equal = TRUE)
armL_t_test <- t.test(x = wt_chrCov_armL_all$x,
                      y = mutant_chrCov_armL_all$x,
                      alternative = "less",
                      paired = FALSE,
                      var.equal = TRUE)
armR_t_test <- t.test(x = wt_chrCov_armR_all$x,
                      y = mutant_chrCov_armR_all$x,
                      alternative = "less",
                      paired = FALSE,
                      var.equal = TRUE)
peri_t_test <- t.test(x = wt_chrCov_peri_all$x,
                      y = mutant_chrCov_peri_all$x,
                      alternative = "greater",
                      paired = FALSE,
                      var.equal = TRUE)
load(
     file = paste0(outDir,
                   wt_libName, "_",
                   mutant_libName,
                   "_chrCov_arm_t_test.RData"))
load(
     file = paste0(outDir,
                   wt_libName, "_",
                   mutant_libName,
                   "_chrCov_armL_t_test.RData"))
load(
     file = paste0(outDir,
                   wt_libName, "_",
                   mutant_libName,
                   "_chrCov_armR_t_test.RData"))
load(
     file = paste0(outDir,
                   wt_libName, "_",
                   mutant_libName,
                   "_chrCov_peri_t_test.RData"))

