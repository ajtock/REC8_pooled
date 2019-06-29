#!/applications/R/R-3.5.0/bin/Rscript

# Calculate wt and kss REC8 mean coverage over qPCR loci:
#Peak 1=Chr5: 11,750,429–11,750,543. 
#Peak 2= Chr1: 15,449,849–15,449,964. 
#Peak 3=Chr2: 660,711–660,826
#Peak 4=chromosome 3: 14,248,393–14,248,483
#Peak 5=chromosome 3: 12,076,700–12,076,835
#Peak 6=chromosome 3: 14,082,421–14,082,549
#Peak 7=chromosome 1: 1,391,551–1,391,664

library(seqmentSeq)
library(genomation)
library(parallel)

chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")

# Create GRanges object of qPCR loci
qPCRloci <- GRanges(seqnames = c("Chr5",
                                 "Chr1",
                                 "Chr2",
                                 "Chr3",
                                 "Chr3",
                                 "Chr3",
                                 "Chr1"),
                    ranges = IRanges(start = c(11750429,
                                               15449849,
                                               660711,
                                               14248393,
                                               12076700,
                                               14082421,
                                               1391551),
                                     end = c(11750543,
                                             15449964,
                                             660826,
                                             14248483,
                                             12076835,
                                             14082549,
                                             1391664)),
                    strand = "*")

## ChIP
ChIP_paths <- c("/home/ajt200/analysis/REC8_pooled/coverage/REC8_HA_Rep1_ChIP_norm_allchrs_coverage_coord_tab.bed",
                "/home/ajt200/analysis/180622_Chris_lambing_ChIP_REC8_HA_Col_kss/WT/coverage/REC8_HA_Rep2_ChIP_norm_allchrs_coverage_coord_tab.bed",
                "/home/ajt200/analysis/REC8_pooled/coverage/REC8_MYC_Rep1_ChIP_norm_allchrs_coverage_coord_tab.bed",
                "/home/ajt200/analysis/180622_Chris_lambing_ChIP_REC8_HA_Col_kss/kss/coverage/kss_REC8_HA_Rep1_ChIP_norm_allchrs_coverage_coord_tab.bed",
                "/home/ajt200/analysis/180830_Chris_lambing_ChIP_rec8HA_kss_rep2/coverage/kss_REC8_HA_Rep2_ChIP_norm_allchrs_coverage_coord_tab.bed")

ChIP_list <- mclapply(seq_along(ChIP_paths), function(x) {
  readGeneric(ChIP_paths[x], meta.col = list(coverage = 4))
}, mc.cores = length(ChIP_paths))

# Obtain indices of overlapping coverage values
ChIP_targetOverlaps_list <- mclapply(seq_along(ChIP_list), function(x) {
   getOverlaps(coordinates = qPCRloci,
               segments = ChIP_list[[x]],
               overlapType = "overlapping",
               whichOverlaps = TRUE,
               ignoreStrand = TRUE)
}, mc.cores = length(ChIP_list))
# Calculate means and SDs
ChIP_targetCovMean_list <- mclapply(seq_along(ChIP_list), function(x) {
  sapply(ChIP_targetOverlaps_list[[x]],
         function(y) mean(ChIP_list[[x]]$coverage[y]))
}, mc.cores = length(ChIP_list))
ChIP_targetCovSD_list <- mclapply(seq_along(ChIP_list), function(x) {
  sapply(ChIP_targetOverlaps_list[[x]],
         function(y) sd(ChIP_list[[x]]$coverage[y]))
}, mc.cores = length(ChIP_list))

# Create and write data.frame
ChIP_df <- data.frame(qPCR_locus = 1:7,
                      chr = as.character(seqnames(qPCRloci)),
                      start = as.integer(start(qPCRloci)),
                      end = as.integer(end(qPCRloci)),
                      wt_REC8_HA_Rep1_ChIP_mean = ChIP_targetCovMean_list[[1]],
                      wt_REC8_HA_Rep1_ChIP_SD = ChIP_targetCovSD_list[[1]],
                      wt_REC8_HA_Rep2_ChIP_mean = ChIP_targetCovMean_list[[2]],
                      wt_REC8_HA_Rep2_ChIP_SD = ChIP_targetCovSD_list[[2]],
                      wt_REC8_Myc_Rep1_ChIP_mean = ChIP_targetCovMean_list[[3]],
                      wt_REC8_Myc_Rep1_ChIP_SD = ChIP_targetCovSD_list[[3]],
                      kss_REC8_HA_Rep1_ChIP_mean = ChIP_targetCovMean_list[[4]],
                      kss_REC8_HA_Rep1_ChIP_SD = ChIP_targetCovSD_list[[4]],
                      kss_REC8_HA_Rep2_ChIP_mean = ChIP_targetCovMean_list[[5]],
                      kss_REC8_HA_Rep2_ChIP_SD = ChIP_targetCovSD_list[[5]])
write.table(ChIP_df,
            file = paste0("REC8_qPCR_loci_wt_and_kss_ChIP_means_and_SDs.tsv"),
            sep = "\t", row.names = F, col.names = T, quote = F)


## log2(ChIP/input)
log2_paths <- c("/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/noZscore/log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep1_input_norm_allchrs_coverage_coord_tab_noZscore.bed",
                "/home/ajt200/analysis/180622_Chris_lambing_ChIP_REC8_HA_Col_kss/WT/coverage/log2ChIPinput/REC8_MYC_Rep1_input/noZscore/log2_REC8_HA_Rep2_ChIP_REC8_MYC_Rep1_input_norm_allchrs_coverage_coord_tab_noZscore.bed",
                "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/noZscore/log2_REC8_MYC_Rep1_ChIP_REC8_MYC_Rep1_input_norm_allchrs_coverage_coord_tab_noZscore.bed",
                "/home/ajt200/analysis/180622_Chris_lambing_ChIP_REC8_HA_Col_kss/kss/coverage/log2ChIPinput/replicate_specific_input/noZscore/log2_kss_REC8_HA_Rep1_ChIP_kss_REC8_HA_Rep1_input_norm_allchrs_coverage_coord_tab_noZscore.bed",
                "/home/ajt200/analysis/180830_Chris_lambing_ChIP_rec8HA_kss_rep2/coverage/log2ChIPinput/kss_REC8_HA_Rep1_input/noZscore/log2_kss_REC8_HA_Rep2_ChIP_kss_REC8_HA_Rep1_input_norm_allchrs_coverage_coord_tab_noZscore.bed")

log2_list <- mclapply(seq_along(log2_paths), function(x) {
  readGeneric(log2_paths[x], meta.col = list(coverage = 4))
}, mc.cores = length(log2_paths))

# Obtain indices of overlapping coverage values
log2_targetOverlaps_list <- mclapply(seq_along(log2_list), function(x) {
   getOverlaps(coordinates = qPCRloci,
               segments = log2_list[[x]],
               overlapType = "overlapping",
               whichOverlaps = TRUE,
               ignoreStrand = TRUE)
}, mc.cores = length(log2_list))
# Calculate means and SDs
log2_targetCovMean_list <- mclapply(seq_along(log2_list), function(x) {
  sapply(log2_targetOverlaps_list[[x]],
         function(y) mean(log2_list[[x]]$coverage[y]))
}, mc.cores = length(log2_list))
log2_targetCovSD_list <- mclapply(seq_along(log2_list), function(x) {
  sapply(log2_targetOverlaps_list[[x]],
         function(y) sd(log2_list[[x]]$coverage[y]))
}, mc.cores = length(log2_list))

log2_df <- data.frame(qPCR_locus = 1:7,
                      chr = as.character(seqnames(qPCRloci)),
                      start = as.integer(start(qPCRloci)),
                      end = as.integer(end(qPCRloci)),
                      wt_REC8_HA_Rep1_log2_mean = log2_targetCovMean_list[[1]],
                      wt_REC8_HA_Rep1_log2_SD = log2_targetCovSD_list[[1]],
                      wt_REC8_HA_Rep2_log2_mean = log2_targetCovMean_list[[2]],
                      wt_REC8_HA_Rep2_log2_SD = log2_targetCovSD_list[[2]],
                      wt_REC8_Myc_Rep1_log2_mean = log2_targetCovMean_list[[3]],
                      wt_REC8_Myc_Rep1_log2_SD = log2_targetCovSD_list[[3]],
                      kss_REC8_HA_Rep1_log2_mean = log2_targetCovMean_list[[4]],
                      kss_REC8_HA_Rep1_log2_SD = log2_targetCovSD_list[[4]],
                      kss_REC8_HA_Rep2_log2_mean = log2_targetCovMean_list[[5]],
                      kss_REC8_HA_Rep2_log2_SD = log2_targetCovSD_list[[5]])
write.table(log2_df,
            file = paste0("REC8_qPCR_loci_wt_and_kss_log2_means_and_SDs.tsv"),
            sep = "\t", row.names = F, col.names = T, quote = F)


## Z-score log2(ChIP/input)
Zlog2_paths <- c("/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep1_input_norm_allchrs_coverage_coord_tab.bed",
                 "/home/ajt200/analysis/180622_Chris_lambing_ChIP_REC8_HA_Col_kss/WT/coverage/log2ChIPinput/REC8_MYC_Rep1_input/log2_REC8_HA_Rep2_ChIP_REC8_MYC_Rep1_input_norm_allchrs_coverage_coord_tab.bed",
                 "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/log2_REC8_MYC_Rep1_ChIP_REC8_MYC_Rep1_input_norm_allchrs_coverage_coord_tab.bed",
                 "/home/ajt200/analysis/180622_Chris_lambing_ChIP_REC8_HA_Col_kss/kss/coverage/log2ChIPinput/replicate_specific_input/log2_kss_REC8_HA_Rep1_ChIP_kss_REC8_HA_Rep1_input_norm_allchrs_coverage_coord_tab.bed",
                 "/home/ajt200/analysis/180830_Chris_lambing_ChIP_rec8HA_kss_rep2/coverage/log2ChIPinput/kss_REC8_HA_Rep1_input/log2_kss_REC8_HA_Rep2_ChIP_kss_REC8_HA_Rep1_input_norm_allchrs_coverage_coord_tab.bed")

Zlog2_list <- mclapply(seq_along(Zlog2_paths), function(x) {
  readGeneric(Zlog2_paths[x], meta.col = list(coverage = 4))
}, mc.cores = length(Zlog2_paths))

# Obtain indices of overlapping coverage values
Zlog2_targetOverlaps_list <- mclapply(seq_along(Zlog2_list), function(x) {
   getOverlaps(coordinates = qPCRloci,
               segments = Zlog2_list[[x]],
               overlapType = "overlapping",
               whichOverlaps = TRUE,
               ignoreStrand = TRUE)
}, mc.cores = length(Zlog2_list))
# Calculate means and SDs
Zlog2_targetCovMean_list <- mclapply(seq_along(Zlog2_list), function(x) {
  sapply(Zlog2_targetOverlaps_list[[x]],
         function(y) mean(Zlog2_list[[x]]$coverage[y]))
}, mc.cores = length(Zlog2_list))
Zlog2_targetCovSD_list <- mclapply(seq_along(Zlog2_list), function(x) {
  sapply(Zlog2_targetOverlaps_list[[x]],
         function(y) sd(Zlog2_list[[x]]$coverage[y]))
}, mc.cores = length(Zlog2_list))

Zlog2_df <- data.frame(qPCR_locus = 1:7,
                       chr = as.character(seqnames(qPCRloci)),
                       start = as.integer(start(qPCRloci)),
                       end = as.integer(end(qPCRloci)),
                       wt_REC8_HA_Rep1_Zlog2_mean = Zlog2_targetCovMean_list[[1]],
                       wt_REC8_HA_Rep1_Zlog2_SD = Zlog2_targetCovSD_list[[1]],
                       wt_REC8_HA_Rep2_Zlog2_mean = Zlog2_targetCovMean_list[[2]],
                       wt_REC8_HA_Rep2_Zlog2_SD = Zlog2_targetCovSD_list[[2]],
                       wt_REC8_Myc_Rep1_Zlog2_mean = Zlog2_targetCovMean_list[[3]],
                       wt_REC8_Myc_Rep1_Zlog2_SD = Zlog2_targetCovSD_list[[3]],
                       kss_REC8_HA_Rep1_Zlog2_mean = Zlog2_targetCovMean_list[[4]],
                       kss_REC8_HA_Rep1_Zlog2_SD = Zlog2_targetCovSD_list[[4]],
                       kss_REC8_HA_Rep2_Zlog2_mean = Zlog2_targetCovMean_list[[5]],
                       kss_REC8_HA_Rep2_Zlog2_SD = Zlog2_targetCovSD_list[[5]])
write.table(Zlog2_df,
            file = paste0("REC8_qPCR_loci_wt_and_kss_Zlog2_means_and_SDs.tsv"),
            sep = "\t", row.names = F, col.names = T, quote = F)

