#############################################################################################
# Calculate average DNA methylation proportion values in windows of defined size            #
# and generate chromosome-scale plots                                                       #
#############################################################################################

library(segmentSeq)
library(genomation)
library(parallel)
library(doParallel)
registerDoParallel(cores = 5)
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

# genomic definitions
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
# pericentromeric regions are as defined in Supplemental Table S26 of Ziolkowski et al. (2017) Genes Dev. 31
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)

winSize <- c(200000)
winName <- c("200kb")

#################################
# make cumulative genomes       #
#################################

sumchr <- cumsum(c(0, chrLens))
print(sumchr)
sumchr_tot <- sumchr[length(sumchr)]
print(sumchr_tot)

centromeres <- sapply(seq_along(centromeres), function(x) {
  centromeres[x] + sumchr[x]
})
print(centromeres)
pericenStart <- sapply(seq_along(pericenStart), function(x) {
  pericenStart[x] + sumchr[x]
})
print(pericenStart)
pericenEnd <- sapply(seq_along(pericenEnd), function(x) {
  pericenEnd[x] + sumchr[x]
})
print(pericenEnd)

# DNA methylation
bedDir <- "/home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/"
outDir <- "/projects/ajt200/REC8_MSH4/data_merged_fastq/coverage/log2ChIPinput/genomeProfiles/"
GSM981060_suvh456_CG <- read.table(file = paste0(bedDir, "GSM981060_suvh456_CG.bed.gr.tab.bed"))
GSM981060_suvh456_CHG <- read.table(file = paste0(bedDir, "GSM981060_suvh456_CHG.bed.gr.tab.bed"))
GSM981060_suvh456_CHH <- read.table(file = paste0(bedDir, "GSM981060_suvh456_CHH.bed.gr.tab.bed"))

#cumMethDat <- NULL
foreach(i = 1:5, .combine = 'c') %dopar% {
  # define windows2 as GRanges object
  windows2 <- seq(1, chrLens[i], by = winSize)
  windows2 <- c(windows2, chrLens[i])
  cum_windows2 <- windows2 + sumchr[i]
  windows2_iranges <- IRanges(start = windows2, width = winSize)
  windows2_granges <- GRanges(seqnames = chrs[i], strand = "+", ranges = windows2_iranges)

  # calculate mean methylation levels (all contexts) within windows2
  kssCG <- GSM981060_suvh456_CG[which(GSM981060_suvh456_CG$V1 == paste0("chr", i)),]
  kssCG_ir_coords <- IRanges(start = kssCG$V2, width = 1)
  kssCG_gr_coords <- GRanges(seqnames = chrs[i], strand = "+", ranges = kssCG_ir_coords)
  overlapskssCG <- getOverlaps(windows2_granges, kssCG_gr_coords, whichOverlaps = TRUE)
  kssCG_win_vals <- sapply(overlapskssCG, function(x) mean(as.numeric(kssCG$V4[x])))

  kssCHG <- GSM981060_suvh456_CHG[which(GSM981060_suvh456_CHG$V1 == paste0("chr", i)),]
  kssCHG_ir_coords <- IRanges(start = kssCHG$V2, width = 1)
  kssCHG_gr_coords <- GRanges(seqnames = chrs[i], strand = "+", ranges = kssCHG_ir_coords)
  overlapskssCHG <- getOverlaps(windows2_granges, kssCHG_gr_coords, whichOverlaps = TRUE)
  kssCHG_win_vals <- sapply(overlapskssCHG, function(x) mean(as.numeric(kssCHG$V4[x])))

  kssCHH <- GSM981060_suvh456_CHH[which(GSM981060_suvh456_CHH$V1 == paste0("chr", i)),]
  kssCHH_ir_coords <- IRanges(start = kssCHH$V2, width = 1)
  kssCHH_gr_coords <- GRanges(seqnames = chrs[i], strand = "+", ranges = kssCHH_ir_coords)
  overlapskssCHH <- getOverlaps(windows2_granges, kssCHH_gr_coords, whichOverlaps = TRUE)
  kssCHH_win_vals <- sapply(overlapskssCHH, function(x) mean(as.numeric(kssCHH$V4[x])))

  methDat <- cbind(cum_windows2, kssCG_win_vals, kssCHG_win_vals, kssCHH_win_vals)
  methDat <- cbind(methDat, rowMeans(methDat[,2:4]))
  colnames(methDat) <- c("cum_windows2", "kssCG_win_vals", "kssCHG_win_vals", "kssCHH_win_vals", "allcntxts_win_vals")
  write.table(methDat, file = paste0(outDir, "kssMeth_chr", i, "_", winName, ".txt"))
  #cumMethDat <- rbind(cumMethDat, methDat)
}
#write.table(cumMethDat, file = paste0(outDir, "kssMeth_genome_200kb.txt"))

chrMethDat <- mclapply(1:5, function(i) {
  read.table(file = paste0(outDir, "kssMeth_chr", i, "_", winName, ".txt"))
}, mc.cores = 5)

# smooth coverage values with MA filter
test <- seq(1, 1000, by = 1)
j = 5
ma <- rep(1, test[j])/test[j]

cumFiltMethDat <- NULL
cumFiltMethDat_noNA <- NULL
for(i in 1:5) {
  filt_kssCG <- stats::filter(chrMethDat[[i]][,2], ma)
  which_na <- which(is.na(filt_kssCG) == TRUE)
  left_na <- which_na[which(which_na < 5)]
  left_val <- filt_kssCG[left_na[length(left_na)]+1]
  filt_kssCG[left_na] <- left_val
  right_na <- which_na[which(which_na > 5)]
  right_val <- filt_kssCG[right_na[1]-1]
  filt_kssCG[right_na] <- right_val
  filt_kssCG_noNA <- filt_kssCG[!is.na(filt_kssCG)]

  filt_kssCHG <- stats::filter(chrMethDat[[i]][,3], ma)
  which_na <- which(is.na(filt_kssCHG) == TRUE)
  left_na <- which_na[which(which_na < 5)]
  left_val <- filt_kssCHG[left_na[length(left_na)]+1]
  filt_kssCHG[left_na] <- left_val
  right_na <- which_na[which(which_na > 5)]
  right_val <- filt_kssCHG[right_na[1]-1]
  filt_kssCHG[right_na] <- right_val
  filt_kssCHG_noNA <- filt_kssCHG[!is.na(filt_kssCHG)]

  filt_kssCHH <- stats::filter(chrMethDat[[i]][,4], ma)
  which_na <- which(is.na(filt_kssCHH) == TRUE)
  left_na <- which_na[which(which_na < 5)]
  left_val <- filt_kssCHH[left_na[length(left_na)]+1]
  filt_kssCHH[left_na] <- left_val
  right_na <- which_na[which(which_na > 5)]
  right_val <- filt_kssCHH[right_na[1]-1]
  filt_kssCHH[right_na] <- right_val
  filt_kssCHH_noNA <- filt_kssCHH[!is.na(filt_kssCHH)]

  filt_allcntxts <- stats::filter(chrMethDat[[i]][,5], ma)
  which_na <- which(is.na(filt_allcntxts) == TRUE)
  left_na <- which_na[which(which_na < 5)]
  left_val <- filt_allcntxts[left_na[length(left_na)]+1]
  filt_allcntxts[left_na] <- left_val
  right_na <- which_na[which(which_na > 5)]
  right_val <- filt_allcntxts[right_na[1]-1]
  filt_allcntxts[right_na] <- right_val
  filt_allcntxts_noNA <- filt_allcntxts[!is.na(filt_allcntxts)]

  filtMethDat <- cbind(chrMethDat[[i]][,1], filt_kssCG, filt_kssCHG, filt_kssCHH, filt_allcntxts)
  filtMethDat_noNA <- cbind(chrMethDat[[i]][,1], filt_kssCG_noNA, filt_kssCHG_noNA, filt_kssCHH_noNA, filt_allcntxts_noNA)
  write.table(filtMethDat, file = paste0(outDir, "filt_kssMeth_chr", i, "_200kb.txt"))
  write.table(filtMethDat_noNA, file = paste0(outDir, "filt_noNA_kssMeth_chr", i, "_200kb.txt"))
  cumFiltMethDat <- rbind(cumFiltMethDat, filtMethDat)
  cumFiltMethDat_noNA <- rbind(cumFiltMethDat_noNA, filtMethDat_noNA)
}
write.table(cumFiltMethDat, file = paste0(outDir, "filt_kssMeth_genome_200kb.txt")) 
write.table(cumFiltMethDat_noNA, file = paste0(outDir, "filt_noNA_kssMeth_genome_200kb.txt"))

