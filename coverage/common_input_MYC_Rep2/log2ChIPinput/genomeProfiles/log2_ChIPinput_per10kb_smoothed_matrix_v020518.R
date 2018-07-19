################################################################################
# Generate dataset-wide matrix of log2 ratio of ChIP and input coverage values #
# and features in 10-kb windows                                                #
################################################################################

# genomic definitions
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
# pericentromeric regions are as defined in Supplemental Table S26 of Ziolkowski et al. (2017) Genes Dev. 31
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)

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

winNames <- c("10kb")

outDir1 <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/genomeProfiles/"
ChIPnames1 <- "REC8_HA_Rep1_ChIP"
names1 <- "REC8-HA_Rep1"

outDir2 <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/genomeProfiles/"
ChIPnames2 <- "REC8_MYC_Rep2_ChIP"
names2 <- "REC8-MYC_Rep2"

outDir3 <- "/home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/genomeProfiles/"
ChIPnames3 <- "REC8_MYC_Rep1_ChIP"
names3 <- "REC8-MYC_Rep1"

outDir4 <- "/projects/ajt200/BAM_masters/nucleosomes/WT/coverage/nakedDNA_untrimmed_input/log2ChIPinput/genomeProfiles/"
ChIPnames4 <- "WT_nuc"
names4 <- "MNase"

outDir5 <- "/projects/ajt200/BAM_masters/SPO11_ChIP/WT/coverage/REC8_MYC_Rep2_input/log2ChIPinput/genomeProfiles/"
ChIPnames5 <- "WT_SPO11_ChIP4"
names5 <- "SPO11-1_ChIP4"

outDir6 <- "/projects/ajt200/BAM_masters/H2A/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames6 <- "H2AZ_ChIP"
names6 <- "H2A.Z"

outDir7 <- "/projects/ajt200/BAM_masters/H3K4me3/replicates/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames7 <- "WT_H3K4me3_ChIP14"
names7 <- "H3K4me3_ChIP14"

outDir8 <- "/projects/ajt200/BAM_masters/H2A/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames8 <- "H2AW_ChIP"
names8 <- "H2A.W"

outDir9 <- "/projects/ajt200/BAM_masters/H3K9me2/WT/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames9 <- "WT_H3K9me2_ChIP"
names9 <- "H3K9me2"

outDir10 <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K4me1/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames10 <- "WT_H3K4me1_Rep1_ChIP"
names10 <- "H3K4me1"

outDir11 <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K4me2/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames11 <- "WT_H3K4me2_Rep1_ChIP"
names11 <- "H3K4me2"

outDir12 <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K27me1/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames12 <- "WT_H3K27me1_Rep1_ChIP"
names12 <- "H3K27me1"

outDir13 <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K27me3/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames13 <- "WT_H3K27me3_Rep1_ChIP"
names13 <- "H3K27me3"

outDir14 <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames14 <- "WT_SPO11_oligo_RPI1"
names14 <- "SPO11-1-oligos_RPI1"

outDir15 <- "/projects/ajt200/BAM_masters/MSH4/WT/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames15 <- "WT_MSH4_ChIP"
names15 <- "MSH4"

outDir16 <- "/projects/ajt200/BAM_masters/PolIV_Law_Jacobsen_2013_Nature/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames16 <- "PolIV_Rep2_ChIP"
names16 <- "Pol_IV"

outDir17 <- "/projects/ajt200/BAM_masters/PolV_Liu_Jacobsen_2018_NatPlants/coverage/log2ChIPinput/genomeProfiles/"
ChIPnames17 <- "PolV_ChIP"
names17 <- "Pol_V"

outDir18 <- "/projects/ajt200/BAM_masters/RNAseq_meiocyte_Walker_Feng_2018_NatGenet/WT/meiocyte/coverage/genomeProfiles/"
ChIPnames18 <- "WT_RNAseq_meiocyte_Rep1"
names18 <- "wt_RNA-seq_meiocyte_Rep1"

outDir19 <- "/projects/ajt200/BAM_masters/RNAseq_meiocyte_Walker_Feng_2018_NatGenet/WT/meiocyte/coverage/genomeProfiles/"
ChIPnames19 <- "WT_RNAseq_meiocyte_Rep2"
names19 <- "wt_RNA-seq_meiocyte_Rep2"

outDir20 <- "/projects/ajt200/BAM_masters/RNAseq_meiocyte_Walker_Feng_2018_NatGenet/WT/meiocyte/coverage/genomeProfiles/"
ChIPnames20 <- "WT_RNAseq_meiocyte_Rep3"
names20 <- "wt_RNA-seq_meiocyte_Rep3"

outDir21 <- "/projects/ajt200/BAM_masters/RNAseq/WT/coverage/genomeProfiles/"
ChIPnames21 <- "WT_RNAseq_Rep1"
names21 <- "wt_RNA-seq_(Kyuha)_Rep1"

outDir22 <- "/home/meiosis/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/WT/coverage/genomeProfiles/"
ChIPnames22 <- "WT_RNAseq_Chris_Rep1"
names22 <- "wt_RNA-seq_(Chris)_Rep1"

outDir23 <- "/home/meiosis/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/WT/coverage/genomeProfiles/"
ChIPnames23 <- "WT_RNAseq_Chris_Rep2"
names23 <- "wt_RNA-seq_(Chris)_Rep2"

outDir24 <- "/home/meiosis/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/kss/coverage/genomeProfiles/"
ChIPnames24 <- "kss_RNAseq_Chris_Rep1"
names24 <- "kss_RNA-seq_(Chris)_Rep1"

outDir25 <- "/home/meiosis/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/kss/coverage/genomeProfiles/"
ChIPnames25 <- "kss_RNAseq_Chris_Rep2"
names25 <- "kss_RNA-seq_(Chris)_Rep2"

outDirCombined <- c(outDir1, outDir2, outDir3, outDir4, outDir5, outDir6, outDir7, outDir8, outDir9,
                    outDir10, outDir11, outDir12, outDir13, outDir14, outDir15, outDir16, outDir17,
                    outDir18, outDir19, outDir20, outDir21, outDir22, outDir23, outDir24, outDir25)
ChIPnamesCombined <- c(ChIPnames1, ChIPnames2, ChIPnames3, ChIPnames4, ChIPnames5, ChIPnames6,
                       ChIPnames7, ChIPnames8, ChIPnames9, ChIPnames10, ChIPnames11, ChIPnames12,
                       ChIPnames13, ChIPnames14, ChIPnames15, ChIPnames16, ChIPnames17,
                       ChIPnames18, ChIPnames19, ChIPnames20, ChIPnames21, ChIPnames22, ChIPnames23, ChIPnames24, ChIPnames25)
namesCombined <- c(names1, names2, names3, names4, names5, names6,
                   names7, names8, names9, names10, names11, names12,
                   names13, names14, names15, names16, names17,
                   names18, names19, names20, names21, names22, names23, names24, names25)

# Load concatenated filtered, log2-transformed coverage datasets
filt_log2trans_list <- lapply(seq_along(ChIPnamesCombined), function(k) {
  read.table(file = paste0(outDirCombined[k], "filt_log2_", ChIPnamesCombined[k], "_genome_norm_coverage_", winNames[1], ".txt"))
})
filt_noNA_log2trans_list <- lapply(seq_along(ChIPnamesCombined), function(k) {
  read.table(file = paste0(outDirCombined[k], "filt_noNA_log2_", ChIPnamesCombined[k], "_genome_norm_coverage_", winNames[1], ".txt"))
})

inDirFeatures <- "/projects/ajt200/REC8_MSH4/data_merged_fastq/coverage/log2ChIPinput/genomeProfiles/"
inDirATGC <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/REC8/coverage/common_input_MYC_Rep2/log2ChIPinput/genomeProfiles/"

# Load concatenated filtered gene dataset
filt_GeneDat <- read.table(file = paste0(inDirFeatures, "filt_gene_density_genome_", winNames[1], ".txt"))
filt_GeneDat_noNA <- read.table(file = paste0(inDirFeatures, "filt_noNA_gene_density_genome_", winNames[1], ".txt"))

# Load concatenated filtered TE dataset
filt_TEDat <- read.table(file = paste0(inDirFeatures, "filt_TE_density_genome_", winNames[1], ".txt"))
filt_TEDat_noNA <- read.table(file = paste0(inDirFeatures, "filt_noNA_TE_density_genome_", winNames[1], ".txt"))

# Load concatenated filtered CO density datasets
filt_CODat <- read.table(file = paste0(inDirFeatures, "filt_WTCO_density_genome_", winNames[1], ".txt"))
filt_CODat_noNA <- read.table(file = paste0(inDirFeatures, "filt_noNA_WTCO_density_genome_", winNames[1], ".txt"))

# Load concatenated filtered AT and GC content datasets
filt_ATcontent <- read.table(file = paste0(inDirATGC, "filt_AT_content_genome_", winNames[1], ".txt"))
filt_ATcontent_noNA <- read.table(file = paste0(inDirATGC, "filt_noNA_AT_content_genome_", winNames[1], ".txt"))
filt_GCcontent <- read.table(file = paste0(inDirATGC, "filt_GC_content_genome_", winNames[1], ".txt"))
filt_GCcontent_noNA <- read.table(file = paste0(inDirATGC, "filt_noNA_GC_content_genome_", winNames[1], ".txt"))

# Load concatenated filtered DNA meth datasets
filt_methDat <- read.table(file = paste0(inDirFeatures, "filt_wtMeth_genome_", winNames[1], ".txt"))
filt_methDat_noNA <- read.table(file = paste0(inDirFeatures, "filt_noNA_wtMeth_genome_", winNames[1], ".txt"))


# Create matrix of filtered values within 10-kb windows
matrix10kbtmp <- sapply(seq_along(filt_log2trans_list), function(x) cbind(filt_log2trans_list[[x]][,2]))
matrix10kb <- cbind(filt_log2trans_list[[1]][,1], matrix10kbtmp, filt_GeneDat[,2], filt_TEDat[,2], filt_CODat[,2], filt_ATcontent[,2], filt_GCcontent[,2], filt_methDat[,2:5])
colnames(matrix10kb) <- c("Windows", namesCombined, "Genes", "TEs", "Crossovers", "%AT", "%GC", "CG_meth", "CHG_meth", "CHH_meth", "DNA_meth_avg")
write.table(matrix10kb, file = paste0(outDir1, "REC8_matrix10kb_smoothed.txt"))

# Create matrix of filtered values (without NAs) within 10-kb windows
matrix10kbtmp_noNA <- sapply(seq_along(filt_noNA_log2trans_list), function(x) cbind(filt_noNA_log2trans_list[[x]][,1]))
matrix10kb_noNA <- cbind(filt_log2trans_list[[1]][,1], matrix10kbtmp_noNA, filt_GeneDat_noNA[,1], filt_TEDat_noNA[,1], filt_CODat_noNA[,1], filt_ATcontent_noNA[,1], filt_GCcontent_noNA[,1], filt_methDat_noNA[,2:5])
colnames(matrix10kb_noNA) <- c("Windows", namesCombined, "Genes", "TEs", "Crossovers", "%AT", "%GC", "CG_meth", "CHG_meth", "CHH_meth", "DNA_meth_avg")
write.table(matrix10kb_noNA, file = paste0(outDir1, "REC8_matrix10kb_smoothed_noNA.txt"))

