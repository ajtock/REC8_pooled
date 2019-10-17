#!/applications/R/R-3.3.2/bin/Rscriot

# Calculate Spearman's rank-order correlation coefficients
# for mean REC8 and mean expression in gene bodies (TSS to TTS),
# for mean REC8 in gene bodies and mean SPO11-1-oligos in gene promoters (TSS-1 to TSS-500), and
# for mean REC8 in gene bodies and mean SPO11-1-oligos in gene terminators (TTS+1 to TTS+500)

library(corrplot)

REC8_HA_Rep2 <- read.table("log2_REC8_HA_Rep2_ChIP_REC8_MYC_Rep1_input_in_gene_bodies/log2_REC8_HA_Rep2_ChIP_REC8_MYC_Rep1_input_in_gene_bodies_meanCov.txt", header = T)
RNAseq_Chris_Rep1 <- read.table("WT_RNAseq_Chris_Rep1_in_gene_bodies/WT_RNAseq_Chris_Rep1_in_gene_bodies_meanCov.txt", header = T)
RNAseq_Chris_Rep2 <- read.table("WT_RNAseq_Chris_Rep2_in_gene_bodies/WT_RNAseq_Chris_Rep2_in_gene_bodies_meanCov.txt", header = T)
RNAseq_meiocyte_Rep1 <- read.table("WT_RNAseq_meiocyte_Rep1_in_gene_bodies/WT_RNAseq_meiocyte_Rep1_in_gene_bodies_meanCov.txt", header = T)
RNAseq_meiocyte_Rep2 <- read.table("WT_RNAseq_meiocyte_Rep2_in_gene_bodies/WT_RNAseq_meiocyte_Rep2_in_gene_bodies_meanCov.txt", header = T)
RNAseq_meiocyte_Rep3 <- read.table("WT_RNAseq_meiocyte_Rep3_in_gene_bodies/WT_RNAseq_meiocyte_Rep3_in_gene_bodies_meanCov.txt", header = T)
SPO11oligos_Rep1_promoters <- read.table("log2wtSPO11oligoRPI1NakedDNA_in_gene_promoters/log2wtSPO11oligoRPI1NakedDNA_in_gene_promoters_meanCov.txt", header = T)
SPO11oligos_Rep2_promoters <- read.table("log2wtSPO11oligoRPI3NakedDNA_in_gene_promoters/log2wtSPO11oligoRPI3NakedDNA_in_gene_promoters_meanCov.txt", header = T)
SPO11oligos_Rep3_promoters <- read.table("log2wtSPO11oligoRPI8NakedDNA_in_gene_promoters/log2wtSPO11oligoRPI8NakedDNA_in_gene_promoters_meanCov.txt", header = T)
SPO11oligos_Rep1_terminators <- read.table("log2wtSPO11oligoRPI1NakedDNA_in_gene_terminators/log2wtSPO11oligoRPI1NakedDNA_in_gene_terminators_meanCov.txt", header = T)
SPO11oligos_Rep2_terminators <- read.table("log2wtSPO11oligoRPI3NakedDNA_in_gene_terminators/log2wtSPO11oligoRPI3NakedDNA_in_gene_terminators_meanCov.txt", header = T)
SPO11oligos_Rep3_terminators <- read.table("log2wtSPO11oligoRPI8NakedDNA_in_gene_terminators/log2wtSPO11oligoRPI8NakedDNA_in_gene_terminators_meanCov.txt", header = T)

allDF <- data.frame(
                    REC8_HA_Rep2$normCov_targetCov,
                    RNAseq_Chris_Rep1$normCov_targetCov,
                    RNAseq_Chris_Rep2$normCov_targetCov,
                    RNAseq_meiocyte_Rep1$normCov_targetCov,
                    RNAseq_meiocyte_Rep2$normCov_targetCov,
                    RNAseq_meiocyte_Rep3$normCov_targetCov,
                    SPO11oligos_Rep1_promoters$normCov_targetCov,
                    SPO11oligos_Rep2_promoters$normCov_targetCov,
                    SPO11oligos_Rep3_promoters$normCov_targetCov,
                    SPO11oligos_Rep1_terminators$normCov_targetCov,
                    SPO11oligos_Rep2_terminators$normCov_targetCov,
                    SPO11oligos_Rep3_terminators$normCov_targetCov
                   )
colnames(allDF) <- c(
                     "REC8-HA",
                     "RNA-seq Rep1 (floral)",
                     "RNA-seq Rep2 (floral)",
                     "RNA-seq Rep1 (meiocytes)",
                     "RNA-seq Rep2 (meiocytes)",
                     "RNA-seq Rep3 (meiocytes)",
                     "SPO11-1-oligos Rep1 (promoters)",
                     "SPO11-1-oligos Rep2 (promoters)",
                     "SPO11-1-oligos Rep3 (promoters)",
                     "SPO11-1-oligos Rep1 (terminators)",
                     "SPO11-1-oligos Rep2 (terminators)",
                     "SPO11-1-oligos Rep3 (terminators)"
                    )
allDF_corrMat <- cor(allDF, method = "spearman", use = "pairwise.complete.obs")
col1 <- colorRampPalette(c("blue", "white", "red"))
pdf("Spearman_corrplot_REC8_HA_Rep2_in_gene_bodies_vs_RNAseq_floral_and_meiocyte_in_gene_bodies_and_SPO11_1_oligos_in_gene_promoters_and_terminators.pdf",
    height = 12, width = 11)
corrplot(allDF_corrMat, method = "color", type = "upper", col = col1(20), tl.col = "black",
         addgrid.col = "white", addCoef.col = "grey90", mar = c(0,0,1,0),
         tl.cex = 1.5, cl.cex = 1.5, number.cex = 1.0,
         title = expression(italic("r"[s]) ~ "for average coverage in gene bodies, promoters and terminators"))
dev.off()


#REC8_RNAseq <- merge(x = REC8,
#                     y = RNAseq,
#                     by.x = "gene_model",
#                     by.y = "gene_model")
#
#print(cor.test(REC8_RNAseq$REC8_norm_geneCov,
#               REC8_RNAseq$RNAseq_norm_geneCov,
#               method = "spearman"))
##
##	Spearman's rank correlation rho
##
##data:  REC8_RNAseq$REC8_norm_geneCov and REC8_RNAseq$RNAseq_norm_geneCov
##S = 4.8305e+12, p-value < 2.2e-16
##alternative hypothesis: true rho is not equal to 0
##sample estimates:
##       rho 
##-0.4396011 
##
##Warning message:
##In cor.test.default(REC8_RNAseq$REC8_norm_geneCov, REC8_RNAseq$RNAseq_norm_geneCov,  :
##  Cannot compute exact p-value with ties

