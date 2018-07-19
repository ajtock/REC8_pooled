#!/applications/R/R-3.3.2/bin/Rscript

# Adjust "signalValue" (region treatment reads; column 7) by region control reads
# signalValue = (region treatment reads+1)/(region control reads+1)
# Use this signalValue with caution - unknown whether these read counts are
# normalised by library size by PeakRanger
# May be preferrable to use qValue (FDR) or pValue in downstream analyses,
# including IDR calculation

# -log10 transform p-values and q-values

# Usage:
# ./REC8_narrowPeak_sigVal_log2TreadsNormCreads_minuslog10PQ_REC8_HA_Rep2_ChIP.R 0.001 0.01

args <- commandArgs(trailingOnly = TRUE)
pval <- args[1]
qval <- args[2]

peakDir <- paste0("/home/meiosis/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p",
                  pval, "_q", qval, "/")

for(i in c("REC8_HA_Rep2")) {
  peaks <- read.table(file = paste0(peakDir, i,"_peaks_peakranger_ranger_p",
                                    pval, "_q", qval,
                                    "_Treads_Creads.narrowPeak.UntransformedPQ"))
  colnames(peaks) <- c("chr", "start0based", "end",
                       "name", "score", "strand",
                       "signalVal", "pValUntrans", "qValUntrans",
                       "summit0based", "treads", "creads")
  peaks <- cbind(peaks[,1:6], log2((peaks[,11]+1)/(peaks[,12]+1)),
                 -log10(peaks[,8]), -log10(peaks[,9]), peaks[,10])
  colnames(peaks) <- c("chr", "start0based", "end",
                       "name", "score", "strand",
                       "TreadsNormCreads", "pVal", "qVal",
                       "summit0based")
  peaks$pVal[which(!is.finite(peaks$pVal))] <- 323
  peaks$qVal[which(!is.finite(peaks$qVal))] <- 323
  head(peaks)
  write.table(peaks, file = paste0(peakDir, i, "_peaks_peakranger_ranger_p",
                                   pval, "_q", qval,
                                   "_log2TreadsNormCreads.narrowPeak"),
              col.names = F, row.names = F, sep = "\t", quote = F)
}

