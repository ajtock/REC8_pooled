#!/bin/bash

# Convert PeakRanger peak loci details file to narrowPeak format

for i in kss_REC8_HA_Rep1
do
  tail -n +24 ${i}_peaks_peakranger_ranger_p0.001_q0.01_details > ${i}_peaks_peakranger_ranger_p0.001_q0.01_noHead
  grep -v 'chloroplast\|mitochondria' ${i}_peaks_peakranger_ranger_p0.001_q0.01_noHead > ${i}_peaks_peakranger_ranger_p0.001_q0.01_noHead_tmp1
  awk 'BEGIN {OFS="\t"}; {$1 = "chr"$1; print}' ${i}_peaks_peakranger_ranger_p0.001_q0.01_noHead_tmp1 > ${i}_peaks_peakranger_ranger_p0.001_q0.01_noHead_tmp2
  awk 'BEGIN {OFS="\t"}; {print $1, $2, $3, $8, $8, $8, $8, $6, $7, $5, $9, $10}' ${i}_peaks_peakranger_ranger_p0.001_q0.01_noHead_tmp2 > ${i}_peaks_peakranger_ranger_p0.001_q0.01_noHead_tmp3
  awk 'BEGIN {OFS="\t"}; {$4 = "."; $5 = "."; $6 = "."; $7 = "."; $10 = $10-$2; print}' ${i}_peaks_peakranger_ranger_p0.001_q0.01_noHead_tmp3 > ${i}_peaks_peakranger_ranger_p0.001_q0.01_Treads_Creads.narrowPeak.UntransformedPQ
  rm ${i}_peaks_peakranger_ranger_p0.001_q0.01_noHead*
done

