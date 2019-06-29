#!/bin/bash

## Use peakranger ranger (v1.18) to call narrow peaks in each REC8 replicate
## PeakRanger manual at http://ranger.sourceforge.net/manual1.18.html

# Usage:
# csmit -m 10G -c 8 "./1_REC8_peaks_peakranger_ranger.sh '/home/ajt200/analysis/REC8_pooled' '/home/ajt200/analysis/REC8_pooled' REC8_MYC_Rep1_ChIP REC8_MYC_Rep1_input 0.001 0.01 200 8" 

ChIP_bamDir=$1
input_bamDir=$2
ChIP_prefix=$3
input_prefix=$4
pval=$5
qval=$6
ext_length=$7
threads=$8

/home/ajt200/tools/PeakRanger-1.18/bin/peakranger ranger \
  --data ${ChIP_bamDir}/${ChIP_prefix}_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam \
  --control ${input_bamDir}/${input_prefix}_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam \
  --format bam \
  --output ${ChIP_prefix}"_rangerPeaks" \
  --pval ${pval} \
  --FDR ${qval} \
  --ext_length ${ext_length} \
  --pad \
  --thread ${threads} \
  --verbose
