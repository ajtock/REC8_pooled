#!/bin/bash

## PeakRanger v1.18
## Use peakranger ranger to call narrow peaks in each REC8 replicate, using input from REC8_MYC_Rep2 as an input control sample
## From the PeakRanger manual at http://ranger.sourceforge.net/manual1.18.html :

# Usage:
# csmit -m 24G -c 24 "./REC8_peaks_peakranger_ranger_REC8_HA_Rep2_ChIP.sh 0.001 0.01" 

### The algorithm

## ranger uses a staged algorithm to discover enriched regions and the summits within them. In the first step, PeakRanger implements a FDR based adapative thresholding algorithm, which was originally proposed by PeakSeq. PeakRanger uses this thresholder to find regions with enriched reads that exceed expects. After that, PeakRanger searches for summits in these regions. The summit-search algorithm first looks for the location with largest number of reads. It then searchs for sub-summits with the sensitivity, the delta -r, specified by the user. Smaller -r will generate more summits.The coverage profiles are smoothed and padded before calling summits. The smoothing grade varies with -b. Higher smoothing bandwidth results less false summits at the cost of degraded summit accuracy .To measure the significance of the enriched regions, PeakRanger uses binormial distribution to model the relative enrichment of sample over control. A p value is generated as a result. Users can thus select highly significant peaks by using a smaller -p.

bamDir=/home/meiosis/ajt200/analysis/REC8_pooled

pval=$1
qval=$2

for i in REC8_HA_Rep2
do
  peakranger ranger -d $bamDir/${i}_ChIP_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam \
                    -c $bamDir/REC8_MYC_Rep2_input_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam \
                    --format bam -o ${i}"_peaks_peakranger_ranger_p"$pval"_q"$qval \
                    -p $pval -q $qval -l 200 --pad -t 24 --verbose
done

