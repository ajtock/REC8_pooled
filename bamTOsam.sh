#!/bin/bash
# use "-h" flag to include SAM header in output

for i in REC8_HA_Rep1_ChIP REC8_MYC_Rep1_ChIP REC8_MYC_Rep2_ChIP 
do
( samtools view -o $samDir/${i}_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.sam ${i}_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam ) &
done
wait

