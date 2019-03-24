#!/bin/bash

# Usage:
# csmit -m 50G -c 32 "bash ./plotFingerprint.sh REC8_HA_Rep1_ChIP REC8_HA_Rep2_ChIP REC8_MYC_Rep1_ChIP kss_REC8_HA_Rep1_ChIP REC8_MYC_Rep2_input REC8_HA_Rep2_input REC8_MYC_Rep1_input kss_REC8_HA_Rep1_input wt_REC8_HA_Rep1_ChIP wt_REC8_HA_Rep2_ChIP wt_REC8_MYC_Rep1_ChIP kss_REC8_HA_Rep1_ChIP wt_REC8_MYC_Rep2_input wt_REC8_HA_Rep2_input wt_REC8_MYC_Rep1_input kss_REC8_HA_Rep1_input 32"

dat1=${1}
dat2=${2}
dat3=${3}
dat4=${4}
dat5=${5}
dat6=${6}
dat7=${7}
dat8=${8}
dat1label=${9}
dat2label=${10}
dat3label=${11}
dat4label=${12}
dat5label=${13}
dat6label=${14}
dat7label=${15}
dat8label=${16}
threads=${17}

plotFingerprint -b ${dat1}_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam \
                   ${dat2}_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam \
                   ${dat3}_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam \
                   ${dat4}_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam \
                   ${dat5}_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam \
                   ${dat6}_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam \
                   ${dat7}_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam \
                   ${dat8}_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam \
                -o REC8_ChIP_input_fingerprints.pdf \
                --extendReads \
                --labels ${dat1label} ${dat2label} ${dat3label} ${dat4label} ${dat5label} ${dat6label} ${dat7label} ${dat8label} \
                --outRawCounts REC8_ChIP_input_fingerprints.tab \
                -p ${threads}
