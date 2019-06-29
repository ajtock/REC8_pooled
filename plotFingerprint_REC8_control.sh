#!/bin/bash

# Usage:
# csmit -m 50G -c 32 "bash ./plotFingerprint_REC8_control.sh REC8_HA_Rep1_ChIP REC8_HA_Rep2_ChIP REC8_MYC_Rep1_ChIP control_HA_Rep1_ChIP control_MYC_Rep1_ChIP REC8_MYC_Rep1_input REC8_MYC_Rep2_input REC8_HA_Rep2_input control_Rep1_input wt_REC8_HA_Rep1_ChIP wt_REC8_HA_Rep2_ChIP wt_REC8_Myc_Rep1_ChIP wt_control_HA_Rep1_ChIP wt_control_Myc_Rep1_ChIP wt_REC8_Myc_Rep1_input wt_REC8_Myc_Rep2_input wt_REC8_HA_Rep2_input wt_control_Rep1_input 32"

dat1=${1}
dat2=${2}
dat3=${3}
dat4=${4}
dat5=${5}
dat6=${6}
dat7=${7}
dat8=${8}
dat9=${9}
dat1label=${10}
dat2label=${11}
dat3label=${12}
dat4label=${13}
dat5label=${14}
dat6label=${15}
dat7label=${16}
dat8label=${17}
dat9label=${18}
threads=${19}

plotFingerprint -b ${dat1}_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam \
                   ${dat2}_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam \
                   ${dat3}_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam \
                   ${dat4}_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam \
                   ${dat5}_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam \
                   ${dat6}_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam \
                   ${dat7}_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam \
                   ${dat8}_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam \
                   ${dat9}_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam \
                -o REC8_control_ChIP_input_fingerprints.pdf \
                --extendReads \
                --labels ${dat1label} ${dat2label} ${dat3label} ${dat4label} ${dat5label} ${dat6label} ${dat7label} ${dat8label} ${dat9label} \
                --outRawCounts REC8_control_ChIP_input_fingerprints.tab \
                -p ${threads}
