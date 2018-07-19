#!/bin/bash

csmit -m 50G -c 1 "Rscript gene_Profiles_commandArgs.R /home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K4me1/coverage/log2ChIPinput/WT_H3K4me1_Rep1_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H3K4me1" & sleep 5;
csmit -m 50G -c 1 "Rscript gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/PolIV_Law_Jacobsen_2013_Nature/coverage/log2ChIPinput/PolIV_Rep2_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed PolIV_Rep2" & sleep 5;
wait

