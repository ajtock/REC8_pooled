#!/bin/bash

csmit -m 20G -c 1 "Rscript peak_Profiles_commandArgs.R /home/meiosis/ajt200/analysis/180622_Chris_lambing_ChIP_REC8_HA_Col_kss/kss/coverage/common_input_MYC_Rep2/log2ChIPinput/log2_kss_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab.bed kss_REC8_HA_Rep1" & sleep 30;
csmit -m 20G -c 1 "Rscript peak_Profiles_commandArgs.R /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/log2wtSPO11oligoRPI1NakedDNA_norm_allchrs_coverage_coord_tab.bed SPO11_1_oligos_RPI1" & sleep 30;
csmit -m 20G -c 1 "Rscript peak_Profiles_commandArgs.R /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/log2wtSPO11oligoRPI3NakedDNA_norm_allchrs_coverage_coord_tab.bed SPO11_1_oligos_RPI3" & sleep 30;
wait

