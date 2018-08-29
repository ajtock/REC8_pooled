#!/bin/bash

csmit -m 20G -c 1 "Rscript TE_DNAfams_Profiles_commandArgs.R /projects/ajt200/BAM_masters/SPO11_ChIP/WT/coverage/REC8_MYC_Rep2_input/log2ChIPinput/log2_WT_SPO11_ChIP4_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab.bed SPO11_ChIP4" & sleep 10;
csmit -m 20G -c 1 "Rscript TE_DNAfams_Profiles_commandArgs.R /projects/ajt200/BAM_masters/SPO11_ChIP/WT/coverage/REC8_MYC_Rep2_input/log2ChIPinput/log2_WT_SPO11_ChIP13_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab.bed SPO11_ChIP13" & sleep 10;

