#!/bin/bash

#csmit -m 50G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/REC8_pooled/coverage/REC8_MYC_Rep1_ChIP_norm_allchrs_coverage_coord_tab.bed REC8_MYC_Rep1_ChIP" & sleep 5;
#csmit -m 50G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/REC8_pooled/coverage/REC8_HA_Rep1_ChIP_norm_allchrs_coverage_coord_tab.bed REC8_HA_Rep1_ChIP" & sleep 5;
#csmit -m 50G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/REC8_pooled/coverage/control_MYC_Rep1_ChIP_norm_allchrs_coverage_coord_tab.bed control_MYC_Rep1_ChIP" & sleep 5;
csmit -m 50G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/REC8_pooled/coverage/control_HA_Rep1_ChIP_norm_allchrs_coverage_coord_tab.bed control_HA_Rep1_ChIP" & sleep 5;
#csmit -m 50G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/REC8_pooled/coverage/REC8_MYC_Rep1_input_norm_allchrs_coverage_coord_tab.bed REC8_MYC_Rep1_input" & sleep 5;
wait

