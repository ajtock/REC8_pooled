#!/bin/bash

csmit -m 30G -c 1 "Rscript peak_Profiles_commandArgs.R /projects/ajt200/BAM_masters/H2A/coverage/log2ChIPinput/H2A_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H2A" & sleep 10;
csmit -m 30G -c 1 "Rscript peak_Profiles_commandArgs.R /projects/ajt200/BAM_masters/H2A/coverage/log2ChIPinput/H2AX_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H2AX" & sleep 10;
wait

