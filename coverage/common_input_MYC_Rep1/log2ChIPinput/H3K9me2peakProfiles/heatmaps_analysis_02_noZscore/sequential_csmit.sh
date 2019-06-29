#!/bin/bash

csmit -m 30G -c 1 "./H3K9me2peak_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/noZscore/log2_REC8_HA_Rep2_ChIP_REC8_MYC_Rep1_input_norm_allchrs_coverage_coord_tab_noZscore.bed REC8_HA_Rep2" & sleep 10;
csmit -m 30G -c 1 "./H3K9me2peak_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/180622_Chris_lambing_ChIP_REC8_HA_Col_kss/kss/coverage/log2ChIPinput/replicate_specific_input/noZscore/log2_kss_REC8_HA_Rep1_ChIP_kss_REC8_HA_Rep1_input_norm_allchrs_coverage_coord_tab_noZscore.bed kss_REC8_HA_Rep1" & sleep 10;
csmit -m 30G -c 1 "./H3K9me2peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/noZscore/log2wtSPO11oligoRPI1NakedDNA_norm_allchrs_coverage_coord_tab_noZscore.bed SPO11_1_oligos_RPI1" & sleep 10;
csmit -m 30G -c 1 "./H3K9me2peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/noZscore/log2wtSPO11oligoRPI3NakedDNA_norm_allchrs_coverage_coord_tab_noZscore.bed SPO11_1_oligos_RPI3" & sleep 10;
csmit -m 30G -c 1 "./H3K9me2peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/noZscore/log2wtSPO11oligoRPI8NakedDNA_norm_allchrs_coverage_coord_tab_noZscore.bed SPO11_1_oligos_RPI8" & sleep 10;
csmit -m 30G -c 1 "./H3K9me2peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/SPO11-oligo/suvh456/coverage/log2ChIPinput/noZscore/log2suvh456SPO11oligoRPI34NakedDNA_norm_allchrs_coverage_coord_tab_noZscore.bed kss_SPO11_1_oligos_RPI34" & sleep 10;
csmit -m 30G -c 1 "./H3K9me2peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/SPO11-oligo/suvh456/coverage/log2ChIPinput/noZscore/log2suvh456SPO11oligoRPI35NakedDNA_norm_allchrs_coverage_coord_tab_noZscore.bed kss_SPO11_1_oligos_RPI35" & sleep 10;
csmit -m 30G -c 1 "./H3K9me2peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/H3K9me2/WT/coverage/log2ChIPinput/noZscore/log2_WT_H3K9me2_ChIP_WT_H3K9me2_input_norm_allchrs_coverage_coord_tab_noZscore.bed H3K9me2" & sleep 10;
csmit -m 30G -c 1 "./H3K9me2peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/H3K9me2/kss/coverage/log2ChIPinput/noZscore/log2_kss_H3K9me2_ChIP_kss_H3K9me2_input_norm_allchrs_coverage_coord_tab_noZscore.bed kss_H3K9me2" & sleep 10;

wait

