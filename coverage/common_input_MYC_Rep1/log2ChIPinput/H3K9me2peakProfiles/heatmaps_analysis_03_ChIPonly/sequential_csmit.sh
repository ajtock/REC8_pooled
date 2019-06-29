#!/bin/bash

csmit -m 30G -c 1 "./H3K9me2peak_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/180622_Chris_lambing_ChIP_REC8_HA_Col_kss/WT/coverage/REC8_HA_Rep2_ChIP_norm_allchrs_coverage_coord_tab.bed REC8_HA_Rep2" & sleep 10;
csmit -m 30G -c 1 "./H3K9me2peak_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/180622_Chris_lambing_ChIP_REC8_HA_Col_kss/kss/coverage/kss_REC8_HA_Rep1_ChIP_norm_allchrs_coverage_coord_tab.bed kss_REC8_HA_Rep1" & sleep 10;
csmit -m 30G -c 1 "./H3K9me2peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/WT_SPO11oligo_RPI1_norm_allchrs_coverage_coord_tab.bed SPO11_1_oligos_RPI1" & sleep 10;
csmit -m 30G -c 1 "./H3K9me2peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/WT_SPO11oligo_RPI3_norm_allchrs_coverage_coord_tab.bed SPO11_1_oligos_RPI3" & sleep 10;
csmit -m 30G -c 1 "./H3K9me2peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/WT_SPO11oligo_RPI8_norm_allchrs_coverage_coord_tab.bed SPO11_1_oligos_RPI8" & sleep 10;
csmit -m 30G -c 1 "./H3K9me2peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/SPO11-oligo/suvh456/coverage/suvh456_SPO11oligo_RPI34_norm_allchrs_coverage_coord_tab.bed kss_SPO11_1_oligos_RPI34" & sleep 10;
csmit -m 30G -c 1 "./H3K9me2peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/SPO11-oligo/suvh456/coverage/suvh456_SPO11oligo_RPI35_norm_allchrs_coverage_coord_tab.bed kss_SPO11_1_oligos_RPI35" & sleep 10;

csmit -m 30G -c 1 "./H3K9me2peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/H3K9me2/WT/coverage/WT_H3K9me2_ChIP_norm_allchrs_coverage_coord_tab.bed H3K9me2" & sleep 10;
csmit -m 30G -c 1 "./H3K9me2peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/H3K9me2/kss/coverage/kss_H3K9me2_ChIP_norm_allchrs_coverage_coord_tab.bed kss_H3K9me2" & sleep 10;

wait

