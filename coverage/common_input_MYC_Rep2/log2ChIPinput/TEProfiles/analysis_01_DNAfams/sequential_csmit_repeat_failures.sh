#!/bin/bash

csmit -m 25G -c 1 "Rscript TE_DNAfams_Profiles_commandArgs.R /projects/ajt200/BAM_masters/RNAseq/WT/bowtie2/coverage/WT_RNAseq_CGATGT_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_Kyuha_Rep2_bowtie2" & sleep 5;
csmit -m 25G -c 1 "Rscript TE_DNAfams_Profiles_commandArgs.R /projects/ajt200/BAM_masters/RNAseq_meiocyte_Walker_Feng_2018_NatGenet/WT/meiocyte/coverage/WT_RNAseq_meiocyte_Rep2_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_meiocyte_Rep2" & sleep 5;

