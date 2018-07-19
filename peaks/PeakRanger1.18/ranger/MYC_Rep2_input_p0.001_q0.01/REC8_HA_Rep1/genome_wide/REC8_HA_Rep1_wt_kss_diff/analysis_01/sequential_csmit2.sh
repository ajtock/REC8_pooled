#!/bin/bash

csmit -m 20G -c 1 "Rscript peak_Profiles_commandArgs.R /projects/ajt200/BAM_masters/RNAseq_meiocyte_Walker_Feng_2018_NatGenet/WT/meiocyte/coverage/WT_RNAseq_meiocyte_Rep2_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_meiocyte_Rep2" & sleep 30;
csmit -m 20G -c 1 "Rscript peak_Profiles_commandArgs.R /projects/ajt200/BAM_masters/H3K4me3/replicates/coverage/log2ChIPinput/WT_H3K4me3_ChIP15_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H3K4me3_ChIP15" & sleep 30;
csmit -m 20G -c 1 "Rscript peak_Profiles_commandArgs.R /projects/ajt200/BAM_masters/H3K9me2/cmt3/coverage/log2ChIPinput/log2_cmt3_H3K9me2_ChIP_cmt3_H3K9me2_input_norm_allchrs_coverage_coord_tab.bed cmt3_H3K9me2" & sleep 30;
wait

