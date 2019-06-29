#!/bin/bash

csmit -m 15G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/log2_REC8_MYC_Rep1_ChIP_REC8_MYC_Rep1_input_norm_allchrs_coverage_coord_tab.bed REC8_MYC_Rep1" & sleep 5;
csmit -m 15G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/log2_REC8_HA_Rep2_ChIP_REC8_MYC_Rep1_input_norm_allchrs_coverage_coord_tab.bed REC8_HA_Rep2" & sleep 5;

csmit -m 15G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/160902_Sasha_ChIP_MSH4_Rep1/coverage/log2ChIPinput/MSH4_Rep1_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed MSH4_Rep1" & sleep 5;

csmit -m 15G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/RNAseq_meiocyte_Walker_Feng_2018_NatGenet/WT/meiocyte/coverage/WT_RNAseq_meiocyte_Rep2_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_meiocyte_Rep2" & sleep 5;
csmit -m 15G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/RNAseq_meiocyte_Walker_Feng_2018_NatGenet/WT/meiocyte/coverage/WT_RNAseq_meiocyte_Rep3_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_meiocyte_Rep3" & sleep 5;

csmit -m 15G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K4me1/coverage/log2ChIPinput/WT_H3K4me1_Rep1_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H3K4me1" & sleep 5;
csmit -m 15G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K4me2/coverage/log2ChIPinput/WT_H3K4me2_Rep1_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H3K4me2" & sleep 5;
csmit -m 15G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/H3K4me3/replicates/coverage/log2ChIPinput/WT_H3K4me3_ChIP12_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H3K4me3_ChIP12" & sleep 5;
csmit -m 15G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/H3K4me3/replicates/coverage/log2ChIPinput/WT_H3K4me3_ChIP14_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H3K4me3_ChIP14" & sleep 5;

csmit -m 15G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/H2A/coverage/log2ChIPinput/H2AX_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H2AX" & sleep 5;

csmit -m 5G -c 1 "./peak_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981060_suvh456_CG.wig.bed.gr.tab.bed kss_mCG" & sleep 5;
csmit -m 5G -c 1 "./peak_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981060_suvh456_CHG.wig.bed.gr.tab.bed kss_mCHG" & sleep 5;

wait

