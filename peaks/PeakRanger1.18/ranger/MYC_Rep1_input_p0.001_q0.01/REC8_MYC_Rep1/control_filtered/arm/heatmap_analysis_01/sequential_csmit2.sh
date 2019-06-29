#!/bin/bash

csmit -m 15G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/RNAseq_meiocyte_Walker_Feng_2018_NatGenet/WT/meiocyte/coverage/WT_RNAseq_meiocyte_Rep3_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_meiocyte_Rep3" & sleep 5;

csmit -m 5G -c 1 "./peak_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981060_suvh456_CG.wig.bed.gr.tab.bed kss_mCG" & sleep 5;
csmit -m 5G -c 1 "./peak_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981060_suvh456_CHG.wig.bed.gr.tab.bed kss_mCHG" & sleep 5;
wait

