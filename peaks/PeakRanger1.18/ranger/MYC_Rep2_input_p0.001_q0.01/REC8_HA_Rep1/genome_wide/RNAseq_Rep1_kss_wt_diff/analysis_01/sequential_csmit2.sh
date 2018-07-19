#!/bin/bash

csmit -m 20G -c 1 "Rscript peak_Profiles_commandArgs.R /home/meiosis/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/WT/coverage/WT_RNAseq_Chris_Rep2_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_Chris_Rep2" & sleep 30;
csmit -m 20G -c 1 "Rscript peak_Profiles_commandArgs.R /projects/ajt200/BAM_masters/RNAseq_meiocyte_Walker_Feng_2018_NatGenet/WT/meiocyte/coverage/WT_RNAseq_meiocyte_Rep1_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_meiocyte_Rep1" & sleep 30;
csmit -m 20G -c 1 "Rscript peak_Profiles_commandArgs.R /projects/ajt200/BAM_masters/RNAseq_meiocyte_Walker_Feng_2018_NatGenet/WT/meiocyte/coverage/WT_RNAseq_meiocyte_Rep2_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_meiocyte_Rep2" & sleep 30;

csmit -m 10G -c 1 "Rscript peak_Profiles_DNAmeth_commandArgs.R /home/meiosis/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM980986_WT_rep2_CG.wig.bed.gr.tab.bed CGmeth" & sleep 30;
wait

