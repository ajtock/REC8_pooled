#!/bin/bash

csmit -m 20G -c 1 "Rscript geneHexile_Profiles_commandArgs.R /projects/ajt200/BAM_masters/RNAseq/WT/coverage/WT_RNAseq_Rep1_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_Kyuha_Rep1" & sleep 5;
csmit -m 20G -c 1 "Rscript geneHexile_Profiles_commandArgs.R /projects/ajt200/BAM_masters/RNAseq/WT/coverage/WT_RNAseq_Rep2_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_Kyuha_Rep2" & sleep 5;
csmit -m 20G -c 1 "Rscript geneHexile_Profiles_commandArgs.R /projects/ajt200/BAM_masters/RNAseq/WT/coverage/WT_RNAseq_Rep3_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_Kyuha_Rep3" & sleep 5;
csmit -m 20G -c 1 "Rscript geneHexile_Profiles_commandArgs.R /projects/ajt200/BAM_masters/RNAseq/WT/bowtie2/coverage/WT_RNAseq_ATCACG_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_Kyuha_Rep1_bowtie2" & sleep 5;
csmit -m 20G -c 1 "Rscript geneHexile_Profiles_commandArgs.R /projects/ajt200/BAM_masters/RNAseq/WT/bowtie2/coverage/WT_RNAseq_CGATGT_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_Kyuha_Rep2_bowtie2" & sleep 5;
csmit -m 20G -c 1 "Rscript geneHexile_Profiles_commandArgs.R /projects/ajt200/BAM_masters/RNAseq/WT/bowtie2/coverage/WT_RNAseq_TTAGGC_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_Kyuha_Rep3_bowtie2" & sleep 5;
wait

