#!/bin/bash

csmit -m 50G -c 1 "Rscript crossover_Profiles_commandArgs.R /projects/ajt200/BAM_masters/SPO11-oligo/arp6/coverage/log2ChIPinput/log2arp6SPO11oligoRPI32NakedDNA_norm_allchrs_coverage_coord_tab.bed arp6_SPO11_1_oligos_RPI32" & sleep 5;
csmit -m 50G -c 1 "Rscript crossover_Profiles_commandArgs.R /projects/ajt200/BAM_masters/SPO11-oligo/arp6/coverage/log2ChIPinput/log2arp6SPO11oligoRPI33NakedDNA_norm_allchrs_coverage_coord_tab.bed arp6_SPO11_1_oligos_RPI33" & sleep 5;
csmit -m 50G -c 1 "Rscript crossover_Profiles_commandArgs.R /projects/ajt200/BAM_masters/nucleosomes/arp6/coverage/nakedDNA_untrimmed_input/log2ChIPinput/log2arp6NucNakedDNAuntrimmed_norm_allchrs_coverage_coord_tab.bed arp6_MNase" & sleep 5;
csmit -m 50G -c 1 "Rscript crossover_Profiles_commandArgs.R /projects/ajt200/BAM_masters/H3K4me3/replicates/coverage/log2ChIPinput/arp6_H3K4me3_ChIP5_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed arp6_H3K4me3_ChIP5" & sleep 5;
csmit -m 50G -c 1 "Rscript crossover_Profiles_commandArgs.R /projects/ajt200/BAM_masters/H3K4me3/replicates/coverage/log2ChIPinput/arp6_H3K4me3_ChIP7_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed arp6_H3K4me3_ChIP7" & sleep 5;
csmit -m 50G -c 1 "Rscript crossover_Profiles_commandArgs.R /projects/ajt200/BAM_masters/H3K4me3/replicates/coverage/log2ChIPinput/arp6_H3K4me3_ChIP16_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed arp6_H3K4me3_ChIP16" & sleep 5;
csmit -m 50G -c 1 "Rscript crossover_Profiles_commandArgs.R /projects/ajt200/BAM_masters/RNAseq/arp6/coverage/arp6_RNAseq_Rep1_norm_allchrs_coverage_coord_tab.bed arp6_RNAseq_Kyuha_Rep1" & sleep 5;
csmit -m 50G -c 1 "Rscript crossover_Profiles_commandArgs.R /projects/ajt200/BAM_masters/RNAseq/arp6/coverage/arp6_RNAseq_Rep2_norm_allchrs_coverage_coord_tab.bed arp6_RNAseq_Kyuha_Rep2" & sleep 5;
csmit -m 50G -c 1 "Rscript crossover_Profiles_commandArgs.R /projects/ajt200/BAM_masters/RNAseq/arp6/coverage/arp6_RNAseq_Rep3_norm_allchrs_coverage_coord_tab.bed arp6_RNAseq_Kyuha_Rep3" & sleep 5;

