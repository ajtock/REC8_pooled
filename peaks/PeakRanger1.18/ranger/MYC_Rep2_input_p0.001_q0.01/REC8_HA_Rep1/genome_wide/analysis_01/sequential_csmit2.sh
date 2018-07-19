#!/bin/bash

csmit -m 30G -c 1 "Rscript peak_Profiles_commandArgs.R /projects/ajt200/BAM_masters/RNAseq_meiocyte_Walker_Feng_2018_NatGenet/WT/meiocyte/coverage/WT_RNAseq_meiocyte_Rep1_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_meiocyte_Rep1" & sleep 10;
csmit -m 30G -c 1 "Rscript peak_Profiles_commandArgs.R /projects/ajt200/BAM_masters/RNAseq_meiocyte_Walker_Feng_2018_NatGenet/WT/meiocyte/coverage/WT_RNAseq_meiocyte_Rep2_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_meiocyte_Rep2" & sleep 10;
csmit -m 30G -c 1 "Rscript peak_Profiles_commandArgs.R /projects/ajt200/BAM_masters/RNAseq_meiocyte_Walker_Feng_2018_NatGenet/WT/meiocyte/coverage/WT_RNAseq_meiocyte_Rep3_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_meiocyte_Rep3" & sleep 10;
csmit -m 30G -c 1 "Rscript peak_Profiles_commandArgs.R /projects/ajt200/BAM_masters/H2A/coverage/log2ChIPinput/H2A_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H2A" & sleep 10;
csmit -m 30G -c 1 "Rscript peak_Profiles_commandArgs.R /projects/ajt200/BAM_masters/H2A/coverage/log2ChIPinput/H2AW_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H2AW" & sleep 10;
csmit -m 30G -c 1 "Rscript peak_Profiles_commandArgs.R /projects/ajt200/BAM_masters/H2A/coverage/log2ChIPinput/H2AX_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H2AX" & sleep 10;
csmit -m 30G -c 1 "Rscript peak_Profiles_commandArgs.R /projects/ajt200/BAM_masters/H2A/coverage/log2ChIPinput/H2AZ_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H2AZ" & sleep 10;
csmit -m 30G -c 1 "Rscript peak_Profiles_commandArgs.R /projects/ajt200/BAM_masters/PolIV_Law_Jacobsen_2013_Nature/coverage/log2ChIPinput/PolIV_Rep2_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed PolIV_Rep2" & sleep 10;
csmit -m 30G -c 1 "Rscript peak_Profiles_commandArgs.R /projects/ajt200/BAM_masters/PolV_Liu_Jacobsen_2018_NatPlants/coverage/log2ChIPinput/PolV_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed PolV" & sleep 10;

csmit -m 10G -c 1 "Rscript peak_Profiles_DNAmeth_commandArgs.R /home/meiosis/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM980986_WT_rep2_CG.wig.bed.gr.tab.bed CGmeth" & sleep 10;
csmit -m 10G -c 1 "Rscript peak_Profiles_DNAmeth_commandArgs.R /home/meiosis/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM980986_WT_rep2_CHG.wig.bed.gr.tab.bed CHGmeth" & sleep 10;
csmit -m 10G -c 1 "Rscript peak_Profiles_DNAmeth_commandArgs.R /home/meiosis/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM980986_WT_rep2_CHH.wig.bed.gr.tab.bed CHHmeth" & sleep 10;
wait

