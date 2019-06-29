#!/bin/bash

csmit -m 20G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/RNAseq_meiocyte_Walker_Feng_2018_NatGenet/WT/meiocyte/coverage/WT_RNAseq_meiocyte_Rep2_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_meiocyte_Rep2" & sleep 10;
csmit -m 20G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/RNAseq_meiocyte_Walker_Feng_2018_NatGenet/WT/meiocyte/coverage/WT_RNAseq_meiocyte_Rep3_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_meiocyte_Rep3" & sleep 10;
csmit -m 20G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K4me1/coverage/log2ChIPinput/WT_H3K4me1_Rep1_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H3K4me1" & sleep 10;
csmit -m 20G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K4me2/coverage/log2ChIPinput/WT_H3K4me2_Rep1_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H3K4me2" & sleep 10;
csmit -m 20G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/H3K4me3/replicates/coverage/log2ChIPinput/WT_H3K4me3_ChIP12_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H3K4me3_ChIP12" & sleep 10;
csmit -m 20G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/H3K4me3/replicates/coverage/log2ChIPinput/WT_H3K4me3_ChIP14_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H3K4me3_ChIP14" & sleep 10;
csmit -m 20G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/H3K4me3/replicates/coverage/log2ChIPinput/WT_H3K4me3_ChIP15_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H3K4me3_ChIP15" & sleep 10;
csmit -m 20G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/H3K9me2/WT/coverage/log2ChIPinput/log2_WT_H3K9me2_ChIP_WT_H3K9me2_input_norm_allchrs_coverage_coord_tab.bed H3K9me2" & sleep 10;
csmit -m 20G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/H3K9me2/kss/coverage/log2ChIPinput/log2_kss_H3K9me2_ChIP_kss_H3K9me2_input_norm_allchrs_coverage_coord_tab.bed kss_H3K9me2" & sleep 10;
csmit -m 20G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/H3K9me2/cmt3/coverage/log2ChIPinput/log2_cmt3_H3K9me2_ChIP_cmt3_H3K9me2_input_norm_allchrs_coverage_coord_tab.bed cmt3_H3K9me2" & sleep 10;
csmit -m 20G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K27me1/coverage/log2ChIPinput/WT_H3K27me1_Rep1_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H3K27me1" & sleep 10;
csmit -m 20G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K27me3/coverage/log2ChIPinput/WT_H3K27me3_Rep1_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H3K27me3" & sleep 10;
csmit -m 20G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/H3K27me3_bud_UWMadison_2015/coverage/log2ChIPinput/log2_H3K27me3_ChIP_SRR1509478_input_norm_allchrs_coverage_coord_tab.bed H3K27me3_SRR1509478" & sleep 10;
csmit -m 20G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/H2A/coverage/log2ChIPinput/H2A_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H2A" & sleep 10;
csmit -m 20G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/H2A/coverage/log2ChIPinput/H2AW_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H2AW" & sleep 10;
csmit -m 20G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/H2A/coverage/log2ChIPinput/H2AX_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H2AX" & sleep 10;
csmit -m 20G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/PolV_Liu_Jacobsen_2018_NatPlants/coverage/log2ChIPinput/PolV_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed PolV" & sleep 10;

csmit -m 10G -c 1 "./peak_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM980986_WT_rep2_CG.wig.bed.gr.tab.bed mCG" & sleep 10;
csmit -m 10G -c 1 "./peak_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM980986_WT_rep2_CHG.wig.bed.gr.tab.bed mCHG" & sleep 10;
csmit -m 10G -c 1 "./peak_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM980986_WT_rep2_CHH.wig.bed.gr.tab.bed mCHH" & sleep 10;

csmit -m 10G -c 1 "./peak_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981060_suvh456_CG.wig.bed.gr.tab.bed kss_mCG" & sleep 10;
csmit -m 10G -c 1 "./peak_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981060_suvh456_CHH.wig.bed.gr.tab.bed kss_mCHH" & sleep 10;

csmit -m 10G -c 1 "./peak_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981003_cmt3_CG.wig.bed.gr.tab.bed cmt3_mCG" & sleep 10;
csmit -m 10G -c 1 "./peak_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981003_cmt3_CHG.wig.bed.gr.tab.bed cmt3_mCHG" & sleep 10;
csmit -m 10G -c 1 "./peak_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981003_cmt3_CHH.wig.bed.gr.tab.bed cmt3_mCHH" & sleep 10;
wait

