#!/bin/bash

csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/log2_REC8_MYC_Rep1_ChIP_REC8_MYC_Rep1_input_norm_allchrs_coverage_coord_tab.bed REC8_MYC_Rep1" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/log2_REC8_HA_Rep2_ChIP_REC8_MYC_Rep1_input_norm_allchrs_coverage_coord_tab.bed REC8_HA_Rep2" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep1_input_norm_allchrs_coverage_coord_tab.bed REC8_HA_Rep1" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/180622_Chris_lambing_ChIP_REC8_HA_Col_kss/kss/coverage/log2ChIPinput/replicate_specific_input/log2_kss_REC8_HA_Rep1_ChIP_kss_REC8_HA_Rep1_input_norm_allchrs_coverage_coord_tab.bed kss_REC8_HA_Rep1" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/180830_Chris_lambing_ChIP_rec8HA_kss_rep2/coverage/log2ChIPinput/kss_REC8_HA_Rep1_input/log2_kss_REC8_HA_Rep2_ChIP_kss_REC8_HA_Rep1_input_norm_allchrs_coverage_coord_tab.bed kss_REC8_HA_Rep2" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/20190320_Chris_Lambing_ChIP_control/fastq_pooled/coverage/log2ChIPinput/REC8_MYC_Rep1_input/log2_control_MYC_Rep1_ChIP_REC8_MYC_Rep1_input_norm_allchrs_coverage_coord_tab.bed control_MYC_Rep1" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/20190320_Chris_Lambing_ChIP_control/fastq_pooled/coverage/log2ChIPinput/REC8_MYC_Rep1_input/log2_control_HA_Rep1_ChIP_REC8_MYC_Rep1_input_norm_allchrs_coverage_coord_tab.bed control_HA_Rep1" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/nucleosomes/WT/coverage/nakedDNA_untrimmed_input/log2ChIPinput/log2wtNucNakedDNAuntrimmed_norm_allchrs_coverage_coord_tab.bed MNase" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/160902_Sasha_ChIP_MSH4_Rep1/coverage/log2ChIPinput/MSH4_Rep1_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed MSH4_Rep1" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/SPO11_ChIP/WT/coverage/REC8_MYC_Rep1_input/log2ChIPinput/log2_WT_SPO11_ChIP4_REC8_MYC_Rep1_input_norm_allchrs_coverage_coord_tab.bed SPO11_1_ChIP4" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/SPO11_ChIP/WT/coverage/REC8_MYC_Rep1_input/log2ChIPinput/log2_WT_SPO11_ChIP13_REC8_MYC_Rep1_input_norm_allchrs_coverage_coord_tab.bed SPO11_1_ChIP13" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/log2wtSPO11oligoRPI1NakedDNA_norm_allchrs_coverage_coord_tab.bed SPO11_1_oligos_RPI1" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/log2wtSPO11oligoRPI3NakedDNA_norm_allchrs_coverage_coord_tab.bed SPO11_1_oligos_RPI3" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/log2wtSPO11oligoRPI8NakedDNA_norm_allchrs_coverage_coord_tab.bed SPO11_1_oligos_RPI8" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/SPO11-oligo/suvh456/coverage/log2ChIPinput/log2suvh456SPO11oligoRPI34NakedDNA_norm_allchrs_coverage_coord_tab.bed kss_SPO11_1_oligos_RPI34" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/SPO11-oligo/suvh456/coverage/log2ChIPinput/log2suvh456SPO11oligoRPI35NakedDNA_norm_allchrs_coverage_coord_tab.bed kss_SPO11_1_oligos_RPI35" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/WT/coverage/WT_RNAseq_Chris_Rep1_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_Chris_Rep1" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/WT/coverage/WT_RNAseq_Chris_Rep2_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_Chris_Rep2" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/kss/coverage/kss_RNAseq_Chris_Rep1_norm_allchrs_coverage_coord_tab.bed kss_RNAseq_Chris_Rep1" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/kss/coverage/kss_RNAseq_Chris_Rep2_norm_allchrs_coverage_coord_tab.bed kss_RNAseq_Chris_Rep2" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/RNAseq/WT/coverage/WT_RNAseq_Rep1_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_Kyuha_Rep1" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/RNAseq/WT/coverage/WT_RNAseq_Rep2_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_Kyuha_Rep2" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/RNAseq/WT/coverage/WT_RNAseq_Rep3_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_Kyuha_Rep3" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/RNAseq_meiocyte_Walker_Feng_2018_NatGenet/WT/meiocyte/coverage/WT_RNAseq_meiocyte_Rep1_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_meiocyte_Rep1" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/RNAseq_meiocyte_Walker_Feng_2018_NatGenet/WT/meiocyte/coverage/WT_RNAseq_meiocyte_Rep2_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_meiocyte_Rep2" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/RNAseq_meiocyte_Walker_Feng_2018_NatGenet/WT/meiocyte/coverage/WT_RNAseq_meiocyte_Rep3_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_meiocyte_Rep3" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K4me1/coverage/log2ChIPinput/WT_H3K4me1_Rep1_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H3K4me1" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K4me2/coverage/log2ChIPinput/WT_H3K4me2_Rep1_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H3K4me2" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/H3K4me3/replicates/coverage/log2ChIPinput/WT_H3K4me3_ChIP12_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H3K4me3_ChIP12" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/H3K4me3/replicates/coverage/log2ChIPinput/WT_H3K4me3_ChIP14_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H3K4me3_ChIP14" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/H3K4me3/replicates/coverage/log2ChIPinput/WT_H3K4me3_ChIP15_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H3K4me3_ChIP15" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/H3K9me2/WT/coverage/log2ChIPinput/log2_WT_H3K9me2_ChIP_WT_H3K9me2_input_norm_allchrs_coverage_coord_tab.bed H3K9me2" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/H3K9me2/kss/coverage/log2ChIPinput/log2_kss_H3K9me2_ChIP_kss_H3K9me2_input_norm_allchrs_coverage_coord_tab.bed kss_H3K9me2" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/H3K9me2/cmt3/coverage/log2ChIPinput/log2_cmt3_H3K9me2_ChIP_cmt3_H3K9me2_input_norm_allchrs_coverage_coord_tab.bed cmt3_H3K9me2" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K27me1/coverage/log2ChIPinput/WT_H3K27me1_Rep1_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H3K27me1" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K27me3/coverage/log2ChIPinput/WT_H3K27me3_Rep1_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H3K27me3" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/H3K27me3_bud_UWMadison_2015/coverage/log2ChIPinput/log2_H3K27me3_ChIP_SRR1509478_input_norm_allchrs_coverage_coord_tab.bed H3K27me3_SRR1509478" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/H2A/coverage/log2ChIPinput/H2A_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H2A" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/H2A/coverage/log2ChIPinput/H2AW_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H2AW" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/H2A/coverage/log2ChIPinput/H2AX_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H2AX" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/H2A/coverage/log2ChIPinput/H2AZ_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed H2AZ" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/PolIV_Law_Jacobsen_2013_Nature/coverage/log2ChIPinput/PolIV_Rep2_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed PolIV_Rep2" & sleep 10;
csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/PolV_Liu_Jacobsen_2018_NatPlants/coverage/log2ChIPinput/PolV_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed PolV" & sleep 10;

csmit -m 10G -c 1 "./gene_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM980986_WT_rep2_CG.wig.bed.gr.tab.bed mCG" & sleep 10;
csmit -m 10G -c 1 "./gene_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM980986_WT_rep2_CHG.wig.bed.gr.tab.bed mCHG" & sleep 10;
csmit -m 10G -c 1 "./gene_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM980986_WT_rep2_CHH.wig.bed.gr.tab.bed mCHH" & sleep 10;

csmit -m 10G -c 1 "./gene_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981060_suvh456_CG.wig.bed.gr.tab.bed kss_mCG" & sleep 10;
csmit -m 10G -c 1 "./gene_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981060_suvh456_CHG.wig.bed.gr.tab.bed kss_mCHG" & sleep 10;
csmit -m 10G -c 1 "./gene_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981060_suvh456_CHH.wig.bed.gr.tab.bed kss_mCHH" & sleep 10;

csmit -m 10G -c 1 "./gene_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981003_cmt3_CG.wig.bed.gr.tab.bed cmt3_mCG" & sleep 10;
csmit -m 10G -c 1 "./gene_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981003_cmt3_CHG.wig.bed.gr.tab.bed cmt3_mCHG" & sleep 10;
csmit -m 10G -c 1 "./gene_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981003_cmt3_CHH.wig.bed.gr.tab.bed cmt3_mCHH" & sleep 10;
wait

