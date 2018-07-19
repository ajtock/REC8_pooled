#!/bin/bash

csmit -m 15G -c 1 "Rscript gene_Profiles_commandArgs.R /home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab.bed REC8_HA_Rep1" & sleep 30;
csmit -m 15G -c 1 "Rscript gene_Profiles_commandArgs.R /home/meiosis/ajt200/analysis/180622_Chris_lambing_ChIP_REC8_HA_Col_kss/WT/coverage/common_input_MYC_Rep2/log2ChIPinput/log2_REC8_HA_Rep2_ChIP_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab.bed REC8_HA_Rep2" & sleep 30;
csmit -m 15G -c 1 "Rscript gene_Profiles_commandArgs.R /home/meiosis/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/log2_REC8_MYC_Rep1_ChIP_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab.bed REC8_MYC_Rep1" & sleep 30;
csmit -m 15G -c 1 "Rscript gene_Profiles_commandArgs.R /home/meiosis/ajt200/analysis/180622_Chris_lambing_ChIP_REC8_HA_Col_kss/kss/coverage/common_input_MYC_Rep2/log2ChIPinput/log2_kss_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab.bed kss_REC8_HA_Rep1" & sleep 30;

csmit -m 15G -c 1 "Rscript gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/nucleosomes/WT/coverage/nakedDNA_untrimmed_input/log2ChIPinput/log2wtNucNakedDNAuntrimmed_norm_allchrs_coverage_coord_tab.bed MNase" & sleep 30;

csmit -m 15G -c 1 "Rscript gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/log2wtSPO11oligoRPI1NakedDNA_norm_allchrs_coverage_coord_tab.bed SPO11_1_oligos_RPI1" & sleep 30;
csmit -m 15G -c 1 "Rscript gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/log2wtSPO11oligoRPI3NakedDNA_norm_allchrs_coverage_coord_tab.bed SPO11_1_oligos_RPI3" & sleep 30;
csmit -m 15G -c 1 "Rscript gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/log2wtSPO11oligoRPI8NakedDNA_norm_allchrs_coverage_coord_tab.bed SPO11_1_oligos_RPI8" & sleep 30;
csmit -m 15G -c 1 "Rscript gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/SPO11-oligo/suvh456/coverage/log2ChIPinput/log2suvh456SPO11oligoRPI34NakedDNA_norm_allchrs_coverage_coord_tab.bed kss_SPO11_1_oligos_RPI34" & sleep 30;
csmit -m 15G -c 1 "Rscript gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/SPO11-oligo/suvh456/coverage/log2ChIPinput/log2suvh456SPO11oligoRPI35NakedDNA_norm_allchrs_coverage_coord_tab.bed kss_SPO11_1_oligos_RPI35" & sleep 30;

csmit -m 15G -c 1 "Rscript gene_Profiles_commandArgs.R /home/meiosis/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/WT/coverage/WT_RNAseq_Chris_Rep1_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_Chris_Rep1" & sleep 30;
csmit -m 15G -c 1 "Rscript gene_Profiles_commandArgs.R /home/meiosis/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/WT/coverage/WT_RNAseq_Chris_Rep2_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_Chris_Rep2" & sleep 30;
csmit -m 15G -c 1 "Rscript gene_Profiles_commandArgs.R /home/meiosis/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/kss/coverage/kss_RNAseq_Chris_Rep1_norm_allchrs_coverage_coord_tab.bed kss_RNAseq_Chris_Rep1" & sleep 30;
csmit -m 15G -c 1 "Rscript gene_Profiles_commandArgs.R /home/meiosis/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/kss/coverage/kss_RNAseq_Chris_Rep2_norm_allchrs_coverage_coord_tab.bed kss_RNAseq_Chris_Rep2" & sleep 30;

csmit -m 15G -c 1 "Rscript gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/H3K9me2/WT/coverage/log2ChIPinput/log2_WT_H3K9me2_ChIP_WT_H3K9me2_input_norm_allchrs_coverage_coord_tab.bed H3K9me2" & sleep 30;
csmit -m 15G -c 1 "Rscript gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/H3K9me2/kss/coverage/log2ChIPinput/log2_kss_H3K9me2_ChIP_kss_H3K9me2_input_norm_allchrs_coverage_coord_tab.bed kss_H3K9me2" & sleep 30;
csmit -m 15G -c 1 "Rscript gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/H3K9me2/cmt3/coverage/log2ChIPinput/log2_cmt3_H3K9me2_ChIP_cmt3_H3K9me2_input_norm_allchrs_coverage_coord_tab.bed cmt3_H3K9me2" & sleep 30;

csmit -m 10G -c 1 "Rscript gene_Profiles_DNAmeth_commandArgs.R /home/meiosis/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM980986_WT_rep2_CG.wig.bed.gr.tab.bed CGmeth" & sleep 30;
csmit -m 10G -c 1 "Rscript gene_Profiles_DNAmeth_commandArgs.R /home/meiosis/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM980986_WT_rep2_CHG.wig.bed.gr.tab.bed CHGmeth" & sleep 30;
csmit -m 10G -c 1 "Rscript gene_Profiles_DNAmeth_commandArgs.R /home/meiosis/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM980986_WT_rep2_CHH.wig.bed.gr.tab.bed CHHmeth" & sleep 30;

csmit -m 10G -c 1 "Rscript gene_Profiles_DNAmeth_commandArgs.R /home/meiosis/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981060_suvh456_CG.bed.gr.tab.bed kss_CGmeth" & sleep 30;
csmit -m 10G -c 1 "Rscript gene_Profiles_DNAmeth_commandArgs.R /home/meiosis/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981060_suvh456_CHG.bed.gr.tab.bed kss_CHGmeth" & sleep 30;
csmit -m 10G -c 1 "Rscript gene_Profiles_DNAmeth_commandArgs.R /home/meiosis/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981060_suvh456_CHH.bed.gr.tab.bed kss_CHHmeth" & sleep 30;

csmit -m 10G -c 1 "Rscript gene_Profiles_DNAmeth_commandArgs.R /home/meiosis/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981003_cmt3_CG.wig.bed.gr.tab.bed cmt3_CGmeth" & sleep 30;
csmit -m 10G -c 1 "Rscript gene_Profiles_DNAmeth_commandArgs.R /home/meiosis/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981003_cmt3_CHG.wig.bed.gr.tab.bed cmt3_CHGmeth" & sleep 30;
csmit -m 10G -c 1 "Rscript gene_Profiles_DNAmeth_commandArgs.R /home/meiosis/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981003_cmt3_CHH.wig.bed.gr.tab.bed cmt3_CHHmeth" & sleep 30;

csmit -m 15G -c 1 "Rscript gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/PolIV_Law_Jacobsen_2013_Nature/coverage/log2ChIPinput/PolIV_Rep2_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed PolIV_Rep2" & sleep 30;
csmit -m 15G -c 1 "Rscript gene_Profiles_commandArgs.R /projects/ajt200/BAM_masters/PolV_Liu_Jacobsen_2018_NatPlants/coverage/log2ChIPinput/PolV_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed PolV" & sleep 30;
wait
