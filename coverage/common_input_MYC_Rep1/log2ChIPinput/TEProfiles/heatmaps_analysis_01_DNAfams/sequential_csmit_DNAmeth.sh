#!/bin/bash

csmit -m 10G -c 1 "./TE_DNAfams_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM980986_WT_rep2_CG.wig.bed.gr.tab.bed mCG" & sleep 10;
csmit -m 10G -c 1 "./TE_DNAfams_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM980986_WT_rep2_CHG.wig.bed.gr.tab.bed mCHG" & sleep 10;
csmit -m 10G -c 1 "./TE_DNAfams_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM980986_WT_rep2_CHH.wig.bed.gr.tab.bed mCHH" & sleep 10;

csmit -m 10G -c 1 "./TE_DNAfams_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981060_suvh456_CG.wig.bed.gr.tab.bed kss_mCG" & sleep 10;
csmit -m 10G -c 1 "./TE_DNAfams_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981060_suvh456_CHG.wig.bed.gr.tab.bed kss_mCHG" & sleep 10;
csmit -m 10G -c 1 "./TE_DNAfams_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981060_suvh456_CHH.wig.bed.gr.tab.bed kss_mCHH" & sleep 10;

csmit -m 10G -c 1 "./TE_DNAfams_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981003_cmt3_CG.wig.bed.gr.tab.bed cmt3_mCG" & sleep 10;
csmit -m 10G -c 1 "./TE_DNAfams_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981003_cmt3_CHG.wig.bed.gr.tab.bed cmt3_mCHG" & sleep 10;
csmit -m 10G -c 1 "./TE_DNAfams_Profiles_DNAmeth_commandArgs.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981003_cmt3_CHH.wig.bed.gr.tab.bed cmt3_mCHH" & sleep 10;
wait

