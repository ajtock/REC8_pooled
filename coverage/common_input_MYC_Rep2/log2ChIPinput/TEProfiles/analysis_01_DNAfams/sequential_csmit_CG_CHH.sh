#!/bin/bash

csmit -m 5G -c 1 "Rscript TE_DNAfams_Profiles_DNAmeth_commandArgs.R /home/meiosis/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM980986_WT_rep2_CG.wig.bed.gr.tab.bed CGmeth" & sleep 5;
csmit -m 5G -c 1 "Rscript TE_DNAfams_Profiles_DNAmeth_commandArgs.R /home/meiosis/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM980986_WT_rep2_CHH.wig.bed.gr.tab.bed CHHmeth" & sleep 5;
wait

