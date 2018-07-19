#!/bin/bash

csmit -m 10G -c 1 "Rscript TE_RNAfams_Profiles_DNAmeth_commandArgs.R /home/meiosis/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM980986_WT_rep2_CG.wig.bed.gr.tab.bed CGmeth" & sleep 30;
csmit -m 10G -c 1 "Rscript TE_RNAfams_Profiles_DNAmeth_commandArgs.R /home/meiosis/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM980986_WT_rep2_CHG.wig.bed.gr.tab.bed CHGmeth" & sleep 30;
csmit -m 10G -c 1 "Rscript TE_RNAfams_Profiles_DNAmeth_commandArgs.R /home/meiosis/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM980986_WT_rep2_CHH.wig.bed.gr.tab.bed CHHmeth" & sleep 30;

csmit -m 10G -c 1 "Rscript TE_RNAfams_Profiles_DNAmeth_commandArgs.R /home/meiosis/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981060_suvh456_CG.bed.gr.tab.bed kss_CGmeth" & sleep 30;
csmit -m 10G -c 1 "Rscript TE_RNAfams_Profiles_DNAmeth_commandArgs.R /home/meiosis/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981060_suvh456_CHG.bed.gr.tab.bed kss_CHGmeth" & sleep 30;
csmit -m 10G -c 1 "Rscript TE_RNAfams_Profiles_DNAmeth_commandArgs.R /home/meiosis/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981060_suvh456_CHH.bed.gr.tab.bed kss_CHHmeth" & sleep 30;

csmit -m 10G -c 1 "Rscript TE_RNAfams_Profiles_DNAmeth_commandArgs.R /home/meiosis/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981003_cmt3_CG.wig.bed.gr.tab.bed cmt3_CGmeth" & sleep 30;
csmit -m 10G -c 1 "Rscript TE_RNAfams_Profiles_DNAmeth_commandArgs.R /home/meiosis/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981003_cmt3_CHG.wig.bed.gr.tab.bed cmt3_CHGmeth" & sleep 30;
csmit -m 10G -c 1 "Rscript TE_RNAfams_Profiles_DNAmeth_commandArgs.R /home/meiosis/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM981003_cmt3_CHH.wig.bed.gr.tab.bed cmt3_CHHmeth" & sleep 30;

wait
