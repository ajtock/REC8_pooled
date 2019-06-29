#!/bin/bash

# wt,kss colours
# REC8: 'red,red4'
# SPO11-1-oligos: 'dodgerblue2,navy'
# H3K9me2: 'green2,darkgreen'
# mC*: 'orange,orange4'

### Custom annotation legends using textGrob()
## top left
#legendPos <- as.numeric(unlist(strsplit("0.02,0.94",
#                                        split = ",")))
## middle left upper
#legendPos <- as.numeric(unlist(strsplit("0.02,0.70",
#                                        split = ",")))
## middle left
#legendPos <- as.numeric(unlist(strsplit("0.02,0.52",
#                                        split = ",")))
## middle left lower
#legendPos <- as.numeric(unlist(strsplit("0.02,0.40",
#                                        split = ",")))
## middle right
#legendPos <- as.numeric(unlist(strsplit("0.70,0.54",
#                                        split = ",")))

./features_avgProfileRibbon_noCorr.R HypoCHG_DMRs_kss 'HypoCHG DMRs' 2000 2kb '2 kb' 20 20bp 'REC8_HA_Rep2,kss_REC8_HA_Rep1' 'wt REC8-HA,kss REC8-HA' 'red,red4' '0.02,0.40'

./features_avgProfileRibbon_noCorr.R HypoCHG_DMRs_kss 'HypoCHG DMRs' 2000 2kb '2 kb' 20 20bp 'SPO11_1_oligos_RPI1,kss_SPO11_1_oligos_RPI34' 'wt SPO11-1,kss SPO11-1' 'dodgerblue2,navy' '0.02,0.40'

./features_avgProfileRibbon_noCorr.R HypoCHG_DMRs_kss 'HypoCHG DMRs' 2000 2kb '2 kb' 20 20bp 'SPO11_1_oligos_RPI3,kss_SPO11_1_oligos_RPI34' 'wt SPO11-1_Rep2,kss SPO11-1_Rep1' 'dodgerblue2,navy' '0.02,0.40'

./features_avgProfileRibbon_noCorr.R HypoCHG_DMRs_kss 'HypoCHG DMRs' 2000 2kb '2 kb' 20 20bp 'SPO11_1_oligos_RPI8,kss_SPO11_1_oligos_RPI35' 'wt SPO11-1_Rep3,kss SPO11-1_Rep2' 'dodgerblue2,navy' '0.02,0.40'

./features_avgProfileRibbon_noCorr.R HypoCHG_DMRs_kss 'HypoCHG DMRs' 2000 2kb '2 kb' 20 20bp 'H3K9me2,kss_H3K9me2' 'wt H3K9me2,kss H3K9me2' 'green2,darkgreen' '0.02,0.40'


