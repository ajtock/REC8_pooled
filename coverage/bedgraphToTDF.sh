#!/bin/bash

#######################################################################################
# Convert concatenated bedgraph files to TDF format for visualisation in IGV          #
#######################################################################################

i=$1

igvtools toTDF ${i}_norm_allchrs_coverage.bedgraph ${i}_norm_allchrs_coverage.bedgraph.tdf /projects/ajt200/TAIR10/tair10.chrom.sizes

