#!/bin/bash

# Concatenate per-chromosome bedgraphs and convert to BED-like format (but with 1-based coordinates)

# Example usage via condor submission system on hydrogen node7:
# csmit -m 10G -c 1 "./2_cat_chr_bg_merged.sh REC8_MYC_Rep2_ChIP"

prefix=$1

# Concatenate chromosome bedgraph files to create genome bedgraph file (0-based)
cat $(find ./ -name ${prefix}"_norm_Chr*_coverage.bedgraph" | sort -V) > ${prefix}_norm_allchrs_coverage.bedgraph
sed -i '1i track type=bedGraph' ${prefix}_norm_allchrs_coverage.bedgraph
# Convert 0-based genome bedgraph file to 1-based bed-like file containing coordinates and library-normalised coverage
awk 'BEGIN {OFS="\t"}; {print $1, $3, $3, $4}' ${prefix}_norm_allchrs_coverage.bedgraph > ${prefix}_norm_allchrs_coverage_coord_tab_tmp.bed
tail -n +2 ${prefix}_norm_allchrs_coverage_coord_tab_tmp.bed > ${prefix}_norm_allchrs_coverage_coord_tab.bed
rm ${prefix}_norm_allchrs_coverage_coord_tab_tmp.bed

