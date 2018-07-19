#!/bin/bash

# Usage:
# zcat_fastq_gz_gzip_best.sh 160902_REC8_MYC_Rep1_ChIP 180427_REC8_MYC_Rep1_ChIP R1 pooled_REC8_MYC_Rep1_ChIP_R1

inName1=$1
inName2=$2
readN=$3
outName=$4

if [ ! -f "$outName.fastq.gz" ]; then
  zcat $inName1"_"$readN".fastq.gz" \
       $inName2"_"$readN".fastq.gz" \
       | gzip -c -k --best > $outName.fastq.gz;
else
  echo "skipping $outName"
fi

