#!/bin/bash

# indir="/Users/ricard/data/NMT-seq_EB/acc/raw/old"
# outdir="/Users/ricard/data/NMT-seq_EB/acc/raw"
indir="/Users/ricard/data/NMT-seq_EB/met/raw/old"
outdir="/Users/ricard/data/NMT-seq_EB/met/raw"

cd $indir
for file in *.gz
do
	plate=$( echo $file | cut -d "_" -f 3 )
	sample=$( echo $file | cut -d "_" -f 4 )
	R=$( echo $file | cut -d "_" -f 6 )
	echo "Processing ${plate}_${sample}_${R}..."
	zmore ${indir}/${file} | awk '{ if ($4>0 || $5>0)  {print $1,$2,$4,$5}}' | sort -k 1,1 -k2,2n | gzip > ${outdir}/${plate}_${sample}_${R}.tsv.gz
done