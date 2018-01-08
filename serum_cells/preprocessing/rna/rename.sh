#!/bin/bash

# indir="/hps/nobackup/stegle/users/ricard/gastrulation/met/original/scMT"
# outdir="/hps/nobackup/stegle/users/ricard/gastrulation/met/raw/scMT"
# for file in ${indir}/*.gz
# do
# 	sample=$( echo $file | cut -d "_" -f 2 )
# 	R=$( echo $file | cut -d "_" -f 7 )
# 	echo "Processing ${sample}_${R}..."
# 	zmore $file | cut -f 1,2,5,6 | gzip > ${outdir}/${sample}_${R}.txt.gz
# done


# indir="/Users/ricard/data/NMT-seq/rna/raw"
# outdir="/Users/ricard/data/NMT-seq/rna/raw"

indir="/Users/ricard/data/NMT-seq_EB/rna/raw/old"
outdir="/Users/ricard/data/NMT-seq_EB/rna/raw"

cd $indir
for file in *bam
do
	sample=$( echo $file | cut -d "_" -f 4,5 )
	# ln -s ${file} ${outdir}/${sample}.bam
	mv ${file} ${outdir}/${sample}.bam
done