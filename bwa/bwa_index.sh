#!/bin/bash

ref=$1

if [ -z $ref ]; then
	echo "Usage: ./bwa_index.sh <ref.fasta>"
	echo "No <ref.fasta> given. Exit."
	exit -1
fi

module load bwa
module load samtools
module load bedtools

if [ ! -f $ref.bwt ]; then
echo "
bwa index $ref"
bwa index $ref
fi

if [ ! -f $ref.fai ]; then
echo "
samtools faidx $ref"
samtools faidx $ref
fi 

# make window for future calcation of depth
if [ ! -f $ref.fai.win10k.bed ]; then
echo "bedtools makewindows -g $ref.fa.fai -w 10000 > $ref.fai.win10k.bed"
bedtools makewindows -g $ref.fai -w 10000 > $ref.fai.win10k.bed
fi
