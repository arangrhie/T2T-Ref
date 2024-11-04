#!/bin/bash

ref=$1

if [ -z $ref ]; then
	echo "Usage: ./bwa_index.sh <ref.fasta>"
	echo "No <ref.fasta> given. Exit."
	exit -1
fi

module load bwa
module load samtools  # v1.21+
module load bedtools

cpu=$SLURM_CPUS_PER_TASK
if [[ -z $cpu ]]; then
  cpu=8
fi

set -x
if [ ! -f $ref.bwt ]; then
  bwa index $ref
fi

if [ ! -f $ref.fai ]; then
  samtools faidx -@$cpu $ref
fi 

set +x
