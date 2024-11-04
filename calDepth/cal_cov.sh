#!/bin/bash
set -o pipefail
set -e

sample=$1
bed=$2
prefix=$3

bam=$sample.dedup.pri.bam
cpu=$SLURM_CPUS_PER_TASK

module load samtools
module load bedtools

set -x
samtools depth -@ $cpu -b $bed $bam | \
  awk '{print $1,$2-1,$2,$3}' OFS='\t' - | \
  bedtools map -a $bed -b - -c 4 -o median > $prefix.coverage_results.bed
set +x
