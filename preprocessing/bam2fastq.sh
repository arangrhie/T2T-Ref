#!/bin/bash 
set -e
set -o pipefail

module load bedtools
module load samtools

# ARGUMENTS
origianl_bam=$1
outFastqName=$2

# sort by read name 
samtools sort -@ $SLURM_CPUS_PER_TASK -n -o ${outFastqName}_nameSort.bam ${origianl_bam}

if [ ! -f fastq/bam2fastq.done ]; then
# bam to fastq
samtools fastq -@ $SLURM_CPUS_PER_TASK ${outFastqName}_nameSort.bam \
                      -1 ${outFastqName}_1.fq \
                      -2 ${outFastqName}_2.fq && rm ${outFastqName}_nameSort.bam
bgzip -f ${outFastqName}_1.fq 
bgzip -f ${outFastqName}_2.fq

if [ -f ${outFastqName}_1.fq.gz ] && [ -f ${outFastqName}_2.fq.gz ];then
	touch fastq/bam2fastq.done
fi

fi

echo -e "${outFastqName}_1.fq.gz\t${outFastqName}_2.fq.gz" > fastq_map
