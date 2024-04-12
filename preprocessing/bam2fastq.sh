#!/bin/bash 
set -e
set -o pipefail

module load bedtools
module load samtools

# ARGUMENTS
origianl_bam=$1
sampleName=$2
cpu=$SLURM_CPUS_PER_TASK
tmp=/lscratch/${SLURM_JOB_ID}

outFastqdir=fastq

# Make output dir
mkdir -p $outFastqdir

# fastq_map.fofn
echo -e "${outFastqdir}/${sampleName}_1.fq.gz\t${outFastqdir}/${sampleName}_2.fq.gz" > fastq_map.fofn

set -x
# sort by read name 
samtools sort -@ $cpu -n -o $tmp/${sampleName}_nameSort.bam ${origianl_bam}

if [ ! -f bam2fastq.done ]; then
# bam to fastq
    samtools fastq -@ $cpu $tmp/${sampleName}_nameSort.bam \
        -1 $tmp/${sampleName}_1.fq \
        -2 $tmp/${sampleName}_2.fq && rm $tmp/${sampleName}_nameSort.bam

    pigz -c $tmp/${sampleName}_1.fq > ${outFastqdir}/${sampleName}_1.fq.gz
    pigz -c $tmp/${sampleName}_2.fq > ${outFastqdir}/${sampleName}_2.fq.gz
fi
rm $tmp/${sampleName}_[12].fq
touch bam2fastq.done
set +x
cd ../