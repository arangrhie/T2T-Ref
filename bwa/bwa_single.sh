#!/bin/bash

ref=$1
fastq_map=$2
sample=$3

if [[ -z $ref || -z $fastq_map || -z $sample ]]; then
	echo "Usage: ./bwa.sh <ref.fasta> <fastq_map> <out>"
	echo "No <ref.fasta> found. Exit."
	exit -1
fi

module load bwa
module load samtools # v1.21+

cpu=$SLURM_CPUS_PER_TASK
idx=1 # fastq_map.fofn has only 1 line
tmp=/lscratch/${SLURM_JOB_ID}
out=$sample
wd=$PWD

line=`sed -n ${idx}p $fastq_map`
r1=`echo $line | awk '{print $1}'`
r2=`echo $line | awk '{print $2}'`

set -o pipefail
set -e
set -x

bwa mem -t $cpu $ref $r1 $r2 > $tmp/$out.sam

samtools fixmate -m -@$cpu $tmp/$out.sam $tmp/$out.fix.bam && rm $tmp/$out.sam

samtools sort -@$cpu -O BAM -o $tmp/$out.bam -T $tmp/$out.tmp $tmp/$out.fix.bam && rm $tmp/$out.fix.bam
samtools index -@$cpu $tmp/$out.bam

# samtools markdup <input.bam> <output.bam>
samtools markdup -r -@$cpu $tmp/$out.bam $tmp/$out.dedup.bam &&
rm $tmp/$out.bam &&

# Filter out secondary alignments
samtools view -@$cpu -F0x100 -hb --write-index -o $tmp/$out.dedup.pri.bam $tmp/$out.dedup.bam &&
mv $tmp/$out.dedup.pri.ba* ${wd}/ || exit -1
rm -rf $tmp/*
touch bwa.done

set +x
