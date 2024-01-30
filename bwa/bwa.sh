#!/bin/bash

ref=$1
fastq_map=$2
wd=$3
sample=$4

if [[ -z $ref || -z $fastq_map || -z $sample ]]; then
	echo "Usage: ./bwa.sh <ref.fasta> <fastq_map> <out>"
	echo "No <ref.fasta> found. Exit."
	exit -1
fi

module load bwa
module load samtools

cpu=$SLURM_CPUS_PER_TASK
export tmp=/lscratch/${SLURM_JOB_ID}
out=$sample

r1=`cut -f 1 $fastq_map`
r2=`cut -f 2 $fastq_map`

set -o pipefail
set -e
set -x

# Align
cd $tmp
cmd="bwa mem -t $cpu $ref $r1 $r2 > $tmp/$out.sam"
echo -e $cmd
eval $cmd

# Fixmate
cmd="samtools fixmate -m -@$cpu $tmp/$out.sam $tmp/$out.fix.bam && rm $tmp/$out.sam"
echo -e $cmd
eval $cmd

# Sorting
cmd="samtools sort -@$cpu -O bam -o $tmp/$out.bam -T $tmp/$out.tmp $tmp/$out.fix.bam"
echo -e $cmd
eval $cmd
cmd="samtools index $tmp/$out.bam && rm $tmp/$out.fix.bam"
echo -e $cmd
eval $cmd

# Mardup
cmd="samtools markdup -r -@$cpu $tmp/$out.bam $tmp/$out.dedup.bam && rm $tmp/$out.bam.*"
echo -e $cmd
eval $cmd
cmd="samtools index $tmp/$out.dedup.bam"
echo -e $cmd
eval $cmd

# Move Output
cmd="cp $tmp/$out.dedup.bam* ${wd}/ && touch ${wd}/bwa.done"
echo -e $cmd
eval $cmd

# Clean up
cmd="rm $tmp/$out.bam*"
echo -e $cmd
eval $cmd
cmd="rm $tmp/$out.dedup.bam*"
echo -e $cmd
eval $cmd

set +x

echo "Done!"
