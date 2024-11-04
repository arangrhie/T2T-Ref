#!/bin/bash

if [[ "$#" -lt 2 ]]; then
  echo "Usage: determine_xy.sh inputBam sampleName"
  exit -1
fi

inputBam=$1
sampleName=$2

cpu=$SLURM_CPUS_PER_TASK
if [[ -z $cpu ]]; then
  cpu=8
fi

tmp=/lscratch/${SLURM_JOB_ID}

module load samtools

set -x

# check sorted by coordinate
if [[ -f $inputBam.bai || -f ${inputBam/.bam/.bai} || -f $inputBam.csi ]]; then
    echo "The BAM file is sorted..."
    ln -s $inputBam $sampleName.sorted.bam
    ln -s $inputBam.bai $sampleName.sorted.bam.bai
else
    echo "The BAM file is NOT sorted..."
    ## Dump the sorted bam in lscratch so we don't need to worry about it
    samtools sort -@$cpu -O BAM -o $tmp/$sampleName.sorted.bam -T $tmp/$sampleName.tmp $inputBam
    ln -s $tmp/$sampleName.sorted.bam $sampleName.sorted.bam
    samtools index $sampleName.sorted.bam
fi

set -e

# Check sex
samtools idxstats -@$cpu $sampleName.sorted.bam > original.bam.idxstats

x_map=$(grep -E "^chrX|^X" original.bam.idxstats | cut -f 3)
x_len=$(grep -E "^chrX|^X" original.bam.idxstats | cut -f 2)
x_cov=$(echo "scale=10; ${x_map}/${x_len}" | bc)

y_map=$(grep -E "^chrY|^Y" original.bam.idxstats | cut -f 3)
y_len=$(grep -E "^chrY|^Y" original.bam.idxstats | cut -f 2)
y_cov=$(echo "scale=10; ${y_map}/${y_len}" | bc)

ratio=$(echo "scale=10; ${x_cov}/${y_cov}" | bc)

if (( $(echo "$ratio > 4.00" | bc -l) )); then
    sex="XX" ## Use XX instead of F
else
    sex="XY" ## Use XY instead of M
fi

# Remove symlinks and index
rm $sampleName.sorted.bam*

# write result
echo -e "${sampleName}\tX:Y_ratio:${ratio}\t${sex}" > sex.determine.txt

touch sexCheck.done