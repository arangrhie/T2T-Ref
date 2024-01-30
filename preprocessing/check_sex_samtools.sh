#!/bin/bash

inputBam=$1
sampleName=$2

module load samtools

# check sorted by coordinate
is_sorted=$(samtools view -H $inputBam | grep "^@HD" | grep -q "SO:coordinate" && echo "true" || echo "false")

if [ "$is_sorted" == "true" ]; then
    echo "The BAM file is sorted by coordinate."
    ln -s $inputBam $sampleName.oriBam.sorted.bam
    samtools index $sampleName.oriBam.sorted.bam
else
    echo "The BAM file is not sorted by coordinate."
    samtools sort $inputBam -o $sampleName.oriBam.sorted.bam
    samtools index $sampleName.oriBam.sorted.bam
fi

# Check sex
samtools idxstats $sampleName.oriBam.sorted.bam > original.bam.idxstats

x_map=$(grep "chrX" original.bam.idxstats | cut -f 3)
x_len=$(grep "chrX" original.bam.idxstats | cut -f 2)
x_cov=$(echo "scale=10; ${x_map}/${x_len}" | bc)

y_map=$(grep "chrY" original.bam.idxstats | cut -f 3)
y_len=$(grep "chrY" original.bam.idxstats | cut -f 2)
y_cov=$(echo "scale=10; ${y_map}/${y_len}" | bc)

ratio=$(echo "scale=10; ${x_cov}/${y_cov}" | bc)

if (( $(echo "$ratio > 4.00" | bc -l) )); then
    sex="F"
else
    sex="M"
fi
rm $sampleName.oriBam.sorted.bam


# write result
echo -e "${sampleName}\tX:Y_ratio:${ratio}\t${sex}" > sex.determine.txt
