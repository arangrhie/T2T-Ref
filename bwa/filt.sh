#!/bin/bash

set -e
set -o pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: filt.sh in.bam"
  echo "Filter out secondary alignments."
  echo "  output: in.pri.bam"
  exit -1
fi

bam=$1
out=`basename $bam`
out=${out/.bam/}

module load samtools

cpu=$SLURM_CPUS_PER_TASK
tmp=/lscratch/${SLURM_JOB_ID}

set -x
samtools view -@$cpu -F0x100 -hb -o $tmp/$out.pri.bam --write-index $bam &&
mv $tmp/$out.pri.ba* . &&
rm $bam && rm $bam.bai || echo "Failed..."
set +x
echo "Done!"
