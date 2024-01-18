#!/bin/sh

if [[ "$#" -lt 2 ]]; then
  echo "Usage: extract_map.sh bam gender [line_num]"
  echo "  bam      input bam file from a different reference (e.g. hg19 or hg38)"
  echo "  gender   XX or XY"
  exit 0
fi

i=$SLURM_ARRAY_TASK_ID
TMP=/lscratch/$SLURM_JOB_ID
if [[ -z $i ]]; then
  i=$1
  TMP="."
fi

if [[ -z $i ]]; then
  echo "No slurm array task ID nor line_num given. Exit"
  exit -1
fi

reads=$1
gender=$2

sample=`echo $bam | sed 's/\.bam$//g'`

if [[ $gender -eq "xy" ]] || [[ $gender -eq "XY" ]]; then
  ref=`cat ref_XY`
  rep=`cat ref_XY_rep`
else
  ref=`cat ref_XX`
  rep=`cat ref_XX_rep`
fi

echo "Aligning $sample ($gender) to $ref"

map=map-ont # map-pb or map-ont

mkdir -p $TMP/out
mkdir -p $sample
cd $sample
out=${sample}_mask

set -o pipefail
set -e

ls $reads

cpus=$SLURM_CPUS_PER_TASK

module load samtools
module load winnowmap

if ! [[ -s $out.bam ]]; then
  set -x
  samtools fastq -@$cpus -T Mm,Ml $reads | \
    winnowmap --MD -W $rep -y -Y -ax $map -t$cpus $ref - > $TMP/$out.sam

  samtools sort -@$cpus -m2G -T $TMP/$out.tmp -O bam -o $out.bam $TMP/$out.sam
  samtools index $out.bam


set +x
fi

cd ../
