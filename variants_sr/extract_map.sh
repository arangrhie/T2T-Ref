#!/bin/bash 

if [[ "$#" -lt 3 ]]; then
  echo "Usage: ./extract_map_calcov_dv.sh sampleInfo.txt refWiY.fa refWoY.fa [idx]"
  echo "  sampleInfo.txt  idx [tab] sampleName [tab] originalBam"
  echo "  refWiY.fa       reference for XY samples. Masking PAR on the Y"
  echo "  refWoY.fo       reference for XX samples. Masking entire Y"
  echo "  idx             idx to proceed"
  echo "                  provide when running this script locally,"
  echo "                  not through the submitter script. OPTIONAL"
  exit -1
fi

set -e
set -o pipefail
set -x

# ARGUMENTS
sampleInfo=$1
refwiY=$2
refwoY=$3
idx=$SLURM_ARRAY_TASK_ID
if [[ -z $idx ]]; then
  idx=$4
fi

PIPELINE='$tools/T2T-Ref'

sample=$(awk -v idx=$idx '$1==idx {print $2}' $sampleInfo)
bam=$(awk -v idx=$idx '$1==idx {print $3}' $sampleInfo)

# Change directory
mkdir -p $sample && cd $sample
mkdir -p logs

# Sex check
echo "# Determine XY - Output: sex.determine.txt"
if [ ! -f sexCheck.done ]; then
  sh $PIPELINE/preprocessing/determine_xy.sh $bam $sample
else
  echo -e "=> Sex check was already done"
fi
echo

sex=`cut -f 3 sex.determine.txt`
echo "== $idx $sample $sex =="
echo

echo "# bam 2 fastq - Output: fastq/${sample}_[12].fq.gz"
if [[ -f bwa.done || ! -f bam2fastq.done ]] ; then
  echo "== extract reads and keep path in fastq_map.fofn =="
  sh $PIPELINE/preprocessing/bam2fastq.sh $bam $sample
  echo "bam2fastq done"
else 
  echo "bam2fastq was already done"
fi
echo

# set ref
if [ "$sex" == "XX" ] ; then
  ref=$refwoY
elif [ "$sex" == "XY" ] ; then
  ref=$refwiY	
fi
echo "Set ref as : $ref"

echo "# BWA - Output: $sample.dedup.pri.bam"
if [ ! -f bwa.done ]; then
  sh $PIPELINE/bwa/bwa_single.sh $ref fastq_map.fofn $sample
  echo "bwa done"
else
  echo "bwa was already done"
fi
echo
