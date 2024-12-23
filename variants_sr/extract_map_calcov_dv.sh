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

# ARGUMENTS
sampleInfo=$1
refWiY=$2
refWoY=$3
idx=$SLURM_ARRAY_TASK_ID
if [[ -z $idx ]]; then
  idx=$4
fi

PIPELINE=$tools/T2T-Ref

sample=$(awk -v idx=$idx '$1==idx {print $2}' $sampleInfo)
bam=$(awk -v idx=$idx '$1==idx {print $3}' $sampleInfo)

# Change directory
mkdir -p $sample && cd $sample
mkdir -p logs

# Sex check
echo "# Determine XY - Output: sex.determine.txt"
if [ ! -f sexCheck.done ]; then
  set -x
  sh $PIPELINE/preprocessing/determine_xy.sh $bam $sample
  set +x
else
  echo -e "=> Sex check was already done"
fi
echo

sex=`cut -f 3 sex.determine.txt`
echo "== $idx $sample $sex =="
echo

echo "# bam 2 fastq - Output: lscratch/.../${sample}_[12].fq.gz"
if [[ ! -f bam2fastq.$sample.done ]] ; then
  echo "== extract reads and keep path in fastq_map.fofn =="
  set -x
  sh $PIPELINE/preprocessing/bam2fastq.sh $bam $sample
  set +x
  echo "bam2fastq done"
else 
  echo "bam2fastq was already done"
fi
echo

# set ref
if [ "$sex" == "XX" ] ; then
  ref=$refWoY
elif [ "$sex" == "XY" ] ; then
  ref=$refWiY
fi
echo "Set ref as : $ref"

echo "# BWA - Output: $sample.dedup.pri.bam"
if [ ! -f bwa.done ]; then
  set -x
  sh $PIPELINE/bwa/bwa_single.sh $ref fastq_map.fofn $sample
  set +x
  echo "bwa done"
else
  echo "bwa was already done"
fi
echo

# Double check the new bam file exists
bam=$sample.dedup.pri.bam
if [[ ! -f $bam ]]; then
  echo "No $bam found. Exit."
  exit -1
fi

### BAM files are ready. Submit DepthCall and DeepVariant jobs
### Note that logs will be present under $sample/logs/
if [ ! -f cal.target.done ] ; then
	echo "# Submit Calculate Coverage"
  set -x
	sh $PIPELINE/calDepth/_submit_cal_depth.sh $sample
  set +x
fi

echo "# Submit DeepVariant"
if [ ! -f deepvariant.step3.done ] ; then
  set -x
	sh $PIPELINE/deepvariant/_submit_deepvariant.sh $ref $bam WGS $sample $sex
  set +x
fi

if [[ -f bwa.done ]]; then
  echo "# Clean up fastq"
  #sh $PIPELINE/variants_sr/_submit_cleanup.sh $sample
  if [[ -d fastq ]]; then
    rm -r fastq
  fi
  rm bam2fastq.$sample.done
fi
echo

