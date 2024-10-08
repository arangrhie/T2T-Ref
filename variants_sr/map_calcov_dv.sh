#!/bin/bash 

if [[ "$#" -lt 3 ]]; then
  echo "Usage: ./map_calcov_dv.sh sampleInfo.txt refWiY.fa refWoY.fa [idx]"
  echo "  sampleInfo.txt  idx [tab] sampleName [tab] XXorXY"
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
refWiY=$2
refWoY=$3
idx=$SLURM_ARRAY_TASK_ID
if [[ -z $idx ]]; then
  idx=$4
fi

PIPELINE=$tools/T2T-Ref

sample=$(awk -v idx=$idx '$1==idx {print $2}' $sampleInfo)
sex=$(awk -v idx=$idx '$1==idx {print $3}' $sampleInfo)

# Change directory
mkdir -p $sample && cd $sample
mkdir -p logs

echo "== $idx $sample $sex =="
echo

# set ref
if [ "$sex" == "XX" ] ; then
  ref=$refWoY
elif [ "$sex" == "XY" ] ; then
  ref=$refWiY
fi
echo "Set ref as : $ref"

echo "# BWA - Output: $sample.dedup.cram"
if [ ! -f bwa.done ]; then
  sh $PIPELINE/bwa/bwa_single.sh $ref fastq_map.fofn $sample
  echo "bwa done"
else
  echo "bwa was already done"
fi
echo

# Double check the new bam file exists
bam=$sample.dedup.bam
if [[ ! -f $bam ]]; then
  echo "No $bam found. Exit."
  exit -1
fi

### BAM files are ready. Submit DepthCall and DeepVariant jobs
### Note that logs will be present under $sample/logs/
if [ ! -f cal.target.done ] ; then
	echo "# Submit Calculate Coverage"
	sh $PIPELINE/calDepth/_submit_cal_depth.sh $sample
fi

echo "# Submit DeepVariant"
if [ ! -f deepvariant.step3.done ] ; then
	sh $PIPELINE/deepvariant/_submit_deepvariant.sh $ref $bam WGS $sample $PWD
fi

# Disable until deepvariant devugging is done
# echo "# Clean up bam"
# sh $PIPELINE/variants_sr/_submit_cleanup.sh $sample
# echo

