#!/bin/bash

if [[ "$#" -lt 1 ]]; then
  echo "Usage: ./_submit_dv_only.sh sample"
  echo "Submit DeepVariant submitter for each sample"
  exit -1
fi

PIPELINE=$tools/T2T-Ref
refWiY=$PIPELINE/ref/chm13v2.0_masked_DJ_5S_rDNA_PHR_PAR_wi_rCRS.fa
refWoY=$PIPELINE/ref/chm13v2.0_masked_DJ_5S_rDNA_PHR_noY_wi_rCRS.fa

sample=$1
cd $sample

if [ ! -f sexCheck.done ]; then
  echo "[[ ERROR ]] :: No sexCheck.done found for $sample"
  exit -1
fi

sex=`cut -f 3 sex.determine.txt`
ref=""
if [ "$sex" == "XX" ] ; then
  ref=$refWoY
elif [ "$sex" == "XY" ] ; then
  ref=$refWiY
fi
echo "Set ref as : $ref"

bam=$sample.dedup.pri.bam

# if [[ ! -s $bam ]]; then
#   echo "No $bam found. Exit."
#   exit -1
# fi

# _submit_deepvariant.sh looks for deepvariant.step*.done
sh $PIPELINE/deepvariant/_submit_deepvariant.sh $ref $bam WGS $sample $sex
