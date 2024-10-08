#!/bin/bash

if [[ "$#" -lt 2 ]]; then
  echo "Usage: ./_submit_map_variant_sr.sh sampleInfo.txt sampleIdx"
  echo "  sampleInfo.txt  Sample information table."
  echo "                  [sampleIdx] [sampleId] [XXorXY]"
  echo "  sampleIdxs      Sample index to submit this pipeline."
  echo "                  e.g. 1-30 for samples from 1 to 30 in sampleInfo.txt"
  exit -1
fi

set -e
set -o pipefail

# Adjust to match local path
PIPELINE=$tools/T2T-Ref
refWiY=$PIPELINE/ref/chm13v2.0_masked_DJ_5S_rDNA_PHR_hardMaskY.fa
refWoY=$PIPELINE/ref/chm13v2.0_masked_DJ_5S_rDNA_PHR_noY.fa

sampleInfo=$1
sampleIndex=$2

array="--array=$2"		# Arguement passed on for --array. e.g. 1-30 for "--array=1-30"

## Submission env variables
mkdir -p logs
name=srWGS_wrapper

script=$PIPELINE/variants_sr/map_calcov_dv.sh

args="$sampleInfo $refWiY $refWoY"
log=logs/$name.%A_%a.log
cpus=24
mem=120g
partition=norm
walltime=24:00:00
extra="--gres=lscratch:800"

set -x
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus \
 -D `pwd` \
 --time=$walltime --error=$log --output=$log \
 $array $extra $script $args | awk '{print $NF}' > ./$name.$sampleIndex.jid
set +x

cat $name.$sampleIndex.jid
