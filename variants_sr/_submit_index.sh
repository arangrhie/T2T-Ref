#!/bin/bash

if [[ "$#" -lt 2 ]]; then
  echo "Usage: ./_submit_index.sh refWiY refWoY"
  exit -1
fi

set -e
set -o pipefail

PIPELINE=$tools/T2T-Ref
refWiY=$1 # $PIPELINE/ref/chm13v2.0_masked_DJ_5S_rDNA_PHR_hardMaskY.fa
refWoY=$2 # $PIPELINE/ref/chm13v2.0_masked_DJ_5S_rDNA_PHR_noY.fa

## Submission env variables
mkdir -p logs
cpus=4
mem=10g
partition=quick
walltime=4:00:00
script=$PIPELINE/bwa/bwa_index.sh

for sexChr in XX XY
do
  if [ "$sexChr" == "XX" ] ; then
    ref=$refWoY
  elif [ "$sexChr" == "XY" ] ; then
    ref=$refWiY
  fi
  if [[ ! -f $ref.bwt ]]; then
	args="$ref"
	name=index.$sexChr
	log=logs/$name.%A.log

	set -x
	sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D `pwd` --time=$walltime --error=$log --output=$log $script $args | awk '{print $NF}' > ./$name.jid
	set +x
	cat $name.jid
  fi
done
