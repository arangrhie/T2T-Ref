#!/bin/bash

sample=$1
target=$2

if [[ "$#" -lt 2 ]]; then
  echo "Usage: ./cal_cov_target.sh sample target.txt"
  echo "  sample      Prefix used in sample.dedup.cram"
  echo "  target.txt  Containing [ region ] [ bed ]"
  echo "*Requires bed present under \$PIPELINE/ref/cal_bed/"
  exit -1
fi

set -e
set -o pipefail

PIPELINE=$tools/T2T-Ref

ln=`wc -l $target | awk '{print $1}'`

# target.txt begins from line 2. line 1 is the background
for i in $(seq 2 $ln)
do
  prefix=`sed -n ${i}p $target | awk '{print $1}'`
  bed=`sed -n ${i}p $target | awk '{print $2}'`
  echo $prefix $bed
  # append path
  bed=$PIPELINE/ref/cal_bed/$bed
  
  if [ ! -f $prefix.coverage_results.bed ]; then
    sh $PIPELINE/calDepth/cal_cov.sh $sample $bed $prefix
  fi
done

touch cal.target.done