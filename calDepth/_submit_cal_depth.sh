#!/bin/bash

sample=$1
extra=$2 # optional. e.g. "--dependency=afterok:[jid]"

if [[ "$#" -lt 1 ]]; then
  echo "Usage: ./_submit_cal_depth.sh sample [extra]"
  exit -1
fi

set -e
set -o pipefail
PIPELINE=$tools/T2T-Ref

mkdir -p logs

# Schedular variables
cpus=20
partition="quick"
walltime=3:59:00
path=$PWD
date=`date | sed -e s'/ /_/g'`


target=$PIPELINE/calDepth/target.txt
bed_prefix=`head -n1 $target`
prefix=`echo $bed_prefix | awk '{print $1}'`
bed=`echo $bed_prefix | awk '{print $2}'`
bed="$PIPELINE/ref/cal_bed/$bed"

echo "## Submit Background Coverage Calculation"

mem=100g
name=${sample}_${prefix}
log=logs/${name}.${date}.%A.log
script=$PIPELINE/calDepth/cal_cov.sh
args="$sample $bed $prefix"

set -x
sbatch -J $name --cpus-per-task=$cpus --mem=$mem \
  --partition=$partition -D $path $extra \
  --time=$walltime --error=$log --output=$log \
  $script $args | tail -n1 | awk '{print $NF}' > cal.background.jid
set +x
echo


echo "## Submit Target Coverage Calculation"
mem=10g
walltime=3:00:00
name=${sample}_target
log=logs/${name}.${date}.%A.log
script=$PIPELINE/calDepth/cal_cov_target.sh
args="$sample $target"

set -x
sbatch -J $name --cpus-per-task=$cpus --mem=$mem \
  --partition=$partition -D $path $extra \
  --time=$walltime --error=$log --output=$log \
  $script $args | tail -n1 | awk '{print $NF}' > cal.target.jid
set +x
echo