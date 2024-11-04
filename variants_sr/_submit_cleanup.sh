#!/bin/bash

sample=$1

if [[ "$#" -lt 1 ]]; then
  echo "Usage: ./_submit_cleanup.sh sample"
  exit -1
fi

set -e
set -o pipefail
PIPELINE=$tools/T2T-Ref

dv_step3_jid=`cat step3.jid`

mkdir -p logs

cpus=2
mem=1g
partition="quick"
walltime=5:00
path=$PWD
name=$sample.cleanup
log=logs/${name}.%A.log
script=$PIPELINE/variants_sr/cleanup.sh
args="$sample"
extra="--dependency=afterok:${dv_step3_jid}"

set -x
sbatch -J $name --cpus-per-task=$cpus --mem=$mem \
  --partition=$partition -D $path $extra \
  --time=$walltime --error=$log --output=$log \
  $script $args | tail -n1 | awk '{print $NF}' > cleanup.jid
set +x

