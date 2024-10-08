#! /bin/sh

if [[ "$#" -lt 1 ]]; then
  echo "Usage: ./_submit_filt.sh in.bam"
  exit -1
fi

bam=$1

mkdir -p logs
PIPELINE=$tools/T2T-Ref

cpus=24
mem=8g
name=filt.${bam/.bam/}
script=$PIPELINE/bwa/filt.sh
args="$bam"
partition=quick
walltime=1:00:00
path=$PWD
log=logs/$name.%A.log
extra="--gres=lscratch:110"

set -x
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args
set +x
