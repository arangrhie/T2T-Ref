#!/bin/sh

if [[ "$#" -lt 1 ]]; then
  echo "Usage: sh submit_filter.sh out-prefix"
  exit -1
fi

out=$1

cpus=16
mem=6g
name=filt
script=filter.sh
partition=norm
walltime=8:00:00
path=`pwd`

args="$out"
# extra="--array=1-23"
extra="--array=1-3,5" # waiting for chrM

mkdir -p logs
log=logs/$name.%A_%a.log

set -x
sbatch -J $name --cpus-per-task=$cpus --mem=$mem --partition=$partition -D $path $extra --time=$walltime --error=$log --output=$log $script $args > filt.jid
set +x
cat filt.jid

