#! /bin/bash

if [[ "$#" -lt 1 ]]; then
  echo "Usage: sh submit_glnexus.sh in.list in.Y.list"
  echo "  in.list    gVCF list, all samples"
  echo "  in.Y.list  gVCF list, for samples having Y (exclude XX samples). OPTIONAL"
  echo "             Will skip calling on Y chromosome if not provided."
  exit -1
fi


cpus=48
mem=260g
name=glnexus
script=glnexus.sh
partition=norm
walltime=24:00:00
path=`pwd`

lst=$1
args="$lst"
# extra="--array=1-23,25" # Failed with --gres=lscratch:3000 on ~1100 samples. Need 4x, configuration unavailable.
extra="--array=25"

mkdir -p logs
log=logs/$name.%A_%a.log

set -x
sbatch -J $name --cpus-per-task=$cpus --mem=$mem --partition=$partition -D $path $extra --time=$walltime --error=$log --output=$log $script $args
set +x

lst=$2
if [[ -z $lst ]]; then
  echo "no Y list given. End here."
  exit 0
fi

args="$lst"
extra="--array=24"

set -x
sbatch -J $name --cpus-per-task=$cpus --mem=$mem --partition=$partition -D $path $extra --time=$walltime --error=$log --output=$log $script $args
set +x
