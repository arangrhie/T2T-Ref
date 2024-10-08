#! /bin/bash

bam=$1
out=$2

if [[ "$#" -lt 2 ]]; then
  echo "Usage: _submit_bam2fastq.sh in.bam out_name"
  echo "  in.bam    input bam file"
  echo "  out_name  output prefix"
  exit -1
fi

cpus=12
mem=32g
name=bam2fq
script=$tools/T2T-Ref/preprocessing/bam2fastq.sh
args="$bam $out"
extra="--gres=lscratch:900"

mkdir -p logs
log=logs/$name.%A.log

echo "\
sbatch --partition=norm -D `pwd` --job-name=$name --time=6:00:00 $extra --cpus-per-task=$cpus --mem=$mem --error=$log --output=$log $script $args"
sbatch --partition=norm -D `pwd` --job-name=$name --time=6:00:00 $extra --cpus-per-task=$cpus --mem=$mem --error=$log --output=$log $script $args

