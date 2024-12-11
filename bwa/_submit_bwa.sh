#! /bin/sh

if [[ "$#" -lt 3 ]]; then
  echo "Usage: ./_submit_bwa.sh <ref.fasta> <fastq.map> <sample>"
  exit -1
fi

ref=$1
fastq_map=$2
sample=$3

mkdir -p logs
PIPELINE=$tools/T2T-Ref

if [ ! -f bwa.done ]; then
  cpus=24
  mem=60g
  name=bwa.$sample
  script=$PIPELINE/bwa/bwa_single.sh
  args="$ref $fastq_map $sample"
  partition=norm
  walltime=30:00:00
  log=logs/$name.%A.log
  extra="--gres=lscratch:1200 $extra"
  
  set -x
  sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $PWD $extra --time=$walltime --error=$log --output=$log $script $args | awk '{print $NF}' > bwa.jid
  set +x
  cat bwa.jid
fi
