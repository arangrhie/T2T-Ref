#! /bin/sh

if [[ "$#" -lt 5 ]]; then
  echo "Usage: ./_submit_bwa.sh <ref.fasta> <fastq.map> <sample> <wd> <sex> <line_num>"
  exit -1
fi

ref=$1
fastq_map=$2
sample=$3
wd=$4
line_num=$6
sexChr=$5

out=$wd/$sample


mkdir -p $wd/logs
PIPELINE=$tools/T2T-Ref

if [ ! -f $wd/bwa.done ]; then
  cpus=24
  mem=120g
  name=map.$sample
  script=$PIPELINE/bwa/bwa.sh
  args="$ref $fastq_map $wd $sample"
  partition=norm
  walltime=2-0
  path=$wd
  log=$wd/logs/$name.%A_%a.log
  extra="--gres=lscratch:1000 $extra"
  
  echo "\
  sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
  sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args | awk '{print $NF}' > $wd/bwa.jid
fi
