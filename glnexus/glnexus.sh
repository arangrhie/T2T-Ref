#!/bin/sh

ulimit -Sn 65536

if [[ $# -lt 1 ]]; then
  echo "Usage: ./glnexus.sh in_gvcf.list [i]"
  echo "Merge gVCFs by chrs"
  echo "in_gvcf.list  List of gvcf files to merge"
  echo "i             ith chromosome to merge (numeric). OPTIONAL"
  echo "              1-25. 23=X, 24=Y, 25=M."
  exit -1
fi

list=$1

arr=$SLURM_ARRAY_TASK_ID
if [[ -z $arr ]]; then
  arr=$2
fi

if [[ -z $arr ]]; then
  echo "no SLURM_ARRAY_ID found. i is empty. exit"
  exit -1
fi

if [[ $arr == "23" ]]; then
  arr=X
elif [[ $arr == "24" ]]; then
  arr=Y
elif [[ $arr == "25" ]]; then
  arr=M
fi

chr="chr$arr"

# Uses too much space. ~10TB / chromosome.
# tmp=/lscratch/${SLURM_JOB_ID}/GLnexus.$chr.DB didn't work out.
tmp=GLnexus.$chr.DB
cpu=24
mem=220

out=glnexus.$chr.bcf
bed=$tools/T2T-Ref/ref/chm13v2.0_$chr.bed

module load glnexus/1.4.1

set -e
set -x
glnexus_cli --config DeepVariantWGS --threads $cpu --mem-gbytes $mem --dir $tmp --bed $bed --list $1 > $out
rm -r $tmp
set +x

