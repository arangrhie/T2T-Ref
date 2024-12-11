#!/bin/sh

if [[ $# -lt 1 ]]; then
  echo "Usage: ./filter.sh out-prefix i"
  echo "Filters for GQ20 and ALT."
  echo "Assumes input glnexus.chrN.bcf exists in the same path."
  echo "  out-prefix  output prefix"
  echo "  i           ith chromosome num. 23=X, 24=Y, 25=M. OPTIONAL."
  exit 0
fi

out=$1

threads=$SLURM_CPUS_PER_TASK
if [[ -z $threads ]]; then
  threads=12
fi

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

module load bcftools # samtools/1.21

set -e
set -x
bcftools view --threads $threads -i 'MIN(FMT/GQ)>=20 & GT[*]="alt"' glnexus.$chr.bcf -W -Oz -o $out.$chr.GQ20_ALT.vcf.gz
if [[ -f parents.list ]]; then
  bcftools view --threads $threads -S parents.list --force-samples $out.$chr.GQ20_ALT.vcf.gz | \
    bcftools +fill-tags --threads $threads -W -Oz -o $out.$chr.GQ20_ALT.AF.parents.vcf.gz
fi
bcftools +fill-tags $out.$chr.GQ20_ALT.vcf.gz | \
  bcftools view --threads $threads --samples "102-00001-02" -W -Oz -o $out.$chr.GQ20_ALT.AF.all.vcf.gz
set +x