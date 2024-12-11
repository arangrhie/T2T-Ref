#! /bin/bash

if [[ "$#" -lt 3 ]]; then
  echo "./step3.sh sample sex par"
  exit -1
fi

if [[ -z $SLURM_CPUS_PER_TASK ]]; then
  CPU=12
else
  CPU=$SLURM_CPUS_PER_TASK
fi

set -o pipefail
set -e
set -x

# Update to DeepVariant v1.6.1 on Sep. 26 2024
module load deepvariant/1.6.1

N_SHARDS=`cat N_SHARD`

SAMPLE=$1
SEX=$2
PAR=$3

REF=`cat REF`
MODE=`cat MODE`
OUT=dv_$MODE
CALL_VARIANTS_OUTPUT="dv_$MODE/call_variants_output.tfrecord.gz"
GVCF_TFRECORDS="${OUT}/examples/examples.gvcf.tfrecord@${N_SHARDS}.gz"

if [[ $SEX = "XY" ]]; then
  haploid="--haploid_contigs chrX,chrY --par_regions_bed ${PAR}"
fi

postprocess_variants \
  --ref "${REF}" $haploid \
  --infile "${CALL_VARIANTS_OUTPUT}" \
  --outfile "${OUT}.$SAMPLE.vcf.gz" \
  --gvcf_outfile "${OUT}.$SAMPLE.gvcf.gz" \
  --nonvariant_site_tfrecord_path $GVCF_TFRECORDS \
  --cpus $CPU \
  --sample_name $SAMPLE  && 
touch deepvariant.step3.done || exit -1
set +x

echo "deepvariant done"

