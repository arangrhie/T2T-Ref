#! /bin/bash

if [[ "$#" -lt 2 ]]; then
  echo "./step3.sh sample sex"
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
REF=`cat REF`
MODE=`cat MODE`
OUT=dv_$MODE
CALL_VARIANTS_OUTPUT="dv_$MODE/call_variants_output.tfrecord.gz"
GVCF_TFRECORDS="${OUT}/examples/examples.gvcf.tfrecord@${N_SHARDS}.gz"
PAR="$tools/T2T-Ref/ref/chm13v2.0_PAR.bed"

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

