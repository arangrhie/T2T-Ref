#!/bin/bash

set -e
set -o pipefail

module load deepvariant/1.5.0 || exit 1
module load parallel

wd=$1
SAMPLE=$2

echo $wd

[[ -d /lscratch/${SLURM_JOB_ID} ]] && cd /lscratch/${SLURM_JOB_ID}

REF=`cat ${wd}/REF`
BAM=`cat ${wd}/BAM`
MODE=`cat ${wd}/MODE`
OUT=dv_$MODE/examples # written in /lscratch by default
N_SHARD=`cat ${wd}/N_SHARD`

echo $REF $BAM $MODE $OUT $N_SHARD

mkdir -p $OUT logs

extra_args=""

if [[ $MODE == "WGS" ]]; then
  extra_args="--channels insert_size"
elif [[ $MODE == "PACBIO" ]]; then
  extra_args="--add_hp_channel true --alt_aligned_pileup diff_channels --max_reads_per_partition 600 --min_mapping_quality 1 --parse_sam_aux_fields True --partition_size 25000 --phase_reads true --pileup_image_width 199 --realign_reads false --sort_by_haplotypes true --track_ref_reads true --vsc_min_fraction_indels 0.12"
elif [[ $MODE == "ONT_R104" ]]; then
  extra_args="--add_hp_channel true --alt_aligned_pileup diff_channels --max_reads_per_partition 600 --min_mapping_quality 5 --parse_sam_aux_fields True --partition_size 25000 --phase_reads true --pileup_image_width 199 --realign_reads false --sort_by_haplotypes true --track_ref_reads true --vsc_min_fraction_indels 0.12 --vsc_min_fraction_snps 0.08"
elif [[ $MODE == "HYBRID_PACBIO_ILLUMINA" ]]; then
  extra_args=""
else
  echo "Unknown $MODE provided. Exit."
  exit -1
fi

echo "make_examples with parallel in $MODE mode"
echo "
make_examples \\
  --mode calling   \\
  --ref "${REF}"   \\
  --reads "${BAM}" \\
  --sample_name "$SAMPLE" \\ 
  --examples $OUT/tfrecord@${N_SHARD}.gz $extra_args \\
  --task {}
"

seq 0 $((N_SHARD-1)) \
  | parallel -j ${SLURM_CPUS_PER_TASK} --eta --halt 2 \
  --joblog "logs/log" --res "logs" \
  make_examples    \
    --mode calling \
    --ref "${REF}" \
    --reads "${BAM}" \
    --sample_name "${SAMPLE}" \
    --examples $OUT/tfrecord@${N_SHARD}.gz $extra_args \
    --task {} \
    || exit 1

mkdir -p $wd/dv_$MODE
cp -r $OUT $wd/dv_$MODE
cp -r logs/* $wd/logs
