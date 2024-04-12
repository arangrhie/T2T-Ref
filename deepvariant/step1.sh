#!/bin/bash
# this is deepvariant_step1.sh

module load deepvariant/1.6.0 || exit 1
module load parallel


set -x 

wd=$(realpath $1)
SAMPLE=$2

echo $wd

[[ -d /lscratch/${SLURM_JOB_ID} ]] && cd /lscratch/${SLURM_JOB_ID} || exit 1

echo $PWD

REF=`cat ${wd}/REF`
BAM=`cat ${wd}/BAM`
MODE=`cat ${wd}/MODE`
OUT=dv_$MODE/examples # written in /lscratch by default
N_SHARDS=`cat ${wd}/N_SHARD`

mkdir -p $OUT logs $wd/logs-parallel-$SLURM_JOB_ID

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

seq 0 $((N_SHARDS-1)) \
    | parallel -P ${N_SHARDS} --halt 2 \
        --joblog "$wd/logs-parallel-$SLURM_JOB_ID/log" --res "$wd/logs-parallel-$SLURM_JOB_ID" \
      make_examples --mode calling \
        --ref "${REF}" \
        --reads "${BAM}" \
        --examples $OUT/examples.tfrecord@${N_SHARDS}.gz \
        --channels insert_size \
        --task {} \
|| exit 1 && touch $wd/deepvariant.step1.done &&


mkdir -p $wd/dv_$MODE &&
cp -r $OUT $wd/dv_$MODE && 
rm -rf "$wd/logs-parallel-$SLURM_JOB_ID"
