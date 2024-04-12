#! /bin/bash

if [[ "$#" -lt 4 ]] ; then
  echo "Usage: ./_submit_deepvariant.sh <ref.fa> <alignment.bam> <mode> <sample> <wd> <jid_to_wait>"
  echo "  mode   WGS | PACBIO | ONT_R104 | HYBRID_PACBIO_ILLUMINA"
  echo "  output file will be generated as dv_mode.vcf.gz"
  exit -1
fi

REF=$(realpath $1)
BAM=$(realpath $2)
MODE=$3
SAMPLE=$4 # name to be appeared in the output VCF SAMPLE field
wd=$(realpath $5)
wait_for=$6 # optional

N_SHARD=16
PIPELINE=$tools/T2T-Ref

echo -e "$REF $BAM $MODE $SAMPLE $wd $wait_for"

# Check $MODE is valid
if [[ $MODE == "WGS" || $MODE == "PACBIO" || $MODE == "ONT_R104" || $MODE == "HYBRID_PACBIO_ILLUMINA" ]]; then
  echo "DeepVariant v1.6.0 in $MODE mode"
else
  echo "Unknown option $MODE"
fi

echo "Check faidx exists"
if [[ ! -e $REF.fai ]]; then
  echo "no .fai found. creating one..."
  samtools faidx $REF
  echo "done"
fi
echo

# Keep MODE and N_SHARD
echo $MODE > ${wd}/MODE
echo $N_SHARD > ${wd}/N_SHARD
echo $REF > ${wd}/REF
echo $BAM > ${wd}/BAM

mkdir -p logs

path=${wd}

if [ ! -f $wd/deepvariant.step1.done ]; then
echo "Step 1. make_examples"
cpus=$N_SHARD
mem=$(($cpus*3))g
gres="lscratch:500" # use /lscratch/ allow up to 500 GB
name=dv_step1
log=${wd}/logs/$name.%A.log
script=$PIPELINE/deepvariant/step1.sh
args="$wd $SAMPLE"
partition=norm
walltime=12:00:00

  if [[ ! -z $wait_for ]]; then
    extra="--dependency=afterok:$wait_for"
  fi
  set -x
  sbatch -J ${SAMPLE}_$name --cpus-per-task=$cpus --mem=$mem --gres=$gres \
    --partition=$partition -D $path \
    $extra --time=$walltime \
    --error=$log --output=$log $script $args > ${wd}/step1.jid
  set +x
fi 


if [ ! -f $wd/deepvariant.step2.done ]; then
  echo "Step 2. variant_call"
  cpus=12 # Can't set mem with gpus
  gres="lscratch:50,gpu:p100:1" # use /lscratch/ allow up to 50GB
  name=dv_step2
  log=${wd}/logs/$name.%A.log
  script=$PIPELINE/deepvariant/step2.sh
  walltime=4:00:00
  args="$wd"
  partition=gpu
  extra="--dependency=afterok:"`cat ${wd}/step1.jid`

  set -x
  sbatch -J ${SAMPLE}_$name --cpus-per-task=$cpus --gres=$gres \
    --partition=$partition \
        -D $path $extra --time=$walltime \
        --error=$log --output=$log $script $args > ${wd}/step2.jid
  set +x
fi

if [ ! -f $wd/deepvariant.step3.done ]; then
  echo "Step 3. postprocess_variants"
  cpus=8
  mem=80g
  name=dv_step3
  log=${wd}/logs/$name.%A.log
  script=$PIPELINE/deepvariant/step3.sh
  args="$SAMPLE $wd"
  partition=quick
  walltime=2:00:00
  extra="--dependency=afterok:"`cat ${wd}/step2.jid`
  set -x
  sbatch -J ${SAMPLE}_$name --cpus-per-task=$cpus --mem=$mem \
    --partition=$partition -D $path \
    $extra --time=$walltime \
    --error=$log --output=$log $script $args > ${wd}/step3.jid
  set +x
fi
