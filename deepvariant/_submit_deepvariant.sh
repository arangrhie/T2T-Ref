#! /bin/bash

if [[ "$#" -lt 5 ]] ; then
  echo "Usage: ./_submit_deepvariant.sh <ref.fa> <alignment.bam> <mode> <sample> <sex> [par] [jid_to_wait]"
  echo "  mode    WGS | PACBIO | ONT_R104 | HYBRID_PACBIO_ILLUMINA"
  echo "  sample  unique identifier for each sample"
  echo "  sex     XX or XY. For XY individuals, X and Y will re-adjust genotypes for hets, except for the PAR"
  echo "  par     PAR bed file for XY individuals. DEFAULT=\$tools/T2T-Ref/ref/chm13v2.0_PAR.bed"
  echo "  output file will be generated as dv_mode.vcf.gz and dv_mode.gvcf.gz"
  exit -1
fi

REF=$(realpath $1)
BAM=$(realpath $2)
MODE=$3
SAMPLE=$4   # name to be appeared in the output VCF SAMPLE field
SEX=$5      # XY? For step3
PAR=$6      # PAR bed file
wait_for=$7 # optional

if [ -f  N_SHARD  ] ; then 
	N_SHARD=$(cat N_SHARD)
else
	N_SHARD=12
fi

PIPELINE=$tools/T2T-Ref
if [[ -z $PAR ]]; then
  PAR=$PIPELINE/ref/chm13v2.0_PAR.bed
fi


echo -e "$REF $BAM $MODE $SAMPLE $SEX $PAR $wait_for"

# Check $MODE is valid
if [[ $MODE == "WGS" || $MODE == "PACBIO" || $MODE == "ONT_R104" || $MODE == "HYBRID_PACBIO_ILLUMINA" ]]; then
  echo "DeepVariant v1.6.1 in $MODE mode"
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
echo $MODE > MODE
echo $N_SHARD > N_SHARD
echo $REF > REF
echo $BAM > BAM

mkdir -p logs

path=$PWD

extra=""

if [ ! -f deepvariant.step1.done ]; then
echo "Step 1. make_examples"
cpus=$N_SHARD
mem=$(($cpus*3))g
gres="lscratch:500" # use /lscratch/ allow up to 500 GB
name=dv_step1
log=logs/$name.%A.log
script=$PIPELINE/deepvariant/step1.sh
args="$PWD $SAMPLE"
partition=norm
walltime=12:00:00

  if [[ ! -z $wait_for ]]; then
    extra="--dependency=afterok:$wait_for"
  fi
  set -x
  sbatch -J ${SAMPLE}_$name --cpus-per-task=$cpus --mem=$mem --gres=$gres \
    --partition=$partition -D $path \
    $extra --time=$walltime \
    --error=$log --output=$log $script $args > step1.jid
  set +x
  extra="--dependency=afterok:"`cat step1.jid`
fi 


if [ ! -f deepvariant.step2.done ]; then
  echo "Step 2. variant_call"
  cpus=16 # Can't set mem with gpus?
  mem=10g
  gres="gpu:1" # p100 takes ~1 hr. k80 ~5 hrs. No more lscratch required. before we used lscratch:50 (50g)
  name=dv_step2
  log=logs/$name.%A.log
  script=$PIPELINE/deepvariant/step2.sh
  walltime=8:00:00 # for k80; took ~5 hrs
  args=""
  partition=gpu

  set -x
  sbatch -J ${SAMPLE}_$name --cpus-per-task=$cpus --mem=$mem --gres=$gres \
    --partition=$partition \
    -D $path $extra --time=$walltime \
    --error=$log --output=$log $script $args > step2.jid
  set +x
  extra="--dependency=afterok:"`cat step2.jid`
fi

if [ ! -f deepvariant.step3.done ]; then
  echo "Step 3. postprocess_variants"
  cpus=12
  mem=72g # 1 sample got killed at 60 G
  name=dv_step3
  log=logs/$name.%A.log
  script=$PIPELINE/deepvariant/step3.sh
  args="$SAMPLE $SEX $PAR"
  partition=quick
  walltime=1:20:00 # 1 sample timed out with 50 min
  
  set -x
  sbatch -J ${SAMPLE}_$name --cpus-per-task=$cpus --mem=$mem \
    --partition=$partition -D $path \
    $extra --time=$walltime \
    --error=$log --output=$log $script $args > step3.jid
  set +x
else
  echo "Nothing to re-run for $SAMPLE. DONE."
fi

echo