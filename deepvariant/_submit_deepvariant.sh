#! /bin/bash

if [[ "$#" -lt 5 ]] ; then
  echo "Usage: ./_submit_deepvariant.sh <ref.fa> <alignment.bam> <mode> <sample> <wd> <jid_to_wait>"
  echo "  mode   WGS | PACBIO | ONT_R104 | HYBRID_PACBIO_ILLUMINA"
  echo "  output file will be generated as dv_mode.vcf.gz"
  exit -1
fi

REF=$1
BAM=$2
MODE=$3
SAMPLE=$4 # name to be appeared in the output VCF SAMPLE field
wd=$5
wait_for=$6 # optional

N_SHARD=48
PIPELINE=/data/Phillippy2/projects/rpc/99.codes

echo -e "$REF $BAM $MODE $SAMPLE $wd $wait_for"

# Check $MODE is valid
if [[ $MODE == "WGS" || $MODE == "PACBIO" || $MODE == "ONT_R104" || $MODE == "HYBRID_PACBIO_ILLUMINA" ]]; then
  echo "DeepVariant v1.5 in $MODE mode"
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

echo "Step 1. make_examples"
cpus=$N_SHARD
mem=$(($cpus*2))g
gres="lscratch:500" # use /lscratch/ allow up to 500 GB
name=dv_step1
script=$PIPELINE/deepvariant/step1.sh
args="$wd $SAMPLE"
partition=norm
walltime=24:00:00
path=${wd}
log=${wd}/logs/$name.%A.log

if [[ ! -z $wait_for ]]; then
	extra="--dependency=afterok:$wait_for"
fi

echo "
sbatch -J $name --cpus-per-task=$cpus --mem=$mem --gres=$gres \\
       --partition=$partition \\
       -D $path \\
       $extra --time=$walltime \\
       --error=$log --output=$log $script $args
"

sbatch -J $name --cpus-per-task=$cpus --mem=$mem --gres=$gres  \
       --partition=$partition \
       -D $path \
       $extra --time=$walltime \
       --error=$log --output=$log $script $args > ${wd}/step1.jid

echo "Step 2. variant_call"

cpus=12
gres="lscratch:50,gpu:p100:1" # use /lscratch/ allow up to 50GB
name=dv_step2
script=$PIPELINE/deepvariant/step2.sh
args="$wd"
partition=gpu
extra="--dependency=afterok:"`cat ${wd}/step1.jid`
log=${wd}/logs/$name.%A.log

echo "
sbatch -J $name --cpus-per-task=$cpus --gres=$gres \\
       --partition=$partition \\
       -D $path \\
       $extra --time=$walltime \\
       --error=$log --output=$log $script $args > step2.jid
"

sbatch -J $name --cpus-per-task=$cpus --gres=$gres  \
       --partition=$partition \
       -D $path \
       $extra --time=$walltime \
       --error=$log --output=$log $script $args > ${wd}/step2.jid

echo "Step 3. postprocess_variants"
cpus=8
mem=48g
name=dv_step3
script=$PIPELINE/deepvariant/step3.sh
args="$SAMPLE $wd"
partition=norm
walltime=12:00:00
extra="--dependency=afterok:"`cat ${wd}/step2.jid`
log=${wd}/logs/$name.%A.log

echo "
sbatch -J $name --cpus-per-task=$cpus --mem=$mem \\
       --partition=$partition \\
       -D $path \\
       $extra --time=$walltime \\
       --error=$log --output=$log $script $args > ${wd}/step3.jid
"
sbatch -J $name --cpus-per-task=$cpus --mem=$mem  \
       --partition=$partition \
       -D $path \
       $extra --time=$walltime \
       --error=$log --output=$log $script $args > ${wd}/step3.jid

