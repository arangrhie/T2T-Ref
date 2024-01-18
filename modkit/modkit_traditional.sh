#!/bin/bash

if [[ "$#" -lt 2 ]]; then
  echo "Usage: modkit_traditional.sh bam ref.fa"
  echo "Run modkit with"
  echo " --preset traditional --force-allow-implicit --suppress-progress"
  echo "And extract 5mC wig and bw files"
fi

bam=$1
ref=$2

cpus=$SLURM_CPUS_PER_TASK

out=`echo $bam | sed 's/\.bam$//g'`
bed=$out.modkit.bed
wig=$out.modkit_5mC.wig
bw=$out.modkit_5mC.bw

mkdir -p logs
ml load modkit
ml load ucsc

set -o pipefail
set -x
modkit pileup $out.pri.bam $bed \
  --log-filepath logs/modkit_pileup.$out.log \
  --ref $ref -t $cpus \
  --preset traditional --force-allow-implicit --suppress-progress

awk '{print $1"\t"$2"\t"$3"\t"$11}' $bed |\
  java -jar -Xmx256m $tools/T2T-Ref/util/bedToWig.jar - "5mc" > $wig
wigToBigWig $wig $ref.fai $bw
set +x

