#!/bin/bash

if [[ "$#" -lt 1 ]]; then
  echo "Usage: ./find_45S.sh asm.fa"
  echo "Outputs:"
  echo "  45S_to_asm.mashmap.out mashmap3 output format"
  echo "  45S_to_asm.mashmap.bed reporting identity in score field as %"
  echo
  exit 0
fi

module load mashmap

asm_fa=$1
asm=`echo $asm_fa | sed 's/\.gz$//g' | sed 's/\.fasta$//g' | sed 's/\.fa$//g'`
cpu=$SLURM_CPUS_PER_TASK

mashmap \
    -t $cpu \
    --noSplit \
    -q $T2T_Ref/rDNA/human_45S.fa \
    -r $asm_fa \
    -s 13332 \
    --pi 85 \
    -f none \
    -o 45S_to_$asm.mashmap.out

cat 45S_to_$asm.mashmap.out |\
  awk -v OFS='\t' '{print $6, $8, $9, $1, $(NF-1), $5}' |\
  awk -F ":" '{print $1"\t"$3}' |\
  awk -v OFS='\t' '{print $1, $2,$3, "45S", (100*$6), $7}' \
   > 45S_to_$asm.mashmap.bed
