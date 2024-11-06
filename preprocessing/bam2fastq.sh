#!/bin/bash 
set -e
set -o pipefail

module load samtools

# ARGUMENTS
inBam=$1
name=$2
cpu=$SLURM_CPUS_PER_TASK
outDir=fastq

if [[ "$#" -lt 2 ]]; then
  echo "Usage: bam2fastq.sh in.bam out_name"
  echo "  in.bam    input bam file"
  echo "  out_name  output prefix"
  exit -1
fi

if [[ -z $SLURM_JOB_ID ]]; then
  echo "Assign lscratch"
  exit -1
fi
tmp=/lscratch/${SLURM_JOB_ID}

if [ ! -f bam2fastq.$name.done ]; then
  set -x
  # sort by read name
  samtools sort -@ $cpu -n -o $tmp/${name}_nameSort.bam ${inBam}

  # bam to fastq
  samtools fastq -@ $cpu $tmp/${name}_nameSort.bam \
      -1 $tmp/${name}_1.fq \
      -2 $tmp/${name}_2.fq && rm $tmp/${name}_nameSort.bam || exit -1

  pigz $tmp/${name}_1.fq
  # mv $tmp/${name}_1.fq.gz ${outDir}/${name}_1.fq.gz

  pigz $tmp/${name}_2.fq
  # mv $tmp/${name}_2.fq.gz ${outDir}/${name}_2.fq.gz

  # fastq_map.fofn
  # echo -e "${PWD}/${outDir}/${name}_1.fq.gz\t${PWD}/${outDir}/${name}_2.fq.gz" > fastq_map.fofn
  echo -e "$tmp/${name}_1.fq.gz\t$tmp/${name}_2.fq.gz" > fastq_map.fofn
  touch bam2fastq.$name.done
  set +x
else
  echo "Found bam2fastq.$name.done. Nothing to do."
fi
