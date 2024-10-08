#!/bin/bash

if [[ "$#" -lt 1 ]]; then
  echo "Usage: ./cleanup.sh sample"
  echo
  echo "Checks for"
  echo "  deepvariant.step1.done"
  echo "  background.coverage_results.bed"
  echo "  cal.target.done"
  echo "and removes"
  echo "  sample.dedup.bam*"
  echo "  fastq/"
  exit -1
fi

sample=$1

echo "## Clean up bam and fastqs"
if [[ -f deepvariant.step3.done
   && -f background.coverage_results.bed
   && -f cal.target.done ]]; then
  set -x
  rm $sample.dedup.bam*
  rm -r fastq
  set +x
fi
