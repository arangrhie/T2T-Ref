#!/bin/bash

if [[ "$#" -lt 1 ]]; then
  echo "Usage: ./cleanup.sh sample"
  echo
  echo "Checks for"
  echo "  deepvariant.step3.done"
  echo "and removes"
  echo "  fastq/ and bam2fastq.sample.done"
  exit -1
fi

sample=$1

echo "## Clean up fastqs"
if [[ -f deepvariant.step3.done
   && -f background.coverage_results.bed
   && -f cal.target.done ]]; then
  set -x
  rm -r fastq
  rm bam2fastq.$sample.done
  set +x
fi
