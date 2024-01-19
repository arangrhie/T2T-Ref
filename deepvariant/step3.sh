#! /bin/bash

if [[ "$#" -lt 1 ]]; then
  echo "./step3.sh sample"
  exit -1
fi

set -o pipefail
set -e

module load deepvariant/1.5.0

SAMPLE=$1
wd=$2
REF=`cat ${wd}/REF`
MODE=`cat ${wd}/MODE`
OUT=dv_$MODE
CALL_VARIANTS_OUTPUT="${wd}/$OUT/call_variants_output.tfrecord.gz"

postprocess_variants \
  --ref "${REF}" \
  --infile "${CALL_VARIANTS_OUTPUT}" \
  --outfile "${wd}/$OUT.vcf.gz" \
  --sample_name $SAMPLE  && touch deepvariant.done
