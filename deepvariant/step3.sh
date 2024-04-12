#! /bin/bash

if [[ "$#" -lt 1 ]]; then
  echo "./step3.sh sample"
  exit -1
fi

set -o pipefail
set -e
set -x 

module load deepvariant/1.6.0

SAMPLE=$1
wd=$2
REF=`cat ${wd}/REF`
MODE=`cat ${wd}/MODE`
OUT=dv_$MODE
CALL_VARIANTS_OUTPUT="${wd}/$OUT/call_variants_output.tfrecord.gz"

postprocess_variants \
  --ref "${REF}" \
  --infile "${CALL_VARIANTS_OUTPUT}" \
  --outfile "${wd}/$OUT.$SAMPLE.vcf.gz" \
  --sample_name $SAMPLE  && touch $wd/deepvariant.step3.done

if [ -f $wd/deepvariant.step3.done ]; then
	rm -rf ${wd}/$OUT
fi
