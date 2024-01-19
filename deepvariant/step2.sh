#! /bin/bash

if [[ "$#" -lt 1 ]]; then
  echo "Usage: ./step2.sh"
  exit -1
fi

module load deepvariant/1.5.0 || exit 1

set -o pipefail
set -e
wd=$1
cd ${wd}


MODE=`cat ${wd}/MODE`
OUT=dv_$MODE
mkdir -p $OUT
mode=`echo $MODE | awk '{print tolower($0)}'`
MODEL="/opt/models/$mode/model.ckpt"
N_SHARD=`cat ${wd}/N_SHARD` # Must match as in step1
CALL_VARIANTS_OUTPUT="${wd}/$OUT/call_variants_output.tfrecord.gz"

set -x

call_variants \
  --outfile "${CALL_VARIANTS_OUTPUT}" \
  --examples "${wd}/$OUT/examples/tfrecord@${N_SHARD}.gz" \
  --checkpoint "${MODEL}" \
  || exit 1
