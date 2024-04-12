#! /bin/bash
set -x
if [[ "$#" -lt 1 ]]; then
  echo "Usage: ./step2.sh"
  exit -1
fi

module load deepvariant/1.6.0 || exit 1

set -o pipefail
set -e
wd=`realpath $1`
cd ${wd}


MODE=`cat ${wd}/MODE`
OUT=dv_$MODE/examples
mkdir -p $OUT
mode=`echo $MODE | awk '{print tolower($0)}'`
MODEL="/opt/models/$mode"
N_SHARD=`cat ${wd}/N_SHARD` # Must match as in step1
CALL_VARIANTS_OUTPUT="${wd}/dv_$MODE/call_variants_output.tfrecord.gz"

set -x

call_variants \
  --outfile "${CALL_VARIANTS_OUTPUT}" \
  --examples $wd/$OUT/examples.tfrecord@${N_SHARD}.gz \
  --checkpoint "${MODEL}" || exit 1 && touch $wd/deepvariant.step2.done &&
  rm -rf ${wd}/$OUT
