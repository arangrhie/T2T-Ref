#! /bin/bash

# ARGUMENTS
sampleInfo=$1
outPrefix=$2
filterOpt_step1=$3 # OPT
filterOpt_step2=$4 # OPT

echo -e "sampleInfo : $sampleInfo"
echo -e "outPrefix : $outPrefix"
echo -e ":filterOpt_step1 : $filterOpt_step1"
echo -e ":filterOpt_step2 : $filterOpt_step2"

# Files
PIPELINE=$tools/T2T-Ref

# RUN WRAPPER
# 01 FILE FIX
set -x
sh $PIPELINE/annovar/filter.deepvariant.headerFixer.sh $sampleInfo $filterOpt_step1
set +x

# 02 MERGE and FILTER
if [ ! -f merge.filter.done ]; then
	set -x
	sh $PIPELINE/annovar/merge.vcfs_GLnexus.sh $sampleInfo $outPrefix $filterOpt_step2
	set +x
fi

# 03 ANNOVAR
if [ ! -f variant_annotation.done ]; then
	set -x
	sh $PIPELINE/annovar/annotation.variants.hs1.sh $outPrefix
	set +x
fi
