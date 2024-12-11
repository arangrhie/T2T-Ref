#! /bin/bash

# ARGUMENS
sampleInfo=$1
outPrefix=$2
filterOpt_step1=$3 # OPT

if [ "$#" -lt 3 ]; then
       filterOpt_step1='INFO/AF > 0.1'
fi
rm -rf GLnexus.DB/

# MODULES
module load glnexus
module load samtools

# RUN CODES
VCFLIST=$(for sample in `cut -f 2 $sampleInfo`; do echo $sample/$sample.dv_WGS.filtered.g.vcf.gz | tr '\n' ' '  ;  done)

glnexus_cli --config DeepVariant $VCFLIST | bcftools view - > ${outPrefix}.vcf
bcftools filter -O z -o ${outPrefix}.filtered.vcf -i  ${outPrefix}.vcf && touch merge.filter.done

if [ -f merge.filter.done ] ; then 
	rm ${outPrefix}.vcf
fi
