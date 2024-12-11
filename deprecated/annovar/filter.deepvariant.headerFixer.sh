#!/bin/bash

# ARGUMENTS
sampleInfo=$1
filterOpt=$2

# VARIANT FILTER
if [ "$#" -lt 2 ]; then
	filterOpt='QUAL>50' # default
fi
echo "filterOpt : $filterOpt"


# MODULES
module load samtools # for bgzip
module load bcftools

# MAIN CODES
for sample in `cut -f 2 $sampleInfo`;
do sex=`cut -f 3 ${sample}/sex.determine.txt`

if [ ! -f ${sample}/vcfFilter.done ]; then
	if [[ -f ${sample}/dv_WGS.vcf.gz ]]; then
		if [[ $sex == "F" ]]; then
			echo -e "$sample is Female"
			gunzip -c ${sample}/dv_WGS.vcf.gz | sed '/contig=<ID=chrX/a ##contig=<ID=chrY,length=62460029>' | bgzip > ${sample}/${sample}.dv_WGS.headerFix.vcf.gz
			bcftools filter -O z -o ${sample}/${sample}.dv_WGS.filtered.g.vcf.gz -i $filterOpt ${sample}/${sample}.dv_WGS.headerFix.vcf.gz && touch ${sample}/vcfFilter.done
			if [ -f ${sample}/vcfFilter.done ]; then
				echo "rm tmp file"
				#	rm ${sample}/${sample}.dv_WGS.headerFix.vcf.gz
			fi
		else
			echo -e "$sample is Male"
			bcftools filter -O z -o ${sample}/${sample}.dv_WGS.filtered.g.vcf.gz -i $filterOpt ${sample}/dv_WGS.vcf.gz && touch ${sample}/vcfFilter.done
		fi

	fi
fi
done
