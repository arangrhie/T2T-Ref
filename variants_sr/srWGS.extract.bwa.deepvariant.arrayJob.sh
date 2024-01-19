#!/bin/bash 
set -e
set -o pipefail

# ARGUMENTS
sampleInfo=$1
refwY=$2 #/data/Phillippy2/projects/rpc/00.rawdata/00.references/chm13v2.rCRS.DJ_PHR_rDNA_masked/masked_DJ_rDNA_PHR.autosome.rCRS.forXY.fa
refwoY=$3 #/data/Phillippy2/projects/rpc/00.rawdata/00.references/chm13v2.rCRS.DJ_PHR_rDNA_masked/masked_DJ_rDNA_PHR.autosome.rCRS.forXX.fa

pipeline='/data/Phillippy2/projects/rpc/99.codes'
wd=$PWD
lineNum=$SLURM_ARRAY_TASK_ID
sampleName=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $sampleInfo)
bam=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $sampleInfo)
wd=`echo $PWD/${sampleName}`

# SET ENV
mkdir -p $PWD/${sampleName}
mkdir -p $PWD/${sampleName}/logs

# change directory
cd $wd
echo $wd

# Sex check
if [ ! -f $wd/sexCheck.done ]; then
	echo -e "Check sex"
	echo -e "sh $pipeline/preprocessing/check_sex_samtools.sh $bam $sampleName 1>> sex.determine.txt 2>> sex.determine.warnings.log && touch $wd/sexCheck.done"
	sh $pipeline/preprocessing/check_sex_samtools.sh $bam $sampleName && touch $wd/sexCheck.done
else
	echo -e "Sex check was already done"
fi

sex=`cut -f 3 sex.determine.txt`

# bam 2 fastq
mkdir -p ${wd}/fastq
if [ ! -f ${wd}/fastq/bam2fastq.done ] ; then
	echo -e "Generating fastqmap\n"
	echo -e "sh $pipeline/preprocessing/bam2fastq.sh $bam ${wd}/fastq/${sampleName}"
	sh $pipeline/preprocessing/bam2fastq.sh $bam ${wd}/fastq/${sampleName}
	echo -e "bam2fastq was done"
else 
	echo -e "bam2fastq was already done"
fi

# run main script bwa -> deepvariant
if [ "$sex" == "F" ] ; then
	ref=$refwoY
	echo "$sampleName $bam $sex $ref"
	cmd="sh ${pipeline}/variants_sr/_submit_bwa_dv.sh ${ref} ${wd}/fastq_map $wd/${sampleName} ${sampleName} ${wd} \"XX\""
	echo -e $cmd && eval $cmd

elif [ "$sex" == "M" ] ; then
	ref=$refwY
	echo "$sampleName $bam $sex $ref"
	echo -e "sh ${pipeline}/variants_sr/_submit_bwa_dv.sh ${ref} $wd/fastq_map ${wd}/${sampleName} ${sampleName} ${wd} XY"
	sh ${pipeline}/variants_sr/_submit_bwa_dv.sh ${ref} ${wd}/fastq_map ${wd}/${sampleName} ${sampleName} ${wd} "XY"
fi

# Calculate Coverage

echo "This is array task ${SLURM_ARRAY_TASK_ID}, the sample name is ${sampleName},the sex is ${sex} and ${ref} for reference" >> array.output.txt
