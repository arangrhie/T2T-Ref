#! /bin/bash

if [[ "$#" -lt 4 ]]; then
  echo "Usage: ./_submit_bwa_dv.sh ref.fasta fastq.map out sample [line_num]"
  echo "  ref.fasta  reference"
  echo "  fastq.map  space or tab teliminted for r1 and r2 in each line"
  echo "             example: /path/to/r1.fq.gz /path/to/r2.fq.gz"
  echo "  out        output prefix. used in bam"
  echo "  sample     sample name to appear in the final VCF"
  echo "  line_num   line_num to launch the mapping OPTIONAL"
  echo
  echo "  output VCF file will be generated as dv_mode.vcf.gz"
  echo
  echo "Not enough parameters given. Exit."
  exit -1
fi

ref=$1
fastq_map=$2
out=$3
sample=$4
wd=$5
sexChr=$6
mode="WGS" # no need to customize this, mode is always WGS

PIPELINE=/data/Phillippy2/projects/rpc/99.codes

if [ ! -f "${wd}/depthCal.done" ] ; then
	echo -e "BWA is starting"
	cmd="sh $PIPELINE/bwa/_submit_bwa.sh $ref $fastq_map $sample $wd $sexChr"
	echo -e $cmd
	eval $cmd	
	wait_for=`tail -n1 ${wd}/bwa.jid`
fi

if [ ! -f ${wd}/deepvariant.done ] ; then
	echo -e "Deepvarint is starting"
	cmd="sh $PIPELINE/deepvariant/_submit_deepvariant.sh $ref ${out}.dedup.bam $mode $sample ${wd} $(tail -n1 ${wd}/bwa.jid)"
	echo -e $cmd
	eval $cmd
fi

