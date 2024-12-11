#! /bin/bash

vcf=$1
ref=$tools/T2T-Ref/ref/chm13v2.0_masked_DJ_5S_rDNA_PHR_PAR_wi_rCRS.fa
gtf=$tools/T2T-Ref/annotation/db/chm13v2.0_RefSeq_Liftoff_v5.2_forVEP.gff3.gz

module load VEP/111

vep --fork 50 \
 -i $vcf \
 --gff ${gtf} \
 -o $vcf.vep \
 --force_overwrite \
 --dir ${VEP_CACHEDIR} \
 --fasta $ref \
 --everything \
 --plugin CSN \
 --plugin Blosum62 \
 --plugin Carol \
 --plugin Condel,$VEP_CACHEDIR/Plugins/config/Condel/config,b \
 --plugin GeneSplicer,$GS/bin/genesplicer,$GS/human,context=200 \
 --plugin Downstream \
 --plugin LoFtool
