#!/bin/bash 
set -e
set -o pipefail

# ARGUMENTS
sampleInfo=$1 # sample_info.txt #1st col=array,  2nd col = sample name , 3rd = bam file
array=$2 # "--array=1-30"

# REFERENCE
#refwY='/data/Phillippy2/projects/rpc/00.rawdata/00.references/chm13v2.rCRS.DJ_PHR_rDNA_masked/masked_DJ_rDNA_PHR.autosome.rCRS.forXY.fa' # soft mask version of chrY
refwY="/data/Phillippy2/projects/rpc/00.rawdata/00.references/chm13v2.rCRS.DJ_PHR_rDNA_masked/masked_DJ_rDNA_PHR.autosome.rCRS.forXY.hardMaskY.fa" # hard mask verion of chrY
refwoY='/data/Phillippy2/projects/rpc/00.rawdata/00.references/chm13v2.rCRS.DJ_PHR_rDNA_masked/masked_DJ_rDNA_PHR.autosome.rCRS.forXX.fa'
pipeline='/data/Phillippy2/projects/rpc/99.codes'

# RUNNING SCRIPTS 
# 01. reference indexing for bwa
sh ~/code/_submit_norm.sh 10 20G bwa_indexing $pipeline/preprocessing/indexing.bwa.sh "$refwY $refwoY"

# 02. running main code
sh ~/code/_submit_norm.sh 1 5G srWGS_wrapper ${pipeline}/variants_sr/srWGS.extract.bwa.deepvariant.arrayJob.sh "$sampleInfo $refwY $refwoY" $array
