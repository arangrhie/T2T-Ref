#! /bin/bash

if [[ $# -lt 2 ]]; then
  echo "Usage: ./annovar.sh in.vcf out_prefix"
  exit 0
fi

# ARGUMENTS
inVCF=$1
outPrefix=$2


thread=12
# annovarDB=$tools/T2T-Ref/annovar/db
annovarDB=/data/Phillippy/tools/T2T-Ref/annotation/db
build='hs1'

# MODULE 
module load annovar

# RUN MAINCODE
convert2annovar.pl -format vcf4 -allsample -withfreq $inVCF > ${outPrefix}.inp

table_annovar.pl \
--thread $thread \
${outPrefix}.inp \
$annovarDB \
-buildver $build \
-out ${outPrefix} \
-remove \
-protocol refGene,curGenev5.2,snp156,1000g2023dec_all,1000g2023dec_afr,1000g2023dec_amr,1000g2023dec_eas,1000g2023dec_eur,1000g2023dec_sas,clinvar_20231217,gwas_20231207,nonSyntenic,hg38_issues,feat,cenSat,sraccess,sraccess_hg38,sraccess_hs1Only \
-operation g,g,f,f,f,f,f,f,f,f,r,r,r,r,r,r,r,r \
-arg ',,,,,,,,,,,,,,,,,' && touch annotate_variants.done
