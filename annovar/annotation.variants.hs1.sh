#! /bin/bash

# ARGUMENTS
outPrefix=$1

thread=10
annovarDB='/data/kimj75/00.Files/references/chm13v2/annovar/filterDB/'
build='hs1'

# MODULE 
module load annovar

# RUN MAINCODE
convert2annovar.pl -format vcf4 ${outPrefix}.vcf > ${outPrefix}.inp

table_annovar.pl \
--thread $thread \
${outPrefix}.inp \
$annovarDB \
-buildver $build \
-out ${outPrefix} \
-remove \
-protocol refGene,curGene,dbsnp156,1000g2023dec_all,1000g2023dec_afr,1000g2023dec_amr,1000g2023dec_eas,1000g2023dec_eur,1000g2023dec_sas,clinvar_20231217,gwas_20231207,nonSyntenic,hg38_issues,feat,cenSat,sraccess,sraccess_hg38,sraccess_hs1Only \
-operation g,g,f,f,f,f,f,f,f,f,r,r,r,r,r,r,r,r \
-arg ',,,,,,,,,,,,,,,,,' && touch variant_annotation.done
