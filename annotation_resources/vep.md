# Variant Effect Predictor (VEP)

There are numerous plugins available for VEP, such as AlphaMissense, CADD, SpliceVault, and others. You can find them on the [VEP plugin homepage](https://plants.ensembl.org/info/docs/tools/vep/script/vep_plugins.html). Here, we introduce how to run VEP using a VCF and GFF file based on the CHM13 reference. You can download the reference genome and gene annotation in GFF format from the [CHM13 GitHub repository](https://github.com/marbl/CHM13?tab=readme-ov-file). The plugins are optional, so you can use others as needed; this is just an example.

## Generating and indexing the GFF file for VEP
```bash
ori_gff=chm13v2.0_RefSeq_Liftoff_v5.1.gff3
grep -v "#" $ori_gff | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > chm13v2.0_RefSeq_Liftoff_v5.1.gff_forVEP.gz
tabix -p gff chm13v2.0_RefSeq_Liftoff_v5.1.gff_forVEP.gz
```

## Running VEP
```bash
vcf=test.vcf.gz
ref=chm13v2.0.fa # reference used for aligning
gff=chm13v2.0_RefSeq_Liftoff_v5.2.gff.gz
VEP_CACHEDIR=./vepCash

mkdir -p $VEP_CACHEDIR

vep --fork 50 \
 -i ${vcf} \
 --gff ${gff} \
 -o ${vcf}.vep \
 --force_overwrite \
 --dir ${VEP_CACHEDIR} \
 --fasta ${ref} \
 --everything
```

