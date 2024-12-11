# ANNOVAR databases for T2T-CHM13v2.0

😊 Welcome! 😊

This repository describes the [ANNOVAR](https://annovar.openbioinformatics.org/) databases 💻 for hs1 (T2T-CHM13v2.0)🧬 and how to use it.

Feel free to leave questions or problems while using this database in the [Issues](T2T-Ref/issues).

## Usage

### 1. Download pre-built databases

Pre-built databases and an example file are available to download from [here](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/annotation/annovar/).

Original resources we used for generating each database and operational details are described in this spreadsheet [T2T-CHM13 Databases and operation](https://docs.google.com/spreadsheets/d/1sgjmGLLbXAZpyNiSUbxiEa1hJEVDxuOqL1yAmpDV5BA/edit?usp=sharing).

The directory should be structured such as:
```bash
- db
  - hs1_${Database_name_1}.txt
  - hs1_${Database_name_2}.txt
  - hs1_${Database_name_3}.txt
  ...
```

After downloading each database, unzip into a `.txt` file:
```bash
bgzip -d filename.txt.gz
```

If you intend to annotate variants using the gene annotation databases such as `hs1_curGene.txt.gz` or `hs1_refGene.txt.gz`, please download `${Database_name}Mrna.fa.gz` file as well.

There are index file(`.txt.idx`) for filter-based databases, such as allele frequency(`hs1_ALL.sites.2023_12.txt`) and dbSNP(`hs1_snp156.txt`).  It is not a mandatory file, but it would be helpful to reduce computing time during variant annotation.

### 2. Make your input compatible with annovar

Make your input (target VCF) file compatible with ANNOVAR. Use the `convert2annovar.pl` script in ANNOVAR to convert your VCF file acceptable to ANNOVAR:

```bash
perl convert2annovar.pl -format vcf4 ${prefix}.vcf > ${prefix}.avinput
```
Adjust parameters based on your specific requirements. More details regarding the input format can be found [here](https://annovar.openbioinformatics.org/en/latest/user-guide/input/).

An example file `example.chr22.inp` used in this description is available [here](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/annotation/annovar/). This file has been converted from VCF to Annovar-compatible format and contains a small subset of variants from Chromosome 22.

### 3. Annotate Variants
Annotate your variants with the following command:

```bash
annovarDB="/path/to/T2T-Ref/annotation/db" # Folder containing all hs1_*.txt files
build="hs1"
outPrefix="OUT" # the prefix you want. The ouput name would be ${OUT}.hs1_multianno.txt

table_annovar.pl \
  --thread 10 \ 
  example.chr22.inp \   # Input file from step 2
  $annovarDB \          # ANNOVAR database direcoty
  -buildver $build \    # hs1 (database prefix)
  -out ${outPrefix} \   # output prefix (example.chr22.out for the example output)
  -remove \
  -protocol refGene,curGene5.2,dbsnp156,1000g2023dec_all,1000g2023dec_afr,1000g2023dec_amr,1000g2023dec_eas,1000g2023dec_eur,1000g2023dec_sas,clinvar_20231217,gwas_20231207,nonSyntenic,hg38_issues,feat,cenSat,sraccess,sraccess_hg38,sraccess_hs1Only \
  -operation g,g,f,f,f,f,f,f,f,f,r,r,r,r,r,r,r,r \
  -arg ',,,,,,,,,,,,,,,,,'
```

Keep in mind that operations such as `g`, `f`, or `r` needs to match each database. You can find the correct operation for each annotation source on [T2T-CHM13 Databases and operation](https://docs.google.com/spreadsheets/d/1sgjmGLLbXAZpyNiSUbxiEa1hJEVDxuOqL1yAmpDV5BA/edit?usp=sharing).

## Advanced: [Building custom ANNOVAR databases](annovar_db.md)

### Written by Juhyun Kim