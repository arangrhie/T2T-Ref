# How to generate ANNOVAR databases

Each ANNOVAR database (file) is paired with a link to its' original source file below.

If the original file was based on GRCh38 or references other than CHM13v2.0, the files were lifted over to CHM13 using crossmap. The chain file used can be downloaded from [here](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/grch38-chm13v2.chain).

## Genome-based DB
* `hs1_refGene.txt`: [ANNOVAR homepage](http://www.openbioinformatics.org/annovar/download/hs1_refGene.txt.gz)
* `hs1_curGene5.2.txt`: [CHM13 github](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz) - This contains curated annotations of the ampliconic genes on the Y chromosome and other protein coding genes, correcting annotation errors in GENCODEv35 CAT/Liftoff and RefSeqv110 annotation.

* If the original file was formatted in GFF, I transformed it to GTF and then used [gtfToGenePred](https://bioconda.github.io/recipes/ucsc-gtftogenepred/README.html) to convert it into GenePred format.
* The gene annotation database in ANNOVAR website used to come with ${prefix}Mrna.fa. This file was generated using [`retrieve_seq_from_fasta.pl`](https://github.com/ronammar/Awesomeomics/raw/master/data/annovar_annotations/annovar/retrieve_seq_from_fasta.pl) script:

  ```shell
  retrieve_seq_from_fasta.pl --format refGene --seqfile chm13v2.0.fa hs1_refGene.txt --out hs1_refGeneMrna.fa
  ```

  * format(hs1_curGene.txt):
    ```
    1  XR_002958507.2  chr1    -       6046    13941   13941   13941   4       6046,12077,13444,13679, 6420,12982,13579,13941,0       LOC124900618    none    none    -1,-1,-1,>
    2  XR_007068557.1_1        chr1    +       15079   21429   21429   21429   2       15079,20565,    15564,21429,    0      LOC124905335_1  none    none    -1,-1,
    3  XM_047436352.1  chr1    -       20528   37628   20949   37628   5       20528,28446,34957,36085,37442,  21087,2862635059,37081,37628,  0       LOC112268260    incmpl  i>
    4  NR_125957.1     chr1    -       52978   54612   54612   54612   3       52978,53559,54521,      53422,53826,54612,     0       LOC101928626    none    none    -1,-1,-1,
    5  NM_001005221.2  chr1    -       111939  112877  111939  112877  1       111939, 112877, 0       OR4F29  incmpl  incmpl 0,
    6  XR_001743907.1_1        chr1    -       115045  117120  117120  117120  2       115045,116798,  115300,117120,  0       LOC107986552_1  none    none    -1,-1,
    ```

## Filter-based DB
* `hs1_dbsnp156.txt`: [NCBI DBsnp ftp server](https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz)
* `hs1_clinvar_20231217.txt`: [CHM13 github](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/liftover/chm13v2.0_ClinVar20220313.vcf.gz)
* `hs1_${population}.sites.2023_12.txt`: [CHM13 github](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/) - 1KGP allele frequency recalled on T2T-CHM13v2.0. Now available for all chromosomes, for the entire 3,202 samples or the unrelated 2504 samples. (popultation : ALL, AFR, AMR, EAS, EUR, and SAS)

  * format (hs1_ALL.sites.2023_12.txt) :
    ```
    chr1    131     A       C       0.00278552
    chr1    131     A       T       0.00278552
    chr1    878     A       AACCCTAACCCTCACCCTCACCCTCACCCTCACCCTCACCCTCACCCTC       0.00019976
    chr1    884     AACCCTCACCCTCACCCTCACCCTCACCCTCACCCTCACCCTCACCCTCACCCTCACCCTCACCCTCACCCTCACCCTCACCCTCACCCTCACCCTCACCCTCACCCTCACCCTCACCCTCACCCTC A       0.00020136
    ```
  
## Region-based DB
* `hs1_gwas_20231207.txt`: [GWAS Catalog](https://www.ebi.ac.uk/gwas/docs/file-downloads) v1.0-associations_e110_r2023-12-07
* `hs1_cenSat.txt`: [CHM13 github](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_censat_v2.1.bed) - A more comprehensive centromere/satellite repeat annotation
* `hs1_nonSyntenic.txt`: [CHM13 github](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-unique_to_hg38.bed) - Regions non-syntenic (unique) compared to GRCh38
* `hs1_hg38_issues.txt`: [CHM13v1 publication](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/publications/Nurk_2021/fig1/) - Issues in GRCh38 reported by GRC in T2T-CHM13 coordinates
* `hs1_sraccess.txt`: [CHM13 github](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/accessibility/combined_mask.bed.gz) - Regions reliably acessible with short reads in CHM13v2.0
* `hs1_sraccess_hg38.txt`: [CHM13 amazon download server](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/annotation/accessibility/) Regions reliably accessible with short reads on GRCh38, liftover to CHM13v2.0
* `hs1_sraccess_hs1Only.txt`: Regions reliably accessible with short reads only in CHM13, not in GRCh38

  * format (hs1_cenSat.txt) :
    ```
    1       chr1    116796047       121405145       ct_1_1(p_arm)   100     .       116796047       121405145       224,224,224
    1       chr1    121405145       121406286       censat_1_1(rnd-6_family-4384)   100     .       121405145       121406286       0,204,204
    1       chr1    121406286       121619169       ct_1_2  100     .       121406286       121619169       224,224,224
    1       chr1    121619169       121625213       hor_1_1(S3C1H2-A,B,C)   100     .       121619169       121625213       255,146,0
    1       chr1    121625213       121667941       hor_1_2(S3C1H2-A,B)     100     .       121625213       121667941       255,146,0
    1       chr1    121667941       121788213       hor_1_3(S3C1H2-B)       100     .       121667941       121788213       255,146,0
    1       chr1    121788213       121790362       ct_1_3  100     .       121788213       121790362       224,224,224
    ```

