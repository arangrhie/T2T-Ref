# Filter, merge VCF

Let's begin by merging 17 samples (4 families).
Collect the list of samples to merge in `test_sample.list`.

Let's filter first, then merge, to speed up the merging.

For each VCF,
```shell
# Input sample list
list=test_sample.list

# Output
mkdir -p analysis

set -e

for sample in $(cat $list)
do
  cd $sample
  VCF=dv_WGS.$sample.vcf.gz
  OUT=dv_WGS.$sample.GQ20
  set -x
  bcftools view -i 'GQ>20' --threads 12 --write-index -Oz $VCF -o $OUT.vcf.gz
  set +x
  cd ../
done

VCFs=`cat $list | awk '{print $1"/dv_WGS."$1".GQ20.vcf.gz"}' | tr '\n' ' '`
bcftools merge --threads 12 --write-index --regions chr22 -Oz -o analysis/test_samples.GQ20.vcf.gz $VCFs
```

Try filter after merging
```
VCFs=`cat $list | awk '{print $1"/dv_WGS."$1".vcf.gz"}' | tr '\n' ' '`
bcftools merge --threads 12 --write-index --regions chr22 -Oz -o analysis/test_samples.vcf.gz $VCFs
bcftools filter -i 'FMT/GQ>20' --threads 12 --write-index -Oz -o test_samples.mrg_GQ20.vcf.gz test_samples.vcf.gz
```
This seems more reasonable.

```shell
cd analysis
python3 count_genotypes.py test_samples.mrg_GQ20.vcf.gz > test_samples.mrg_GQ20.gt_count

# Include missings
cat test_samples.mrg_GQ20.gt_count | awk '{ref=$1+$3; alt=$2; print 100*alt/(ref+alt)}' | awk '$1>0' > test_samples.mrg_GQ20.gt_count.af

Rscript draw_dist.R test_samples.mrg_GQ20.gt_count.af
```

Let's try it on chr22, all parental samples. First, merge by 100 samples to check memory footprint. It used  ~2g. Trial with 500 samples used ~12g but runs much more slower. Trial with 200 samples worked in ~2.5 hrs, with 4g mem and 2~4 cpus. All split files have been merged per chromosomes for speed-up.


```shell
cd ../
list=parents.100.list
VCFs=`cat $list | awk '{print $1"/dv_WGS."$1".vcf.gz"}' | tr '\n' ' '`
bcftools merge --threads 24 --no-version --write-index --regions chr22 -Oz -o analysis/parents_chr22.500.vcf.gz $VCFs
## This takes a lot of memory, we need to split and merge
```

Let's use `parents_chr22.500.vcf.gz`; this merged version has the ALTs split by line. Can't use it.

```shell
cd analysis
IN_VCF=parents_chr22.500.vcf.gz
OUT_VCF=${IN_VCF/.vcf.gz/.GQ20.vcf.gz}
OUT=${OUT_VCF/.vcf.gz/}

bcftools filter -i 'FMT/GQ>20' --threads 24 --write-index -Oz -o $OUT_VCF $IN_VCF

# With all
python3 count_genotypes.py $OUT_VCF > $OUT.gt_count
cat $OUT.gt_count |\
 awk '{ref=$1+$3; alt=$2; print 100*alt/(ref+alt)}' |\
 awk '$1>0' > $OUT.gt_count.af_all

Rscript draw_dist.R $OUT.gt_count.af_all

# No missing
cat $OUT.gt_count |\
 awk '{ref=$1; alt=$2; print 100*alt/(ref+alt)}' |\
 awk '$1>0' > $OUT.gt_count.af_nom

Rscript draw_dist.R $OUT.gt_count.af_nom
```