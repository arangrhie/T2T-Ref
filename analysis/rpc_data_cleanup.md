# RPC Data Cleanup

## Rename sample dir and SAMPLE name in VCF
First, collect sample names that has to be replaced from the spreadsheet; `INOVA_sample_List_Fixed.xlsx`.

Save as `sample_names_to_replace.map`.

Copy one sample dir as `test`, and add `test2 test` in the first line of the `map`.
Test it works with the below as expected; using `for i in $(seq 1 1)`.

Run the following for all samples:
```
nohup ./fix_sample_names.sh > logs/fix_sample_names.log &
```

Below is how `fix_sample_names.sh` has been scripted:

```shell
#!/bin/sh

module load samtools

map=sample_names_to_replace.map

set -e

for i in $(seq 2 3031)
do
  echo "-- $i --"
  fixedName=`sed -n ${i}p $map | awk '{print $1}'`
  falseName=`sed -n ${i}p $map | awk '{print $2}'`
  set -x
  if ! [[ -d $fixedName ]]; then
  	mv $falseName $fixedName
  else
    echo "$fixedName already exists. Abort."
	exit -1
  fi
  cd $fixedName
  rename -v $falseName $fixedName $falseName.copynum_dip.txt
  sed -i "s/$falseName/$fixedName/g" sex.determine.txt
  sed -i "s/$falseName/$fixedName/g" $fixedName.copynum_dip.txt
  echo "default $fixedName" > dv_fix_sample_name.txt
  bcftools reheader -s dv_fix_sample_name.txt -o dv_WGS.$fixedName.vcf.gz dv_WGS.$falseName.vcf.gz
  bcftools index -t --threads 8 dv_WGS.$fixedName.vcf.gz
  rm dv_WGS.$falseName.vcf.gz dv_WGS.$falseName.vcf.gz.tbi
  cd ../
  set +x
  echo
done
```

## Rename SAMPLE in VCF
The rest didn't have to change the sample IDs. Only the `default` name that has been set as SAMPLE in DeepVariant has to be fixed.

These sample names have been saved in `sample_names_to_replace_in_vcf.list`.

The following has been run:
```
nohup ./fix_sample_names_in_vcf.sh >> logs/fix_sample_names_in_vcf.log
```

Below is how `fix_sample_names_in_vcf.sh` has been scripted:
```shell
#!/bin/sh

module load samtools

list=sample_names_to_replace_in_vcf.list

set -e

for i in $(seq 3 1157)
do
  echo "-- $i --"
  fixedName=`sed -n ${i}p $list`
  set -x
  cd $fixedName
  rename -v $fixedName default dv_WGS.$fixedName.vcf.gz*
  echo "default $fixedName" > dv_fix_sample_name.txt
  bcftools reheader -s dv_fix_sample_name.txt -o dv_WGS.$fixedName.vcf.gz dv_WGS.default.vcf.gz
  bcftools index -t --threads 8 dv_WGS.$fixedName.vcf.gz
  rm dv_WGS.default.vcf.gz dv_WGS.default.vcf.gz.tbi
  cd ../
  set +x
  echo
done
```
