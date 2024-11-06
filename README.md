# T2T-Ref
Pipeline for using T2T-CHM13v2.0 as a reference

## Mapping and variant calling for short-reads
This is an experimental pipeline for testing short-read based analysis variant calling methods using T2T-CHM13v2 as the reference. Acknowledging the limitations of difficult-to-map regions with short-reads, we developed a pipeline utilizing a masked reference. Masking has been designed to 1) complement sex-chromosome differences between XX and XY individuals and 2) detect robertsonian translocations, often resulting in copy loss in the acrocentric p-arms. The mitochondrial sequence has been replaced with the Cambridge Reference Sequence (rCRS, [NC_012920](https://www.ncbi.nlm.nih.gov/nuccore/251831106)) to abide with the rich annotation and resources.

The masked references described below is available to download [here](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/analysis_set/masked_DJ_rDNA_PHR_5S_wi_rCRS/).

### Sex chromosome complement masking
The pipeline takes a mapped bam file on hg19 (or hg38), and determins the presence of Y chromosome. For XX samples, the Y chromosome has been entirely hard-masked in the reference. For XY samples, pseudo-autosomal region (PAR) on the Y has been hard-masked to call variants in diploid mode on the PAR region. This approach has been tested to reduce false positives and recover more variants in the PAR. More details are described in the [T2T-HG002Y paper](https://doi.org/10.1038/s41586-023-06457-y).

Masked PAR regions on the Y (includes telomere):
```
chrY    0       2458320       PAR1
chrY    62122809        62460029       PAR2
```

Masked Y for XX individuals (entire Y chromosome has been hard-masked):
```
chrY	0	62460029
```

### Acrocentric chromosomal arm masking
Each sex chromosome complement references are masked on the multi-copy regions of the acrocentric short arms in order to leave one representative copy present in the entire genome. This approach provides opportunities to further investigate highly homologous regions, which has been discovered to contain regions with high similarity and recombination between heterologous chromosomes, called pseudo-homologous region (PHR, [Guarracino et al.](https://doi.org/10.1038/s41586-023-05976-y)). In particular, our pipeline is designed to capture copy number differences and allelic variation within the distal junctions (DJs), rDNA repeat units, and the PHR surrounding the SST1 satellite.

Each of the five acrocentric short-arms contains a highly similar DJ, defined in earlier work ([van Sluis et al.](https://doi.org/10.1101/gad.331892.119)) as beginning before a CER block with an HERV-K integration, stretching through a large inverted repeat, and ending at the start of the rDNA repeat array. For our purposes, we defined the DJ to be the region from the end of the CER block to the beginning of the rDNA array, and masked this region on chromosomes 14, 15, 21, and 22, leaving the chromosome 13 copy as an unmasked representative copy.

The masked DJ locations are:
```
chr14   1755893 2099537 chr14_DJ
chr15   2164544 2506442 chr15_DJ
chr21   2762174 3108298 chr21_DJ
chr22   4448073 4793794 chr22_DJ
```

Kept:
```
chr13   5424523 5770548 chr13_DJ
```

Each of the five acrocentric short-arms contain a highly similar rDNA repeat array, made up of head-to-tail repeats of the ~43kb rDNA unit, which contains the 45S gene that is translated to rRNA and forms part of the final ribosome. In all but chromosome 13, we masked this entire array. In chromosome 13, we masked the entire array, aside from the first repeat unit, left to be the representative copy. The coordinates of the rDNA arrays come from the initial [T2T-CHM13 assembly paper (Nurk et al.)](https://doi.org/10.1126/science.abj6987) and the first repeat unit were defined using a self-mummer plot of chromosome 13.

Masked rDNAs:
```
chr13   5817416 9348041 chr13_rDNA
chr14   2099537 2817811 chr14_rDNA
chr15   2506442 4707485 chr15_rDNA
chr21   3108298 5612715 chr21_rDNA
chr22   4793794 5720650 chr22_rDNA
```

Kept:
```
chr13   5770548 5817416 chr13_first_rDNA
```

Chromosomes 13, 14, and 21 contain a region identified in [Guarracino et al.](https://doi.org/10.1038/s41586-023-05976-y) as residing within the proximal pseudo homologous region, consisting of a large inverted repeat surrounding some single-copy region containing an SST1 satellite block. This single-copy region containing the SST1 array was identified in the relevant chromosomes using self-mummer plots, and masked in all but chromosome 13. In chromosome 13, the single copy region surrounding the SST1 was defined separately for estimating coverage.

Masked PHR:
```
chr14   6582691 7707094 chr14_PHR
chr21   8706141 9832742 chr21_PHR
```

Kept:
```
chr13   11644595        12819139        chr13_PHR
```

Used for coverage estimates:
```
chr13   11644595        12301366        chr13_PHR_arm1
chr13   12440010        12819139        chr13_PHR_arm2
```

### Chromosome 1 5S rDNA Masking
Chromosome 1 contains a highly similar 5S rDNA repeat array, made up of head-to-tail repeats of the ~2.2kb 5S rDNA unit on the (-) strand, which contains the 5S gene that is translated to rRNA and forms part of the final ribosome. We extracted a single unit of this array using a self-mummer plot of chromosome 1, then aligned it back to chromosome 1 to define the array boundaries. We masked the entire array, aside from the first repeat unit, left to be the representative copy.

Masked 5S array:
```
chr1    227743723       228022744       chr1_5S_array
```

Kept:
```
chr1    228022744       228024984       chr1_first_5S_unit
```

## Dependency
For the short-read mapping and variant calling pipeline. These tools are assumed to be loadable as a module.
* bwa
* samtools
* bcftools
* bedtools
* parallel
* glnexus
* DeepVariant v1.6.0
* annovar 2020-06-08

## Quick start
1. Set env variable tools, which this git repo has been cloned under
	```
	export tools=/path/to/tools
	```

2. Download the reference files and ANNOVAR dbs from aws (assuming aws cli is available):
	```
	cd $tools/T2T-Ref/ref/
	aws s3 cp --no-sign-request --recursive s3://human-pangenomics/T2T/CHM13/assemblies/analysis_set/masked_DJ_rDNA_PHR_5S_wi_rCRS/ .
	cd $tools/T2T-Ref/annovar/db/
	aws s3 cp --no-sign-request --recursive s3://human-pangenomics/T2T/CHM13/assemblies/annotation/annovar/ .
	```
	Otherwise, download the [reference files](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/analysis_set/masked_DJ_rDNA_PHR_5S_wi_rCRS/) and [ANNOVAR dbs](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/annotation/annovar/), using https links.

3. Submit the pipeline
	Currently, the pipeline is built to work on Slurm environments, set up at hpc.nih.gov.

	The pipeline assumes read mappings (bam) to older reference (e.g. hg19) exists, and uses it to determine the presence of Y in each sample.

	Prepare an input file containing `sampleIndex <tab> sampleName <tab> path/to/old.bam`.

	For example, a typical sample_info file may look like this:
	```
	head sample_info.Inova.txt -n2
	1       102-00001-01    /data/path/to/Illumina/LP6005637-DNA_A01/Assembly/LP6005637-DNA_A01.bam
	2       102-00001-02    /data/path/to/Illumina/LP6005637-DNA_B01/Assembly/LP6005637-DNA_B01.bam
	```

	Submit to the schedular:
	```
	cd /path/to/produce/output/
	# Launch samples in lines 1-100 and 205
	$tools/T2T-Ref/variants_sr/_submit_extract_map_variant_sr.sh sample_info.Inova.txt 1-100,205
	```
	This will create a directory named after the `sampleName`, and place the output files there.

4. Post processing - copy number collection
	Once the mapping and copy number calculation has been done, it's possible to collect copy numbers of the DJ, rDNA (18S and 5S), and PHR. 
	Under `/path/to/produce/output/`, where all the output directories are available, run this script:
	```
	$tools/T2T-Ref/variants_sr/collect_copynum_dip.sh
	```
	This script will output the entire result in `copynum_dip.txt`, and prints potential Robertsonian translocation carriers information.
	It is currently relying on the DJ count, defined as less than 8.5 copies. Copy numbers are estimated in diploid (e.g. 10 acrocentric haplotypes in total for 5 acrocentric chromosomal arms).

### Inputs
* BAM file aligned to any reference.
* The table of array, sample and the location of Bam file. 

### Outputs
* Estimated coverage files:`${region}.coverage_results.bed`:
    * `$region` can be:
      * background
      * DJ
      * PHR_forCov
      * PHR_keep
      * rDNA
      * 5S
      * rDNA_18S
* VCF files generated with DeepVariant.	
