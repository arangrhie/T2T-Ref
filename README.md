# T2T-Ref
Pipeline for using T2T-CHM13v2.0 as a reference

## Mapping and variant calling for short-reads
This is an experimental pipeline for testing short-read based analysis variant calling methods using T2T-CHM13v2 as the reference. Acknowledging the limitations of difficult-to-map regions with short-reads, we developed a pipeline utilizing a masked reference. Masking has been designed to 1) complement sex-chromosome differences between XX and XY individuals and 2) detect robertsonian translocations, often resulting in copy loss in the acrocentric p-arms. The mitochondrial sequence has been replaced with the Cambridge Reference Sequence (rCRS, [NC_012920](https://www.ncbi.nlm.nih.gov/nuccore/251831106)) to abide with the rich annotation and resources.

The masked references described below is available to download [here](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/analysis_set/masked_DJ_rDNA_PHR/).

### Sex chromosome complement masking
The pipeline takes a mapped bam file on hg19, and determins the presence of Y chromosome. For XX samples, the Y chromosome has been entirely hard-masked in the reference. For XY samples, pseudo-autosomal region (PAR) on the Y has been hard-masked to call variants in diploid mode on the PAR region. This approach has been tested to reduce false positives and recover more variants in the PAR. More details are described in the [T2T-HG002Y paper](https://doi.org/10.1038/s41586-023-06457-y).

Masked PAR regions on the Y (includes telomere):
```
chrY    0       2458320       PAR1
chrY    62122809        62460029       PAR2
```

Masked Y for XX individuals (entire Y chromosome has been masked):
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
