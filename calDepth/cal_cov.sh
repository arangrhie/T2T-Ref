#! /bin/bash
sample=$1

module load samtools
module load bedtools

# Background
# bed="$refDir/masked_DJ_rDNA_PHR.autosome.rCRS.forXY.fa.fai.win10k.autosoneONLY.woAcrocentric.bed"
bed="/data/NHGRIgeno/juhyun_work/test_cnv/chm13.v2.autosome.wo_acro.bed"
prefix="background"
samtools depth -@ $SLURM_CPUS_PER_TASK  -b $bed $sample.dedup.bam > $sample.$prefix.coverage_results.depth
awk '{print $1,$2,$2,$3}' OFS='\t' $sample.$prefix.coverage_results.depth > $sample.$prefix.coverage_results.depth.tmp
bedtools map -a $bed -b $sample.$prefix.coverage_results.depth.tmp -c 4 -o median > $sample.$prefix.coverage_results.txt && rm $sample.$prefix.coverage_results.depth.tmp

# rDNA
bed="/data/Phillippy2/projects/rpc/01.mappoed/00.chm13v2.v0.3/first_rDNA_chm13.bed"
prefix="rDNA"
samtools depth -@ $SLURM_CPUS_PER_TASK -b $bed $sample.dedup.bam > $sample.$prefix.coverage_results.depth
awk '{print $1,$2,$2,$3}' OFS='\t' $sample.$prefix.coverage_results.depth > $sample.$prefix.coverage_results.depth.tmp
bedtools map -a $bed -b $sample.$prefix.coverage_results.depth.tmp -c 4 -o median > $sample.$prefix.coverage_results.txt && rm $sample.$prefix.coverage_results.depth.tmp

# DJ
bed='/data/NHGRIgeno/juhyun_work/test_cnv/DJ_on_chr13.bed'
prefix="DJ"
samtools depth -@ SLURM_CPUS_PER_TASK -b $bed $sample.dedup.bam > $sample.$prefix.coverage_results.depth
awk '{print $1,$2,$2,$3}' OFS='\t' $sample.$prefix.coverage_results.depth > $sample.$prefix.coverage_results.depth.tmp
bedtools map -a $bed -b $sample.$prefix.coverage_results.depth.tmp -c 4 -o median > $sample.$prefix.coverage_results.txt && rm $sample.$prefix.coverage_results.depth.tmp

# PHR
bed="/data/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/beds_for_rpc/mask_DJ_rDNA_PHR/PHR_keep.bed"
prefix="PHR_keep"
samtools depth -@ $SLURM_CPUS_PER_TASK  -b $bed $sample.dedup.bam > $sample.$prefix.coverage_results.depth
awk '{print $1,$2,$2,$3}' OFS='\t' $sample.$prefix.coverage_results.depth > $sample.$prefix.coverage_results.depth.tmp
bedtools map -a $bed -b $sample.$prefix.coverage_results.depth.tmp -c 4 -o median > $sample.$prefix.coverage_results.txt && rm $sample.$prefix.coverage_results.depth.tmp

bed='/data/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/beds_for_rpc/mask_DJ_rDNA_PHR/PHR_keep_for_cov.bed'
prefix="PHR_forCov"
samtools depth -@ $SLURM_CPUS_PER_TASK -b $bed $sample.dedup.bam > $sample.$prefix.coverage_results.depth
awk '{print $1,$2,$2,$3}' OFS='\t' $sample.$prefix.coverage_results.depth > $sample.$prefix.coverage_results.depth.tmp
bedtools map -a $bed -b $sample.$prefix.coverage_results.depth.tmp -c 4 -o median > $sample.$prefix.coverage_results.txt && rm $sample.$prefix.coverage_results.depth.tmp

if [ $(ls -l ./ | wc -l) -eq 5 ]; then touch depthCal.done ; fi
