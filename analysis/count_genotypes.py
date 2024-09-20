#!/bin/python3

import sys
import gzip

def count_genotypes(vcf_file):
    ref_count = 0
    alt_count = 0
    missing = 0

    # Open the VCF file, handling both plain and gzipped files
    if vcf_file.endswith('.gz'):
        open_func = gzip.open
        mode = 'rt'
    else:
        open_func = open
        mode = 'r'

    with open_func(vcf_file, mode) as file:
        for line in file:
            if line.startswith('#'):
                continue  # Skip header lines
            columns = line.strip().split('\t')
            format_column = columns[8].split(':')
            gt_index = format_column.index('GT')
            ad_index = format_column.index('AD')

            # Initialize
            ref_count = 0
            alt_count = 0
            missing = 0

            for sample in columns[9:]:
                gt = sample.split(':')[gt_index]
                alleles = gt.replace('|', '/').split('/')
                for allele in alleles:
                    if allele == '0':
                        ref_count += 1
                    elif allele != '.':
                        alt_count += 1
                    elif allele == '.':
                        # For missings, if there are AD re-assign based on variant allele frequency
                        ad = sample.split(':')[ad_index]
                        if ad == '.':
                            missing += 1
                        else:
                            ads = ad.split(',')
                            ref_ad = int(ads[0])
                            alt_ad = 0
                            for allele_d in ads[1:]:
                                if allele_d != '.':
                                    alt_ad += int(allele_d)
                            vaf = alt_ad * 100 / (ref_ad + alt_ad)
                            if (vaf > 70):   alt_count += 2
                            elif (vaf > 30):
                                ref_count += 1
                                alt_count += 1
                            else: ref_count += 2

            print(f"{ref_count}\t{alt_count}\t{missing}")

    return

if __name__ == "__main__":
    vcf_file = sys.argv[1]
    count_genotypes(vcf_file)
