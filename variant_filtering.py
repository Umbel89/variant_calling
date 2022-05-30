#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  variant_filtering.py
#  
#  Copyright 2020 Petros Skiadas <p.skiadas@uu.nl>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

import sys, os, glob, subprocess
from subprocess import check_call

#filtering options
filter_cmd = f'--filter-name "QD4" --filter-expression "QD < 4.0" \
               --filter-name "FS60" --filter-expression "FS > 60.0" \
               --filter-name "MQ20" --filter-expression "MQ < 20.0"'
               # ~ --filter-name "MAF01" --filter-expression "MAF < 0.1" \
               # ~ --filter-name "SOR3" --filter-expression "SOR > 3.0"\
               # ~ --filter-name "MQRS-3" --filter-expression "MQRankSum < -3.0" \
               # ~ --filter-name "ReadPosRS" --filter-expression "ReadPosRankSum < -1.0 || ReadPosRankSum > 3.5"'


def init(input_vcf, reference_fasta, output_dir, sample):
    
    #filter sample
    filtered_vcf = filter_vcf (sample, reference_fasta, input_vcf, output_dir)
    #measure after filtering
    measure_variants (f'{sample}_filtered', reference_fasta, filtered_vcf, output_dir)
    #remove filtered and unused variants
    removed_vcf = remove_filtered (sample, reference_fasta, filtered_vcf, output_dir)
    #remove repetative regions
    unique_vcf = remove_repeats (sample, reference_fasta, removed_vcf, output_dir)
    #measure after removing repeats
    measure_variants (f'{sample}_unique', reference_fasta, unique_vcf, output_dir)


def filter_vcf (sample, reference_fasta, sample_vcf, output_dir):
    
    filtered_vcf = f'{output_dir}/{sample}_filtered.vcf.gz'
    
    if not os.path.isfile(filtered_vcf) or not os.path.getsize(filtered_vcf) > 0:
        cmd = f'gatk VariantFiltration -R {reference_fasta} -V {sample_vcf} -O {filtered_vcf} {filter_cmd}'
        check_call(cmd, shell=True)
    
    return filtered_vcf


def remove_filtered (sample, reference_fasta, filtered_vcf, output_dir):
    
    removed_vcf = f'{output_dir}/{sample}_filtered_removed.vcf.gz'
    #remove lines that didn't pass the filters, unused alternative alleles
    if not os.path.isfile(removed_vcf) or not os.path.getsize(removed_vcf) > 0:
        cmd = f'gatk SelectVariants -R {reference_fasta} -V {filtered_vcf} -O {removed_vcf} \
                --exclude-non-variants True --exclude-filtered True --remove-unused-alternates True' #--max-indel-size 0' 
        check_call(cmd, shell=True)
    
    return removed_vcf


def remove_repeats (sample, reference_fasta, removed_vcf, output_dir):
    
    reference = reference_fasta.split('/')[-1].split('_')[0]
    repeat_gff = f"/net/virus/linuxhome/petros/Pe_genome_assemblies/{reference}_nanopore/final/repeats/{reference}_assembly.fasta.out.gff"

    unique_vcf = f'{output_dir}/{sample}_unique.vcf'
    #remove lines that didn't pass the filters, unused alternative alleles
    if not os.path.isfile(unique_vcf) or not os.path.getsize(unique_vcf) > 0:
        cmd = f'bedtools intersect -a {removed_vcf} -b {repeat_gff} -header -wa -v > {unique_vcf}'
        check_call(cmd, shell=True)
    
    unique_vcf = compress_input (unique_vcf)
    
    return unique_vcf


def compress_input (input_fn):
    
    if not input_fn.endswith('.gz'):
        cmd = f'bgzip {input_fn}'
        check_call(cmd, shell=True)
        input_fn = input_fn+'.gz'
    
    #if not indexed, create index for vcf files
    if not os.path.isfile(input_fn+'.tbi'):
        cmd = f'tabix -p vcf {input_fn}'
        check_call(cmd, shell=True)
    
    return input_fn


def measure_variants (sample, reference_fasta, filtered_vcf, output_dir):
    
    comp_fn = f'{output_dir}/{sample}.vchk'
    plot_dir = f'{output_dir}/{sample}/'
    
    cmd = f'bcftools stats -f .,PASS -F {reference_fasta} -s- {filtered_vcf} > {comp_fn}'
    check_call(cmd, shell=True)
    
    cmd = f'plot-vcfstats -p {plot_dir} -s {comp_fn}'
    check_call(cmd, shell=True)


if __name__ == '__main__':
    
    sample, reference_fasta, output_dir, input_vcf = sys.argv[1:]
    
    init(input_vcf, reference_fasta, output_dir, sample)
