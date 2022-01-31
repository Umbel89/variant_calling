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


def parse_sample_list ():
    
    #parse sample list from mutli vcf file
    cmd = f'bcftools query -l {input_vcf}'
    vcf_list = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    vcf_list = vcf_list.communicate()[0].strip().decode('UTF-8')

    #make output into a list
    sample_list = vcf_list.split('\n')
    
    return sample_list


def parse_sample (sample):
    
    sample_vcf = f'{output_dir}/{sample}.vcf.gz'
    
    #parse sample from multi vcf file | remove all the lines that have no information or are identical to the reference
    if not os.path.isfile(sample_vcf) or not os.path.getsize(sample_vcf) > 0:
        cmd = f'bcftools view -s {sample} {input_vcf} | bcftools view -e \'GT~"\."\' - | bcftools view -e \'GT="ref"\' - -o {sample_vcf} -O z'
        check_call(cmd, shell=True)
        
    #index bgzip vcf file
    if not os.path.isfile(sample_vcf+'.tbi') or not os.path.getsize(sample_vcf+'.tbi') > 0:
        cmd = f'tabix -p vcf {sample_vcf}'
        check_call(cmd, shell=True)
    
    return sample_vcf


def filter_vcf (sample, sample_vcf):
    
    filtered_vcf = f'{output_dir}/{sample}_filtered.vcf.gz'
    
    if not os.path.isfile(filtered_vcf) or not os.path.getsize(filtered_vcf) > 0:
        cmd = f'gatk VariantFiltration -R {reference_fasta} -V {sample_vcf} -O {filtered_vcf} {filter_cmd}'
        check_call(cmd, shell=True)
    
    return filtered_vcf


def remove_filtered (sample, filtered_vcf):
    
    removed_vcf = f'{output_dir}/{sample}_filtered_removed.vcf.gz'
    #remove lines that didn't pass the filters, unused alternative alleles, indels larger than 1
    if not os.path.isfile(removed_vcf) or not os.path.getsize(removed_vcf) > 0:
        cmd = f'gatk SelectVariants -R {reference_fasta} -V {filtered_vcf} -O {removed_vcf} \
                --exclude-non-variants True --exclude-filtered True --remove-unused-alternates True' #--max-indel-size 0
        check_call(cmd, shell=True)
        
    hom_vcf = f'{output_dir}/{sample}_filtered_removed_hom.vcf.gz'
    
    if not os.path.isfile(hom_vcf) or not os.path.getsize(hom_vcf) > 0:
        cmd = f'bcftools view -i \'GT="hom"\' {removed_vcf} -o {hom_vcf} -O z'
        check_call(cmd, shell=True)
        
    #index bgzip vcf file
    if not os.path.isfile(hom_vcf+'.tbi') or not os.path.getsize(hom_vcf+'.tbi') > 0:
        cmd = f'tabix -p vcf {hom_vcf}'
        check_call(cmd, shell=True)
    
    return removed_vcf


def measure_variants (sample, filtered_vcf):
    
    comp_fn = f'{output_dir}/{sample}_filtered.vchk'
    plot_dir = f'{output_dir}/{sample}_filtered/'
    
    cmd = f'bcftools stats -f .,PASS -F {reference_fasta} -s- {filtered_vcf} > {comp_fn}'
    check_call(cmd, shell=True)
    
    cmd = f'plot-vcfstats -p {plot_dir} -s {comp_fn}'
    check_call(cmd, shell=True)


if __name__ == '__main__':
    
    reference_fasta = "/net/virus/linuxhome/petros/variant_calling/Pe14/reference/Pe14_assembly.fasta"
    input_vcf = "/net/virus/linuxhome/petros/variant_calling/Pe14/genotyped/Pe14.vcf.gz"
    output_dir = "/net/virus/linuxhome/petros/variant_calling/Pe14/genotyped/"
    if not os.path.exists(os.path.dirname(output_dir+'/')):
        os.makedirs(os.path.dirname(output_dir+'/'))
    #split multi vcf file into each sample
    split_multi_vcf = False
    #name or list of names of the output file ('all' to run for all samples in vcf file)
    input_sample = 'Pe14'
    #filtering options
    filter_cmd = f'--filter-name "QD4" --filter-expression "QD < 4.0" \
                   --filter-name "FS60" --filter-expression "FS > 60.0" \
                   --filter-name "MQ20" --filter-expression "MQ < 20.0" \
                   --filter-name "SOR3" --filter-expression "SOR > 3.0"\
                   --filter-name "MQRS-3" --filter-expression "MQRankSum < -3.0" \
                   --filter-name "ReadPosRS" --filter-expression "ReadPosRankSum < -1.0 || ReadPosRankSum > 3.5"'
                   
    #if input is all, parse sample list from multi vcf file
    if input_sample == 'all':
        sample_list = parse_sample_list ()
    else:
        sample_list = input_sample.split(',')
    
    for i, sample in enumerate(sample_list, 1):
        print (i, sample)
        #if the input is multi vcf file, parse one sample
        if split_multi_vcf:
            sample_vcf = parse_sample (sample)
        else:
            sample_vcf = input_vcf
        #filter sample
        filtered_vcf = filter_vcf (sample, sample_vcf)
        removed_vcf = remove_filtered (sample, filtered_vcf)
        measure_variants (sample, filtered_vcf)
