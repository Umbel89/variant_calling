#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  compare_vcf.py
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

import sys, os, glob
from subprocess import check_call


def zip_index_vcf (vcf_file):
    
    cmd = f'bgzip {vcf_file}'
    check_call(cmd, shell=True)
    
    vcf_file = vcf_file+'.gz'
    
    cmd = f'tabix -p vcf {vcf_file}'
    check_call(cmd, shell=True)
    
    return vcf_file


def compare_vcf (ref_vcf, query_vcf, sample):
    
    comp_fn = f'{output_dir}/{sample}vs{query}.vchk'
    plot_dir = f'{output_dir}/{sample}vs{query}/'
    
    cmd = f'bcftools stats -F {reference_fasta} {ref_vcf} {query_vcf} > {comp_fn}'
    check_call(cmd, shell=True)
    
    cmd = f'plot-vcfstats -p {plot_dir} -s {comp_fn}'
    check_call(cmd, shell=True)


if __name__ == '__main__':
    
    query = 'JK'
    ref_vcf_dir = '/net/virus/linuxhome/petros/variant_calling/variants/sample_vcf/'
    query_vcf_dir = '/net/virus/linuxhome/petros/variant_calling/variants/sample_vcf_JK/'
    reference_fasta = '/net/virus/linuxhome/petros/variant_calling/reference/pfs1_scaffolds_pb.fasta'
    output_dir = '/net/virus/linuxhome/petros/variant_calling/variants/comparison/'
    if not os.path.exists(os.path.dirname(output_dir)):
        os.makedirs(os.path.dirname(output_dir))
    
    #find reference vcf files (gz or not)
    ref_vcf_list = glob.glob(ref_vcf_dir+'*_filtered_removed.vcf')
    ref_vcf_list.extend(glob.glob(ref_vcf_dir+'*_filtered_removed.vcf.gz'))
    #for each reference file, 
    for ref_vcf in ref_vcf_list:
        #exlude g.vcf files
        if not ref_vcf.endswith('.g.vcf'):
            sample = ref_vcf.split('/')[-1].split('_')[0]
            #locate the query file
            query_vcf = glob.glob(f'{query_vcf_dir}/{sample}_filtered_removed.vcf')
            query_vcf.extend(glob.glob(f'{query_vcf_dir}/{sample}_filtered_removed.vcf.gz'))
            #if query exists
            if query_vcf:
                query_vcf = query_vcf[0]
                print (sample)
                #check if vcf files are gz
                if not ref_vcf.endswith('.gz'):
                    ref_vcf= zip_index_vcf (ref_vcf)
                if not query_vcf.endswith('.gz'):
                    query_vcf = zip_index_vcf (query_vcf)
                #compare vcf files
                compare_vcf (ref_vcf, query_vcf, sample)
