#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  haplotype_caller.py
#  
#  Copyright 2020 Petros <petros@pskiadas.gr>
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

import glob, sys, os
import numpy as np
from subprocess import check_call
from multiprocessing import Process


def init (sample_list, reference_fasta, output_dir, input_type, genotyping, ploidy, parallelise, input_intervals, all_sites):
    
    reference_dictionary (reference_fasta)
    #make list into an array
    sample_array = np.array (sample_list)
    split_dict = {}
    #split sample list into multiple (number=parallelise) sublists and run them in parallel
    for i, sample_subarray in enumerate(np.array_split(sample_array, parallelise)):
        sample_sublist = list (sample_subarray)
        #run process in parallel
        split_dict[f'genotype_pr{i}'] = Process(target=genotype_bam, args=(sample_sublist, reference_fasta, output_dir, input_type, genotyping, ploidy, input_intervals, all_sites))
        split_dict[f'genotype_pr{i}'].start()
    #join processes when they are finished
    for split_pr in split_dict.keys():
        split_dict[split_pr].join()


def genotype_bam (sample_list, reference_fasta, output_dir, input_type, genotyping, ploidy, input_intervals, all_sites):
    
    for sample in sample_list:
        input_bam = f'{output_dir}/mapped_filtered/{sample}.filtered.bam'
        run_haplotypecaller (input_bam, sample, reference_fasta, output_dir, input_type, genotyping, ploidy, input_intervals, all_sites)


def reference_dictionary (reference_fasta):
    
    #create reference dictionary
    reference_dict = '.'.join(reference_fasta.split('.')[:-1])+'.dict'
    if not os.path.isfile(reference_dict):
        cmd = f'gatk CreateSequenceDictionary -R {reference_fasta} -O {reference_dict}'
        check_call(cmd, shell=True)


def run_haplotypecaller (input_bam, sample, reference_fasta, output_dir, input_type, genotyping, ploidy, input_intervals, all_sites):
    
    if genotyping == 'joint':
        erc_mode = 'GVCF'
        output_fn = f'{output_dir}/genotyped/{sample}.g.vcf.gz'
    else:
        erc_mode = 'NONE'
        output_fn = f'{output_dir}/genotyped/{sample}.vcf.gz'
    
    if input_type == 'transcriptomic':
        min_thr = '20'
    elif input_type == 'genomic':
        min_thr = '30'
        
    if all_sites:
        output_mode = 'EMIT_ALL_ACTIVE_SITES'
    else:
        output_mode = 'EMIT_VARIANTS_ONLY'
    
    if input_intervals == 'all':
        intervals = ''
    else:
        intervals = f'-L {input_intervals}'
    
    if not os.path.isfile(output_fn):
        cmd = f'gatk --java-options "-Xmx8g -Djava.io.tmpdir={output_dir}/genotyped/" HaplotypeCaller -ERC {erc_mode} \
                --verbosity ERROR -VS LENIENT --native-pair-hmm-threads 8 \
                -ploidy {ploidy} -stand-call-conf {min_thr} -I {input_bam} -O {output_fn} -R {reference_fasta} \
                --output-mode {output_mode} {intervals}'
        check_call(cmd, shell=True)


if __name__ == '__main__':
    
    sample_list, reference_fasta, output_dir, input_type, genotyping = sys.argv[1:]
    
    init (list(sample_list.split(',')), reference_fasta, output_dir, input_type, genotyping, ploidy, parallelise, input_intervals, all_sites)
