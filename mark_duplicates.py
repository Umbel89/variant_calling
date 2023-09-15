#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  untitled.py
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

import os, sys
import numpy as np
import subprocess
from multiprocessing import Process


def init (sample_list, sample_dict, output_dir, threads, parallelise):
    
    #make list into an array
    sample_array = np.array (sample_list)
    split_dict = {}
    #split sample list into multiple (number=parallelise) sublists and run them in parallel
    for i, sample_subarray in enumerate(np.array_split(sample_array, int(parallelise/2))):
        sample_sublist = list (sample_subarray)
        #run process in parallel
        split_dict[f'filter_pr{i}'] = Process(target=filter_bam, args=(sample_sublist, sample_dict, output_dir, threads))
        split_dict[f'filter_pr{i}'].start()
    #join processes when they are finished
    for split_pr in split_dict.keys():
        split_dict[split_pr].join()


def filter_bam (sample_list, sample_dict, output_dir, threads):
    
    for sample in sample_list:
        #create a list with all the bam files per sample
        input_list = []
        subsample_list = sample_dict[sample]
        for subsample in subsample_list:
            sample_join = ''.join([sample, subsample[:-1]])
            input_bam = f'{output_dir}/mapped/{sample_join}.bam'
            input_list.append(input_bam)
        #mark duplicates
        filtered_bam = mark_duplicates (sample, input_list, output_dir)
        #bam stats
        samtools_stats (sample, filtered_bam, output_dir, threads)


def mark_duplicates (sample, input_list, output_dir):
    
    filtered_bam = f'{output_dir}/mapped_filtered/{sample}.filtered.bam'
    
    #for multiple bam files per sample, build input command
    input_cmd = ''
    for bam_file in input_list:
        input_cmd += f' -I {bam_file}'
    
    if os.path.isfile(filtered_bam):
        print (f"Filtered bam file of sample {sample} already exists.")
    else:
        print (f"Running mark duplicates for sample {sample}.")
        # ~ cmd = f'gatk MergeSamFiles{input_cmd} -O {filtered_bam} -AS'
        cmd = f'gatk MarkDuplicates{input_cmd} -O {filtered_bam} -M {filtered_bam}-metrics.txt -AS\
                --REMOVE_DUPLICATES true --VERBOSITY ERROR --CREATE_INDEX true --TMP_DIR {output_dir}'
        subprocess.run(cmd, shell=True, check=True)
    
    if not os.path.isfile(filtered_bam+'.bai'):
        cmd = f'samtools index -@ 4 -b {filtered_bam}'
        subprocess.run(cmd, shell=True, check=True)
        
    return filtered_bam


def samtools_stats (sample, bam_file, output_dir, threads):
    
    output_fn = f'{output_dir}/mapped_filtered/{sample}'
    
    if not os.path.isfile(f'{output_fn}_stats.txt'):
        cmd = f'samtools stats -@ {threads} -d {bam_file} > {output_fn}_stats.txt'
        subprocess.run(cmd, shell=True, check=True)
        
        cmd = f'plot-bamstats -p {output_fn}/ {output_fn}_stats.txt'
        subprocess.run(cmd, shell=True, check=True)


if __name__ == '__main__':
    
    sample_list, output_dir, threads = sys.argv[1:]
    
    init (list(sample_list.split(',')), output_dir, threads)
