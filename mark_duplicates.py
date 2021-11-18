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

import glob, os, sys
from subprocess import check_call


def init (sample_list, output_dir, threads):
    
    for sample in sample_list:
        input_bam = glob.glob(f'{output_dir}/mapped/{sample}_*.bam')
        #mark duplicates
        filtered_bam = mark_duplicates (sample, input_bam, output_dir)
        #bam stats
        samtools_stats (sample, filtered_bam, output_dir, threads)


def mark_duplicates (sample, input_bam, output_dir):
    
    filtered_bam = f'{output_dir}/mapped_filtered/{sample}.filtered.bam'
    
    #for multiple bam files per sample, build input command
    input_cmd = ''
    for bam_file in input_bam:
        input_cmd += f' -I {bam_file}'
    
    if os.path.isfile(filtered_bam):
        print (f"Filtered bam file of sample {sample} already exists.")
    else:
        print (f"Running mark duplicates for sample {sample}.")
        cmd = f'gatk MarkDuplicates{input_cmd} -O {filtered_bam} -M {filtered_bam}-metrics.txt \
                --REMOVE_DUPLICATES true --VERBOSITY ERROR --CREATE_INDEX true --TMP_DIR {output_dir}'
        check_call(cmd, shell=True)
    
    # ~ if not os.path.isfile(filtered_bam+'.bai'):
        # ~ cmd = f'{samtools} index -@ 4 -b {filtered_bam}'
        # ~ check_call(cmd, shell=True)
        
    return filtered_bam


def samtools_stats (sample, bam_file, output_dir, threads):
    
    output_fn = f'{output_dir}/mapped_filtered/{sample}'
    
    if not os.path.isfile(f'{output_fn}_stats.txt'):
        cmd = f'samtools stats -@ {threads} -d {bam_file} > {output_fn}_stats.txt'
        check_call(cmd, shell=True)
        
        cmd = f'plot-bamstats -p {output_fn}/ {output_fn}_stats.txt'
        check_call(cmd, shell=True)


if __name__ == '__main__':
    
    sample_list, output_dir, threads = sys.argv[1:]
    
    init (list(sample_list.split(',')), output_dir, threads)
