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

import glob, os
import subprocess


def parse_sample_list ():
    """Parse each line of the file and return it into a list.
    """
    sample_list = []
    
    #check if input is file
    with open (sample_txt) as the_file:
        for line in the_file:
            sample_name = line.strip()
            sample_list.append(sample_name)
    
    return sample_list


def samtools_stats (bam_file, sample):
    
    output_fn = f'{input_dir}/{sample}'
    
    if not os.path.isfile(f'{output_fn}_stats.txt'):
        cmd = f'samtools stats -@ 16 -d {bam_file} > {output_fn}_stats.txt'
        subprocess.run(cmd, shell=True, check=True)
        
        cmd = f'plot-bamstats -p {output_fn}/ {output_fn}_stats.txt'
        subprocess.run(cmd, shell=True, check=True)


if __name__ == '__main__':
    
    input_dir = '/net/virus/linuxhome/petros/variant_calling/uu_mapped_filtered/Pfs15/'
    
    sample_list = parse_sample_list ()
    bam_list = glob.glob(input_dir+'/*.bam')
    
    for i, bam_file in enumerate(bam_list, 1):
        sample = bam_file.split('/')[-1].split('.')[0]
        print (i, sample)
        samtools_stats (bam_file, sample)
