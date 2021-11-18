#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  kmer_plots.py
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

import argparse, os, glob
from subprocess import check_call, call
#get current directory
cwd = os.getcwd()
#set parameters
kmer_size = 21
cov_min = 5
cov_max = 15000


def init (sample_list, reads_dir, reads_fn, output_dir, ploidy, threads):
    
    for sample_name in sample_list:
        print (f"Starting k-mer analysis for sample {sample_name}.")
        #make output directories
        plot_dir = f'{output_dir}/kmer-plots/{sample_name}_{cov_max}_{kmer_size}/'
        if not os.path.exists(os.path.dirname(plot_dir)):
            os.makedirs(os.path.dirname(plot_dir))
        if not os.path.exists(os.path.dirname(plot_dir+'/tmp/')):
            os.makedirs(os.path.dirname(plot_dir+'/tmp/'))
        #change to the output directory
        os.chdir(plot_dir)
        #add input files
        input_fastq = glob.glob(f'{reads_dir}/{sample_name}/{sample_name}*{reads_fn}')
        with open ('FILES', 'w') as input_file:
            for fastq in input_fastq:
                input_file.write(fastq+'\n')
        #run kmc
        run_kmc (threads)
        #build smudgeplot
        build_smudgeplot (sample_name)
        #build genoscope graph
        run_genoscope (sample_name, ploidy)
        
        #change back to the original directory
        os.chdir(cwd)


def run_kmc (threads):
    
    # ~ #activate conda enviroment
    # ~ cmd = 'conda activate k-mer'
    # ~ call(["activate", 'k-mer'])
    
    cmd = f'kmc -k{kmer_size} -t{threads} -m64 -ci{cov_min} -cs{cov_max} -cx{cov_max} @FILES kmer_counts tmp'
    check_call(cmd, shell=True)
    
    cmd = f'kmc_dump kmer_counts kmer_k{kmer_size}.dump'# -ci{cov_min} -cx{cov_max}
    check_call(cmd, shell=True)
    
    #create histrogram for genoscope graphs
    cmd = f'kmc_tools transform kmer_counts histogram kmer_k{kmer_size}.hist -cx{cov_max}'
    check_call(cmd, shell=True)


def build_smudgeplot (sample_name):
    
    cmd = f'smudgeplot.py hetkmers -o kmer_pairs < kmer_k{kmer_size}.dump'
    check_call(cmd, shell=True)
    
    cmd = f'smudgeplot.py plot -o {sample_name} -t {sample_name} -q 0.99 kmer_pairs_coverages.tsv'
    check_call(cmd, shell=True)


def run_genoscope (sample_name, ploidy):
    
    cmd = f'genomescope2 -i kmer_k{kmer_size}.hist -k {kmer_size} -p {ploidy} -o . -n {sample_name}_genomescope' #-l {kmer_cov}
    check_call(cmd, shell=True)

    # ~ cmd = 'conda deactivate'
    # ~ check_call(cmd, shell=True)

if __name__ == '__main__':
    
    #input
    init (list(sample_list.split(',')), reads_dir, reads_fn, output_dir, ploidy, threads)
