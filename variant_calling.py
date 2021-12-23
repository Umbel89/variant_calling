#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  variant_calling.py
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

import numpy as np
import argparse, os, subprocess, re, logging, sys, time, datetime
from multiprocessing import Process
from shutil import copy2

#import python scripts from the directory of the current script
sdir = f'{os.path.dirname(os.path.realpath(__file__))}/'
sys.path.append('sdir')
import bwa_mapping, mark_duplicates, haplotype_caller, joint_genotyping, kmer_plots


def init (reference_fasta, output_dir, threads, projectname):
    
    if not os.path.exists(os.path.dirname(output_dir+'/')):
        os.makedirs(os.path.dirname(output_dir+'/'))
    if not os.path.exists(os.path.dirname(output_dir+'/mapped/')):
        os.makedirs(os.path.dirname(output_dir+'/mapped/'))
    if not os.path.exists(os.path.dirname(output_dir+'/mapped_filtered/')):
        os.makedirs(os.path.dirname(output_dir+'/mapped_filtered/'))
    if not os.path.exists(os.path.dirname(output_dir+'/logs/')):
        os.makedirs(os.path.dirname(output_dir+'/logs/'))
    if not os.path.exists(os.path.dirname(output_dir+'/reference/')):
        os.makedirs(os.path.dirname(output_dir+'/reference/'))
    if not os.path.exists(os.path.dirname(output_dir+'/genotyped/')):
        os.makedirs(os.path.dirname(output_dir+'/genotyped/'))
    if not os.path.exists(os.path.dirname(output_dir+'/variants/')) and args.genotyping == 'joint':
        os.makedirs(os.path.dirname(output_dir+'/variants/'))
    if k_mers and not os.path.exists(os.path.dirname(output_dir+'/kmer-plots/')):
        os.makedirs(os.path.dirname(output_dir+'/kmer-plots/'))
    
    #write arguments in a txt file in log directory
    arg_fn = f'{output_dir}/logs/{projectname}_arguments.txt'
    with open (arg_fn, 'a') as output_arg:
        output_arg.write (f'## {start_time}\n')
        for argument in vars(args):
            output_arg.write (f'{argument}: {getattr(args, argument)}\n')
    
    #copy reference to its directory
    reference = reference_fasta.split('/')[-1]
    new_reference = output_dir+'/reference/'+reference
    if not os.path.isfile(new_reference):
        copy2 (reference_fasta, new_reference)
    
    #create log file
    logging.basicConfig(filename=f'{output_dir}/logs/{projectname}_variant-calling.log', format='%(levelname)s: %(message)s', level=logging.DEBUG)
    logging.info (f'Threads: {threads}')
    
    return new_reference


def check_sample_input (reads_dir, input_samples, reads_fn, output_dir, projectname):
    
    sample_dict = {}
    
    #parse list of subdirectories and files
    sample_lib = list(os.walk (reads_dir, followlinks=True))
    
    #if there is an input sample list
    if input_samples:
        sample_list = parse_sample_list (input_samples)
    else:
        sample_list = sample_lib[0][1]
    sample_fn = f'{output_dir}/logs/{projectname}_samples.txt'
    
    with open (sample_fn, 'w') as output_list:
        for sample_info in sample_lib:
            sample_dir = os.path.realpath(sample_info[0])
            sample_name = sample_dir.split('/')[-1]
            if sample_name in sample_list:
                #find upload file in the sample directory
                sample_reads = [fn for fn in sample_info[2] if fn.endswith(reads_fn)]
                #if there are no sample read files, print error and exit
                if not sample_reads:
                    print (f"The sample directory {sample_dir} contains no fastq file(s) with the suffix {reads_fn}.")
                    sys.exit()
                #if two or more files per sample, parse their names
                elif len(sample_reads) >= 2:
                    #check if input is file
                    for sample_fn in sample_reads:
                        os.path.isfile (os.path.join(sample_dir,sample_fn))
                    #parse names and add in a dictionary
                    sample_dict[sample_name] = []
                    for fastq in sample_reads:
                        subsample = re.search(sample_name+r'(.*)R.*$', fastq).group(1)
                        sample_dict[sample_name] += [subsample]
                #if there is not the right amount of input data per sample, print error and exit
                else:
                    print (f"The sample directory {sample_dir} should contain at least two fastq files.")
                    sys.exit()
                #write sample name to the txt file
                output_list.write (sample_name+'\n')
    
    return sample_list, sample_dict


def parse_sample_list (input_samples):
    """Check if input is a file, then parse each line of the file and 
    return it into a list.
    """
    sample_list = []
    
    with open (input_samples) as the_file:
        for line in the_file:
            sample_name = line.strip()
            sample_list.append(sample_name)
    
    return sample_list


def infasta (string):
    """Check if (input) file exists and is fasta."""
    
    fasta_list = ('.fasta', '.fa', '.fna')
    
    string = inFile (string)
    if not string.endswith(fasta_list):
        print (f"The input reference has to end with {fasta_list}.")
        raise IOError("not a fasta file [%s]" % string)
    
    return string


def inFile (string):
    """Check if (input) file exists."""
    
    if string != "" and string is not None:
        string = os.path.abspath(string)
        if not os.path.isfile(string):
            print (f"The input is not a file.")
            raise IOError("not a file [%s]" % string)
    
    return string


def check_name (input_name):
    """Check for illegal characters in input."""
    
    illegal_characters = ["[", "[", "{", "}", ":", ";", ">" "<", "@", "!", "#", "$", "^", "&", "*", "(", ")", "=", "+", "|", "~", '"', "'", ".", "/", "\\" ]
    
    if any (character in input_name for character in illegal_characters):
        print ("Project name can contain only lowercase letters, numeric characters, dashes (-), and underscores (_).")
        sys.exit()
    elif len(input_name) > 38:
        print ("Project name can be up to 38 characters.")
        sys.exit()
    
    return input_name


def create_switches (steps):
    
    switches = []
    
    if steps == '0':
        switches = [False, False, False, False, False]
    else:
        for num in range (1, 6):
            if str(num) in steps:
                switches.append(True)
            else:
                switches.append(False)
    
    return switches


def setting ():
    """Create options."""
    
    #range of options
    options_input = ['genomic', 'transcriptomic']
    options_genotyping = ['joint', 'individual']
    #default options
    def_readfn = '_filtered_rmMin30G.fastq.gz'
    def_switches = '2345'
    def_sample_list = ''
    def_intervals = ''
    def_ploidy = 2
    def_threads = 24
    def_memory = 32
    def_parallelise = int(def_threads/6)
    #input parser
    parser = argparse.ArgumentParser(description="Run the mapping and variant calling pipeline.")
    parser.add_argument('-r', '--reference', help="Reference fasta file.", metavar="FILE", type=infasta, required=True)
    parser.add_argument('-s', '--samples', help="Directory with one subdirectory per sample with the sample reads.", metavar="STR", type=str, required=True)
    parser.add_argument('-o', '--output', help="Directory where output will be written.", metavar="STR", type=str, required=True)
    parser.add_argument('-n', '--projectname', help="Name of the Project. Can contain only lowercase or upercase letters, numeric characters, and dashes (-). Maximum length: 38 characters", type=check_name, metavar="STR", required=True)
    #optional arguments
    parser.add_argument('-f', '--readfn', default=def_readfn, help= f"Extension of sample reads filename. [default= {def_readfn}]", type=str, metavar="STR")
    parser.add_argument('-p', '--ploidy', default=def_ploidy, help=f"Ploidy of the samples. [default={def_ploidy}]", type=int, metavar="INT")
    parser.add_argument('-e', '--sample_list', default=def_intervals, help="A file with a list of sample names, one per line. Only the sample names in this file will be used for variant calling.", type=inFile, metavar="FILE")
    parser.add_argument('-i', '--input', default=options_input[0], help=f"Type of Input Data. [default={options_input[0]}]", choices=options_input)
    parser.add_argument('-g', '--genotyping', default=options_genotyping[0], choices=options_genotyping, help=f"Type of Genotyping. [default={options_genotyping[0]}]")
    parser.add_argument('-a', '--all_sites', action='store_true', help="The output vcf from joint genotyping includes variant and invariant sites.")
    parser.add_argument('-t', '--threads', default=def_threads, help=f"Number of threads to use on the mapping", type=int, metavar="INT")
    parser.add_argument('-m', '--memory', default=def_memory, help=f"Gibabytes of memory to use during joint genotyping - dbimport.", type=int, metavar="INT")
    parser.add_argument('-c', '--parallelise', default=def_parallelise, choices=range(1,13), help=f"Number of parallel run of haplotype caller, between 1 and 12. Recomended 8 cpus, 8 Gb of memory per run. [default={def_parallelise}]", type=int, metavar="INT")
    parser.add_argument('-l', '--intervals', default=def_intervals , help=f"Intervals of bam files where Haplotype Caller, and Joint Genotyping will run. Input can be either one interval or a file [list or bed format] with one interval per line. [list=<chr>:<start>-<stop>] [example=chr2:1-2000]", type=inFile)
    parser.add_argument('-sw', '--switches', default=def_switches, help=f"Disable steps of the pipeline. [1-5=Run, 0=Disable] [default={def_switches}] steps: [1 kmer analysis, 2 mapping, 3 mark duplicates, 4 haplotype caller, 5 joint genotyping]", metavar="STR", type=create_switches)
    args = parser.parse_args()
    
    return args


if __name__ == '__main__':
    
    start = time.time()
    start_time = time.strftime("%Y-%m-%d %H:%M:%S")
    
    #import user input
    args = setting()
    #step switches
    k_mers = args.switches[0]
    mapping = args.switches[1]
    markduplicates = args.switches[2]
    haplotyper = args.switches[3]
    genotyping = args.switches[4]
    
    reference_fasta = init (args.reference, args.output, args.threads, args.projectname)
    #parse samples
    sample_list, sample_dict = check_sample_input (args.samples, args.sample_list, args.readfn, args.output, args.projectname)
    
    #fastq k-mer analysis
    if k_mers:
        kmer_plots.init (sample_list, args.samples, args.readsfn, args.output, args.ploidy, args.ncpus)
    #run mapping
    if mapping:
        bwa_mapping.init (sample_dict, args.samples, args.output, reference_fasta, args.threads, args.projectname, args.genotyping, args.readfn)
    if markduplicates:
        mark_duplicates.init (sample_list, sample_dict, args.output, args.threads, args.parallelise)
    #run variant calling
    if haplotyper:
        haplotype_caller.init (sample_list, reference_fasta, args.output, args.input, args.genotyping, args.ploidy, args.parallelise, args.intervals, args.all_sites)
    if args.genotyping == 'joint' and genotyping:
        joint_genotyping.init (sample_list, reference_fasta, args.output, args.projectname, args.threads, args.intervals, args.all_sites, args.memory)
