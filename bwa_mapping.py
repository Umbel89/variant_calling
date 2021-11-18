#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  bwa_mapping.py
#  
#  Copyright 2018 Petros Skiadas <p.skiadas@genetwister.nl>
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
import time, datetime, logging
from subprocess import check_call, Popen, PIPE, STDOUT

start = time.time()


def init (sample_dict, reads_dir, output_dir, reference_fasta, threads, projectname, genotyping):
    
    #create log file
    log_fn = f'{output_dir}/logs/{projectname}_bwa-mapping.log'
    logging.basicConfig(filename=log_fn, format='%(levelname)s: %(message)s', level=logging.DEBUG)
    logging.info (f'Threads: {threads}')
    
    #index reference fasta if it's not already indexed
    index_fasta(reference_fasta)
    
    #return pairs
    for sample, subsample, R1, R2 in return_pairs (sample_dict, reads_dir):
        sample_start = time.time()
        #create read groups for bam files
        read_group = find_read_group (R1, sample, subsample, genotyping, projectname)
        #map reads to reference
        run_bwa (read_group, reference_fasta, R1, R2, sample, subsample, output_dir, threads)
        mapping = time.time()
        #write times to log
        log_time (sample, sample_start, mapping)
    
    end = datetime.timedelta(seconds=round(time.time() - start))
    logging.info (f'Total Time: {end}')


def return_pairs (sample_dict, reads_dir):
    
    for sample, subsample_list in sample_dict.items():
        sample_dir = reads_dir+sample
        for subsample in subsample_list:
            sample_join = '_'.join([sample, subsample])
            R1 = f'{reads_dir}/{sample}/{sample_join}_R1_filtered_rmMin30G.fastq.gz'
            R2 = f'{reads_dir}/{sample}/{sample_join}_R2_filtered_rmMin30G.fastq.gz'
            yield sample, subsample, R1, R2


def find_read_group (R1, sample, subsample, genotyping, projectname):
    
    if subsample:
        group_id = subsample
    else:
        cmd = f'zcat {R1} | head -n 1'
        
        result = Popen(cmd, stdout=PIPE, shell=True, stderr=STDOUT)
        fasta_header = result.communicate()[0].strip()[1:].decode('UTF-8')
        #remove any tabs or spaces
        fasta_header = ''.join(fasta_header.split())
        
        read_id = '.'.join(fasta_header.split(':')[:2])
        flowcell = '.'.join(fasta_header.split(':')[2:4])
        group_id = read_id+flowcell
        
    #create read groups for the gvcf or the individual mode, read group strings need to be connected with \\t and not \t otherwise gatk will return an error
    if genotyping == 'joint':
        read_group = f'@RG\\tID:{group_id}\\tSM:{sample}\\tLB:{projectname}\\tPL:ILLUMINA'
    elif genotyping == 'individual':
        read_group = f'@RG\\tID:{group_id}_{sample}\\tSM:{sample}\\tLB:{projectname}\\tPL:ILLUMINA'
    
    return read_group


def index_fasta (reference_fasta):
    
    #index reference
    if not os.path.isfile(reference_fasta+'.bwt') or not os.path.getsize(reference_fasta+'.bwt') > 0:
        print ("Indexing reference with bwa.")
        cmd = f'bwa index {reference_fasta}' 
        check_call(cmd, shell=True)
        
        index_time = datetime.timedelta(seconds=round(time.time() - start))
        logging.info (f'Indexing Time: {index_time}')
    
    if not os.path.isfile(reference_fasta+'.fai') or not os.path.getsize(reference_fasta+'.fai') > 0:
        cmd = f'samtools faidx {reference_fasta}'
        check_call(cmd, shell=True)


def run_bwa (read_group, reference_fasta, R1, R2, sample, subsample, output_dir, threads):
    
    output_name = '_'.join([sample, subsample])
    output_fn = f'{output_dir}/mapped/{output_name}.bam'
    
    if not os.path.isfile(output_fn) or not os.path.getsize(output_fn) > 0:
        #map with bwa and convert to bam and sort with samtools
        print (f"Mapping with bwa mem for sample {sample}.")
        cmd = f'bwa mem -M -R "{read_group}" -t {threads} {reference_fasta} {R1} {R2} \
                | samtools sort -@ {threads} -m 1G - -o {output_fn}'
        check_call(cmd, shell=True)
    else:
        print (f"bam file of sample {sample} already exists.")
    
    index_bam (output_fn)


def index_bam (bam_file):
    
    if not os.path.isfile(f'{bam_file}.bai'):
        cmd = f'samtools index -@ 4 -b {bam_file}'
        check_call(cmd, shell=True)


def log_time (sample, sample_start, mapping):
    
    #calculate times
    sample_time = datetime.timedelta(seconds=round(mapping - sample_start))
    #write info to log file
    logging.info (f'{sample}: {sample_time}')


if __name__ == '__main__':
    
    #input
    sample_dict, reads_dir, output_dir, reference_fasta = sys.argv[1:3]
    threads, projectname, genotyping = sys.argv[3:]
    
    init (sample_dict, reads_dir, output_dir, reference_fasta, threads, projectname, genotyping)
