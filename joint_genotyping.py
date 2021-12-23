#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  joint_genotyping.py
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

import sys, os, glob, re
import time, datetime, logging, subprocess
from subprocess import check_call


def init (sample_list, reference_fasta, output_dir, projectname, threads, input_intervals, all_sites, memory):
    
    start = time.time()
    #create log file
    log_fn = f'{output_dir}/logs/{projectname}_joint-genotyping.log'
    logging.basicConfig(filename=log_fn, format='%(levelname)s: %(message)s', level=logging.DEBUG)
    
    input_list = glob.glob(f'{output_dir}/genotyped/*.g.vcf.gz')
    
    #preprocessing of the gvcf files
    check_input (sample_list, input_list, output_dir)
    #if there are no intervals make a file for with all the regions
    if not input_intervals:
        intervals_fn = create_intervals (reference_fasta, output_dir)
    else:
        intervals_fn = input_intervals
    preprocessing = time.time ()
    
    print ("Joint Genotyping")
    genotyping_start = time.time()
    database = run_dbimport (intervals_fn, sample_list, output_dir, threads,  input_intervals, memory)
    database_time = time.time()
    output_vcf = run_genotype (database, sample_list, reference_fasta, output_dir, projectname, all_sites, intervals_fn)
    genotyping_end = time.time()
    #calculate vcf statistics
    bcftools_stats (output_vcf, reference_fasta, output_dir, projectname)
    
    log_time (start, genotyping_start, preprocessing, database_time, genotyping_end)


def check_input (sample_list, input_list, output_dir):
    
    print ("Preprocessing")
    
    error_list = []
    #check if input gvcf exists for all samples
    for sample in sample_list:
        input_gvcf = f'{output_dir}/genotyped/{sample}.g.vcf.gz'
        if not input_gvcf in input_list:
            error_list.append(sample)
        else:
            index_gvcf (input_gvcf)
    #if any gvcf is missing, print error and exit
    if error_list:
        print (f"The following list of samples were not genotyped succesfully:\n{', '.join(error_list)}.")
        sys.exit(1)


def index_gvcf (input_gvcf):
    
    #if not indexed, create index for gvcf files
    if not os.path.isfile(input_gvcf+'.tbi'):
        print (f"Creating gvcf index for sample {sample}")
        cmd = f'tabix -p vcf {input_gvcf}'
        check_call(cmd, shell=True)


def create_intervals (reference_fasta, output_dir):
    
    ref_fai = reference_fasta+'.fai'
    intervals_fn = f'{output_dir}/reference/intervals.list'
    
    with open (intervals_fn, 'w') as output_fn:
        with open (ref_fai) as the_file:
            for line in the_file:
                contig, size = line.split('\t')[:2]
                output_fn.write (f'{contig}:{1}-{size}\n')
            
    return intervals_fn


def run_dbimport (intervals_fn, sample_list, output_dir, threads, input_intervals, memory):
    
    workspace_dir = f'{output_dir}/variants/db_workspace'
    db_input = ''
    merge_intervals = ''
    reader_threads = 1
    
    for sample in sample_list:
        input_gvcf = f'{output_dir}/genotyped/{sample}.g.vcf.gz'
        db_input+=(' -V '+input_gvcf)
    
    if not input_intervals:
        merge_intervals = '--merge-input-intervals'
    else:
        merge_intervals = ''
    
    #check if database already exists
    complete=True
    with open (intervals_fn) as the_file:
        for line in the_file:
            contig, start_end = line.strip().split(':')#re.split('[:-]', line.strip())
            start, end = start_end.split('-')
            if not os.path.isdir(f'{workspace_dir}/{contig}${start}${end}'):
                complete=False
    
    if complete:
        print ("Database already exists.")
    else:
        print ("Creating database of the gvcf files.")
        
        cmd = f'gatk --java-options "-Xmx{memory}g -Xms8g -Djava.io.tmpdir={output_dir}/variants/" GenomicsDBImport {db_input} \
                -L {intervals_fn} --tmp-dir {output_dir}/variants/ \
                --genomicsdb-workspace-path {workspace_dir} --batch-size 70 \
                --seconds-between-progress-updates 120 --reader-threads {threads} {merge_intervals}'
        check_call(cmd, shell=True)
    
    return workspace_dir


def run_genotype (database, sample_list, reference_fasta, output_dir, projectname, all_sites, intervals_fn):
    
    output_vcf = f'{output_dir}/variants/{projectname}.vcf.gz'
    
    #output vcf with all the reference sites
    if all_sites:
        vcf_type = '--all-sites --annotate-with-num-discovered-alleles -stand-call-conf 0'
    else:
        vcf_type = ''
    
    if not os.path.isfile(output_vcf) or not os.path.getsize(output_vcf) > 0:
        print (f"Joint genotyping")
        cmd = f'gatk --java-options "-Xmx32g -Djava.io.tmpdir={output_dir}/variants/" GenotypeGVCFs -V gendb://{database} \
                -R {reference_fasta} -O {output_vcf} --tmp-dir {output_dir}/variants/ -L {intervals_fn}\
                -G StandardAnnotation --seconds-between-progress-updates 120 {vcf_type}'
        check_call(cmd, shell=True)
    else:
        print (f"vcf file for project {projectname} already exists.")
    
    #parse sample list from vcf file
    cmd = f'bcftools query -l {output_vcf}'
    vcf_list = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    vcf_list = vcf_list.communicate()[0].strip().decode('UTF-8')
    
    #check if sample list in vcf is the same with input sample list
    if vcf_list.split('\n').sort() == sample_list.sort():
        print ("Joint Genotyping for the input samples finished successfully.")
    else:
        print ("The samples of the output vcf do not much the input samples.")
    
    return output_vcf


def bcftools_stats (output_vcf, reference_fasta, output_dir, projectname):
    
    comp_fn = f'/{output_dir}/variants/{projectname}.vchk'
    plot_dir = f'/{output_dir}/variants/{projectname}/'
    
    if not os.path.isfile(comp_fn):
        cmd = f'bcftools stats -F {reference_fasta} -s- {output_vcf} > {comp_fn}'
        check_call(cmd, shell=True)
    
    if not os.path.exists(os.path.dirname(plot_dir+'/')):
        cmd = f'plot-vcfstats -p {plot_dir} -s {comp_fn}'
        check_call(cmd, shell=True)


def log_time (start, genotyping_start, preprocessing, database_time, genotyping_end):
    
    #calculate times
    preprocessing_time = datetime.timedelta(seconds=round(preprocessing - start))
    database_real = datetime.timedelta(seconds=round(database_time - preprocessing))
    genotyping_time = datetime.timedelta(seconds=round(genotyping_end - database_time))
    total_time = datetime.timedelta(seconds=round(genotyping_end - start))
    #write info to log file
    logging.info (f'Preprocessing Time: {preprocessing_time}')
    logging.info (f'Database Time: {database_real}')
    logging.info (f'Joint Genotyping Time: {genotyping_time}')
    logging.info (f'Total Time: {total_time}')


if __name__ == '__main__':
    
    #import directories and project name
    sample_list, reference_fasta, output_dir, projectname = sys.argv[1:]
    
    init (list(sample_list.split(',')), reference_fasta, output_dir, projectname, threads, input_intervals, all_sites)
