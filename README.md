# Variant Calling
Variant calling pipeline for short-read sequences using a single reference genome with GATK.

## Prerequisites
- [The Genome Analysis Toolkit (GATK) v4.3.0.0](https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip)
- [bwa-mem2 v.2.2.1](https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2)

## Installation
```bash
git clone https://github.com/Umbel89/variant_calling.git
```

## Usage
```
python variant_calling.py --help

usage: variant_calling.py [-h] -r FILE -s STR -o STR -n STR [-f STR] [-p INT] [-e STR] [-i {genomic,transcriptomic}]
                          [-g {joint,individual}] [-a] [-t INT] [-m INT] [-c INT] [-l INTERVALS] [-sw STR]

Run the mapping and variant calling pipeline.

optional arguments:
  -h, --help            show this help message and exit
  -r FILE, --reference FILE
                        Reference fasta file.
  -s STR, --samples STR
                        Directory with one subdirectory per sample with the sample reads.
  -o STR, --output STR  Directory where output will be written.
  -n STR, --projectname STR
                        Name of the Project. Can contain only lowercase or upercase letters, numeric characters, and
                        dashes (-). Maximum length: 38 characters
  -f STR, --readfn STR  Extension of sample reads filename. [default= _filtered_rmMin30G.fastq.gz]
  -p INT, --ploidy INT  Ploidy of the samples. [default=2]
  -e STR, --sample_list STR
                        A file or a comma seperated list of sample names, one per line. Only the sample names in this
                        file will be used for variant calling.
  -i {genomic,transcriptomic}, --input {genomic,transcriptomic}
                        Type of Input Data. [default=genomic]
  -g {joint,individual}, --genotyping {joint,individual}
                        Type of Genotyping. [default=joint]
  -a, --all_sites       The output vcf from joint genotyping includes variant and invariant sites.
  -t INT, --threads INT
                        Number of threads to use on the mapping
  -m INT, --memory INT  Gibabytes of memory to use during joint genotyping - dbimport.
  -c INT, --parallelise INT
                        Number of parallel run of haplotype caller, between 1 and 12. Recomended 8 cpus, 8 Gb of
                        memory per run. [default=4]
  -l INTERVALS, --intervals INTERVALS
                        Intervals of bam files where Haplotype Caller, and Joint Genotyping will run. Input can be
                        either one interval or a file [list or bed format] with one interval per line.
                        [list=<chr>:<start>-<stop>] [example=chr2:1-2000]
  -sw STR, --switches STR
                        Disable steps of the pipeline. [1-5=Run, 0=Disable] [default=2345] steps: [1 kmer analysis, 2
                        mapping, 3 mark duplicates, 4 haplotype caller, 5 joint genotyping]
```

## Cite
[Skiadas, P. *et al*. (2022), Environmental Microbiology: Sexual reproduction contributes to the evolution of resistance-breaking isolates of the spinach pathogen *Peronospora effusa*](https://enviromicro-journals-onlinelibrary-wiley-com.proxy.library.uu.nl/doi/10.1111/1462-2920.15944)
