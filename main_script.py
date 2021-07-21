#!/usr/bin/env python
# coding: utf-8

"""
1. Align metaT to reference genome with bwa
2. Filter SAM file by percent identity cutoff
3. Sort paired-end entries by name
4. Count reads with htseq-count

   Python >= 3.6
"""
import os
import subprocess
import metaT2genome.filter_by_identity as fi
import metaT2genome.helper_functions as hf

__author__ = 'Semidán Robaina Estévez'
__email__ = 'srobaina@ull.edu.es'

work_dir = os.getcwd()
fasta_file = ''
gtf_file = ''
data_dir = '/usr/gonzalez/metagenomes/salazar2019/download/pe'
identity_cutoff = 95
n_threads = 1

if not hf.has_bwa_index(work_dir):
    hf.makeBWAindex(fasta_file)

# Iterate over conditions
paired_fastqs = hf.getFastqPairedFiles(data_dir) 
for condition, (fastq_1_file, fastq_2_file) in paired_fastqs.items():
    
    print(f'Processing condition: {condition}')
    hf.terminalExecute('conda activate samtools')
    
    print('\t1.Aligning metaT to reference genome')
    hf.bwaAlign(fasta_file, fastq_1_file, fastq_2_file, 
                n_threads=n_threads, output_dir=f'{condition}.sam',
                only_mapped=True, additional_params=None)
    
    print(f'\t2.Filtering SAM at {identity_cutoff}% identity')
    fi.filterSAMbyIdentity(f'{condition}.sam', identity_cutoff=identity_cutoff,
                           output_path=f'{condition}_filtered_at_{identity_cutoff}.sam')
    
    print('\t3.Sorting SAM by name')
    # required by htseq-count
    hf.sortSAMbyName(f'{condition}_filtered_at_{identity_cutoff}.sam',
                     output_dir=f'{condition}_filtered_at_{identity_cutoff}_sorted.sam')

    terminalExecute('conda deactivate samtools')
    terminalExecute('conda activate htseq')
    
    print('\t4.Counting reads\n')
    hf.htseqCount(f'{condition}_filtered_at_{identity_cutoff}_sorted.sam',
                  gtf_file, feature_type='gene',
                  feature_id='gene_id', output_dir=f'{condition}_counts.tsv',
                  additional_params=None)

    terminalExecute('conda deactivate htseq')