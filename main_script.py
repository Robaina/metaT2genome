#!/usr/bin/env python3
# coding: utf-8

"""
1. Align metaT to reference genome with bwa
2. Filter SAM file by percent identity cutoff
3. Sort paired-end entries by name
4. Count reads with htseq-count

   Semidán Robaina Estévez (srobaina@ull.edu.es)
   Python >= 3.6
"""
import os
from metaT2genome.src.align import has_bwa_index, makeBWAindex, bwaAlign
from metaT2genome.src.count import htseqCount
from metaT2genome.src.filtersam import filterSAMbyIdentity
from metaT2genome.src.utils import (terminalExecute, sortSAMbyName,
                                    deleteTemporaryFiles, getFastqPairedFiles)

work_dir = os.getcwd()
fasta_file = 'MIT9301/Prochlorococcus_marinus_str_MIT_9301.fasta'
gtf_file = 'MIT9301/MIT9301.gtf'
data_dir = '/usr/gonzalez/metagenomes/salazar2019/download/pe'
identity_cutoff = 95
n_threads = 20

if not has_bwa_index(work_dir):
    print('Building bwa index\n')
    makeBWAindex(fasta_file)

# Iterate over conditions
paired_fastqs = getFastqPairedFiles(data_dir) 
paired_fastqs = {k: v for k,v in paired_fastqs.items() if k == 'ERS488299'} # PRUEBA
for condition, (fastq_1_file, fastq_2_file) in paired_fastqs.items():
    
    print(f'Processing condition: {condition}')
    
    print('\t1.Aligning metaT to reference genome')
    bwaAlign(fasta_file, os.path.join(data_dir, fastq_1_file),
             os.path.join(data_dir, fastq_2_file), 
             n_threads=n_threads, output_dir=f'temp/{condition}.sam',
             only_mapped=True, additional_params=None)
    
    print(f'\t2.Filtering SAM at {identity_cutoff}% identity')
    filterSAMbyIdentity(f'temp/{condition}.sam', identity_cutoff=identity_cutoff,
                        output_path=f'temp/{condition}_filtered_at_{identity_cutoff}.sam')
    
    print('\t3.Sorting SAM by name')
    # required by htseq-count
    sortSAMbyName(f'temp/{condition}_filtered_at_{identity_cutoff}.sam',
                  output_dir=f'sam_files/{condition}_filtered_at_{identity_cutoff}_sorted.sam')
    
    print('\t4.Counting reads\n')
    htseqCount(f'sam_files/{condition}_filtered_at_{identity_cutoff}_sorted.sam',
               gtf_file, feature_type='gene',
               feature_id='gene_id', output_dir=f'counts/{condition}_counts.tsv',
               additional_params=None)

    deleteTemporaryFiles('temp')
    