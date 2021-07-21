#!/usr/bin/env python
# coding: utf-8

"""
1. Align metaT to reference genome with bwa
2. Filter SAM file by percent identity cutoff
3. Count reads with htseq-count

python >= 3.6
"""
import os
import subprocess
import metaT2genome.filter_by_identity as fi
import metaT2genome.helper_functions as hf

__author__ = 'Semidán Robaina Estévez'
__email__ = 'srobaina@ull.edu.es'

work_dir = os.getcwd()
hf.terminalExecute('conda activate samtools')

# Make index if not present
if not hf.has_bwa_index(work_dir):
    hf.makeBWAindex(fasta_file)

# Align 
hf.bwaAlign(fasta_file, fastq_1_file, fastq_2_file, 
            n_threads=1, output_dir=None)

# Filter by identity
fi.filterSAMbyIdentity(input_path, identity_cutoff=95, output_path=None)

# Sort-by-name filtered SAM (required by htseq-count)
hf.sortSAMbyName(sam_file, output_dir=None)

terminalExecute('conda deactivate samtools')
terminalExecute('conda activate htseq')

# Call htseq-count (using default parameters)
hf.htseqCount(sorted_sam, gtf_file, feature_type='gene',
              feature_id='gene_id', output_dir=None)

terminalExecute('conda deactivate htseq')