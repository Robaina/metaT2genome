#!/usr/bin/env python3
# coding: utf-8

"""
1. Filter SAM file by percent identity cutoff
2. Sort paired-end entries by name
3. Count reads with htseq-count
4. TPM-normalize counts

   Semidán Robaina Estévez (srobaina@ull.edu.es)
   Python >= 3.6
"""
import os
from metaT2genome.src.count import (htseqCount, tpmNormalizeHtseqOutput,
                                    aggregateTPMresults)
from metaT2genome.src.filtersam import filterSAMbyIdentity, extractSegmentsWithMDtag
from metaT2genome.src.utils import sortSAMbyName, deleteTemporaryFiles

work_dir = os.getcwd()
gtf_file = 'MIT9301/MIT9301.gtf'
gbk_file = 'MIT9301/MIT9301.gb'
samfiles_dir = 'sam_files'
identity_cutoff = 100

filtered_dir = f'filtered_{identity_cutoff}_sam_files'
if not os.path.exists(filtered_dir):
        os.makedirs(filtered_dir)
        
counts_dir = f'counts_{identity_cutoff}'
if not os.path.exists(counts_dir):
        os.makedirs(counts_dir)
        
tpm_dir = f'tpm_{identity_cutoff}'
if not os.path.exists(tpm_dir):
        os.makedirs(tpm_dir)


# Iterate over conditions
failed_conditions = []
samfiles = os.listdir(samfiles_dir)
n_conds = len(samfiles)


for n, sam_file in enumerate(samfiles):
    
    condition = os.path.basename(samfile).split('.sam')[0]
    print(f'Processing condition ({n + 1}/{n_conds}): {condition}')
    
    try:
        print(f'\t1.Filtering SAM at {identity_cutoff}% identity')
        extractSegmentsWithMDtag(os.path.join(samfiles_dir, f'{condition}.sam'),
                                 output_dir=f'sam_files_with_md/{condition}.sam')
        
        filterSAMbyIdentity(f'sam_files_with_md/{condition}.sam',
                            identity_cutoff=identity_cutoff,
                            output_path=f'temp/{condition}_filtered_at_{identity_cutoff}.sam')

        print('\t2.Sorting SAM by name')
        # required by htseq-count
        sortSAMbyName(f'temp/{condition}_filtered_at_{identity_cutoff}.sam',
                      output_dir=f'{filtered_dir}/{condition}_filtered_at_{identity_cutoff}_sorted.sam')

        print('\t3.Counting reads')
        htseqCount(f'{filtered_dir}/{condition}_filtered_at_{identity_cutoff}_sorted.sam',
                   gtf_file, feature_type='gene',
                   feature_id='gene_id', output_dir=os.path.join(counts_dir, f'{condition}_counts.tsv'))

        print('\t4.TPM normalizing counts\n')
        tpmNormalizeHtseqOutput(os.path.join(counts_dir, f'{condition}_counts.tsv'), gbk_file,
                                output_dir=os.path.join(tpm_dir, f'{condition}_tpm.tsv')

        deleteTemporaryFiles('temp')
        
    except Exception:
        failed_conditions.append(condition)
        print(f'Failed condition: {condition} with exception {Exception}')


aggregateTPMresults(tpm_dir=tpm_dir, output_dir=f'tpm_{identity_cutoff}.tsv')
    