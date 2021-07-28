#!/usr/bin/env python3
# coding: utf-8

"""
1. Filter SAM file by percent matched cutoff
2. Count reads with htseq-count
3. TPM-normalize counts

   Semidán Robaina Estévez (srobaina@ull.edu.es)
   Python >= 3.6
"""
import os
from metaT2genome.src.count import (htseqCount, tpmNormalizeHtseqOutput,
                                    aggregateTPMresults, aggregateCountsresults)
from metaT2genome.src.filtersam import filterSAMbyPercentMatched
from metaT2genome.src.utils import deleteTemporaryFiles

work_dir = os.getcwd()
gtf_file = 'MIT9301/MIT9301.gtf'
gbk_file = 'MIT9301/MIT9301.gb'

matched_cutoff = 50


# Iterate over identity cutoff
for identity_cutoff in [100, 98, 95]:
    
    counts_dir = f'counts_{identity_cutoff}'
    if not os.path.exists(counts_dir):
            os.makedirs(counts_dir)

    tpm_dir = f'tpm_{identity_cutoff}'
    if not os.path.exists(tpm_dir):
            os.makedirs(tpm_dir)
    
    filtered_dir = f'filtered_{identity_cutoff}_sam_files'
    samfiles = os.listdir(filtered_dir)
    n_conds = len(samfiles)

    # Iterate over conditions
    failed_conditions = []
    for n, samfile in enumerate(samfiles):

        condition = os.path.basename(samfile).split('.sam')[0]
        print(f'Processing condition ({n + 1}/{n_conds}): {condition}')

        try:
            print(f'\t1.Filtering SAM at {matched_cutoff}% percent matched')
            filterSAMbyPercentMatched(f'{filtered_dir}/{samfile}',
                                      matched_cutoff=matched_cutoff,
                                      output_path=f'temp/{condition}_PM_{matched_cutoff}.sam')

            
            print('\t2.Counting reads')
            htseqCount(f'temp/{condition}_PM_{matched_cutoff}.sam',
                       gtf_file, feature_type='gene',
                       feature_id='gene_id', output_dir=os.path.join(counts_dir, f'{condition}_counts.tsv'))

            print('\t3.TPM normalizing counts\n')
            tpmNormalizeHtseqOutput(os.path.join(counts_dir, f'{condition}_counts.tsv'), gbk_file,
                                    output_dir=os.path.join(tpm_dir, f'{condition}_tpm.tsv'))

            deleteTemporaryFiles('temp')

        except Exception:
            failed_conditions.append(condition)
            print(f'Failed condition: {condition} with exception {Exception}')


    aggregateTPMresults(tpm_dir=tpm_dir, output_dir=f'tpm_{identity_cutoff}.tsv')
    aggregateCountsresults(counts_dir=counts_dir, output_dir=f'counts_{identity_cutoff}.tsv')
    