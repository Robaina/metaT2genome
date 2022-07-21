"""
Functions to call mapped read counters and process count data
"""
import os
import numpy as np
import pandas as pd
from .utils import terminalExecute
from .gbk import getGeneLengthsFromGBK


def htseqCount(sorted_sam: str, gtf_file: str, feature_type='gene',
               feature_id='gene_id', output_dir=None,
               additional_params: str=None, suppress_output=False) -> None:
    """
    Count reads in SAM file through hseq-count
    Additional hseq-count params can be passed as a string to the argument: 
    additional_params
    """
    if output_dir is None:
        output_dir = f'{os.path.splitext(sorted_sam)[0]}_counts.tsv'
    if additional_params is None:
        add_args = ''
    else:
        add_args = additional_params
    htseq_command = (f'htseq-count --order name --stranded yes --type {feature_type} '
                     f'--idattr {feature_id} --quiet {sorted_sam} '
                     f'{gtf_file} {add_args} > {output_dir}')
    terminalExecute(htseq_command, suppress_output=suppress_output)
    
def getCountDataframeFromHtseqOutput(counts_tsv: str):
    """
    Get counts pandas dataframe from htseq-count tsv output.
    Suitable for input to tpmNormalizeCounts
    """
    counts = pd.read_csv(counts_tsv, delimiter='\t', header=None,
                     names=['gene_id', 'counts'])
    counts['gene_id'] = counts['gene_id'].apply(lambda s: s.strip())
    counts = counts.set_index('gene_id')
    for i, row in counts.iterrows():  # Remove summary data entries
        if row.name.startswith('__'):
            counts = counts.drop(row.name)
    return counts
        
def tpmNormalizeCounts(counts, gbk_file: str, feature_type: str = 'gene'):
    """
    Apply TPM normalization to read counts (pandas dataframe)
    tpm.to_csv('tpm.tsv', sep='\t') to save as tsv
    """
    gene_lengths = getGeneLengthsFromGBK(gbk_file, feature_type=feature_type)
    counts_gene_lengths = {}
    for gene_id in counts.index:
        if gene_id in gene_lengths.keys():
            counts_gene_lengths[gene_id] = gene_lengths[gene_id]
        else:
            counts_gene_lengths[gene_id] = np.nan
            
    rpk = counts.divide(list(counts_gene_lengths.values()), axis=0)
    tpm = rpk.divide(rpk.sum(axis=0).values / 1e6, axis=1)
    tpm.rename(columns={'counts': 'tpm'})
    return tpm

def tpmNormalizeHtseqOutput(counts_tsv: str, gbk_file: str,
                            output_dir: str=None,
                            feature_type: str = 'gene') -> None:
    """
    TPM-normalize hsetq tsv output and export to tsv
    """
    if output_dir is None:
        output_dir = f'{os.path.splitext(counts_tsv)[0]}_tpm.tsv'
    
    counts = getCountDataframeFromHtseqOutput(counts_tsv)
    tpm = tpmNormalizeCounts(counts, gbk_file, feature_type=feature_type)
    tpm.to_csv(output_dir, sep='\t')
    
def aggregateTPMresults(tpm_dir: str, sep_type: str='\t',
                        output_dir: str=None,
                        name_split_pattern: str = None) -> None:
    """
    Aggregate sample tpm files into a single file
    tpm_dir: path to directory containing tpm files
    """
    if output_dir is None:
        output_dir = os.path.join(tpm_dir, 'aggregated_tpm.tsv')
    tpm_files = os.listdir(tpm_dir)
    li = []
    conditions = []
    for fname in tpm_files:
        if name_split_pattern is not None:
            sample_id = os.path.basename(fname).split(name_split_pattern)[0]
        else:
            sample_id = os.path.basename(fname)
        conditions.append(sample_id)
        df = pd.read_csv(os.path.join(tpm_dir, fname), index_col='gene_id',
                         header=0, sep=sep_type)
        li.append(df)

    df = pd.concat(li, axis=1, ignore_index=True)
    df.columns = conditions
    df.to_csv(output_dir, sep=sep_type)
    
def aggregateCountsresults(counts_dir: str, sep_type: str='\t',
                           output_dir: str=None,
                           name_split_pattern: str = None) -> None:
    """
    Aggregate sample counts files into a single file
    counts_dir: path to directory containing tpm files
    """
    if output_dir is None:
        output_dir = os.path.join(counts_dir, 'aggregated_counts.tsv')
    counts_files = os.listdir(counts_dir)
    li = []
    conditions = []
    for fname in counts_files:
        if name_split_pattern is not None:
            sample_id = os.path.basename(fname).split(name_split_pattern)[0]
        else:
            sample_id = os.path.basename(fname)
        conditions.append(sample_id)
        df = pd.read_csv(os.path.join(counts_dir, fname),
                         header=None, sep=sep_type)
        df.columns = ['gene_id', 'counts']
        df = df.set_index('gene_id')
        # Remove summary data entries
        for i, row in df.iterrows():
            if row.name.startswith('__'):
                df = df.drop(row.name)
        li.append(df)
    
    df = pd.concat(li, axis=1, ignore_index=True)
    df.columns = conditions
    df.to_csv(output_dir, sep=sep_type)

    