"""
Functions to call mapped read counters and process count data
"""
import numpy as np
import pandas as pd
from .utils import terminalExecute
from .gbk import getGeneLengthsFromGBK


def htseqCount(sorted_sam: str, gtf_file: str, feature_type='gene',
               feature_id='gene_id', output_dir=None,
               additional_params: str=None) -> None:
    """
    Count reads in SAM file through hseq-count
    Additional hseq-count params can be passed as a string to the argument: 
    additional_params
    """
    if output_dir is None:
        output_dir = f'{os.path.splitext(sorted_sam)[0]}_counts.tsv'
    htseq_command = (f'htseq-count --order name --stranded yes --type {feature_type} '
                     f'--idattr {feature_id} --quiet {sorted_sam} '
                     f'{gtf_file} > {output_dir}')
    terminalExecute(htseq_command)
    
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
        
def tpmNormalizeCounts(counts, gbk_file: str):
    """
    Apply TPM normalization to read counts (pandas dataframe)
    tpm.to_csv('tpm.tsv', sep='\t') to save as tsv
    """
    gene_lengths = getGeneLengthsFromGBK('Data/MIT9301.gb')
    counts_gene_lengths = {}
    for gene_id in counts.index:
        if gene_id in gene_lengths.keys():
            counts_gene_lengths[gene_id] = gene_lengths[gene_id]
        else:
            counts_gene_lengths[gene_id] = np.nan
            
    rpk = counts.divide(list(counts_gene_lengths.values()), axis=0)
    tpm = rpk.divide(rpk.sum(axis=0).values / 1e6, axis=1)
    return tpm

def tpmNormalizeHtseqOutput(counts_tsv: str, gbk_file: str,
                            output_dir: str=None) -> None:
    """
    TPM-normalize hsetq tsv output and export to tsv
    """
    if output_dir is None:
        output_dir = f'{os.path.splitext(counts_tsv)[0]}_tpm.tsv'
    
    counts = getCountDataframeFromHtseqOutput(counts_tsv)
    tpm = tpmNormalizeCounts(counts, gbk_file)
    tpm.to_csv(output_dir, sep='\t')