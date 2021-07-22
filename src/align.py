"""
Functions to call aligners
"""

import os
from .utils import terminalExecute


def has_bwa_index(fasta_file: str, work_dir=None) -> bool:
    """
    Check if bwa index files exist in work_dir
    """
    if work_dir is None:
        work_dir = os.getcwd()
    dir_files = os.listdir(work_dir)
    index_extensions = ['amb', 'ann', 'bwt', 'pac', 'sa']
    return all([
        any([ext in fname for fname in dir_files])
        for ext in index_extensions]) == True

def makeBWAindex(fasta_file: str, db_prefix=None, suppress_output=False) -> None:
    """
    Make BWA index from fasta file
    """
    if db_prefix is None:
        db_prefix = os.path.basename(fasta_file).split('.')[0]
    bwa_index_command = f'bwa index -p {db_prefix} {fasta_file}'
    terminalExecute(bwa_index_command,suppress_output=suppress_output)
    
def bwaAlign(fasta_file: str, fastq_1_file: str, fastq_2_file: str=None, 
             n_threads: int=1, output_dir: str=None, only_mapped: bool=False,
             db_prefix=None, additional_params: str=None, suppress_output=False) -> None:
    """
    Align sequences to reference genome through BWA-mem
    only_mapped: return only primarily aligned fragments.
    SAM flag 4: Read unmapped (-F 4: filter out unmapped reads)
    """
    if db_prefix is None:
        db_prefix = os.path.basename(fasta_file).split('.')[0]
    if output_dir is None:
        output_dir = f'{os.path.basename(fastq_1_file).split("_1.")[0]}.sam'
    if only_mapped:
        output_str = f'| samtools view -S -h -F 4 - > {output_dir}'
    else:
        output_str = f'> {output_dir}'
    bwa_command = (f'bwa mem -M -t {n_threads} {db_prefix} '
                   f'{fastq_1_file} {fastq_2_file} {output_str}')
    terminalExecute(bwa_command, suppress_output=suppress_output)