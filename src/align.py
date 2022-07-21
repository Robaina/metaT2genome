"""
Functions to call aligners
"""

import os
from .utils import terminalExecute


def has_bwa_index(work_dir: str = None, genome_prefix: str = None) -> bool:
    """
    Check if bwa index files exist in work_dir
    """
    if work_dir is None:
        work_dir = os.getcwd()
    if genome_prefix is None:
        genome_prefix = ''
    dir_files = os.listdir(work_dir)
    index_extensions = ['.amb', '.ann', '.bwt', '.pac', '.sa']
    index_files = [genome_prefix + ext for ext in index_extensions]
    return all([
        any([ifname in fname for fname in dir_files])
        for ifname in index_files]) == True

def makeBWAindex(fasta_file: str, db_prefix=None, suppress_output=False) -> None:
    """
    Make BWA index from fasta file
    """
    if db_prefix is None:
        db_prefix = os.path.basename(fasta_file).split('.')[0]
    bwa_index_command = f'bwa index -p {db_prefix} {fasta_file}'
    terminalExecute(bwa_index_command,suppress_output=suppress_output)
    
def bwaAlign(fasta_file: str, fastq_1_file: str, fastq_2_file: str=None, 
             n_threads: int=1, output_dir: str=None, db_prefix=None,
             additional_params: str=None, suppress_output=False) -> None:
    """
    Align sequences to reference genome through BWA-mem
    only_mapped: return only primarily aligned fragments.
    SAM flag 4: Read unmapped (-F 4: filter out unmapped reads)
    """
    if db_prefix is None:
        db_prefix = os.path.basename(fasta_file).split('.')[0]
    if output_dir is None:
        output_dir = f'{os.path.basename(fastq_1_file).split("_1.")[0]}.sam'
    if additional_params is None:
        add_args = ''
    else:
        add_args = additional_params
    bwa_command = (f'bwa mem -M -t {n_threads} {db_prefix} '
                   f'{fastq_1_file} {fastq_2_file} {add_args} > {output_dir}')
    terminalExecute(bwa_command, suppress_output=suppress_output)