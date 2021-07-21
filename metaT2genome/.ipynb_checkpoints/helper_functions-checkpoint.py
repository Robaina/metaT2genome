import os
import subprocess


def terminalExecute(command_str: str) -> None:
    """
    Execute given command in terminal through Python's subprocess
    """
    subprocess.run(command_str.split(' '))
    
def deleteTemporaryFiles() -> None:
    """
    Remove temporary SAM files from directory
    """

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

def makeBWAindex(fasta_file: str) -> None:
    """
    Make BWA index from fasta file
    """
    bwa_index_command = f'bwa index {fasta_file}'
    terminalExecute(bwa_index_command)

def bwaAlign(fasta_file: str, fastq_1_file: str, fastq_2_file=None, 
             n_threads=1, output_dir=None, 
             additional_params=None) -> None:
    """
    Align sequences to reference genome through BWA-mem
    """
    bwa_command = (f'bwa mem -M -t {n_threads} {additional_params} {fasta_file} '
                   f' {fastq_1_file} {fastq_2_file} > {out_sam_file}')
    terminalExecute(bwa_command)
    
def sortSAMbyName(sam_file: str, output_dir=None) -> None:
    """
    Sort SAM entries by name. Required by htseq-count to process paired-end data.
    """
    if output_dir is None:
        output_dir = f'{sam_file.split(".sam")[0]}_sorted.sam'
    samtools_command = f'samtools sort -n -O sam {sam_file} > {output_dir}'
    terminalExecute(samtools_command)  
    
def htseqCount(sorted_sam: str, gtf_file: str, feature_type='gene',
               feature_id='gene_id', output_dir=None, additional_params: str=None) -> None:
    """
    Count reads in SAM file through hseq-count
    Additional hseq-count params can be passed as a string to the argument: additional_params
    """
    if output_dir is None:
        output_dir = f'{sorted_sam.split(".sam")[0]}_counts.tsv'
    htseq_command = (f'htseq-count --order name --stranded yes --type {feature_type} '
                     f'--idattr {feature_id} {additional_params} {sorted_sam} '
                     f'{gtf_file} > {output_dir}')
    terminalExecute(htseq_command)