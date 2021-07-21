import os
import subprocess
import numpy as np


def terminalExecute(command_str: str) -> None:
    """
    Execute given command in terminal through Python's subprocess
    """
    subprocess.run(command_str.split(' '), shell=True, capture_output=True)
    
def deleteTemporaryFiles(dir_path: str) -> None:
    """
    Remove temporary SAM files from directory
    """
    for fname in os.listdir(die_path):
        os.remove(os.path.join(dir, fname))

def getFastqPairedFiles(data_dir: str, pattern: tuple=('_1.', '_2.')) -> dict:
    """
    Group paired-end fastq files by condition
    """
    fnames = os.listdir(data_dir)
    conditions = np.unique([fname.split(pattern[0])[0].split(pattern[1])[0] 
                            for fname in fnames]).tolist()
    return {
        condition: [fname for fname in fnames if condition in fname] 
        for condition in conditions
    }
    
def sortSAMbyName(sam_file: str, output_dir=None) -> None:
    """
    Sort SAM entries by name. Required by htseq-count to process paired-end data.
    """
    if output_dir is None:
        output_dir = f'{sam_file.split(".sam")[0]}_sorted.sam'
    samtools_command = f'samtools sort -n -O sam {sam_file} > {output_dir}'
    terminalExecute(samtools_command)  
    