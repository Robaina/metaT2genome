"""
Functions for general purposes
"""
import os
import shutil
import random
import string
import numpy as np


def terminalExecute(command_str: str, suppress_output=False) -> None:
    """
    Execute given command in terminal through Python
    """
    if suppress_output:
        suppress_code = '>/dev/null 2>&1'
        command_str = f'{command_str} {suppress_code}'
    os.system(command_str)
    
def deleteTemporaryFiles(dir_path: str) -> None:
    """
    Remove temporary SAM files from directory
    """
    for fname in os.listdir(dir_path):
        os.remove(os.path.join(dir_path, fname))

def getFastqPairedFiles(data_dir: str, pattern: tuple=('_1.', '_2.')) -> dict:
    """
    Group paired-end fastq files by condition
    """
    fnames = os.listdir(data_dir)
    conditions = np.unique([fname.split(pattern[0])[0].split(pattern[1])[0] 
                            for fname in fnames]).tolist()
    return {
        condition: np.sort(
            [fname for fname in fnames if condition in fname]
        ).tolist() 
        for condition in conditions
    }

def filterFastqByReadNames(path_to_fastq: str, path_to_names: str, output_file: str = None) -> None:
    """
    Filter (paired) fastq files by list of read names (using seqkit)
    """
    if output_file is None:
        output_file = setDefaultOutputPath(path_to_fastq, tag="_filtered")
    cmd_str = f"seqkit grep -f {path_to_names} -o {output_file} {path_to_fastq}"
    terminalExecute(cmd_str)

def sortSAMbyName(sam_file: str, output_dir=None, suppress_output=False) -> None:
    """
    Sort SAM entries by name. Required by htseq-count to process paired-end data.
    """
    if output_dir is None:
        output_dir = f'{sam_file.split(".sam")[0]}_sorted.sam'
    samtools_command = f'samtools sort -n -O sam {sam_file} > {output_dir}'
    terminalExecute(samtools_command, suppress_output=suppress_output)

def setDefaultOutputPath(input_path: str, tag: str = None,
                         extension: str = None,
                         only_filename: bool = False,
                         only_dirname: bool = False) -> str:
    """
    Get default path to output file or directory
    """
    basename = os.path.basename(input_path)
    dirname = os.path.abspath(os.path.dirname(input_path))
    fname, ext = os.path.splitext(basename)
    if extension is None:
        extension = ext
    if tag is None:
        tag = ''
    default_file = f'{fname}{tag}{extension}'
    if only_filename:
        return default_file
    if only_dirname:
        return os.path.abspath(dirname)
    else:
        return os.path.abspath(os.path.join(dirname, default_file))

class TemporaryFilePath:
    """
    Custom context manager to create a temporary file
    which is removed when exiting context manager
    """
    def __init__(self,
                 work_dir: str = None,
                 extension: str = None,
                 create_file: bool = False):
        self.work_dir = work_dir or ''
        self.extension = extension or ''
        self.create_file = create_file

    def __enter__(self):
        temp_id = ''.join(
        random.choice(string.ascii_lowercase) for i in range(10)
        )
        self.file_path = os.path.join(
            self.work_dir, f'temp_{temp_id}{self.extension}'
            )
        if self.create_file:
            os.mkdir(self.file_path)
        return self.file_path

    def __exit__(self, exc_type, exc_val, exc_tb):
        if os.path.exists(self.file_path):
            os.remove(self.file_path)
            
class TemporaryDirectoryPath:
    """
    Custom context manager to create a temporary directory
    which is removed when exiting context manager
    """
    def __init__(self, work_dir: str = None):
        self.work_dir = work_dir or ''

    def __enter__(self):
        temp_id = ''.join(
        random.choice(string.ascii_lowercase) for i in range(10)
        )
        self.dir_path = os.path.join(
            self.work_dir, f'temp_{temp_id}/'
            )
        os.mkdir(self.dir_path)
        return self.dir_path

    def __exit__(self, exc_type, exc_val, exc_tb):
        if os.path.exists(self.dir_path):
            shutil.rmtree(self.dir_path)
    