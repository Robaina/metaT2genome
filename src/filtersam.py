#!/usr/bin/env python3
# coding: utf-8

import os
import re
import pysam
from .utils import terminalExecute


def extractSegmentsWithMDtag(sam_dir: str, output_dir,
                             suppress_output=False) -> None:
    """
    Use samtools to filter out segments that do not have an MD tag
    """
    # if output_dir is None:
    #     output_dir = f'{name}_only_md.sam'
    samtools_command = f'samtools view -h -d MD {sam_dir} > {output_dir}'
    terminalExecute(samtools_command, suppress_output=suppress_output)
    
def sumMatchesAndMismatches(segment):
    """
    Get total matches/mismatches from CIGAR string (M field)
    Code dictionary:
    M	BAM_CMATCH	0
    I	BAM_CINS	1
    D	BAM_CDEL	2
    N	BAM_CREF_SKIP	3
    S	BAM_CSOFT_CLIP	4
    H	BAM_CHARD_CLIP	5
    P	BAM_CPAD	6
    =	BAM_CEQUAL	7
    X	BAM_CDIFF	8
    B	BAM_CBACK	9
    """
    return sum(
        [value for (code, value) in segment.cigartuples if code == 0]
    )

def getNumberOfMatches(segment):
    """
    Get numnber of matches from alignment
    Do not consider insertion/deletion as mismatches
    """
    parsed_MD = segment.get_aligned_pairs(with_seq=True)
    return len([
        base for (read_pos, ref_pos, base) in parsed_MD 
        if ((base is not None and base.isupper()) and read_pos is not None)
    ])

def getQueryLength(segment):
    """
    Compute query length from CIGAR field corresponding to
    query sequence. 
    Following: https://stackoverflow.com/questions/39710796/infer-the-length-of-a-sequence-using-the-cigar
    
    Cigar fields which 'consume sequence': M, I, S, =, X
    
    Code dictionary:
    M	BAM_CMATCH	0
    I	BAM_CINS	1
    D	BAM_CDEL	2
    N	BAM_CREF_SKIP	3
    S	BAM_CSOFT_CLIP	4
    H	BAM_CHARD_CLIP	5
    P	BAM_CPAD	6
    =	BAM_CEQUAL	7
    X	BAM_CDIFF	8
    B	BAM_CBACK	9
    """
    codes = [0, 1, 4, 7, 8]
    return sum(
        [value for (code, value) in segment.cigartuples if code in codes]
    )

def percent_matched(segment):
    """
    Compute percentage of sequence that has been matched to reference
    """
    seq_length = getQueryLength(segment)
    n_matches = getNumberOfMatches(segment)
    return 100 * (n_matches / seq_length)
    
def percent_identity(segment):
    """
    Compute percent identity from MD tag of aligned segment.
    segment: pysam AlignedSegment object.
    """
    return 100 * (getNumberOfMatches(segment) / sumMatchesAndMismatches(segment))

def has_MD_tag(segment):
    return 'MD' in [tag for (tag, _) in segment.get_tags()]

def filterSAMbyIdentity(input_path, identity_cutoff=95, output_path=None):
    """
    Filter aligned segments in BAM or SAM file with percent identity
    equal or above identity_cutoff value.
    """
    file_ext = re.search('.(s|b)am', input_path).group()
    if output_path is None:
        output_path = (f'{input_path.split(file_ext)[0]}'
                       f'.identity_filtered_at_{identity_cutoff}{file_ext}') 
    
    samfile = pysam.AlignmentFile(input_path, 'r')
    filtered_sam = pysam.AlignmentFile(output_path, 'w', template=samfile)
    for segment in samfile:
        if (has_MD_tag(segment) and percent_identity(segment) >= identity_cutoff):
            filtered_sam.write(segment)
            
    filtered_sam.close()
    samfile.close()
    
def filterSAMbyPercentMatched(input_path, matched_cutoff=50,
                              output_path=None):
    """
    Filter aligned segments in BAM or SAM file with percent of matched
    based equal or higher than matched_cutoff. 
    
    Percent of matched bases is computed as the fraction of matches in
    the total query length.
    """
    file_ext = re.search('.(s|b)am', input_path).group()
    if output_path is None:
        output_path = (f'{input_path.split(file_ext)[0]}'
                       f'_matched_filtered_at_{matched_cutoff}{file_ext}') 
    
    samfile = pysam.AlignmentFile(input_path, 'r')
    filtered_sam = pysam.AlignmentFile(output_path, 'w', template=samfile)
    for segment in samfile:
        if (has_MD_tag(segment) and percent_matched(segment) >= matched_cutoff):
            filtered_sam.write(segment)
            
    filtered_sam.close()
    samfile.close()

def filterSAMbyReadLength(input_path: str, min_length: int = 100,
                          output_path: str = None) -> None:
    """
    Filter aligned segments by original read length (per CIGAR string) 
    """
    file_ext = re.search('.(s|b)am', input_path).group()
    if output_path is None:
        output_path = (f'{input_path.split(file_ext)[0]}'
                       f'_length_filtered_at_{min_length}{file_ext}') 
    
    samfile = pysam.AlignmentFile(input_path, 'r')
    filtered_sam = pysam.AlignmentFile(output_path, 'w', template=samfile)
    for segment in samfile:
        seq_length = getQueryLength(segment)
        if seq_length >= min_length:
            filtered_sam.write(segment)
            
    filtered_sam.close()
    samfile.close()

def filterSAMbyReadNames(input_path: str,
                         output_path: str = None,
                         query_name_txt: str = None) -> None:
    """
    Extract aligned segments whose name containedd in query_name_txt:
    samtools view -N qnames_list.txt -o filtered_output.bam input.bam
    from here: https://bioinformatics.stackexchange.com/questions/3380/how-to-subset-a-bam-by-a-list-of-qnames

    query_name_txt: text file containing a list of read names (matching names in SAM file)
    """
    if output_path is None:
        basename, ext = os.path.splitext(input_path)
        output_path = (f'{basename}'
                       f'_filtered_by_name{ext}') 
    cmd_str= (
        f"samtools view -h -N {query_name_txt} -o {output_path} {input_path}"
    )
    terminalExecute(cmd_str)