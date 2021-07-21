#!/usr/bin/env python
# coding: utf-8
# python >= 3.6
# Semidán Robaina Estévez (srobaina@ull.edu.es)

import pysam
import sys
import re
import time
import numpy as np


def sumMatchesAndMismatches(segment):
    """
    Get total matches/mismatches from CIGAR string (M field)
    """
    return sum(
        [value for (code, value) in segment.cigartuples if code == 0]
    )

def getNumberOfMatches(segment):
    """
    Get numnber of matches from alignedment
    Do not consider insertion/deletion as mismatches
    """
    parsed_MD = segment.get_aligned_pairs(with_seq=True)
    return len([
        base for (read_pos, ref_pos, base) in parsed_MD 
        if ((base is not None and base.isupper()) and read_pos is not None)
    ])
    
def percent_identity(segment):
    """
    Compute percent identity from MD tag of aligned segment.
    segment: pysam AlignedSegment object.
    """
    return 100 * (getNumberOfMatches(segment) / sumMatchesAndMismatches(segment))

def has_MD_tag(segment):
    return 'MD' in [tag for (tag, _) in segment.get_tags()]

def file_has_MD_tags(samfile):
    """
    Check if SAM file has MD tags (at least in the first few segments)
    """
    return any([has_MD_tag(seg) for seg in samfile.head(100)])

def filterSAMbyIdentity(input_path, identity_cutoff=95, output_path=None):
    """
    Filter aligned segments in BAM or SAM file with percent identity
    equal or above identity_cutoff value.
    """
    print(f'Filtering {input_path}')
    file_ext = re.search('.(s|b)am', input_path).group()
    if output_path is None:
        output_path = (f'{input_path.split(file_ext)[0]}'
                       f'.identity_filtered_at_{identity_cutoff}{file_ext}') 
    
    samfile = pysam.AlignmentFile(input_path, 'r')
    if not file_has_MD_tags(samfile):
        raise ValueError(f'MD tags not found in {file_ext[1:].upper()} file')
    filtered_sam = pysam.AlignmentFile(output_path, 'wb', template=samfile)
    for segment in samfile:
        if (has_MD_tag(segment) and percent_identity(segment) >= identity_cutoff):
            filtered_sam.write(segment)
            
    filtered_sam.close()
    samfile.close()
    
def computeSAMstatistics(input_path, identity_cutoff=95):
    """
    Some basic statistics about percent identity
    """
    file_ext = re.search('.(s|b)am', input_path).group()
    samfile = pysam.AlignmentFile(input_path, 'r')
    if not file_has_MD_tags(samfile):
        raise ValueError(f'MD tags not found in {file_ext[1:].upper()} file')
    
    total_segments, segments_without_md, segments_above_cutoff = 0, 0, 0
    identity_distribution = []
    for segment in samfile:
        total_segments += 1
        if has_MD_tag(segment):
            I = percent_identity(segment)
            if I >= identity_cutoff:
                segments_above_cutoff += 1
                identity_distribution.append(I)
        else:
            segments_without_md += 1
    samfile.close()
    out = {'n_total': total_segments, 'n_without_md': segments_without_md,
           'n_above_cutoff': segments_above_cutoff, 'Idist': identity_distribution}
    
    print(f'Identity stats for: {input_file}')
    print(f'Total segments: {out["n_above_cutoff"]}')
    print(f'% without MD: {100 * out["n_without_md"] / out["n_total"]}')
    print(f'% above cutoff: {100 * out["n_above_cutoff"] / out["n_total"]}')
    print(f'Average Identity: {np.mean(out["Idist"])}%')
    
    return out
    


# Run script: python3 filter_by_identity.py input.bam identity_cutoff [output_path]
if __name__ == "__main__":
    """
    Filter records in SAM/BAM file by given percent identity.
    Usage: 
    python3 filter_by_identity.py input.{bam|sam} identity_cutoff [output_path]
    """
    
    input_file = sys.argv[1]
    if (len(sys.argv) > 2):
        identity_cutoff = int(sys.argv[2])
    else:
        identity_cutoff = 95
    if (len(sys.argv) > 3):
        output_path = int(sys.argv[3])
    else:
        output_path = None
    
    start = time.time()
    filterSAMbyIdentity(input_file, identity_cutoff, output_path)
    # out = computeSAMstatistics(input_file, identity_cutoff)
    end = time.time()
    
    print(f'Execution time: {end - start}')