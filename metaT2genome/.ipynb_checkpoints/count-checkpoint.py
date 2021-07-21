from .utils import terminalExecute

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