"""
Functions and classes to deal with GeneBank files
"""

from Bio import SeqIO
import os


class GenomeGBK:

    def __init__(self, path_to_gbk):
        self._gbk = next(SeqIO.parse(path_to_gbk, 'genbank'))

    @property
    def meta(self):
        return dict(self._gbk.features[0].qualifiers)

    @property
    def features(self):
        return [f for f in self._gbk.features[1:]]

    def getGeneInfo(self, gene_id: str):
        try:
            gene, cds = [f for f in self._gbk.features[1:]
                         if gene_id.lower() in f.qualifiers['locus_tag'][0].lower()]
            res = dict(cds.qualifiers)
            res.update({'location': gene.location})
            return res
        except Exception:
            raise ValueError(f'Gene {gene_id} not found in GBK')

    def has_EC_number(self, gene_id: str):
        return 'EC_number' in self.getGeneInfo(gene_id).keys()

    def getEnzymeGene(self, ec_number: str):
        try:
            return [f.qualifiers['locus_tag'][0] for f in self._gbk.features[1:]
                    if ('EC_number' in f.qualifiers.keys()
                        and ec_number in f.qualifiers['EC_number'][0])]
        except Exception:
            raise ValueError(f'Enzyme {ec_number} not found in GBK')

def convertGBKtoGTFfile(gbk_file: str, chr_id: str=None,
                        output_dir: str=None) -> None:
    """
    Convert gbk file format to simpler GTF format,
    addressed to single chromosome organism (specified with
    chr_id or taken as species name in gbk)
    """
    gbk = GenomeGBK(gbk_file) 
    if chr_id is None:
        chr_id = gbk.meta['organism'][0]
    if output_dir is None:
        output_dir = f'{os.path.splitext(gbk_file)[0]}.gtf'
        
    with open(output_dir, 'w') as gtf:
        for feature in gbk.features:
            if feature.type == 'gene':
                locus_tag = feature.qualifiers['locus_tag'][0]
                start_pos =  feature.location.start.position
                end_pos =  feature.location.end.position
                strand = '+' if feature.location.strand == 1 else '-'

                line = (f'{chr_id}\tGenebank\tgene\t{start_pos}\t{end_pos}\t.'
                        f'\t{strand}\t.\tgene_id {locus_tag} ; transcript_id {locus_tag}')
                gtf.write(f'{line}\n')
                
def getGeneLengthsFromGBK(gbk_file: str, feature_type: str='gene') -> dict:
    """
    Obtain gene lengths in bp for genes in gbk file
    """
    gbk = GenomeGBK(gbk_file)            
    return {
        feature.qualifiers['locus_tag'][0]: (feature.location.end.position - 
                                             feature.location.start.position)
        for feature in gbk.features if feature.type.lower() in feature_type.lower()
    }

def extractLocusTagInfoFromGBK(path_to_gbk: str,
                               fields: list = ['gene', 'product', 'EC_number'],
                               feature_type: str = 'cds') -> dict:
    """
    Assign info to locus tags
    """
    def assign_field_info(field):
        return feature.qualifiers[field][0] if field in feature.qualifiers else 'unspecified'

    gbk = GenomeGBK(path_to_gbk)
    locus_tag_info = {}
    for feature in gbk.features:
        if feature.type.lower() in feature_type.lower():
            locus_tag_info[feature.qualifiers['locus_tag'][0]] = {
                field: assign_field_info(field)
                for field in fields
            }
    return locus_tag_info