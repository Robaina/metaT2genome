#!/usr/bin/env python
# coding: utf-8
# python >= 3.6
# Semidán Robaina Estévez (srobaina@ull.edu.es)

from Bio import SeqIO
import sys


class GenomeGBK:

    def __init__(self, path_to_gbk):
        self._gbk = list(SeqIO.parse(path_to_gbk, 'genbank'))[0]

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



if __name__ == '__main__':
    """
    Convert GenBank file to GTF. Script addressed to prokaryotes
    
    call: python3 gbk2gtf.py file.gbk [chromosome_id] > file.gtf
    Optional [chromosome_id], chromosome or genome name represented in SAM file.
    """
    input_file = sys.argv[1]
    gbk = GenomeGBK(input_file)
    
    if (len(sys.argv) > 2):
        chr_id = sys.argv[2]
    else:
        chr_id = gbk.meta['organism'][0]

    for feature in gbk.features:
    #     if feature.type == 'cds':
    #         data = feature.qualifiers
    #         locus_tag = data['locus_tag'][0]
    #         gene_id = data['gene'][0]
        if feature.type == 'gene':
            locus_tag = feature.qualifiers['locus_tag'][0]
            start_pos =  feature.location.start.position
            end_pos =  feature.location.end.position
            strand = '+' if feature.location.strand == 1 else '-'

            print((f'{chr_id}\tGenebank\tgene\t{start_pos}\t{end_pos}\t.\t{strand}\t.'
                   f'\tgene_id {locus_tag} ; transcript_id {locus_tag}'))
    
    
    
    