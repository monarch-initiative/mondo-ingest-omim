"""The one place a CURIE prefix is given its IRI."""

from curies import Converter

OBO = 'http://purl.obolibrary.org/obo/'

CONVERTER = Converter.from_prefix_map(
    {
        'OMIM': 'https://omim.org/entry/',
        'OMIMPS': 'https://omim.org/phenotypicSeries/PS',
        'MONDO': OBO + 'MONDO_',
        'RO': OBO + 'RO_',
        'CHR': OBO + 'CHR_',
        'HGNC': 'http://identifiers.org/hgnc/',
        'HGNC.symbol': 'https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/',
        'NCBIGene': 'https://www.ncbi.nlm.nih.gov/gene/',
        'PMID': 'http://www.ncbi.nlm.nih.gov/pubmed/',
        'UMLS': 'http://linkedlifedata.com/resource/umls/id/',
        'Orphanet': 'http://www.orpha.net/ORDO/Orphanet_',
        'skos': 'http://www.w3.org/2004/02/skos/core#',
        'semapv': 'https://w3id.org/semapv/vocab/',
    }
)

expand = CONVERTER.expand_strict
