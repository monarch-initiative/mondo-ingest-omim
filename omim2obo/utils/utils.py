"""Misc utilities"""
from typing import Dict, List, Optional, Tuple, Union

import pandas as pd
from rdflib import URIRef

from omim2obo.config import DISEASE_GENE_EXCLUSIONS_PATH, DISEASE_GENE_PROTECTED_PATH
from omim2obo.namespaces import ORCID


# todo: also in mondo-ingest. Refactor into mondolib: https://github.com/monarch-initiative/mondolib/issues/13
def remove_angle_brackets(uris: Union[str, List[str]]) -> Union[str, List[str]]:
    """Remove angle brackets from URIs, e.g.:
    <https://omim.org/entry/100050> --> https://omim.org/entry/100050"""
    str_input = isinstance(uris, str)
    uris = [uris] if str_input else uris
    uris2 = []
    for x in uris:
        x = x[1:] if x.startswith('<') else x
        x = x[:-1] if x.endswith('>') else x
        uris2.append(x)
    return uris2[0] if str_input else uris2


def get_d2g_protected(
    path=DISEASE_GENE_PROTECTED_PATH
) -> Dict[Tuple[str, str], Tuple[str, Optional[URIRef]]]:
    """Get information for manually curated disease-gene association protections.

    Protected gene associations are gene associations that we want to add even if they (no longer) appear in the OMIM source data.
    This function handles all protected types (digenic, included entry, manually reviewed, etc.), not just digenic.

    :return: Dictionary with (phenotype MIM, gene MIM) keys and (HGNC id, curator ORCID) values.
    """
    df = pd.read_csv(path, sep='\t').fillna('')
    for col in ['phenotype_mim', 'gene_mim', 'hgnc_id']:
        df[col] = df[col].apply(lambda x: x.split(':')[1])
    
    return {
        (x['phenotype_mim'], x['gene_mim']): (
            x['hgnc_id'],
            ORCID[x['orcid'].replace('https://orcid.org/', '')] if x['orcid'] else None
        )
        for x in df.to_dict(orient='records')
    }


def get_d2g_protected_rows(path=DISEASE_GENE_PROTECTED_PATH) -> List[Dict[str, str]]:
    """Get manually curated disease-gene protections, one dict per row.

    Unlike get_d2g_protected(), this keeps each row intact. A row's `mondo_id` says which Mondo class
    the protected association belongs to, and that pairing has to survive: one OMIM phenotype entry
    can cover several Mondo classes (e.g. OMIM 'included terms' that Mondo splits into gene-specific
    subtypes), and then each gene belongs to one of them rather than to all of them.

    Rows missing any of the fields needed to place the association are skipped.

    :return: List of dicts with bare 'phenotype_mim', 'gene_mim' and 'mondo_id' ids.
    """
    df = pd.read_csv(path, sep='\t').fillna('')
    cols = ['phenotype_mim', 'gene_mim', 'mondo_id']
    rows = [{col: str(x[col]).strip() for col in cols} for x in df.to_dict(orient='records')]
    return [
        {col: val.split(':')[1] for col, val in row.items()}
        for row in rows if all(row.values())
    ]


def get_d2g_exclusions_by_curator(path=DISEASE_GENE_EXCLUSIONS_PATH) -> Dict[str, Optional[URIRef]]:
    """Get disease-gene exclusions

    Situations where the pipeline logic would otherwise keep a disease-gene association, we want to exclude it.

    :return: Dict[str, str]: Phenotype MIM as keys, ORCID of curator as values
    """
    df = pd.read_csv(path, sep='\t').fillna('')
    df['omim_id'] = df['omim_id'].apply(lambda x: x.split(':')[1])
    phenotype_mim_orcid_map = {x['omim_id']: x['orcid'] for x in df.to_dict(orient='records')}
    return {k: ORCID[v] if v else None for k, v in phenotype_mim_orcid_map.items()}
