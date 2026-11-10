"""Reference data from outside OMIM's downloads: mappings to other vocabularies, citations, and HGNC.

Mappings are SSSOM's own `Mapping`, read and in time written by the sssom package.
"""
import csv
from collections.abc import Iterable
from dataclasses import dataclass
from typing import TextIO

from sssom.parsers import parse_sssom_table
from sssom_schema import Mapping


def parse_sssom(file: TextIO) -> list[Mapping]:
    """An SSSOM TSV, such as the OMIM mappings Mondo publishes."""
    return parse_sssom_table(file).to_mappings()


def parse_entry_mappings(lines: Iterable[str]) -> list[Mapping]:
    """data/mappings.tsv: the UMLS and Orphanet identifiers OMIM's API gives for each entry.

        mim     is_phenotype  date_fetched  umls_ids           orphanet_ids
        900001  True          2026-01-01    C9900001|C9900002  99001

    OMIM does not say how it arrived at these, hence the justification.
    """
    return [
        Mapping(
            subject_id=f'OMIM:{row["mim"]}',
            predicate_id='skos:exactMatch',
            object_id=f'{prefix}:{local_id}',
            mapping_justification='semapv:UnspecifiedMatching',
        )
        for row in csv.DictReader(lines, delimiter='\t')
        for prefix, column in (('UMLS', 'umls_ids'), ('Orphanet', 'orphanet_ids'))
        for local_id in _pipe_separated(row[column])
    ]


@dataclass(frozen=True, slots=True)
class PubmedReference:
    """A publication cited by an OMIM entry."""
    mim: str
    pmid: str


def parse_pubmed_refs(lines: Iterable[str]) -> list[PubmedReference]:
    """data/pubmed-refs.tsv: the references OMIM's API gives for each entry, `|`-separated."""
    return [
        PubmedReference(row['mim'], pmid)
        for row in csv.DictReader(lines, delimiter='\t')
        for pmid in _pipe_separated(row['pmid_refs'])
    ]


@dataclass(frozen=True, slots=True)
class HgncGene:
    hgnc_id: str
    """As a CURIE, e.g. `HGNC:5`."""
    symbol: str


def parse_hgnc(lines: Iterable[str]) -> list[HgncGene]:
    """HGNC's complete set, reduced to the two of its columns that the pipeline reads."""
    genes = []
    table = csv.DictReader(lines, delimiter='\t')
    for row in table:
        if not (row['hgnc_id'] and row['symbol']):
            raise ValueError(f'hgnc_complete_set.txt, line {table.line_num}: missing hgnc_id or symbol')
        genes.append(HgncGene(row['hgnc_id'], row['symbol']))
    return genes


def _pipe_separated(field: str) -> list[str]:
    return field.split('|') if field else []
