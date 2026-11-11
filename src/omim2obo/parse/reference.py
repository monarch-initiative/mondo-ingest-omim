"""Reference data from outside OMIM's downloads: mappings to other vocabularies, citations, and HGNC.

SSSOM files are read by the sssom package, which knows the format; what comes out is our own record.
"""
import csv
from collections.abc import Iterable
from dataclasses import dataclass
from typing import TextIO

from sssom.parsers import parse_sssom_table


@dataclass(frozen=True, slots=True)
class Mapping:
    """A statement that two identifiers, given as CURIEs, denote the same thing. Fields are SSSOM's."""
    subject_id: str
    predicate_id: str
    object_id: str
    mapping_justification: str
    subject_label: str | None = None
    object_label: str | None = None


def parse_sssom(file: TextIO) -> list[Mapping]:
    """An SSSOM TSV, such as the OMIM mappings Mondo publishes."""
    return [
        Mapping(
            subject_id=str(m.subject_id),
            predicate_id=str(m.predicate_id),
            object_id=str(m.object_id),
            mapping_justification=str(m.mapping_justification),
            subject_label=m.subject_label and str(m.subject_label),
            object_label=m.object_label and str(m.object_label),
        )
        for m in parse_sssom_table(file).to_mappings()
    ]


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


def parse_known_capitalizations(lines: Iterable[str]) -> dict[str, str]:
    """data/known_capitalizations.tsv: words whose capitals survive a title being lowercased.

    Only the first two columns are read, the lowercase word and its replacement.
    """
    _header, *rows = csv.reader(lines, delimiter='\t')
    return {row[0]: row[1] for row in rows}


def _pipe_separated(field: str) -> list[str]:
    return field.split('|') if field else []
