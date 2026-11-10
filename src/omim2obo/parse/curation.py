"""The files Mondo's curators maintain in data/, saying where OMIM's associations should be overridden.

These are edited by hand, so unlike OMIM's own files they are checked field by field, and an error
names the row a curator needs to look at.
"""
import csv
import re
from collections.abc import Iterable, Iterator
from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class ProtectedRow:
    """One line of protected-disease-gene.tsv: a disease-gene association a curator has asserted,
    to be kept whatever OMIM's files say."""
    phenotype_mim: str
    gene_mim: str
    hgnc_id: str
    mondo_id: str
    """The Mondo class the association is meant for. Not always the one the OMIM entry maps to: for
    an entry with included titles, it may be the class for one of those."""
    mondo_label: str
    protection_type: str
    """The curator's reason, from a loose vocabulary: `digenic`, `included entry`, `manually
    reviewed`, …"""
    orcid: str | None
    comment: str


@dataclass(frozen=True, slots=True)
class ExclusionRow:
    """One line of exclusions-disease-gene.tsv: a phenotype for which a curator has rejected the
    disease-gene associations the pipeline would otherwise assert."""
    phenotype_mim: str
    mondo_id: str
    mondo_label: str
    orcid: str | None
    reason: str


def parse_protected(lines: Iterable[str]) -> list[ProtectedRow]:
    return [
        ProtectedRow(
            phenotype_mim=_mim(row, 'phenotype_mim'),
            gene_mim=_mim(row, 'gene_mim'),
            hgnc_id=_curie(row, 'hgnc_id', _HGNC),
            mondo_id=_curie(row, 'mondo_id', _MONDO),
            mondo_label=row['mondo_label'],
            protection_type=row['type'],
            orcid=row['orcid'] or None,
            comment=row['comment'],
        )
        for row in _rows('protected-disease-gene.tsv', lines)
    ]


def parse_exclusions(lines: Iterable[str]) -> list[ExclusionRow]:
    return [
        ExclusionRow(
            phenotype_mim=_mim(row, 'omim_id'),
            mondo_id=_curie(row, 'mondo_id', _MONDO),
            mondo_label=row['mondo_label'],
            orcid=row['orcid'] or None,
            reason=row['exclusion_reason_comment'],
        )
        for row in _rows('exclusions-disease-gene.tsv', lines)
    ]


_OMIM = re.compile(r'OMIM:(\d{6})')
_MONDO = re.compile(r'MONDO:\d{7}')
_HGNC = re.compile(r'HGNC:\d+')


class _Row(dict[str, str]):
    """A row that remembers where it came from, for error messages."""
    where: str


def _rows(file_name: str, lines: Iterable[str]) -> Iterator[_Row]:
    reader = csv.DictReader(lines, delimiter='\t', strict=True)
    for fields in reader:
        if None in fields or None in fields.values():
            raise ValueError(f'{file_name}, line {reader.line_num}: wrong number of columns')
        row = _Row(fields)
        row.where = f'{file_name}, line {reader.line_num}'
        yield row


def _curie(row: _Row, column: str, pattern: re.Pattern[str]) -> str:
    value = row[column]
    if not pattern.fullmatch(value):
        raise ValueError(f'{row.where}: {column} is {value!r}, expected something like {pattern.pattern!r}')
    return value


def _mim(row: _Row, column: str) -> str:
    """The MIM number of an `OMIM:` CURIE."""
    return _curie(row, column, _OMIM).removeprefix('OMIM:')
