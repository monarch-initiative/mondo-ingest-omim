"""OMIM's downloadable files. See https://omim.org/downloads"""

import re
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from enum import Enum

from omim2obo.entry_prefix import EntryPrefix


def data_lines(lines: Iterable[str]) -> Iterator[str]:
    """The lines of an OMIM file that carry data.

    Every file opens with `#` comment lines, the last of which is the column header, and several
    close with more of them documenting the fields.
    """
    for line in lines:
        line = line.rstrip('\n')
        if line and not line.startswith('#'):
            yield line


class MappingKey(Enum):
    """The phenotype mapping key: how the phenotype came to be associated with the locus."""

    UNKNOWN_DEFECT = 1
    LINKAGE = 2
    MOLECULAR_BASIS_KNOWN = 3
    CONTIGUOUS_GENE_SYNDROME = 4


@dataclass(frozen=True, slots=True, kw_only=True)
class MorbidmapRow:
    """One line of morbidmap.txt: a phenotype, and a gene or locus it has been mapped to."""

    phenotype_label: str
    """As OMIM writes it, including the markers that qualify the association: `{braces}` for a
    susceptibility, `[brackets]` for a nondisease, and a leading `?` for a provisional one."""

    phenotype_mim: str | None
    """None when the phenotype has no entry of its own and is described within the gene's."""

    mapping_key: MappingKey
    gene_symbols: tuple[str, ...]
    gene_mim: str
    cyto_location: str


# "Each Phenotype is followed by its MIM number, if different from that of the locus/gene, and then
# followed by its phenotype mapping key in parentheses."
#
#     17,20-lyase deficiency, isolated, 202110 (3)
#     Blood group, Colton (3)
#
# The commas are OMIM's convention, not a guarantee: a handful of rows lack the one before the MIM
# number, and one without a MIM number keeps the comma that would have preceded it.
_PHENOTYPE = re.compile(r'(?P<label>.*?)(?:,? (?P<mim>\d{6}))?,? \((?P<key>[1-4])\)')


def parse_morbidmap(lines: Iterable[str]) -> list[MorbidmapRow]:
    rows = []
    for line in data_lines(lines):
        phenotype, gene_symbols, gene_mim, cyto_location = line.split('\t')
        match = _PHENOTYPE.fullmatch(phenotype)
        if match is None:
            raise ValueError(f'morbidmap.txt: cannot read the phenotype in {line!r}')
        rows.append(
            MorbidmapRow(
                phenotype_label=match['label'],
                phenotype_mim=match['mim'],
                mapping_key=MappingKey(int(match['key'])),
                gene_symbols=tuple(gene_symbols.split(', ')),
                gene_mim=gene_mim,
                cyto_location=cyto_location,
            )
        )
    return rows


@dataclass(frozen=True, slots=True, kw_only=True)
class TitleRow:
    """One line of mimTitles.txt, for an entry that exists.

    Each of the three fields holds `;;`-separated titles, each optionally followed by `; SYMBOL`s.
    They are kept as OMIM writes them.
    """

    prefix: EntryPrefix
    mim: str
    preferred: str
    alternative: str
    included: str
    """Titles of other phenotypes described within this entry: 'Other entities represented in this
    entry' on omim.org."""


@dataclass(frozen=True, slots=True, kw_only=True)
class RemovedRow:
    """One line of mimTitles.txt, for an entry that no longer exists: prefix `Caret`."""

    mim: str
    moved_to: tuple[str, ...]
    """The entries that replaced it. Empty if it was removed from the database outright."""


def parse_mim_titles(lines: Iterable[str]) -> list[TitleRow | RemovedRow]:
    rows = []
    for line in data_lines(lines):
        prefix, mim, preferred, alternative, included = line.split('\t')
        if prefix == 'Caret':
            # MOVED TO 900001 | MOVED TO 900001 AND 900002 | REMOVED FROM DATABASE
            if not preferred.startswith(('MOVED TO ', 'REMOVED FROM DATABASE')):
                raise ValueError(
                    f'mimTitles.txt: cannot read what became of removed entry in {line!r}'
                )
            rows.append(RemovedRow(mim=mim, moved_to=tuple(re.findall(r'\d{6}', preferred))))
        else:
            rows.append(
                TitleRow(
                    prefix=EntryPrefix(prefix),
                    mim=mim,
                    preferred=preferred,
                    alternative=alternative,
                    included=included,
                )
            )
    return rows


class Mim2GeneEntryType(Enum):
    GENE = 'gene'
    GENE_PHENOTYPE = 'gene/phenotype'
    PHENOTYPE = 'phenotype'
    PREDOMINANTLY_PHENOTYPES = 'predominantly phenotypes'
    MOVED_REMOVED = 'moved/removed'


@dataclass(frozen=True, slots=True, kw_only=True)
class Mim2GeneRow:
    """One line of mim2gene.txt: an entry's identifiers in other gene databases."""

    mim: str
    entry_type: Mim2GeneEntryType
    entrez_id: str | None
    hgnc_symbol: str | None
    ensembl_id: str | None


def parse_mim2gene(lines: Iterable[str]) -> list[Mim2GeneRow]:
    rows = []
    for line in data_lines(lines):
        mim, entry_type, entrez_id, hgnc_symbol, ensembl_id = line.split('\t')
        rows.append(
            Mim2GeneRow(
                mim=mim,
                entry_type=Mim2GeneEntryType(entry_type),
                entrez_id=entrez_id or None,
                hgnc_symbol=hgnc_symbol or None,
                ensembl_id=ensembl_id or None,
            )
        )
    return rows


@dataclass(frozen=True, slots=True, kw_only=True)
class Genemap2Row:
    """One line of genemap2.txt, reduced to the two of its fourteen columns the pipeline reads."""

    mim: str
    approved_symbol: str | None


def parse_genemap2(lines: Iterable[str]) -> list[Genemap2Row]:
    rows = []
    for line in data_lines(lines):
        fields = line.split('\t')
        rows.append(Genemap2Row(mim=fields[5], approved_symbol=fields[8] or None))
    return rows


@dataclass(frozen=True, slots=True, kw_only=True)
class PhenotypicSeries:
    """A group of entries for clinically similar phenotypes, from phenotypicSeries.txt."""

    number: str
    """With its prefix, e.g. `PS900001`."""

    title: str
    members: tuple[str, ...]


def parse_phenotypic_series(lines: Iterable[str]) -> list[PhenotypicSeries]:
    """A series is spread over several lines: one giving its title, then one for each member.

    PS900001    Wobbly knee syndrome
    PS900001    900001    Wobbly knee syndrome 1
    PS900001    900002    Wobbly knee syndrome 2
    """
    titles: dict[str, str] = {}
    members: dict[str, list[str]] = {}
    for line in data_lines(lines):
        match line.split('\t'):
            case [number, title]:
                titles[number] = title
                members[number] = []
            case [number, mim, _phenotype]:
                members[number].append(mim)
            case _:
                raise ValueError(f'phenotypicSeries.txt: expected two or three fields in {line!r}')
    return [
        PhenotypicSeries(number=number, title=title, members=tuple(members[number]))
        for number, title in titles.items()
    ]
