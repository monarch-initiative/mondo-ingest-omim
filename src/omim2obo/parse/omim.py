"""OMIM's downloadable files. See https://omim.org/downloads"""
import re
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from enum import Enum


def data_lines(lines: Iterable[str]) -> Iterator[str]:
    """The lines of an OMIM file that carry data.

    Every file opens with `#` comment lines, the last of which is the column header, and several close
    with more of them documenting the fields.
    """
    for line in lines:
        line = line.rstrip('\n')
        if line and not line.startswith('#'):
            yield line


class MappingKey(Enum):
    """The phenotype mapping key: how the phenotype came to be associated with the locus.

    Descriptions are OMIM's own, from the footer of morbidmap.txt.
    """
    UNKNOWN_DEFECT = 1
    """The disorder is placed on the map based on its association with a gene, but the underlying
    defect is not known."""
    LINKAGE = 2
    """The disorder has been placed on the map by linkage or other statistical method; no mutation has
    been found."""
    MOLECULAR_BASIS_KNOWN = 3
    """The molecular basis for the disorder is known; a mutation has been found in the gene."""
    CONTIGUOUS_GENE_SYNDROME = 4
    """A contiguous gene deletion or duplication syndrome, multiple genes are deleted or duplicated
    causing the phenotype."""


@dataclass(frozen=True, slots=True)
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
        rows.append(MorbidmapRow(
            phenotype_label=match['label'],
            phenotype_mim=match['mim'],
            mapping_key=MappingKey(int(match['key'])),
            gene_symbols=tuple(gene_symbols.split(', ')),
            gene_mim=gene_mim,
            cyto_location=cyto_location,
        ))
    return rows
