"""Associations that were made, but look odd enough for a curator to check. See README.md."""

from collections.abc import Iterable
from dataclasses import dataclass

from omim2obo.associations import Association, OmimMorbidmap
from omim2obo.entry_prefix import EntryPrefix
from omim2obo.parse.omim import MorbidmapRow, RemovedRow
from omim2obo.predicates import Predicate
from omim2obo.titles import EntryTitles

REVIEW_CASE_NAMES = {
    1: 'D2G: digenic',
    2: 'D2G: self-referential',
    3: 'D2G: somatic',
    4: 'D2G: Phenotype is gene',
    5: 'D2G: Phenotype type error',
}
# What review.tsv has always called the entry types that a phenotype should not have.
UNEXPECTED_PHENOTYPE_TYPES = {
    EntryPrefix.GENE: (4, 'GENE'),
    EntryPrefix.GENE_AND_PHENOTYPE: (5, 'HAS_AFFECTED_FEATURE'),
    EntryPrefix.OTHER: (5, 'SUSPECTED'),
    None: (5, 'OBSOLETE'),
}


@dataclass(frozen=True, slots=True, kw_only=True)
class ReviewCase:
    code: int
    value: str


def review_cases(
    associations: Iterable[Association],
    morbidmap: Iterable[MorbidmapRow],
    entries: Iterable[EntryTitles],
    removed: Iterable[RemovedRow],
) -> list[ReviewCase]:
    """Only associations OMIM alone accounts for, and only the disease-defining ones."""
    morbidmap = list(morbidmap)
    prefix_of: dict[str, EntryPrefix | None] = {entry.mim: entry.prefix for entry in entries}
    prefix_of |= {row.mim: None for row in removed}

    # A phenotype's own entry sometimes appears in morbidmap as a gene, for a phenotype with no
    # MIM number: the entry is being described as its own cause.
    self_references: dict[str, list[MorbidmapRow]] = {}
    for row in morbidmap:
        if not row.phenotype_mim:
            self_references.setdefault(row.gene_mim, []).append(row)

    cases = []
    self_referential = 0
    for association in _in_legacy_order(associations, morbidmap):
        if association.predicate is not Predicate.GERMLINE_MUTATION:
            continue
        if not isinstance(association.provenance, OmimMorbidmap):
            continue
        phenotype_mim, label = association.phenotype_mim, association.provenance.phenotype_label
        key = association.provenance.mapping_key.value
        # QUIRK: the ", " that separated the label from its MIM number in morbidmap is kept.
        info = (
            f'(Phenotype: {phenotype_mim} {label}, ), (Map key: {key}), '
            f'(Gene: {association.gene_mim})'
        )

        if 'digenic' in label.lower():
            cases.append(ReviewCase(code=1, value=info))
        if 'somatic' in label.lower():
            cases.append(ReviewCase(code=3, value=info))
        if rows := self_references.get(phenotype_mim):
            self_referential += 1
            cases.append(ReviewCase(code=2, value=f'{self_referential}: {info}'))
            cases += [
                ReviewCase(
                    code=2,
                    value=f'{self_referential}: (Phenotype: {row.phenotype_label}), '
                    f'(Map key: {row.mapping_key.value}), (Gene: {phenotype_mim})',
                )
                for row in rows
            ]
        if phenotype_mim in prefix_of and prefix_of[phenotype_mim] in UNEXPECTED_PHENOTYPE_TYPES:
            code, type_name = UNEXPECTED_PHENOTYPE_TYPES[prefix_of[phenotype_mim]]
            cases.append(ReviewCase(code=code, value=f'(Phenotype MIM type {type_name}), {info}'))
    return sorted(cases, key=lambda case: (case.code, case.value))


def _in_legacy_order(
    associations: Iterable[Association], morbidmap: list[MorbidmapRow]
) -> list[Association]:
    """QUIRK: self-referential cases are numbered in the order the legacy build met them, which
    was gene by gene (in order of each gene's first row), not row by row."""
    first_row_of_gene: dict[str, int] = {}
    row_of_pair: dict[tuple[str | None, str], int] = {}
    for index, row in enumerate(morbidmap):
        first_row_of_gene.setdefault(row.gene_mim, index)
        row_of_pair.setdefault((row.phenotype_mim, row.gene_mim), index)
    return sorted(
        (a for a in associations if (a.phenotype_mim, a.gene_mim) in row_of_pair),
        key=lambda a: (first_row_of_gene[a.gene_mim], row_of_pair[a.phenotype_mim, a.gene_mim]),
    )
