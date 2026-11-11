"""Projecting OMIM's entries onto Mondo's classes."""

from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass

from omim2obo.associations import Association, CuratorProtected
from omim2obo.parse.omim import MorbidmapRow
from omim2obo.parse.reference import Mapping
from omim2obo.predicates import Predicate


@dataclass(frozen=True, slots=True)
class MondoGeneRow:
    """A Mondo class has its material basis in germline mutation in an HGNC gene."""

    mondo_id: str
    hgnc_id: str
    phenotype_mim: str
    gene_mim: str


def project(
    associations: Iterable[Association],
    mondo_mappings: Iterable[Mapping],
    hgnc_mappings: Iterable[Mapping],
) -> set[MondoGeneRow]:
    mondo_ids = _objects_by_mim(_flipped(mondo_mappings), 'MONDO:')
    hgnc_ids = _objects_by_mim(hgnc_mappings, 'HGNC:')

    rows = set()
    for association in associations:
        if association.predicate is not Predicate.GERMLINE_MUTATION:
            continue
        phenotype_mim, gene_mim, provenance = (
            association.phenotype_mim,
            association.gene_mim,
            association.provenance,
        )
        # A curator pairs one class with one gene. That is as well as, not instead of, whatever
        # the entry and its gene are mapped to.
        if isinstance(provenance, CuratorProtected):
            rows.add(MondoGeneRow(provenance.mondo_id, provenance.hgnc_id, phenotype_mim, gene_mim))
        rows |= {
            MondoGeneRow(mondo_id, hgnc_id, phenotype_mim, gene_mim)
            for mondo_id in mondo_ids[phenotype_mim]
            for hgnc_id in hgnc_ids[gene_mim]
        }
    return rows


def susceptibility_subset(
    morbidmap: Iterable[MorbidmapRow], mondo_mappings: Iterable[Mapping]
) -> set[tuple[str, str]]:
    """(Mondo class, phenotype MIM) wherever OMIM puts the phenotype in braces, its notation for
    susceptibility to a multifactorial disorder or to infection (https://omim.org/help/faq#1_6)."""
    mondo_ids = _objects_by_mim(_flipped(mondo_mappings), 'MONDO:')
    return {
        (mondo_id, row.phenotype_mim)
        for row in morbidmap
        if row.phenotype_mim and row.phenotype_label.startswith('{')
        for mondo_id in mondo_ids[row.phenotype_mim]
    }


def _flipped(mappings: Iterable[Mapping]) -> list[Mapping]:
    """Mondo publishes its mappings with the Mondo class as subject."""
    return [
        Mapping(m.object_id, m.predicate_id, m.subject_id, m.mapping_justification)
        for m in mappings
    ]


def _objects_by_mim(mappings: Iterable[Mapping], object_prefix: str) -> dict[str, set[str]]:
    objects = defaultdict(set)
    for mapping in mappings:
        if (
            mapping.predicate_id == 'skos:exactMatch'
            and mapping.subject_id.startswith('OMIM:')
            and mapping.object_id.startswith(object_prefix)
        ):
            objects[mapping.subject_id.removeprefix('OMIM:')].add(mapping.object_id)
    return objects
