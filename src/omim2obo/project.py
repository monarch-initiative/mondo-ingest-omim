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
    associations = list(associations)
    mondo_ids = _objects_by_mim(_flipped(mondo_mappings), 'MONDO:')
    hgnc_ids = _objects_by_mim(hgnc_mappings, 'HGNC:')

    # QUIRK (https://github.com/monarch-initiative/omim/pull/332): a protected row's Mondo class
    # and HGNC gene are pooled with everything else known about its OMIM entries, so every gene of
    # the phenotype reaches every one of its classes, whichever the curator paired.
    for association in associations:
        if isinstance(association.provenance, CuratorProtected):
            mondo_ids[association.phenotype_mim].add(association.provenance.mondo_id)
            hgnc_ids[association.gene_mim].add(association.provenance.hgnc_id)

    return {
        MondoGeneRow(mondo_id, hgnc_id, association.phenotype_mim, association.gene_mim)
        for association in associations
        if association.predicate is Predicate.GERMLINE_MUTATION
        for mondo_id in mondo_ids[association.phenotype_mim]
        for hgnc_id in hgnc_ids[association.gene_mim]
    }


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
