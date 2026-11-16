"""Disease-gene associations: what OMIM's morbidmap and Mondo's curators together assert."""

from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass

from omim2obo.parse.curation import ExclusionRow, ProtectedRow
from omim2obo.parse.omim import MappingKey, MorbidmapRow
from omim2obo.predicates import (
    MORBIDMAP_PHENOTYPE_MAPPING_KEY_PREDICATES,
    Predicate,
    p2g_is_definitive,
)


@dataclass(frozen=True, slots=True, kw_only=True)
class OmimMorbidmap:
    """OMIM says so, in a row of morbidmap."""

    mapping_key: MappingKey
    phenotype_label: str


@dataclass(frozen=True, slots=True, kw_only=True)
class CuratorProtected:
    """A curator asserted it, in protected-disease-gene.tsv."""

    orcid: str | None
    protection_type: str
    mondo_id: str
    hgnc_id: str
    omim: OmimMorbidmap | None
    """The morbidmap row for the same phenotype and gene, if OMIM has one."""


@dataclass(frozen=True, slots=True, kw_only=True)
class CuratorExcluded:
    """OMIM says so, but a curator has excluded the phenotype, which weakens it to the hedge."""

    orcid: str | None
    reason: str
    omim: OmimMorbidmap


type Provenance = OmimMorbidmap | CuratorProtected | CuratorExcluded


@dataclass(frozen=True, slots=True, kw_only=True)
class Association:
    phenotype_mim: str
    gene_mim: str
    predicate: Predicate
    provenance: Provenance


def classify(
    morbidmap: Iterable[MorbidmapRow],
    protected: Iterable[ProtectedRow],
    exclusions: Iterable[ExclusionRow],
) -> list[Association]:
    rows_by_phenotype: dict[str, list[MorbidmapRow]] = defaultdict(list)
    for row in morbidmap:
        if row.phenotype_mim:  # https://github.com/monarch-initiative/omim/issues/78
            rows_by_phenotype[row.phenotype_mim].append(row)

    protections: dict[tuple[str, str], list[ProtectedRow]] = defaultdict(list)
    for protection in protected:
        protections[protection.phenotype_mim, protection.gene_mim].append(protection)

    exclusion_of = {exclusion.phenotype_mim: exclusion for exclusion in exclusions}

    associations = []
    in_morbidmap = set()
    for phenotype_mim, rows in rows_by_phenotype.items():
        # The one decision that looks beyond its own row: a phenotype with several rows has several
        # genes or mechanisms, and none of them alone is taken to define it.
        sole_row = len(rows) == 1

        for row in rows:
            pair = phenotype_mim, row.gene_mim
            omim = OmimMorbidmap(mapping_key=row.mapping_key, phenotype_label=row.phenotype_label)
            in_morbidmap.add(pair)

            if pair in protections:
                associations += [
                    _association(pair, Predicate.GERMLINE_MUTATION, _protected(protection, omim))
                    for protection in protections[pair]
                ]
            elif row.mapping_key is MappingKey.UNKNOWN_DEFECT:
                continue
            elif exclusion := exclusion_of.get(phenotype_mim):
                associations.append(
                    _association(
                        pair,
                        Predicate.CAUSES_OR_CONTRIBUTES,
                        CuratorExcluded(orcid=exclusion.orcid, reason=exclusion.reason, omim=omim),
                    )
                )
            elif row.mapping_key is not MappingKey.MOLECULAR_BASIS_KNOWN:
                predicate = (
                    MORBIDMAP_PHENOTYPE_MAPPING_KEY_PREDICATES[row.mapping_key]
                    if sole_row
                    else Predicate.CAUSES_OR_CONTRIBUTES
                )
                associations.append(_association(pair, predicate, omim))
            elif sole_row and p2g_is_definitive(row.phenotype_label):
                associations.append(_association(pair, Predicate.GERMLINE_MUTATION, omim))

    associations += [
        _association(pair, Predicate.GERMLINE_MUTATION, _protected(protection, omim=None))
        for pair, pair_protections in protections.items()
        if pair not in in_morbidmap
        for protection in pair_protections
    ]
    return associations


def _association(
    pair: tuple[str, str], predicate: Predicate, provenance: Provenance
) -> Association:
    phenotype_mim, gene_mim = pair
    return Association(
        phenotype_mim=phenotype_mim, gene_mim=gene_mim, predicate=predicate, provenance=provenance
    )


def _protected(protection: ProtectedRow, omim: OmimMorbidmap | None) -> CuratorProtected:
    return CuratorProtected(
        orcid=protection.orcid,
        protection_type=protection.protection_type,
        mondo_id=protection.mondo_id,
        hgnc_id=protection.hgnc_id,
        omim=omim,
    )
