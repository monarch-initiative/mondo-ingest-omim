"""Disease-gene associations: what OMIM's morbidmap and Mondo's curators together assert."""
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass
from enum import Enum

from omim2obo.parse.curation import ExclusionRow, ProtectedRow
from omim2obo.parse.omim import MappingKey, MorbidmapRow


class Predicate(Enum):
    GERMLINE_MUTATION = 'RO:0004003'
    """has material basis in germline mutation in: disease to gene. The only predicate that says the
    gene defines the disease. Rendered together with its inverse, RO:0004013."""
    CAUSES = 'RO:0003303'
    """causes condition: gene to disease."""
    CONTRIBUTES = 'RO:0003304'
    """contributes to condition: gene to disease."""
    CAUSES_OR_CONTRIBUTES = 'RO:0003302'
    """causes or contributes to condition: gene to disease. The hedge, for when more than one gene or
    mechanism is in play. https://github.com/monarch-initiative/omim/issues/75 and /81"""


# Why each key maps as it does: https://github.com/monarch-initiative/omim/issues/79#issuecomment-1319408780
# Key 1, an association whose underlying defect is not known, asserts nothing.
NON_CAUSAL_PREDICATE = {
    MappingKey.LINKAGE: Predicate.CAUSES,
    MappingKey.CONTIGUOUS_GENE_SYNDROME: Predicate.CONTRIBUTES,
}


@dataclass(frozen=True, slots=True)
class OmimMorbidmap:
    """OMIM says so, in a row of morbidmap."""
    mapping_key: MappingKey
    phenotype_label: str


@dataclass(frozen=True, slots=True)
class CuratorProtected:
    """A curator asserted it, in protected-disease-gene.tsv."""
    orcid: str | None
    protection_type: str
    mondo_id: str
    omim: OmimMorbidmap | None
    """The morbidmap row for the same phenotype and gene, if OMIM has one."""


@dataclass(frozen=True, slots=True)
class CuratorExcluded:
    """OMIM says so, but a curator has excluded the phenotype, which weakens it to the hedge."""
    orcid: str | None
    reason: str
    omim: OmimMorbidmap


type Provenance = OmimMorbidmap | CuratorProtected | CuratorExcluded


@dataclass(frozen=True, slots=True)
class Association:
    phenotype_mim: str
    gene_mim: str
    predicate: Predicate
    provenance: Provenance


def is_definitive(phenotype_label: str) -> bool:
    """Whether OMIM presents the phenotype without qualification: not a `[nondisease]`, not a
    `{susceptibility}`, and not `?provisional`. https://omim.org/help/faq#1_6"""
    return not phenotype_label.startswith(('[', '{', '?'))


def classify(
    morbidmap: Iterable[MorbidmapRow], protected: Iterable[ProtectedRow], exclusions: Iterable[ExclusionRow],
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
            omim = OmimMorbidmap(row.mapping_key, row.phenotype_label)
            in_morbidmap.add(pair)

            if pair in protections:
                associations += [
                    Association(*pair, Predicate.GERMLINE_MUTATION, CuratorProtected(
                        protection.orcid, protection.protection_type, protection.mondo_id, omim))
                    for protection in protections[pair]
                ]
            elif row.mapping_key is MappingKey.UNKNOWN_DEFECT:
                continue
            elif exclusion := exclusion_of.get(phenotype_mim):
                associations.append(Association(*pair, Predicate.CAUSES_OR_CONTRIBUTES, CuratorExcluded(
                    exclusion.orcid, exclusion.reason, omim)))
            elif row.mapping_key is not MappingKey.MOLECULAR_BASIS_KNOWN:
                predicate = NON_CAUSAL_PREDICATE[row.mapping_key] if sole_row else Predicate.CAUSES_OR_CONTRIBUTES
                associations.append(Association(*pair, predicate, omim))
            elif sole_row and is_definitive(row.phenotype_label):
                associations.append(Association(*pair, Predicate.GERMLINE_MUTATION, omim))

    associations += [
        Association(*pair, Predicate.GERMLINE_MUTATION, CuratorProtected(
            protection.orcid, protection.protection_type, protection.mondo_id, omim=None))
        for pair, pair_protections in protections.items() if pair not in in_morbidmap
        for protection in pair_protections
    ]
    return associations
