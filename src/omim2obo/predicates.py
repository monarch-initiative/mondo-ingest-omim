"""The predicates a disease-gene association can take, and how morbidmap's mapping keys choose them."""
from enum import Enum

from omim2obo.parse.omim import MappingKey


class Predicate(Enum):
    GERMLINE_MUTATION = 'RO:0004003'  # disease to gene; rendered together with its inverse, RO:0004013
    CAUSES = 'RO:0003303'  # gene to disease, as are the rest
    CONTRIBUTES = 'RO:0003304'
    CAUSES_OR_CONTRIBUTES = 'RO:0003302'


# MORBIDMAP_PHENOTYPE_MAPPING_KEY_MEANINGS: The Monarch Ingest and Kevin S have identified some ECO properties that
# ...probably apply to each of these. See:
# https://github.com/monarch-initiative/monarch-ingest/blob/main/monarch_ingest/ingests/omim/omim-translation.yaml
# "1": "inference from background scientific knowledge used in manual assertion"
# "2": "genomic context evidence"
# "3": "sequencing assay evidence"
# "4": "sequencing assay evidence"
MORBIDMAP_PHENOTYPE_MAPPING_KEY_MEANINGS = {
    MappingKey.UNKNOWN_DEFECT:
         'The disorder is placed on the map based on its association with a gene, but the underlying defect is '
         'not known.',
    MappingKey.LINKAGE:
         'The disorder has been placed on the map by linkage or other statistical method; no mutation has '
         'been found.',
    MappingKey.MOLECULAR_BASIS_KNOWN: 'The molecular basis for the disorder is known; a mutation has been found in the gene.',
    MappingKey.CONTIGUOUS_GENE_SYNDROME:
         'A contiguous gene deletion or duplication syndrome, multiple genes are deleted or duplicated causing '
         'the phenotype.',
}
# todo: double check / handle edge cases if/as needed for any concerns:
#  - https://github.com/monarch-initiative/omim/issues/79#issuecomment-1319408780
# MORBIDMAP_PHENOTYPE_MAPPING_KEY_PREDICATES
# - Gene-to-Disease predicates
# - Provenance https://github.com/monarch-initiative/omim/issues/79#issuecomment-1319408780
MORBIDMAP_PHENOTYPE_MAPPING_KEY_PREDICATES = {
    MappingKey.UNKNOWN_DEFECT: None,  # association with unknown defect
    # RO:0003303 (causes condition)
    # A relationship between an entity (e.g. a genotype, genetic variation, chemical, or environmental exposure) and a
    # condition (a phenotype or disease), where the entity has some causal role for the condition.
    # https://www.ebi.ac.uk/ols/ontologies/ro/properties?iri=http://purl.obolibrary.org/obo/RO_0003303
    MappingKey.LINKAGE: Predicate.CAUSES,
    # RO:0004013 (is causal germline mutation in)
    # Relates a gene to condition, such that a mutation in this gene is sufficient to produce the condition and that
    # can be passed on to offspring[modified from orphanet].
    # https://www.ebi.ac.uk/ols/ontologies/ro/properties?iri=http://purl.obolibrary.org/obo/RO_0004013
    MappingKey.MOLECULAR_BASIS_KNOWN: Predicate.GERMLINE_MUTATION,  # RO:0004003, of which this is the inverse
    # RO:0003304 (contributes to condition)
    # A relationship between an entity (e.g. a genotype, genetic variation, chemical, or environmental exposure) and a
    # condition (a phenotype or disease), where the entity has some contributing role that influences the condition.
    # https://www.ebi.ac.uk/ols/ontologies/ro/properties?iri=http://purl.obolibrary.org/obo/RO_0003304
    MappingKey.CONTIGUOUS_GENE_SYNDROME: Predicate.CONTRIBUTES,
}

# RO:0003302 (causes or contributes to condition)
# - Note: This is used in main.py, but collecting documentation for it here.
# A relationship between an entity (e.g. a genotype, genetic variation, chemical, or environmental exposure)
# and a condition (a phenotype or disease), where the entity has some causal or contributing role that
# influences the condition.
# https://www.ebi.ac.uk/ols/ontologies/ro/properties?iri=http://purl.obolibrary.org/obo/RO_0003302
# Provenance for this decision:
# - Multiple rows, same mapping key: https://github.com/monarch-initiative/omim/issues/75
# - Multiple rows, diff mapping keys: https://github.com/monarch-initiative/omim/issues/81


def p2g_is_definitive(label: str) -> bool:
    """Is phenotype to gene association definitive?

    AKA doesn't have the special syntax asdescribed in: https://www.omim.org/help/faq#1_6

    1. Brackets, "[ ]", Non-diseases: indicate "nondiseases," mainly genetic variations that lead to apparently
    abnormal laboratory test values (e.g., dysalbuminemic euthyroidal hyperthyroxinemia).
    2. Braces, "{ }", Susceptibility: indicate mutations that contribute to susceptibility to multifactorial
    disorders (e.g., diabetes, asthma) or to susceptibility to infection (e.g., malaria).
    3. A question mark, "?", Provisionality: before the phenotype name indicates that the relationship between
    the phenotype and gene is provisional. More details about this relationship are provided in the comment field
    of the map and in the gene and phenotype OMIM entries.
    """
    return not any(label.startswith(x) for x in ['[', '{', '?'])
