"""Fixtures are synthetic: MIM numbers in the 9xxxxx range, phenotypes 9000xx and genes 9001xx."""

import pytest

from omim2obo.associations import (
    Association,
    CuratorExcluded,
    CuratorProtected,
    OmimMorbidmap,
    classify,
)
from omim2obo.parse.curation import ExclusionRow, ProtectedRow
from omim2obo.parse.omim import MappingKey, MorbidmapRow
from omim2obo.predicates import Predicate

KEY_1, KEY_2, KEY_3, KEY_4 = MappingKey
ORCID = 'https://orcid.org/0000-0000-0000-0001'


def row(phenotype='900001', gene='900100', key=KEY_3, label='Wobbly knee syndrome') -> MorbidmapRow:
    return MorbidmapRow(label, phenotype, key, ('FAKE',), gene, '1p36')


def protect(
    phenotype='900001', gene='900100', mondo='MONDO:9900001', type_='manually reviewed'
) -> ProtectedRow:
    return ProtectedRow(
        phenotype, gene, 'HGNC:99901', mondo, 'wobbly knee syndrome', type_, ORCID, ''
    )


def exclude(phenotype='900001') -> ExclusionRow:
    return ExclusionRow(phenotype, 'MONDO:9900001', 'wobbly knee syndrome', ORCID, 'somatic')


def predicates(associations) -> list[tuple[str, Predicate]]:
    return [(a.gene_mim, a.predicate) for a in associations]


def test_one_gene_with_known_molecular_basis_defines_the_disease():
    """The overwhelmingly common case."""
    assert classify([row()], [], []) == [
        Association(
            '900001',
            '900100',
            Predicate.GERMLINE_MUTATION,
            OmimMorbidmap(KEY_3, 'Wobbly knee syndrome'),
        )
    ]


@pytest.mark.parametrize(
    'key, predicate', [(KEY_2, Predicate.CAUSES), (KEY_4, Predicate.CONTRIBUTES)]
)
def test_one_gene_by_linkage_or_contiguous_deletion(key, predicate):
    assert predicates(classify([row(key=key)], [], [])) == [('900100', predicate)]


def test_unknown_defect_asserts_nothing():
    assert classify([row(key=KEY_1)], [], []) == []


def test_phenotype_without_its_own_entry_asserts_nothing():
    assert classify([row(phenotype=None)], [], []) == []


@pytest.mark.parametrize(
    'label', ['{Wobbly knee, susceptibility to}', '[Knee wobble QTL]', '?Wobbly knee syndrome']
)
def test_qualified_phenotype_does_not_define_a_disease(label):
    assert classify([row(label=label)], [], []) == []


def test_qualified_phenotype_still_gets_a_non_causal_predicate():
    assert predicates(
        classify([row(key=KEY_2, label='{Wobbly knee, susceptibility to}')], [], [])
    ) == [('900100', Predicate.CAUSES)]


def test_several_genes_are_hedged_and_none_defines_the_disease():
    """Key 3 rows of a phenotype with several rows produce nothing at all."""
    associations = classify(
        [row(gene='900100', key=KEY_2), row(gene='900101', key=KEY_4), row(gene='900102')], [], []
    )
    assert predicates(associations) == [
        ('900100', Predicate.CAUSES_OR_CONTRIBUTES),
        ('900101', Predicate.CAUSES_OR_CONTRIBUTES),
    ]


def test_unknown_defect_row_still_counts_towards_the_hedge():
    assert predicates(
        classify([row(gene='900100', key=KEY_1), row(gene='900101', key=KEY_2)], [], [])
    ) == [('900101', Predicate.CAUSES_OR_CONTRIBUTES)]


def test_same_gene_under_two_mapping_keys():
    """Modelled on OMIM:613443 and OMIM:608232: an included title for a contiguous deletion (key 4)
    alongside the entry's own phenotype (key 3), both on one gene. Two rows, so hedged."""
    associations = classify(
        [
            row(key=KEY_4, label='Chromosome 99q deletion syndrome'),
            row(key=KEY_3, label='Neurodevelopmental disorder with wobbling'),
        ],
        [],
        [],
    )
    assert associations == [
        Association(
            '900001',
            '900100',
            Predicate.CAUSES_OR_CONTRIBUTES,
            OmimMorbidmap(KEY_4, 'Chromosome 99q deletion syndrome'),
        )
    ]


def test_same_gene_under_two_labels_counts_as_two_rows():
    """Modelled on OMIM:202110. Nothing is asserted, though only one gene is involved."""
    assert (
        classify(
            [
                row(label='17,20-wobblase deficiency, isolated'),
                row(label='17-alpha-wobblase deficiency'),
            ],
            [],
            [],
        )
        == []
    )


def test_exclusion_weakens_every_association_of_the_phenotype():
    """Modelled on OMIM:130060 (Ehlers-Danlos syndrome, arthrochalasia type)."""
    assert classify([row()], [], [exclude()]) == [
        Association(
            '900001',
            '900100',
            Predicate.CAUSES_OR_CONTRIBUTES,
            CuratorExcluded(ORCID, 'somatic', OmimMorbidmap(KEY_3, 'Wobbly knee syndrome')),
        )
    ]


def test_exclusion_does_not_resurrect_an_unknown_defect():
    assert classify([row(key=KEY_1)], [], [exclude()]) == []


def test_digenic_disease_keeps_both_genes():
    """Modelled on OMIM:619151 (AMED syndrome). Without protection, two rows would assert nothing."""
    associations = classify(
        [row(gene='900100'), row(gene='900101')],
        [protect(gene='900100', type_='digenic'), protect(gene='900101', type_='digenic')],
        [],
    )
    assert predicates(associations) == [
        ('900100', Predicate.GERMLINE_MUTATION),
        ('900101', Predicate.GERMLINE_MUTATION),
    ]
    assert all(isinstance(a.provenance, CuratorProtected) for a in associations)


def test_included_terms_each_keep_their_own_mondo_class():
    """Modelled on OMIM:231200 (Bernard-Soulier): one entry, three genes, each belonging to a
    different Mondo class. The association records which, so nothing downstream has to guess (#332)."""
    genes_and_classes = [
        ('900100', 'MONDO:9900002'),
        ('900101', 'MONDO:9900003'),
        ('900102', 'MONDO:9900004'),
    ]
    associations = classify(
        [row(gene=gene) for gene, _ in genes_and_classes],
        [
            protect(gene=gene, mondo=mondo, type_='included entry')
            for gene, mondo in genes_and_classes
        ],
        [],
    )
    assert [(a.gene_mim, a.provenance.mondo_id) for a in associations] == genes_and_classes


def test_protection_carries_the_morbidmap_row_it_coincides_with():
    [association] = classify([row(key=KEY_2)], [protect()], [])
    assert association.provenance == CuratorProtected(
        ORCID, 'manually reviewed', 'MONDO:9900001', OmimMorbidmap(KEY_2, 'Wobbly knee syndrome')
    )


def test_protection_names_a_mondo_class_for_a_gene_omim_already_lists():
    """Modelled on OMIM:137215 (CDH1 / DGLBC): the gene is the entry's own, under key 3. What the
    protected row adds is the Mondo class for the entry's included title."""
    assert classify([row()], [protect(mondo='MONDO:9900002', type_='included entry')], []) == [
        Association(
            '900001',
            '900100',
            Predicate.GERMLINE_MUTATION,
            CuratorProtected(
                ORCID,
                'included entry',
                'MONDO:9900002',
                OmimMorbidmap(KEY_3, 'Wobbly knee syndrome'),
            ),
        )
    ]


def test_protection_of_a_gene_morbidmap_does_not_list_for_the_phenotype():
    associations = classify([row(gene='900100')], [protect(gene='900101')], [])
    assert associations == [
        Association(
            '900001',
            '900100',
            Predicate.GERMLINE_MUTATION,
            OmimMorbidmap(KEY_3, 'Wobbly knee syndrome'),
        ),
        Association(
            '900001',
            '900101',
            Predicate.GERMLINE_MUTATION,
            CuratorProtected(ORCID, 'manually reviewed', 'MONDO:9900001', omim=None),
        ),
    ]


def test_protection_overrides_exclusion_and_unknown_defect():
    associations = classify([row(key=KEY_1)], [protect()], [exclude()])
    assert predicates(associations) == [('900100', Predicate.GERMLINE_MUTATION)]


def test_pair_protected_for_two_mondo_classes_gives_two_associations():
    associations = classify(
        [row()], [protect(mondo='MONDO:9900001'), protect(mondo='MONDO:9900002')], []
    )
    assert [a.provenance.mondo_id for a in associations] == ['MONDO:9900001', 'MONDO:9900002']


def test_associations_are_hashable():
    assert len(set(classify([row(), row()], [protect()], []))) == 1
