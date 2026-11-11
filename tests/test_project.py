from omim2obo.associations import Association, CuratorProtected, OmimMorbidmap
from omim2obo.parse.omim import MappingKey, MorbidmapRow
from omim2obo.parse.reference import Mapping
from omim2obo.predicates import Predicate
from omim2obo.project import MondoGeneRow, project, susceptibility_subset

OMIM_SAYS = OmimMorbidmap(MappingKey.MOLECULAR_BASIS_KNOWN, 'Wobbly knee syndrome')


def mondo(mondo_id, mim) -> Mapping:
    return Mapping(mondo_id, 'skos:exactMatch', f'OMIM:{mim}', 'semapv:ManualMappingCuration')


def hgnc(mim, hgnc_id) -> Mapping:
    return Mapping(f'OMIM:{mim}', 'skos:exactMatch', hgnc_id, 'semapv:UnspecifiedMatching')


def germline(gene='900100', provenance=OMIM_SAYS) -> Association:
    return Association('900001', gene, Predicate.GERMLINE_MUTATION, provenance)


def protected(mondo_id, hgnc_id) -> CuratorProtected:
    return CuratorProtected(None, 'included entry', mondo_id, hgnc_id, omim=None)


MONDO = [mondo('MONDO:9900001', '900001')]
HGNC = [hgnc('900100', 'HGNC:99901')]


def test_association_lands_on_the_mapped_class_and_gene():
    assert project([germline()], MONDO, HGNC) == {
        MondoGeneRow('MONDO:9900001', 'HGNC:99901', '900001', '900100')
    }


def test_only_germline_mutation_is_projected():
    hedge = Association('900001', '900100', Predicate.CAUSES_OR_CONTRIBUTES, OMIM_SAYS)
    assert project([hedge], MONDO, HGNC) == set()


def test_unmapped_phenotype_or_gene_is_dropped():
    assert project([germline()], [], HGNC) == set()
    assert project([germline()], MONDO, []) == set()


def test_protected_row_supplies_its_own_class_and_gene():
    association = germline(provenance=protected('MONDO:9900002', 'HGNC:99901'))
    assert project([association], [], []) == {
        MondoGeneRow('MONDO:9900002', 'HGNC:99901', '900001', '900100')
    }


def test_QUIRK_protected_class_receives_every_gene_of_the_phenotype():
    """#332, modelled on OMIM:231200 (Bernard-Soulier): each included-entry class was paired with
    one gene, and gets them all."""
    associations = [
        germline('900100', protected('MONDO:9900002', 'HGNC:99901')),
        germline('900101', protected('MONDO:9900003', 'HGNC:99902')),
    ]
    assert {(row.mondo_id, row.hgnc_id) for row in project(associations, [], [])} == {
        ('MONDO:9900002', 'HGNC:99901'),
        ('MONDO:9900002', 'HGNC:99902'),
        ('MONDO:9900003', 'HGNC:99901'),
        ('MONDO:9900003', 'HGNC:99902'),
    }


def morbidmap(label, phenotype='900001') -> MorbidmapRow:
    return MorbidmapRow(
        label, phenotype, MappingKey.MOLECULAR_BASIS_KNOWN, ('WBK1',), '900100', '1p36'
    )


def test_susceptibility_is_a_label_in_braces():
    rows = [morbidmap('{Wobbly knee, susceptibility to}'), morbidmap('Wobbly knee', '900002')]
    mappings = [*MONDO, mondo('MONDO:9900002', '900002')]
    assert susceptibility_subset(rows, mappings) == {('MONDO:9900001', '900001')}
