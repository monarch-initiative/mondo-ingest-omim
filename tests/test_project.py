from omim2obo.associations import Association, CuratorProtected, OmimMorbidmap
from omim2obo.parse.omim import MappingKey, MorbidmapRow
from omim2obo.parse.reference import Mapping
from omim2obo.predicates import Predicate
from omim2obo.project import MondoGeneRow, project, susceptibility_subset

OMIM_SAYS = OmimMorbidmap(
    mapping_key=MappingKey.MOLECULAR_BASIS_KNOWN, phenotype_label='Wobbly knee syndrome'
)


def mondo(mondo_id, mim) -> Mapping:
    return Mapping(
        subject_id=mondo_id,
        predicate_id='skos:exactMatch',
        object_id=f'OMIM:{mim}',
        mapping_justification='semapv:ManualMappingCuration',
    )


def hgnc(mim, hgnc_id) -> Mapping:
    return Mapping(
        subject_id=f'OMIM:{mim}',
        predicate_id='skos:exactMatch',
        object_id=hgnc_id,
        mapping_justification='semapv:UnspecifiedMatching',
    )


def germline(gene='900100', provenance=OMIM_SAYS) -> Association:
    return Association(
        phenotype_mim='900001',
        gene_mim=gene,
        predicate=Predicate.GERMLINE_MUTATION,
        provenance=provenance,
    )


def protected(mondo_id, hgnc_id) -> CuratorProtected:
    return CuratorProtected(
        orcid=None, protection_type='included entry', mondo_id=mondo_id, hgnc_id=hgnc_id, omim=None
    )


MONDO = [mondo('MONDO:9900001', '900001')]
HGNC = [hgnc('900100', 'HGNC:99901')]


def test_association_lands_on_the_mapped_class_and_gene():
    assert project([germline()], MONDO, HGNC) == {
        MondoGeneRow(
            mondo_id='MONDO:9900001',
            hgnc_id='HGNC:99901',
            phenotype_mim='900001',
            gene_mim='900100',
        )
    }


def test_only_germline_mutation_is_projected():
    hedge = Association(
        phenotype_mim='900001',
        gene_mim='900100',
        predicate=Predicate.CAUSES_OR_CONTRIBUTES,
        provenance=OMIM_SAYS,
    )
    assert project([hedge], MONDO, HGNC) == set()


def test_unmapped_phenotype_or_gene_is_dropped():
    assert project([germline()], [], HGNC) == set()
    assert project([germline()], MONDO, []) == set()


def test_protected_row_supplies_its_own_class_and_gene():
    association = germline(provenance=protected('MONDO:9900002', 'HGNC:99901'))
    assert project([association], [], []) == {
        MondoGeneRow(
            mondo_id='MONDO:9900002',
            hgnc_id='HGNC:99901',
            phenotype_mim='900001',
            gene_mim='900100',
        )
    }


def test_protected_class_receives_only_the_gene_it_was_paired_with():
    """#332, modelled on OMIM:231200 (Bernard-Soulier): three included-entry classes, one gene
    each."""
    associations = [
        germline('900100', protected('MONDO:9900002', 'HGNC:99901')),
        germline('900101', protected('MONDO:9900003', 'HGNC:99902')),
    ]
    assert {(row.mondo_id, row.hgnc_id) for row in project(associations, [], [])} == {
        ('MONDO:9900002', 'HGNC:99901'),
        ('MONDO:9900003', 'HGNC:99902'),
    }


def test_protected_association_still_reaches_the_mapped_class():
    """Modelled on OMIM:137215 (CDH1): the entry's own class keeps the gene."""
    association = germline(provenance=protected('MONDO:9900002', 'HGNC:99901'))
    assert {row.mondo_id for row in project([association], MONDO, HGNC)} == {
        'MONDO:9900001',
        'MONDO:9900002',
    }


def morbidmap(label, phenotype='900001') -> MorbidmapRow:
    return MorbidmapRow(
        phenotype_label=label,
        phenotype_mim=phenotype,
        mapping_key=MappingKey.MOLECULAR_BASIS_KNOWN,
        gene_symbols=('WBK1',),
        gene_mim='900100',
        cyto_location='1p36',
    )


def test_susceptibility_is_a_label_in_braces():
    rows = [morbidmap('{Wobbly knee, susceptibility to}'), morbidmap('Wobbly knee', '900002')]
    mappings = [*MONDO, mondo('MONDO:9900002', '900002')]
    assert susceptibility_subset(rows, mappings) == {('MONDO:9900001', '900001')}
