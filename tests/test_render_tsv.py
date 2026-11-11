from omim2obo.associations import Association, OmimMorbidmap
from omim2obo.entry_prefix import EntryPrefix
from omim2obo.parse.omim import MappingKey
from omim2obo.predicates import Predicate
from omim2obo.project import MondoGeneRow
from omim2obo.render.tsv import render_mondo_genes, render_qc, render_susceptibility
from omim2obo.titles import EntryTitles, TitleGroup

NONE = TitleGroup((), (), (), ())
KNEE = EntryTitles('900001', EntryPrefix.PHENOTYPE, 'wobbly "knee" syndrome', ('WKS',), NONE, NONE)
TRAIT = EntryTitles('900002', EntryPrefix.OTHER, 'knee wobbliness', (), NONE, NONE)
WBK1 = EntryTitles('900100', EntryPrefix.GENE, 'wobbly knee gene 1', ('WBK1',), NONE, NONE)


def association(phenotype, predicate) -> Association:
    provenance = OmimMorbidmap(MappingKey.MOLECULAR_BASIS_KNOWN, 'Wobbly knee syndrome')
    return Association(phenotype, '900100', predicate, provenance)


def test_mondo_genes_template():
    rows = [
        MondoGeneRow('MONDO:9900002', 'HGNC:99901', '900001', '900100'),
        MondoGeneRow('MONDO:9900001', 'HGNC:99901', '900001', '900100'),
    ]
    assert render_mondo_genes(rows) == (
        'mondo_id\thgnc_id\tomim_disease_xref\tomim_gene\n'
        "ID\tSC 'has material basis in germline mutation in' some %\t>A oboInOwl:source\t\n"
        'http://purl.obolibrary.org/obo/MONDO_9900001\thttp://identifiers.org/hgnc/99901'
        '\tOMIM:900001\thttps://omim.org/entry/900100\n'
        'http://purl.obolibrary.org/obo/MONDO_9900002\thttp://identifiers.org/hgnc/99901'
        '\tOMIM:900001\thttps://omim.org/entry/900100\n'
    )


def test_susceptibility_template():
    assert render_susceptibility({('MONDO:9900001', '900001')}) == (
        'mondo_id\tsubset\tomim_id\r\n'
        'ID\tAI oboInOwl:inSubset\t>A oboInOwl:source\r\n'
        'MONDO:9900001\thttp://purl.obolibrary.org/obo/mondo#omim_susceptibility\tOMIM:900001\r\n'
    )


def test_qc_states_germline_mutation_from_both_ends():
    qc = render_qc([association('900001', Predicate.GERMLINE_MUTATION)], [KNEE, WBK1])
    assert qc.splitlines() == [
        '?SubId\t?SubLabel\t?SubBiolink\t?Pred\t?PredDirection\t?ObjId\t?ObjLabel\t?ObjBiolink',
        '"OMIM:900001"\t"wobbly \\"knee\\" syndrome"\t"Disease"\t"RO:0004003"\t"Disease2Gene"'
        '\t"OMIM:900100"\t"WBK1"\t"Gene"',
        '"OMIM:900100"\t"WBK1"\t"Gene"\t"RO:0004013"\t"Gene2Disease"'
        '\t"OMIM:900001"\t"wobbly \\"knee\\" syndrome"\t"Disease"',
    ]


def test_qc_states_other_predicates_from_the_gene_only():
    qc = render_qc([association('900002', Predicate.CAUSES)], [TRAIT, WBK1])
    assert qc.splitlines()[1:] == [
        '"OMIM:900100"\t"WBK1"\t"Gene"\t"RO:0003303"\t"Gene2Disease"'
        '\t"OMIM:900002"\t"knee wobbliness"\t'
    ]


def test_qc_omits_relationships_to_entries_without_titles():
    assert len(render_qc([association('900001', Predicate.CAUSES)], [WBK1]).splitlines()) == 1
