from omim2obo.genes import approved_symbols, hgnc_mappings
from omim2obo.parse.omim import Genemap2Row, Mim2GeneEntryType, Mim2GeneRow
from omim2obo.parse.reference import HgncGene, Mapping


def mim2gene(mim, symbol) -> Mim2GeneRow:
    return Mim2GeneRow(mim, Mim2GeneEntryType.GENE, None, symbol, None)


def test_symbol_comes_from_either_file():
    symbols = approved_symbols(
        [mim2gene('900100', 'WBK1'), mim2gene('900101', None)],
        [Genemap2Row('900100', 'WBK1'), Genemap2Row('900101', 'WBK2')],
    )
    assert symbols == {'900100': 'WBK1', '900101': 'WBK2'}


def test_entry_the_files_disagree_on_has_no_symbol():
    assert approved_symbols([mim2gene('900100', 'WBK1')], [Genemap2Row('900100', 'WBK1A')]) == {}


def test_symbol_is_resolved_to_its_hgnc_id():
    mappings = hgnc_mappings(
        {'900100': 'WBK1', '900101': 'NOT-IN-HGNC'}, [HgncGene('HGNC:99901', 'WBK1')]
    )
    assert mappings == [
        Mapping('OMIM:900100', 'skos:exactMatch', 'HGNC:99901', 'semapv:UnspecifiedMatching')
    ]
