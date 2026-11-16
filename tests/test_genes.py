from omim2obo.genes import approved_symbols, hgnc_mappings
from omim2obo.parse.omim import Genemap2Row, Mim2GeneEntryType, Mim2GeneRow
from omim2obo.parse.reference import HgncGene, Mapping


def mim2gene(mim, symbol) -> Mim2GeneRow:
    return Mim2GeneRow(
        mim=mim,
        entry_type=Mim2GeneEntryType.GENE,
        entrez_id=None,
        hgnc_symbol=symbol,
        ensembl_id=None,
    )


def test_symbol_comes_from_either_file():
    symbols = approved_symbols(
        [mim2gene('900100', 'WBK1'), mim2gene('900101', None)],
        [
            Genemap2Row(mim='900100', approved_symbol='WBK1'),
            Genemap2Row(mim='900101', approved_symbol='WBK2'),
        ],
    )
    assert symbols == {'900100': 'WBK1', '900101': 'WBK2'}


def test_entry_the_files_disagree_on_has_no_symbol():
    assert (
        approved_symbols(
            [mim2gene('900100', 'WBK1')], [Genemap2Row(mim='900100', approved_symbol='WBK1A')]
        )
        == {}
    )


def test_symbol_is_resolved_to_its_hgnc_id():
    mappings = hgnc_mappings(
        {'900100': 'WBK1', '900101': 'NOT-IN-HGNC'}, [HgncGene(hgnc_id='HGNC:99901', symbol='WBK1')]
    )
    assert mappings == [
        Mapping(
            subject_id='OMIM:900100',
            predicate_id='skos:exactMatch',
            object_id='HGNC:99901',
            mapping_justification='semapv:UnspecifiedMatching',
        )
    ]
