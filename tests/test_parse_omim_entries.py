"""Fixtures are synthetic: MIM numbers in the 9xxxxx range and invented names."""

import pytest

from omim2obo.parse.omim import (
    EntryPrefix,
    Genemap2Row,
    Mim2GeneEntryType,
    Mim2GeneRow,
    PhenotypicSeries,
    RemovedRow,
    TitleRow,
    parse_genemap2,
    parse_mim2gene,
    parse_mim_titles,
    parse_phenotypic_series,
)


def titles(prefix, mim, preferred, alternative='', included='') -> str:
    return f'{prefix}\t{mim}\t{preferred}\t{alternative}\t{included}\n'


def test_title_fields_are_kept_as_written():
    """Modelled on OMIM:231200 (Bernard-Soulier): an entry whose included titles name the
    gene-specific subtypes that Mondo models as separate classes."""
    line = titles(
        'Number Sign',
        '900001',
        'WOBBLY KNEE SYNDROME; WKS',
        alternative='KNEE WOBBLE;; GENU INSTABILE; GI',
        included='WOBBLY KNEE SYNDROME, TYPE A1, INCLUDED;; WOBBLY KNEE SYNDROME, TYPE B, INCLUDED',
    )
    assert parse_mim_titles([line]) == [
        TitleRow(
            prefix=EntryPrefix.PHENOTYPE,
            mim='900001',
            preferred='WOBBLY KNEE SYNDROME; WKS',
            alternative='KNEE WOBBLE;; GENU INSTABILE; GI',
            included=(
                'WOBBLY KNEE SYNDROME, TYPE A1, INCLUDED;; WOBBLY KNEE SYNDROME, TYPE B, INCLUDED'
            ),
        )
    ]


@pytest.mark.parametrize('prefix', EntryPrefix)
def test_every_prefix_of_a_live_entry(prefix):
    [row] = parse_mim_titles([titles(prefix.value, '900001', 'SOMETHING')])
    assert row.prefix is prefix


@pytest.mark.parametrize(
    'preferred, moved_to',
    [
        ('MOVED TO 900002', ('900002',)),
        ('MOVED TO 900002 AND 900003', ('900002', '900003')),
        ('MOVED TO 900002, 900003, AND 900004', ('900002', '900003', '900004')),
        ('MOVED TO {900002}', ('900002',)),
        ('REMOVED FROM DATABASE', ()),
    ],
)
def test_removed_entry_says_where_it_went(preferred, moved_to):
    """Each shape here occurs in the real file."""
    assert parse_mim_titles([titles('Caret', '900001', preferred)]) == [
        RemovedRow(mim='900001', moved_to=moved_to)
    ]


def test_removed_entry_with_an_unknown_fate_is_an_error():
    with pytest.raises(ValueError, match='what became of'):
        parse_mim_titles([titles('Caret', '900001', 'GONE FISHING')])


def test_unknown_prefix_is_an_error():
    with pytest.raises(ValueError):
        parse_mim_titles([titles('Ampersand', '900001', 'SOMETHING')])


def test_mim2gene_absent_identifiers_are_none():
    assert parse_mim2gene(
        [
            '# MIM Number\tMIM Entry Type\tEntrez Gene ID\tApproved Gene Symbol\tEnsembl Gene ID\n',
            '900100\tgene\t99901\tFAKE1\tENSG00000999001\n',
            '900001\tphenotype\t\t\t\n',
        ]
    ) == [
        Mim2GeneRow(
            mim='900100',
            entry_type=Mim2GeneEntryType.GENE,
            entrez_id='99901',
            hgnc_symbol='FAKE1',
            ensembl_id='ENSG00000999001',
        ),
        Mim2GeneRow(
            mim='900001',
            entry_type=Mim2GeneEntryType.PHENOTYPE,
            entrez_id=None,
            hgnc_symbol=None,
            ensembl_id=None,
        ),
    ]


def test_genemap2_keeps_mim_and_approved_symbol():
    def line(mim, approved_symbol):
        return (
            '\t'.join(
                [
                    'chr1',
                    '1',
                    '1000',
                    '1p36',
                    '',
                    mim,
                    'FAKE1, FK1',
                    'Fake gene 1',
                    approved_symbol,
                    '99901',
                    'ENSG00000999001',
                    '',
                    'Wobbly knee syndrome, 900001 (3), Autosomal dominant',
                    '',
                ]
            )
            + '\n'
        )

    assert parse_genemap2([line('900100', 'FAKE1'), line('900101', '')]) == [
        Genemap2Row(mim='900100', approved_symbol='FAKE1'),
        Genemap2Row(mim='900101', approved_symbol=None),
    ]


def test_phenotypic_series_gathers_its_members():
    assert parse_phenotypic_series(
        [
            '# Phenotypic Series Number\tMIM Number\tPhenotype\n',
            'PS900001\tWobbly knee syndrome\n',
            'PS900001\t900001\tWobbly knee syndrome 1\n',
            'PS900001\t900002\tWobbly knee syndrome 2\n',
            'PS900010\tKnee wobble, isolated\n',
            'PS900010\t900010\tKnee wobble, isolated, 1\n',
        ]
    ) == [
        PhenotypicSeries(
            number='PS900001', title='Wobbly knee syndrome', members=('900001', '900002')
        ),
        PhenotypicSeries(number='PS900010', title='Knee wobble, isolated', members=('900010',)),
    ]


def test_phenotypic_series_member_is_kept_once_however_often_it_is_listed():
    """As in the real file, where an entry appears once for each of its phenotypes in the series."""
    [series] = parse_phenotypic_series(
        [
            'PS900001\tWobbly knee syndrome\n',
            'PS900001\t900001\tWobbly knee syndrome 1\n',
            'PS900001\t900001\t{Wobbly knee, susceptibility to}\n',
        ]
    )
    assert series.members == ('900001',)


def test_phenotypic_series_member_may_come_before_its_title():
    [series] = parse_phenotypic_series(
        ['PS900001\t900001\tWobbly knee syndrome 1\n', 'PS900001\tWobbly knee syndrome\n']
    )
    assert series.members == ('900001',)


@pytest.mark.parametrize(
    'lines',
    [
        ['PS900001\t900001\tWobbly knee syndrome 1\n'],
        ['PS900001\t900001\tWobbly knee syndrome 1\textra\n'],
    ],
    ids=['no title', 'four fields'],
)
def test_phenotypic_series_malformed_file_is_an_error(lines):
    with pytest.raises(ValueError, match='phenotypicSeries.txt'):
        parse_phenotypic_series(lines)
