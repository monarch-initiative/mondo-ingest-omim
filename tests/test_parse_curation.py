"""Fixtures are synthetic: MIM numbers in the 9xxxxx range, Mondo and HGNC ids in ranges those
namespaces do not use, and invented names."""

import pytest

from omim2obo.parse.curation import ExclusionRow, ProtectedRow, parse_exclusions, parse_protected

ORCID = 'https://orcid.org/0000-0000-0000-0001'
PROTECTED_HEADER = 'phenotype_mim\tmondo_id\tmondo_label\ttype\tgene_mim\thgnc_id\torcid\tcomment\n'
EXCLUSIONS_HEADER = 'omim_id\tmondo_id\tmondo_label\torcid\texclusion_reason_comment\n'


def protected(
    phenotype='OMIM:900001',
    mondo='MONDO:9900001',
    label='wobbly knee syndrome',
    type_='digenic',
    gene='OMIM:900100',
    hgnc='HGNC:99901',
    orcid=ORCID,
    comment='',
) -> str:
    return '\t'.join([phenotype, mondo, label, type_, gene, hgnc, orcid, comment]) + '\n'


def test_digenic_disease_has_a_row_for_each_gene():
    """Modelled on OMIM:619151 (AMED syndrome): two genes, one disease, both required."""
    assert parse_protected(
        [
            PROTECTED_HEADER,
            protected(label='"WOBBLY KNEE SYNDROME, DIGENIC"'),
            protected(
                label='"WOBBLY KNEE SYNDROME, DIGENIC"', gene='OMIM:900101', hgnc='HGNC:99902'
            ),
        ]
    ) == [
        ProtectedRow(
            '900001',
            '900100',
            'HGNC:99901',
            'MONDO:9900001',
            'WOBBLY KNEE SYNDROME, DIGENIC',
            'digenic',
            ORCID,
            '',
        ),
        ProtectedRow(
            '900001',
            '900101',
            'HGNC:99902',
            'MONDO:9900001',
            'WOBBLY KNEE SYNDROME, DIGENIC',
            'digenic',
            ORCID,
            '',
        ),
    ]


def test_included_entry_rows_name_a_mondo_class_other_than_the_entrys_own():
    """Modelled on OMIM:231200 (Bernard-Soulier): one entry, whose included titles Mondo models as
    gene-specific subtypes. Each row pairs the same phenotype MIM with another gene and class."""
    rows = parse_protected(
        [
            PROTECTED_HEADER,
            protected(
                mondo='MONDO:9900002', type_='included entry', gene='OMIM:900100', hgnc='HGNC:99901'
            ),
            protected(
                mondo='MONDO:9900003', type_='included entry', gene='OMIM:900101', hgnc='HGNC:99902'
            ),
        ]
    )
    assert [(r.phenotype_mim, r.gene_mim, r.mondo_id) for r in rows] == [
        ('900001', '900100', 'MONDO:9900002'),
        ('900001', '900101', 'MONDO:9900003'),
    ]


def test_blank_orcid_is_none():
    [row] = parse_protected([PROTECTED_HEADER, protected(orcid='')])
    assert row.orcid is None


def test_exclusion():
    """Modelled on OMIM:130060 (Ehlers-Danlos syndrome, arthrochalasia type)."""
    assert parse_exclusions(
        [
            EXCLUSIONS_HEADER,
            f'OMIM:900001\tMONDO:9900001\twobbly knee syndrome\t{ORCID}'
            '\t"caused by a somatic mutation, so not germline"\n',
        ]
    ) == [
        ExclusionRow(
            '900001',
            'MONDO:9900001',
            'wobbly knee syndrome',
            ORCID,
            'caused by a somatic mutation, so not germline',
        ),
    ]


@pytest.mark.parametrize(
    'bad_row, complaint',
    [
        (protected(phenotype='900001'), "line 2: phenotype_mim is '900001'"),
        (protected(gene='OMIM:9001'), "line 2: gene_mim is 'OMIM:9001'"),
        (protected(hgnc='FAKE1'), "line 2: hgnc_id is 'FAKE1'"),
        (protected(mondo='MONDO_9900001'), "line 2: mondo_id is 'MONDO_9900001'"),
        (protected().replace('\t', ' ', 1), 'line 2: wrong number of columns'),
        (protected().rstrip('\n') + '\textra\n', 'line 2: wrong number of columns'),
    ],
)
def test_malformed_row_names_the_line_to_fix(bad_row, complaint):
    with pytest.raises(ValueError, match=complaint):
        parse_protected([PROTECTED_HEADER, bad_row])
