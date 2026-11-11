"""Fixtures are synthetic: MIM numbers in the 9xxxxx range and invented names. OMIM assigns nothing
above the 6xxxxx range, and its real files are licensed and cannot be committed."""

import pytest

from omim2obo.parse.omim import MappingKey, MorbidmapRow, parse_morbidmap

HEADER = [
    '# Copyright (c) Example. Not OMIM data.\n',
    '# Phenotype\tGene/Locus And Other Related Symbols\tMIM Number\tCyto Location\n',
]


def row(phenotype: str, symbols='FAKE1, FK1', gene_mim='900100', cyto='1p36.33') -> str:
    return f'{phenotype}\t{symbols}\t{gene_mim}\t{cyto}\n'


def test_phenotype_with_its_own_entry():
    """The overwhelmingly common case."""
    assert parse_morbidmap(HEADER + [row('Wobbly knee syndrome, 900001 (3)')]) == [
        MorbidmapRow(
            phenotype_label='Wobbly knee syndrome',
            phenotype_mim='900001',
            mapping_key=MappingKey.MOLECULAR_BASIS_KNOWN,
            gene_symbols=('FAKE1', 'FK1'),
            gene_mim='900100',
            cyto_location='1p36.33',
        )
    ]


def test_phenotype_described_within_the_gene_entry_has_no_mim():
    """Modelled on the blood groups, e.g. `Blood group, Colton (3)` on OMIM:107776."""
    [parsed] = parse_morbidmap([row('Blood group, Invented (3)')])
    assert parsed.phenotype_label == 'Blood group, Invented'
    assert parsed.phenotype_mim is None


def test_commas_and_digits_inside_the_label_are_not_mistaken_for_the_mim():
    """Modelled on OMIM:202110, `17,20-lyase deficiency, isolated, 202110 (3)`."""
    [parsed] = parse_morbidmap([row('17,20-wobblase deficiency, isolated, 900002 (3)')])
    assert parsed.phenotype_label == '17,20-wobblase deficiency, isolated'
    assert parsed.phenotype_mim == '900002'


def test_mim_without_its_comma_is_still_the_mim():
    """Modelled on OMIM:260350, which morbidmap writes as `Pancreatic cancer, somatic 260350 (3)`."""
    [parsed] = parse_morbidmap([row('Wobbly knee, somatic 900003 (3)')])
    assert parsed.phenotype_label == 'Wobbly knee, somatic'
    assert parsed.phenotype_mim == '900003'


def test_dangling_comma_is_not_part_of_the_label():
    """Modelled on OMIM:188830, which morbidmap writes as `Adrenocortical tumor, somatic, (3)`."""
    [parsed] = parse_morbidmap([row('Wobbly knee tumor, somatic, (3)')])
    assert parsed.phenotype_label == 'Wobbly knee tumor, somatic'
    assert parsed.phenotype_mim is None


@pytest.mark.parametrize(
    'label',
    [
        '{Wobbly knee, susceptibility to}',
        '[Knee wobble level QTL 1]',
        '?Wobbly knee syndrome 2',
        '{?Wobbly knee, susceptibility to, 3}',
    ],
)
def test_markers_qualifying_the_association_stay_on_the_label(label):
    """https://omim.org/help/faq#1_6. What they mean for the association is decided downstream; the
    parser's job is to not lose them."""
    [parsed] = parse_morbidmap([row(f'{label}, 900004 (3)')])
    assert parsed.phenotype_label == label


@pytest.mark.parametrize('key', MappingKey)
def test_every_mapping_key(key):
    [parsed] = parse_morbidmap([row(f'Wobbly knee syndrome, 900001 ({key.value})')])
    assert parsed.mapping_key is key


def test_same_phenotype_and_gene_under_two_labels_gives_two_rows():
    """Modelled on OMIM:613443, an entry with an included title: morbidmap names the entry's own
    phenotype and the included one separately, each with its own mapping key."""
    parsed = parse_morbidmap(
        [
            row('Chromosome 99q deletion syndrome, 900005 (4)'),
            row('Neurodevelopmental disorder with wobbling, 900005 (3)'),
        ]
    )
    assert [(r.phenotype_label, r.mapping_key) for r in parsed] == [
        ('Chromosome 99q deletion syndrome', MappingKey.CONTIGUOUS_GENE_SYNDROME),
        ('Neurodevelopmental disorder with wobbling', MappingKey.MOLECULAR_BASIS_KNOWN),
    ]


def test_comments_before_and_after_the_data_are_skipped():
    footer = ['#\n', '# Phenotype Mapping key - Appears in parentheses after a disorder :\n']
    assert len(parse_morbidmap(HEADER + [row('Wobbly knee syndrome, 900001 (3)')] + footer)) == 1


def test_unreadable_phenotype_is_an_error():
    with pytest.raises(ValueError, match='cannot read the phenotype'):
        parse_morbidmap([row('Wobbly knee syndrome, 900001')])
