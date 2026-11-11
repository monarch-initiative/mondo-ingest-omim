"""Fixtures are synthetic: identifiers in ranges their namespaces do not use, and invented names."""

import io

import pytest

from omim2obo.parse.reference import (
    HgncGene,
    Mapping,
    PubmedReference,
    parse_entry_mappings,
    parse_hgnc,
    parse_known_capitalizations,
    parse_pubmed_refs,
    parse_sssom,
)


def mapping(subject_id, object_id, justification, **labels) -> Mapping:
    return Mapping(subject_id, 'skos:exactMatch', object_id, justification, **labels)


def test_sssom_metadata_header_is_skipped_and_labels_kept():
    assert parse_sssom(
        io.StringIO(
            ''.join(
                [
                    '# curie_map:\n',
                    '#   MONDO: http://purl.obolibrary.org/obo/MONDO_\n',
                    '#   OMIM: https://omim.org/entry/\n',
                    '#   OMIMPS: https://omim.org/phenotypicSeries/PS\n',
                    '# mapping_set_id: http://example.org/mappings.sssom.tsv\n',
                    'subject_id\tsubject_label\tpredicate_id\tobject_id\tobject_label\tmapping_justification\n',
                    'MONDO:9900001\twobbly knee syndrome\tskos:exactMatch\tOMIM:900001\tWobbly knee syndrome'
                    '\tsemapv:ManualMappingCuration\n',
                    'MONDO:9900002\t\tskos:exactMatch\tOMIMPS:900001\t\tsemapv:ManualMappingCuration\n',
                ]
            )
        )
    ) == [
        mapping(
            'MONDO:9900001',
            'OMIM:900001',
            'semapv:ManualMappingCuration',
            subject_label='wobbly knee syndrome',
            object_label='Wobbly knee syndrome',
        ),
        mapping('MONDO:9900002', 'OMIMPS:900001', 'semapv:ManualMappingCuration'),
    ]


def test_mappings_are_hashable():
    assert (
        len(
            {
                mapping('OMIM:900001', 'UMLS:C9900001', 'semapv:UnspecifiedMatching')
                for _ in range(2)
            }
        )
        == 1
    )


def test_entry_mappings_are_one_per_identifier():
    assert parse_entry_mappings(
        [
            'mim\tis_phenotype\tdate_fetched\tumls_ids\torphanet_ids\n',
            '900001\tTrue\t2026-01-01\tC9900001|C9900002\t99001\n',
            '900002\tTrue\t2026-01-01\t\t\n',
            '900100\tFalse\t2026-01-01\t\t99002\n',
        ]
    ) == [
        mapping('OMIM:900001', 'UMLS:C9900001', 'semapv:UnspecifiedMatching'),
        mapping('OMIM:900001', 'UMLS:C9900002', 'semapv:UnspecifiedMatching'),
        mapping('OMIM:900001', 'Orphanet:99001', 'semapv:UnspecifiedMatching'),
        mapping('OMIM:900100', 'Orphanet:99002', 'semapv:UnspecifiedMatching'),
    ]


def test_pubmed_refs_are_one_per_publication():
    assert parse_pubmed_refs(
        [
            'mim\tis_phenotype\tdate_fetched\tpmid_refs\n',
            '900001\tTrue\t2026-01-01\t99000001|99000002\n',
            '900002\tTrue\t2026-01-01\t\n',
        ]
    ) == [PubmedReference('900001', '99000001'), PubmedReference('900001', '99000002')]


HGNC_HEADER = 'hgnc_id\tsymbol\tname\tlocus_group\n'


def test_hgnc_keeps_id_and_symbol():
    assert parse_hgnc(
        [
            HGNC_HEADER,
            'HGNC:99901\tFAKE1\tfake gene 1\tprotein-coding gene\n',
            'HGNC:99902\tFAKE2\t"fake gene 2, with a comma"\tprotein-coding gene\n',
        ]
    ) == [HgncGene('HGNC:99901', 'FAKE1'), HgncGene('HGNC:99902', 'FAKE2')]


def test_hgnc_gene_without_a_symbol_is_an_error():
    with pytest.raises(ValueError, match='line 2: missing hgnc_id or symbol'):
        parse_hgnc([HGNC_HEADER, 'HGNC:99901\t\tfake gene 1\tprotein-coding gene\n'])


def test_known_capitalizations():
    assert parse_known_capitalizations(
        [
            'lower_name\tcap_name\tpattern\n',
            'wobbleton\tWobbleton\t%[ ]wobbleton[ ]%\n',
        ]
    ) == {'wobbleton': 'Wobbleton'}


def test_known_capitalizations_QUIRK_stray_tab_shifts_the_columns():
    """Modelled on the `de lange` row, which is why OMIM:122470 is "cornelia lange lange syndrome"."""
    assert parse_known_capitalizations(
        [
            'lower_name\tcap_name\tpattern\n',
            'von\twobbleton\tvon Wobbleton %[ ]von wobbleton[ ]%\n',
        ]
    ) == {'von': 'wobbleton'}
