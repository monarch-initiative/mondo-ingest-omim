"""Fixtures are synthetic: identifiers in ranges their namespaces do not use, and invented names."""
import pytest

from omim2obo.parse.reference import (
    HgncGene, Mapping, PubmedReference, parse_entry_mappings, parse_hgnc, parse_pubmed_refs, parse_sssom)


def test_sssom_metadata_header_is_skipped_and_labels_kept():
    assert parse_sssom([
        '# curie_map:\n',
        '#   MONDO: http://purl.obolibrary.org/obo/MONDO_\n',
        '# mapping_set_id: http://example.org/mappings.sssom.tsv\n',
        'subject_id\tsubject_label\tpredicate_id\tobject_id\tobject_label\tmapping_justification\n',
        'MONDO:9900001\twobbly knee syndrome\tskos:exactMatch\tOMIM:900001\tWobbly knee syndrome'
        '\tsemapv:ManualMappingCuration\n',
        'MONDO:9900002\t\tskos:exactMatch\tOMIMPS:900001\t\tsemapv:ManualMappingCuration\n',
    ]) == [
        Mapping('MONDO:9900001', 'skos:exactMatch', 'OMIM:900001', 'wobbly knee syndrome',
                'Wobbly knee syndrome', 'semapv:ManualMappingCuration'),
        Mapping('MONDO:9900002', 'skos:exactMatch', 'OMIMPS:900001', None, None,
                'semapv:ManualMappingCuration'),
    ]


def test_entry_mappings_are_one_per_identifier():
    assert parse_entry_mappings([
        'mim\tis_phenotype\tdate_fetched\tumls_ids\torphanet_ids\n',
        '900001\tTrue\t2026-01-01\tC9900001|C9900002\t99001\n',
        '900002\tTrue\t2026-01-01\t\t\n',
        '900100\tFalse\t2026-01-01\t\t99002\n',
    ]) == [
        Mapping('OMIM:900001', 'skos:exactMatch', 'UMLS:C9900001'),
        Mapping('OMIM:900001', 'skos:exactMatch', 'UMLS:C9900002'),
        Mapping('OMIM:900001', 'skos:exactMatch', 'Orphanet:99001'),
        Mapping('OMIM:900100', 'skos:exactMatch', 'Orphanet:99002'),
    ]


def test_pubmed_refs_are_one_per_publication():
    assert parse_pubmed_refs([
        'mim\tis_phenotype\tdate_fetched\tpmid_refs\n',
        '900001\tTrue\t2026-01-01\t99000001|99000002\n',
        '900002\tTrue\t2026-01-01\t\n',
    ]) == [PubmedReference('900001', '99000001'), PubmedReference('900001', '99000002')]


HGNC_HEADER = 'hgnc_id\tsymbol\tname\tlocus_group\n'


def test_hgnc_keeps_id_and_symbol():
    assert parse_hgnc([
        HGNC_HEADER,
        'HGNC:99901\tFAKE1\tfake gene 1\tprotein-coding gene\n',
        'HGNC:99902\tFAKE2\t"fake gene 2, with a comma"\tprotein-coding gene\n',
    ]) == [HgncGene('HGNC:99901', 'FAKE1'), HgncGene('HGNC:99902', 'FAKE2')]


def test_hgnc_gene_without_a_symbol_is_an_error():
    with pytest.raises(ValueError, match='line 2: missing hgnc_id or symbol'):
        parse_hgnc([HGNC_HEADER, 'HGNC:99901\t\tfake gene 1\tprotein-coding gene\n'])
