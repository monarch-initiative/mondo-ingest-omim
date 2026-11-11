"""The whole build, on the smallest inputs that exercise every file it reads."""

from click.testing import CliRunner

from omim2obo.__main__ import main

OMIM = {
    'mimTitles.txt': (
        'Number Sign\t900001\tWOBBLY KNEE SYNDROME; WKS\t\t\n'
        'Number Sign\t900002\tWOBBLY KNEE, SUSCEPTIBILITY TO\t\t\n'
        'Asterisk\t900100\tWOBBLY KNEE GENE 1; WBK1\t\t\n'
    ),
    'mim2gene.txt': '900001\tphenotype\t\t\t\n900100\tgene\t99\tWBK1\tENSG00000000099\n',
    'genemap2.txt': 'chr1\t1\t2\t1p36\t\t900100\tWBK1\tWobbly knee gene 1\tWBK1\t99\t\t\t\t\n',
    'morbidmap.txt': (
        'Wobbly knee syndrome, 900001 (3)\tWBK1\t900100\t1p36\n'
        '{Wobbly knee, susceptibility to}, 900002 (3)\tWBK1\t900100\t1p36\n'
    ),
    'phenotypicSeries.txt': 'PS900001\tWobbly knee\nPS900001\t900001\tWobbly knee syndrome\n',
}
DATA = {
    'known_capitalizations.tsv': 'lower_name\tcap_name\tpattern\n',
    'hgnc/hgnc_complete_set.txt': 'hgnc_id\tsymbol\nHGNC:99901\tWBK1\n',
    'protected-disease-gene.tsv': (
        'phenotype_mim\tmondo_id\tmondo_label\ttype\tgene_mim\thgnc_id\torcid\tcomment\n'
    ),
    'exclusions-disease-gene.tsv': (
        'omim_id\tmondo_id\tmondo_label\torcid\texclusion_reason_comment\n'
    ),
    'pubmed-refs.tsv': 'mim\tis_phenotype\tdate_fetched\tpmid_refs\n900001\tTrue\t2026-01-01\t1\n',
    'mappings.tsv': (
        'mim\tis_phenotype\tdate_fetched\tumls_ids\torphanet_ids\n'
        '900001\tTrue\t2026-01-01\tC9900001\t\n'
    ),
}
MONDO = (
    'subject_id\tpredicate_id\tobject_id\tmapping_justification\n'
    'MONDO:9900001\tskos:exactMatch\tOMIM:900001\tsemapv:ManualMappingCuration\n'
    'MONDO:9900002\tskos:exactMatch\tOMIM:900002\tsemapv:ManualMappingCuration\n'
)


def test_build(tmp_path):
    for directory, files in (('omim', OMIM), ('data', DATA), ('.', {'mondo.sssom.tsv': MONDO})):
        for name, text in files.items():
            path = tmp_path / directory / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text)
    options = {
        '--mim-titles': 'omim/mimTitles.txt',
        '--morbidmap': 'omim/morbidmap.txt',
        '--mim2gene': 'omim/mim2gene.txt',
        '--genemap2': 'omim/genemap2.txt',
        '--phenotypic-series': 'omim/phenotypicSeries.txt',
        '--protected': 'data/protected-disease-gene.tsv',
        '--exclusions': 'data/exclusions-disease-gene.tsv',
        '--known-capitalizations': 'data/known_capitalizations.tsv',
        '--hgnc': 'data/hgnc/hgnc_complete_set.txt',
        '--pubmed-refs': 'data/pubmed-refs.tsv',
        '--entry-mappings': 'data/mappings.tsv',
        '--mondo-mappings': 'mondo.sssom.tsv',
        '--out': 'out',
        '--version': '2026-01-01',
    }
    arguments = []
    for flag, value in options.items():
        arguments += [flag, value if flag == '--version' else str(tmp_path / value)]

    result = CliRunner().invoke(main, ['build', *arguments], catch_exceptions=False)

    assert result.exit_code == 0, result.output
    assert (tmp_path / 'out/mondo-omim-genes.robot.tsv').read_text().splitlines()[2:] == [
        'http://purl.obolibrary.org/obo/MONDO_9900001\thttp://identifiers.org/hgnc/99901'
        '\tOMIM:900001\thttps://omim.org/entry/900100'
    ]
    susceptibility = tmp_path / 'out/mondo-omim-susceptibility-subset.robot.tsv'
    assert susceptibility.read_text().count('MONDO:99') == 1
    assert len((tmp_path / 'out/disease-gene-relationships-qc.tsv').read_text().splitlines()) == 3
    triples = (tmp_path / 'out/omim.nt').read_text()
    assert '<https://omim.org/entry/900001> <http://www.w3.org/2000/01/rdf-schema#label>' in triples
