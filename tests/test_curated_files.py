"""The curated files as committed. A row the build would reject fails here, on the pull request
that adds it."""

from pathlib import Path

from omim2obo.parse.curation import parse_exclusions, parse_protected

DATA = Path(__file__).parent.parent / 'data'


def test_protected_file_parses():
    with open(DATA / 'protected-disease-gene.tsv') as file:
        assert parse_protected(file)


def test_exclusions_file_parses():
    with open(DATA / 'exclusions-disease-gene.tsv') as file:
        assert parse_exclusions(file)
