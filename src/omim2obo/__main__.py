"""Build OMIM's release artefacts for Mondo."""

from datetime import date
from pathlib import Path

import click

from omim2obo.project import project, susceptibility_subset
from omim2obo.release import Inputs, read_release
from omim2obo.render.owl import render_owl
from omim2obo.render.tsv import (
    render_mondo_genes,
    render_qc,
    render_review,
    render_susceptibility,
)
from omim2obo.review import review_cases


@click.group(help=__doc__)
def main() -> None:
    pass


@main.command()
@click.argument('what', type=click.Choice(['files', 'cache']), required=False)
def fetch(what: str | None) -> None:
    """Download OMIM's files and refresh the API cache, into data/. Needs the keys in .env.

    Give WHAT to do only one of the two.
    """
    # Imported here because it reads .env on import
    from omim2obo.fetch import download_mim_files, update_cache__pubmed_refs_and_mappings

    if what != 'cache':
        download_mim_files()
    if what != 'files':
        update_cache__pubmed_refs_and_mappings()


def input_file(flag: str, default: str, help: str, **kwargs):
    path = click.Path(exists=True, dir_okay=False, path_type=Path)
    return click.option(flag, type=path, default=default, show_default=True, help=help, **kwargs)


@main.command()
@input_file('--mim-titles', 'data/mimTitles.txt', 'OMIM: every entry, with its titles')
@input_file('--morbidmap', 'data/morbidmap.txt', 'OMIM: phenotypes and their genes')
@input_file('--mim2gene', 'data/mim2gene.txt', "OMIM: entries' gene identifiers elsewhere")
@input_file('--genemap2', 'data/genemap2.txt', 'OMIM: the gene map, for HGNC symbols')
@input_file('--phenotypic-series', 'data/phenotypicSeries.txt', 'OMIM: series and their members')
@input_file(
    '--protected',
    'data/protected-disease-gene.tsv',
    'curated: associations to assert whatever OMIM says',
)
@input_file('--exclusions', 'data/exclusions-disease-gene.tsv', 'curated: phenotypes to hold back')
@input_file(
    '--known-capitalizations',
    'data/known_capitalizations.tsv',
    'curated: words that keep their capitals in titles',
)
@input_file('--hgnc', 'data/hgnc/hgnc_complete_set.txt', "HGNC's complete gene set")
@input_file('--pubmed-refs', 'data/pubmed-refs.tsv', "cache of OMIM's API: references")
@input_file(
    '--entry-mappings', 'data/mappings.tsv', "cache of OMIM's API: UMLS and Orphanet mappings"
)
@input_file(
    '--mondo-mappings',
    ('mondo_exactmatch_omim.sssom.tsv', 'mondo_exactmatch_omimps.sssom.tsv'),
    "Mondo's SSSOM mappings to OMIM entries and to phenotypic series; give once for each file",
    multiple=True,
)
@click.option(
    '--version', default=lambda: date.today().isoformat(), help='release date  [default: today]'
)
@click.option(
    '--out',
    type=click.Path(file_okay=False, path_type=Path),
    default='.',
    show_default=True,
    help='directory to write to',
)
def build(version: str, out: Path, **inputs) -> None:
    """Build from files already on disk: omim.nt (for `robot convert`) and four TSVs."""
    release = read_release(version, Inputs(**inputs))
    gene_rows = project(release.associations, release.mondo_mappings, release.hgnc_mappings)
    susceptibility = susceptibility_subset(release.morbidmap, release.mondo_mappings)
    review = review_cases(release.associations, release.morbidmap, release.entries, release.removed)

    out.mkdir(parents=True, exist_ok=True)
    # N-Triples, not Turtle: prefixes bound in Turtle carry through robot into omim.owl's text.
    render_owl(release).serialize(out / 'omim.nt', format='nt', encoding='utf-8')
    for name, text in (
        ('mondo-omim-genes.robot.tsv', render_mondo_genes(gene_rows)),
        ('disease-gene-relationships-qc.tsv', render_qc(release.associations, release.entries)),
        ('mondo-omim-susceptibility-subset.robot.tsv', render_susceptibility(susceptibility)),
        ('review.tsv', render_review(review)),
    ):
        (out / name).write_text(text, newline='')


if __name__ == '__main__':
    main()
