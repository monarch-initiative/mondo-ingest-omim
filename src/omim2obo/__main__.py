"""Build the release artefacts from files already on disk.

    python -m omim2obo --omim omim-data --mondo-mappings mondo_exactmatch_omim.sssom.tsv ...

Writes omim.nt, for `robot convert` to turn into omim.owl, and the four TSVs.
"""

from argparse import ArgumentParser
from datetime import date
from pathlib import Path

from omim2obo.project import project, susceptibility_subset
from omim2obo.release import read_release
from omim2obo.render.owl import render_owl
from omim2obo.render.tsv import (
    render_mondo_genes,
    render_qc,
    render_review,
    render_susceptibility,
)
from omim2obo.review import review_cases


def main() -> None:
    parser = ArgumentParser(prog='omim2obo', description=__doc__.splitlines()[0])
    parser.add_argument('--omim', type=Path, default='data', help="OMIM's .txt downloads")
    parser.add_argument('--data', type=Path, default='data', help='curated files and caches')
    parser.add_argument(
        '--mondo-mappings',
        type=Path,
        nargs='+',
        default=[
            Path('mondo_exactmatch_omim.sssom.tsv'),
            Path('mondo_exactmatch_omimps.sssom.tsv'),
        ],
        help="Mondo's SSSOM mappings to OMIM entries and phenotypic series",
    )
    parser.add_argument('--version', default=date.today().isoformat(), help='release date')
    parser.add_argument('--out', type=Path, default='.', help='directory to write to')
    args = parser.parse_args()

    release = read_release(
        version=args.version, omim=args.omim, data=args.data, mondo_mappings=args.mondo_mappings
    )
    gene_rows = project(release.associations, release.mondo_mappings, release.hgnc_mappings)
    susceptibility = susceptibility_subset(release.morbidmap, release.mondo_mappings)

    review = review_cases(release.associations, release.morbidmap, release.entries, release.removed)

    args.out.mkdir(parents=True, exist_ok=True)
    # N-Triples, not Turtle: prefixes bound in Turtle carry through robot into omim.owl's text.
    render_owl(release).serialize(args.out / 'omim.nt', format='nt', encoding='utf-8')
    for name, text in (
        ('mondo-omim-genes.robot.tsv', render_mondo_genes(gene_rows)),
        ('disease-gene-relationships-qc.tsv', render_qc(release.associations, release.entries)),
        ('mondo-omim-susceptibility-subset.robot.tsv', render_susceptibility(susceptibility)),
        ('review.tsv', render_review(review)),
    ):
        (args.out / name).write_text(text, newline='')


if __name__ == '__main__':
    main()
