"""The released TSVs, byte for byte as the SPARQL-and-pandas route used to write them."""

import csv
import io
from collections.abc import Iterable

from omim2obo.associations import Association
from omim2obo.entry_prefix import BIOLINK_CATEGORY
from omim2obo.predicates import Predicate
from omim2obo.project import MondoGeneRow
from omim2obo.titles import EntryTitles

CAUSAL_GERMLINE_MUTATION_IN = 'RO:0004013'
"""Inverse of `Predicate.GERMLINE_MUTATION`, stated from the gene."""


def render_mondo_genes(rows: Iterable[MondoGeneRow]) -> str:
    """mondo-omim-genes.robot.tsv, a ROBOT template."""
    lines = sorted(
        (
            row.mondo_id.replace('MONDO:', 'http://purl.obolibrary.org/obo/MONDO_'),
            row.hgnc_id.replace('HGNC:', 'http://identifiers.org/hgnc/'),
            f'https://omim.org/entry/{row.gene_mim}',
            f'OMIM:{row.phenotype_mim}',
        )
        for row in rows
    )
    return _tsv(
        [
            ('mondo_id', 'hgnc_id', 'omim_disease_xref', 'omim_gene'),
            (
                'ID',
                "SC 'has material basis in germline mutation in' some %",
                '>A oboInOwl:source',
                '',
            ),
            *((mondo, hgnc, phenotype, gene) for mondo, hgnc, gene, phenotype in lines),
        ],
        lineterminator='\n',
    )


def render_susceptibility(subset: Iterable[tuple[str, str]]) -> str:
    """mondo-omim-susceptibility-subset.robot.tsv, a ROBOT template."""
    return _tsv(
        [
            ('mondo_id', 'subset', 'omim_id'),
            ('ID', 'AI oboInOwl:inSubset', '>A oboInOwl:source'),
            *(
                (
                    mondo_id,
                    'http://purl.obolibrary.org/obo/mondo#omim_susceptibility',
                    f'OMIM:{mim}',
                )
                for mondo_id, mim in sorted(subset)
            ),
        ],
        lineterminator='\r\n',
    )


def render_qc(associations: Iterable[Association], entries: Iterable[EntryTitles]) -> str:
    """disease-gene-relationships-qc.tsv: every relationship between two labelled entries, in both
    directions where the predicate has an inverse."""
    entry_of = {entry.mim: entry for entry in entries}
    relationships = set()
    for a in associations:
        if a.phenotype_mim not in entry_of or a.gene_mim not in entry_of:
            continue
        phenotype, gene = entry_of[a.phenotype_mim], entry_of[a.gene_mim]
        if a.predicate is Predicate.GERMLINE_MUTATION:
            relationships.add((phenotype, a.predicate.value, 'Disease2Gene', gene))
            relationships.add((gene, CAUSAL_GERMLINE_MUTATION_IN, 'Gene2Disease', phenotype))
        else:
            relationships.add((gene, a.predicate.value, 'Gene2Disease', phenotype))

    def category(entry: EntryTitles) -> str:
        return BIOLINK_CATEGORY.get(entry.prefix, '')

    rows = sorted(
        (direction, category(subject), predicate, subject.mim, category(object_), object_.mim)
        + (subject.label, object_.label)
        for subject, predicate, direction, object_ in relationships
    )
    header = '?SubId ?SubLabel ?SubBiolink ?Pred ?PredDirection ?ObjId ?ObjLabel ?ObjBiolink'
    lines = [header.split()] + [
        [
            _quoted(f'OMIM:{subject_mim}'),
            _quoted(subject_label),
            _quoted(subject_category),
            _quoted(predicate),
            _quoted(direction),
            _quoted(f'OMIM:{object_mim}'),
            _quoted(object_label),
            _quoted(object_category),
        ]
        for (
            direction,
            subject_category,
            predicate,
            subject_mim,
            object_category,
            object_mim,
            subject_label,
            object_label,
        ) in rows
    ]
    return ''.join('\t'.join(line) + '\n' for line in lines)


def _quoted(value: str) -> str:
    """A string as SPARQL's TSV results format writes it; an unbound variable is left empty."""
    if not value:
        return ''
    escaped = value.replace('\\', '\\\\').replace('"', '\\"')
    return f'"{escaped}"'


def _tsv(rows: Iterable[Iterable[str]], lineterminator: str) -> str:
    out = io.StringIO()
    csv.writer(out, delimiter='\t', lineterminator=lineterminator).writerows(rows)
    return out.getvalue()
