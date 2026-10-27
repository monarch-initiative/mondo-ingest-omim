"""Tests for sparql/mondo-omim-genes.sparql: projecting OMIM gene relations onto Mondo classes

Of interest is the pairing in data/protected-disease-gene.tsv. One OMIM phenotype entry can cover
several Mondo classes, so joining every mapped Mondo class against every gene of the entry gives each
of them all of the genes.
"""
from typing import Dict, List, Tuple

from rdflib import BNode, Graph, Namespace, RDF, RDFS, OWL, SKOS

from omim2obo.config import ROOT_DIR
from omim2obo.namespaces import HGNC, MONDO, MONDONS, OMIM, RO


QUERY = (ROOT_DIR / 'sparql' / 'mondo-omim-genes.sparql').read_text()
BERNARD_SOULIER_GENES = {'606672': '4439', '138720': '4440', '173515': '4444'}  # gene MIM: HGNC id


def omim_entry_graph(
    entry: str,
    genes: Dict[str, str],
    exact_mondo: List[str],
    protected: List[Tuple[str, str]],
) -> Graph:
    """Build the part of omim.owl that the query reads.

    :param genes: {gene MIM: HGNC id} of the entry's disease-gene associations.
    :param exact_mondo: Mondo ids mapped to the entry by mondo_exactmatch_omim.sssom.tsv.
    :param protected: (Mondo id, gene MIM) pairs as omim2obo.main asserts them for protected rows.
    """
    graph = Graph()
    graph.add((OMIM[entry], RDF.type, OWL.Class))
    for gene_mim, hgnc_id in genes.items():
        graph.add((OMIM[gene_mim], RDF.type, OWL.Class))
        graph.add((OMIM[gene_mim], SKOS.exactMatch, HGNC[hgnc_id]))
        restriction = BNode()
        graph.add((restriction, RDF.type, OWL.Restriction))
        graph.add((restriction, OWL.onProperty, RO['0004003']))
        graph.add((restriction, OWL.someValuesFrom, OMIM[gene_mim]))
        graph.add((OMIM[entry], RDFS.subClassOf, restriction))
    for mondo_id in exact_mondo:
        graph.add((OMIM[entry], SKOS.exactMatch, MONDO[mondo_id]))
    for mondo_id, gene_mim in protected:
        graph.add((OMIM[entry], SKOS.exactMatch, MONDO[mondo_id]))
        graph.add((MONDO[mondo_id], MONDONS.protected_gene, OMIM[gene_mim]))
    return graph


def mondo_gene_pairs(graph: Graph) -> List[Tuple[str, str]]:
    """Run the query, returning (Mondo id, HGNC id) pairs."""
    return sorted(
        (str(row.mondo_id).split('_')[-1], str(row.hgnc_id).rsplit('/', 1)[-1])
        for row in graph.query(QUERY))


def test_included_terms_keep_only_their_own_gene():
    """OMIM:231200 (Bernard-Soulier syndrome) covers 3 'included terms', which Mondo splits into 3
    gene-specific subtypes. Each subtype gets its own gene, and the entry's own Mondo class gets no
    gene association at all."""
    graph = omim_entry_graph(
        entry='231200', genes=BERNARD_SOULIER_GENES, exact_mondo=['0009276'],
        protected=[('1060238', '606672'), ('1060239', '138720'), ('1060237', '173515')])

    assert mondo_gene_pairs(graph) == [
        ('1060237', '4444'), ('1060238', '4439'), ('1060239', '4440')]


def test_unprotected_entry_is_untouched():
    """An entry with no protected rows keeps whatever the join produces, multiple genes included."""
    graph = omim_entry_graph(
        entry='608133', genes={'179605': '9942', '600000': '10254'},
        exact_mondo=['0011974'], protected=[])

    assert mondo_gene_pairs(graph) == [('0011974', '10254'), ('0011974', '9942')]


def test_one_gene_may_be_curated_onto_several_classes():
    """OMIM:217090, where curators scoped the single gene to the entry's class and 2 others."""
    graph = omim_entry_graph(
        entry='217090', genes={'173350': '9071'}, exact_mondo=['0009009'],
        protected=[('0009009', '173350'), ('0100537', '173350'), ('0100538', '173350')])

    assert mondo_gene_pairs(graph) == [
        ('0009009', '9071'), ('0100537', '9071'), ('0100538', '9071')]


def test_curated_class_need_not_be_an_exact_match():
    """A curated Mondo class that is not an exact match of the entry still receives the gene, and the
    exact match does not."""
    graph = omim_entry_graph(
        entry='619805', genes={'600000': '14438'}, exact_mondo=['0030721'],
        protected=[('0030515', '600000')])

    assert mondo_gene_pairs(graph) == [('0030515', '14438')]


def test_unprotected_gene_does_not_leak_onto_a_curated_class():
    """Curating one gene of an entry must neither restrict that entry's other genes nor hand them to
    the curated class."""
    graph = omim_entry_graph(
        entry='608133', genes={'179605': '9942', '600000': '10254'},
        exact_mondo=['0011974'], protected=[('1060145', '179605')])

    assert mondo_gene_pairs(graph) == [('0011974', '10254'), ('1060145', '9942')]
