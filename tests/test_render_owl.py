from rdflib import OWL, RDF, RDFS, SKOS, Graph, Literal, URIRef

from omim2obo.associations import (
    Association,
    CuratorExcluded,
    CuratorProtected,
    OmimMorbidmap,
)
from omim2obo.entry_prefix import EntryPrefix
from omim2obo.parse.omim import MappingKey, MorbidmapRow, RemovedRow
from omim2obo.parse.reference import Mapping
from omim2obo.predicates import Predicate
from omim2obo.release import Release
from omim2obo.render.owl import (
    BIOLINK,
    GENE,
    OBOINOWL,
    SSSOM,
    TERM_REPLACED_BY,
    iri,
    omim,
    render_owl,
)
from omim2obo.titles import EntryTitles, TitleGroup

KNEE, WBK1, WBK2 = omim('900001'), omim('900100'), omim('900101')
ORCID = 'https://orcid.org/0000-0000-0000-0001'
OTHER_ORCID = 'https://orcid.org/0000-0000-0000-0002'
KEY_2 = OmimMorbidmap(mapping_key=MappingKey.LINKAGE, phenotype_label='Wobbly knee syndrome')
KEY_3 = OmimMorbidmap(
    mapping_key=MappingKey.MOLECULAR_BASIS_KNOWN, phenotype_label='Wobbly knee syndrome'
)


def render(**inputs) -> Graph:
    return render_owl(Release(**{'version': '2026-01-01', **inputs}))


def restrictions(graph: Graph) -> set[tuple]:
    """(subject, property, filler, evidence, source) for each `SubClassOf … some …` axiom."""
    found = set()
    for subject, restriction in graph.subject_objects(RDFS.subClassOf):
        if (restriction, RDF.type, OWL.Restriction) not in graph:
            continue
        axiom = graph.value(predicate=OWL.annotatedTarget, object=restriction)
        evidence = graph.value(axiom, BIOLINK.has_evidence) if axiom else None
        found.add(
            (
                subject,
                graph.value(restriction, OWL.onProperty),
                graph.value(restriction, OWL.someValuesFrom),
                evidence and str(evidence)[:13],
                graph.value(axiom, OBOINOWL.source) if axiom else None,
            )
        )
    return found


def association(provenance, predicate=Predicate.GERMLINE_MUTATION, gene='900100') -> Association:
    return Association(
        phenotype_mim='900001', gene_mim=gene, predicate=predicate, provenance=provenance
    )


def protected(orcid=ORCID, omim=KEY_3) -> CuratorProtected:
    return CuratorProtected(
        orcid=orcid,
        protection_type='manually reviewed',
        mondo_id='MONDO:9900001',
        hgnc_id='HGNC:99901',
        omim=omim,
    )


def test_version():
    graph = render(version='2026-01-01')
    ontology = URIRef('http://purl.obolibrary.org/obo/mondo/omim.owl')
    assert graph.value(ontology, OWL.versionIRI) == URIRef(
        'http://purl.obolibrary.org/obo/mondo/releases/2026-01-01/omim.owl'
    )


def test_gene_is_labelled_by_symbol_and_phenotype_by_title():
    none = TitleGroup(titles=(), symbols=(), former_titles=(), former_symbols=())
    graph = render(
        entries=[
            EntryTitles(
                mim='900100',
                prefix=EntryPrefix.GENE,
                preferred_title='wobbly knee gene 1',
                preferred_symbols=('WBK1',),
                alternative=none,
                included=none,
            ),
            EntryTitles(
                mim='900001',
                prefix=EntryPrefix.PHENOTYPE,
                preferred_title='wobbly knee',
                preferred_symbols=('WKS',),
                alternative=none,
                included=none,
            ),
        ]
    )
    assert graph.value(WBK1, RDFS.label) == Literal('WBK1')
    assert (WBK1, RDFS.subClassOf, GENE) in graph
    assert graph.value(KNEE, RDFS.label) == Literal('wobbly knee')
    assert graph.value(KNEE, BIOLINK.category) == BIOLINK.Disease


def test_removed_entry_points_to_its_replacements():
    graph = render(
        removed=[
            RemovedRow(mim='900001', moved_to=('900002',)),
            RemovedRow(mim='900003', moved_to=('900004', '900005')),
        ]
    )
    assert graph.value(KNEE, OWL.deprecated) == Literal(True)
    assert graph.value(KNEE, TERM_REPLACED_BY) == omim('900002')
    assert set(graph.objects(omim('900003'), OBOINOWL.consider)) == {
        omim('900004'),
        omim('900005'),
    }
    assert graph.value(KNEE, RDFS.label) is None


def test_QUIRK_entry_removed_without_replacement_is_labelled():
    graph = render(removed=[RemovedRow(mim='900001', moved_to=())])
    assert graph.value(KNEE, RDFS.label) == Literal('removed from database')


def test_germline_mutation_is_stated_from_both_ends():
    assert restrictions(render(associations=[association(KEY_3)])) == {
        (KNEE, iri('RO:0004003'), WBK1, 'Evidence: (3)', None),
        (WBK1, iri('RO:0004013'), KNEE, 'Evidence: (3)', None),
    }


def test_other_predicates_are_stated_from_the_gene():
    graph = render(associations=[association(KEY_2, Predicate.CAUSES)])
    assert restrictions(graph) == {(WBK1, iri('RO:0003303'), KNEE, 'Evidence: (2)', None)}


def test_exclusion_names_its_curator():
    excluded = CuratorExcluded(orcid=ORCID, reason='somatic', omim=KEY_3)
    graph = render(associations=[association(excluded, Predicate.CAUSES_OR_CONTRIBUTES)])
    assert restrictions(graph) == {(WBK1, iri('RO:0003302'), KNEE, 'Evidence: (3)', URIRef(ORCID))}


def test_protection_maps_its_phenotype_and_gene():
    graph = render(associations=[association(protected())])
    assert (KNEE, SKOS.exactMatch, iri('MONDO:9900001')) in graph
    assert (WBK1, SKOS.exactMatch, iri('HGNC:99901')) in graph


def test_QUIRK_protection_claims_key_3_evidence_whatever_omim_says():
    for omim_says in (KEY_2, None):
        graph = render(associations=[association(protected(omim=omim_says))])
        assert (KNEE, iri('RO:0004003'), WBK1, 'Evidence: (3)', URIRef(ORCID)) in restrictions(
            graph
        )


def test_QUIRK_only_the_last_curator_to_protect_a_pair_is_named():
    graph = render(
        associations=[association(protected(ORCID)), association(protected(OTHER_ORCID))]
    )
    assert {source for *_, source in restrictions(graph)} == {URIRef(OTHER_ORCID), None}


def test_basis_in_an_omim_gene_is_also_stated_of_its_hgnc_gene():
    graph = render(
        associations=[association(KEY_3)],
        hgnc_mappings=[
            Mapping(
                subject_id='OMIM:900100',
                predicate_id='skos:exactMatch',
                object_id='HGNC:99901',
                mapping_justification='semapv:UnspecifiedMatching',
            )
        ],
    )
    assert (KNEE, iri('RO:0004003'), iri('HGNC:99901'), None, None) in restrictions(graph)


def test_mondo_mapping_is_stated_both_ways_round():
    mapping = Mapping(
        subject_id='MONDO:9900001',
        predicate_id='skos:exactMatch',
        object_id='OMIM:900001',
        mapping_justification='semapv:ManualMappingCuration',
        subject_label='wks',
    )
    graph = render(mondo_mappings=[mapping])
    mondo = iri('MONDO:9900001')
    assert (mondo, SKOS.exactMatch, KNEE) in graph
    assert (KNEE, SKOS.exactMatch, mondo) in graph
    axiom = graph.value(predicate=OWL.annotatedSource, object=mondo)
    assert graph.value(axiom, SSSOM.subject_label) == Literal('wks')
    assert graph.value(axiom, SSSOM.object_label) is None


def test_gene_location_comes_from_its_first_morbidmap_row():
    def row(location):
        return MorbidmapRow(
            phenotype_label='Wobbly knee',
            phenotype_mim='900001',
            mapping_key=MappingKey.LINKAGE,
            gene_symbols=(),
            gene_mim='900100',
            cyto_location=location,
        )

    graph = render(morbidmap=[row('1p36'), row('2q11')])
    assert {filler for _, _, filler, _, _ in restrictions(graph)} == {iri('CHR:9606chr1p36')}
