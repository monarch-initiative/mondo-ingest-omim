"""omim.owl, as an RDF graph for `robot convert` to turn into functional syntax."""

from collections.abc import Iterable

from rdflib import OWL, RDF, RDFS, SKOS, BNode, Graph, Literal, Namespace, URIRef
from rdflib.term import Node

from omim2obo.associations import Association, CuratorExcluded, CuratorProtected
from omim2obo.entry_prefix import BIOLINK_CATEGORY, EntryPrefix
from omim2obo.parse.omim import (
    MappingKey,
    Mim2GeneEntryType,
    Mim2GeneRow,
    MorbidmapRow,
    PhenotypicSeries,
    RemovedRow,
)
from omim2obo.parse.reference import Mapping, PubmedReference
from omim2obo.predicates import MORBIDMAP_PHENOTYPE_MAPPING_KEY_MEANINGS, Predicate
from omim2obo.titles import EntryTitles

OBO = 'http://purl.obolibrary.org/obo/'
BIOLINK = Namespace('https://w3id.org/biolink/vocab/')
OBOINOWL = Namespace('http://www.geneontology.org/formats/oboInOwl#')
SSSOM = Namespace('https://w3id.org/sssom/')
MONDO_EXCLUSION_REASON = URIRef(OBO + 'MONDO_exclusionReason')
OMIM_INCLUDED = URIRef(OBO + 'mondo#omim_included')
ABBREVIATION = URIRef(OBO + 'OMO_0003000')
TERM_REPLACED_BY = URIRef(OBO + 'IAO_0100001')
MENTIONS = URIRef(OBO + 'IAO_0000142')
GENE = URIRef(OBO + 'SO_0000704')
HUMAN = URIRef(OBO + 'NCBITaxon_9606')
HAS_PHENOTYPE = URIRef(OBO + 'RO_0002200')
IS_SUBSEQUENCE_OF = URIRef(OBO + 'RO_0002525')
CAUSAL_GERMLINE_MUTATION_IN = URIRef(OBO + 'RO_0004013')

IRI_PREFIXES = {
    'OMIM': 'https://omim.org/entry/',
    'OMIMPS': 'https://omim.org/phenotypicSeries/PS',
    'MONDO': OBO + 'MONDO_',
    'RO': OBO + 'RO_',
    'CHR': OBO + 'CHR_',
    'HGNC': 'http://identifiers.org/hgnc/',
    'HGNC.symbol': 'https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/',
    'NCBIGene': 'https://www.ncbi.nlm.nih.gov/gene/',
    'PMID': 'http://www.ncbi.nlm.nih.gov/pubmed/',
    'UMLS': 'http://linkedlifedata.com/resource/umls/id/',
    'Orphanet': 'http://www.orpha.net/ORDO/Orphanet_',
    'skos': str(SKOS),
    'semapv': 'https://w3id.org/semapv/vocab/',
}

type Annotations = Iterable[tuple[URIRef, Node]]


def iri(curie: str) -> URIRef:
    prefix, local_id = curie.split(':', 1)
    return URIRef(IRI_PREFIXES[prefix] + local_id)


def omim(mim: str) -> URIRef:
    return iri(f'OMIM:{mim}')


def render_owl(
    *,
    version: str,
    entries: Iterable[EntryTitles],
    removed: Iterable[RemovedRow],
    mim2gene: Iterable[Mim2GeneRow],
    symbols: dict[str, str],
    hgnc_mappings: Iterable[Mapping],
    series: Iterable[PhenotypicSeries],
    morbidmap: Iterable[MorbidmapRow],
    associations: Iterable[Association],
    references: Iterable[PubmedReference],
    entry_mappings: Iterable[Mapping],
    mondo_mappings: Iterable[Mapping],
) -> Graph:
    """`version` is the release date, as YYYY-MM-DD."""
    graph = Graph()
    _add_ontology(graph, version)
    for entry in entries:
        _add_entry(graph, entry)
    for row in removed:
        _add_removed(graph, row)
    _add_gene_identifiers(graph, mim2gene, symbols, hgnc_mappings)
    for a_series in series:
        _add_series(graph, a_series)
    _add_cyto_locations(graph, morbidmap)
    associations = list(associations)
    # QUIRK: where a pair is protected by more than one row, only the last row's curator is named.
    last_curator = {
        (a.phenotype_mim, a.gene_mim): a.provenance.orcid
        for a in associations
        if isinstance(a.provenance, CuratorProtected)
    }
    for association in associations:
        _add_association(graph, association, last_curator)
    for reference in references:
        graph.add((omim(reference.mim), MENTIONS, iri(f'PMID:{reference.pmid}')))
    for mapping in entry_mappings:
        graph.add((iri(mapping.subject_id), iri(mapping.predicate_id), iri(mapping.object_id)))
    for mapping in mondo_mappings:
        _add_mondo_mapping(graph, mapping)
    _add_hgnc_links(graph)
    return graph


def _add_ontology(graph: Graph, version: str) -> None:
    ontology = URIRef(OBO + 'mondo/omim.owl')
    graph.add((ontology, RDF.type, OWL.Ontology))
    graph.add((ontology, OWL.versionIRI, URIRef(f'{OBO}mondo/releases/{version}/omim.owl')))
    graph.add((ontology, OWL.versionInfo, Literal(version)))
    for annotation_property in (
        OBOINOWL.hasSynonymType,
        OBOINOWL.source,
        OMIM_INCLUDED,
        ABBREVIATION,
        BIOLINK.has_evidence,
        SSSOM.mapping_justification,
        SSSOM.subject_label,
        SSSOM.object_label,
    ):
        graph.add((annotation_property, RDF.type, OWL.AnnotationProperty))
    graph.add((HUMAN, RDF.type, OWL.Class))
    graph.add((HUMAN, RDFS.label, Literal('Homo sapiens')))


def _add_entry(graph: Graph, entry: EntryTitles) -> None:
    uri = omim(entry.mim)
    graph.add((uri, RDF.type, OWL.Class))
    graph.add((uri, RDFS.label, Literal(entry.label)))

    if category := BIOLINK_CATEGORY.get(entry.prefix):
        graph.add((uri, BIOLINK.category, BIOLINK[category]))
    if category == 'Gene':
        graph.add((uri, RDFS.subClassOf, GENE))
        graph.add((uri, MONDO_EXCLUSION_REASON, URIRef(OBO + 'MONDO_nonDisease')))
    elif entry.prefix is EntryPrefix.OTHER:
        graph.add((uri, MONDO_EXCLUSION_REASON, URIRef(OBO + 'MONDO_excludeTrait')))

    deprecated = (OWL.deprecated, Literal(True))
    abbreviation = (OBOINOWL.hasSynonymType, ABBREVIATION)
    exact, related = OBOINOWL.hasExactSynonym, OBOINOWL.hasRelatedSynonym

    for title in (entry.preferred_title, *entry.alternative.titles):
        graph.add((uri, exact, Literal(title)))
    for symbol in (*entry.preferred_symbols, *entry.alternative.symbols):
        _add_annotated(graph, uri, exact, Literal(symbol), [abbreviation])
    for title in entry.alternative.former_titles:
        _add_annotated(graph, uri, related, Literal(title), [deprecated])
    for symbol in entry.alternative.former_symbols:
        _add_annotated(graph, uri, related, Literal(symbol), [deprecated, abbreviation])

    if entry.included.titles or entry.included.symbols:
        comment = "This term has one or more labels that end with ', INCLUDED'."
        graph.add((uri, RDFS.comment, Literal(comment)))
    # Included symbols are abbreviations too, but omim_included is not a synonym property, so they
    # cannot be given a synonym type.
    for name in (*entry.included.titles, *entry.included.symbols):
        graph.add((uri, OMIM_INCLUDED, Literal(name)))
    for name in (*entry.included.former_titles, *entry.included.former_symbols):
        _add_annotated(graph, uri, OMIM_INCLUDED, Literal(name), [deprecated])


def _add_removed(graph: Graph, row: RemovedRow) -> None:
    uri = omim(row.mim)
    graph.add((uri, RDF.type, OWL.Class))
    graph.add((uri, OWL.deprecated, Literal(True)))
    match row.moved_to:
        case []:
            # QUIRK: an entry removed outright is given OMIM's placeholder title as its name.
            graph.add((uri, RDFS.label, Literal('removed from database')))
            graph.add((uri, OBOINOWL.hasExactSynonym, Literal('removed from database')))
        case [replacement]:
            graph.add((uri, TERM_REPLACED_BY, omim(replacement)))
        case replacements:
            for replacement in replacements:
                graph.add((uri, OBOINOWL.consider, omim(replacement)))


def _add_gene_identifiers(
    graph: Graph,
    mim2gene: Iterable[Mim2GeneRow],
    symbols: dict[str, str],
    hgnc_mappings: Iterable[Mapping],
) -> None:
    for row in mim2gene:
        if not row.entrez_id:
            continue
        ncbi_gene = iri(f'NCBIGene:{row.entrez_id}')
        if row.entry_type in (Mim2GeneEntryType.GENE, Mim2GeneEntryType.GENE_PHENOTYPE):
            graph.add((omim(row.mim), SKOS.exactMatch, ncbi_gene))
        elif row.entry_type in (
            Mim2GeneEntryType.PHENOTYPE,
            Mim2GeneEntryType.PREDOMINANTLY_PHENOTYPES,
        ):
            _add_restriction(graph, ncbi_gene, HAS_PHENOTYPE, omim(row.mim))
    for mim, symbol in symbols.items():
        graph.add((omim(mim), SKOS.exactMatch, iri(f'HGNC.symbol:{symbol}')))
    for mapping in hgnc_mappings:
        graph.add((iri(mapping.subject_id), SKOS.exactMatch, iri(mapping.object_id)))


def _add_series(graph: Graph, series: PhenotypicSeries) -> None:
    uri = iri(f'OMIMPS:{series.number.removeprefix("PS")}')
    graph.add((uri, RDF.type, OWL.Class))
    graph.add((uri, RDFS.label, Literal(series.title)))
    graph.add((uri, BIOLINK.category, BIOLINK.Disease))
    for mim in series.members:
        graph.add((omim(mim), RDFS.subClassOf, uri))


def _add_cyto_locations(graph: Graph, morbidmap: Iterable[MorbidmapRow]) -> None:
    location_of: dict[str, str] = {}
    for row in morbidmap:
        location_of.setdefault(row.gene_mim, row.cyto_location)
    for gene_mim, location in location_of.items():
        if location:
            # 9606 is the NCBI Taxonomy ID of Homo sapiens
            _add_restriction(
                graph, omim(gene_mim), IS_SUBSEQUENCE_OF, iri(f'CHR:9606chr{location}')
            )


def _add_association(
    graph: Graph, association: Association, last_curator: dict[tuple[str, str], str | None]
) -> None:
    phenotype, gene = omim(association.phenotype_mim), omim(association.gene_mim)
    provenance = association.provenance
    source, mapping_keys = None, set()
    match provenance:
        case CuratorProtected(omim=omim_says):
            source = last_curator[association.phenotype_mim, association.gene_mim]
            # QUIRK: a protected association always claims OMIM's strongest evidence, whatever
            # OMIM gives it, or whether OMIM lists it at all.
            mapping_keys = {MappingKey.MOLECULAR_BASIS_KNOWN}
            if omim_says:
                mapping_keys.add(omim_says.mapping_key)
            graph.add((phenotype, SKOS.exactMatch, iri(provenance.mondo_id)))
            graph.add((gene, SKOS.exactMatch, iri(provenance.hgnc_id)))
        case CuratorExcluded(orcid=orcid, omim=omim_says):
            source = orcid
            mapping_keys = {omim_says.mapping_key}
        case omim_says:
            mapping_keys = {omim_says.mapping_key}

    for mapping_key in mapping_keys:
        meaning = MORBIDMAP_PHENOTYPE_MAPPING_KEY_MEANINGS[mapping_key]
        evidence = Literal(f'Evidence: ({mapping_key.value}) {meaning}')
        annotations = [(BIOLINK.has_evidence, evidence), (RDFS.comment, evidence)]
        if source:
            annotations.append((OBOINOWL.source, URIRef(source)))
        predicate = iri(association.predicate.value)
        if association.predicate is Predicate.GERMLINE_MUTATION:
            _add_restriction(graph, phenotype, predicate, gene, annotations)
            _add_restriction(graph, gene, CAUSAL_GERMLINE_MUTATION_IN, phenotype, annotations)
        else:
            _add_restriction(graph, gene, predicate, phenotype, annotations)


def _add_mondo_mapping(graph: Graph, mapping: Mapping) -> None:
    """Both ways round: as Mondo publishes it, with its provenance, and from the OMIM entry."""
    subject, predicate, object_ = (
        iri(mapping.subject_id),
        iri(mapping.predicate_id),
        iri(mapping.object_id),
    )
    annotations = [(SSSOM.mapping_justification, iri(mapping.mapping_justification))]
    if mapping.subject_label:
        annotations.append((SSSOM.subject_label, Literal(mapping.subject_label)))
    if mapping.object_label:
        annotations.append((SSSOM.object_label, Literal(mapping.object_label)))
    _add_annotated(graph, subject, predicate, object_, annotations)
    graph.add((object_, predicate, subject))


def _add_hgnc_links(graph: Graph) -> None:
    """Wherever a class has its basis in an OMIM gene, say so of the HGNC gene as well."""
    germline_mutation = iri(Predicate.GERMLINE_MUTATION.value)
    links = set()
    for restriction in graph.subjects(OWL.onProperty, germline_mutation):
        gene = graph.value(restriction, OWL.someValuesFrom)
        for subject in graph.subjects(RDFS.subClassOf, restriction):
            for match in graph.objects(gene, SKOS.exactMatch):
                if match.startswith(IRI_PREFIXES['HGNC']):
                    links.add((subject, match))
    for subject, hgnc_gene in links:
        graph.add((hgnc_gene, RDF.type, OWL.Class))
        _add_restriction(graph, subject, germline_mutation, hgnc_gene)


def _add_restriction(
    graph: Graph, subject: URIRef, predicate: URIRef, filler: URIRef, annotations: Annotations = ()
) -> None:
    """`subject SubClassOf predicate some filler`"""
    restriction = BNode()
    graph.add((restriction, RDF.type, OWL.Restriction))
    graph.add((restriction, OWL.onProperty, predicate))
    graph.add((restriction, OWL.someValuesFrom, filler))
    _add_annotated(graph, subject, RDFS.subClassOf, restriction, annotations)


def _add_annotated(
    graph: Graph, subject: URIRef, predicate: URIRef, object_: Node, annotations: Annotations
) -> None:
    graph.add((subject, predicate, object_))
    annotations = list(annotations)
    if annotations:
        axiom = BNode()
        graph.add((axiom, RDF.type, OWL.Axiom))
        graph.add((axiom, OWL.annotatedSource, subject))
        graph.add((axiom, OWL.annotatedProperty, predicate))
        graph.add((axiom, OWL.annotatedTarget, object_))
        for annotation in annotations:
            graph.add((axiom, *annotation))
