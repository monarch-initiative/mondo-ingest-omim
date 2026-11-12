from omim2obo.associations import Association, CuratorProtected, OmimMorbidmap
from omim2obo.entry_prefix import EntryPrefix
from omim2obo.parse.omim import MappingKey, MorbidmapRow, RemovedRow
from omim2obo.predicates import Predicate
from omim2obo.render.tsv import render_review
from omim2obo.review import ReviewCase, review_cases
from omim2obo.titles import EntryTitles, TitleGroup

KEY_3 = MappingKey.MOLECULAR_BASIS_KNOWN
NONE = TitleGroup(titles=(), symbols=(), former_titles=(), former_symbols=())


def row(label, phenotype='900001', gene='900100', key=KEY_3) -> MorbidmapRow:
    return MorbidmapRow(
        phenotype_label=label,
        phenotype_mim=phenotype,
        mapping_key=key,
        gene_symbols=(),
        gene_mim=gene,
        cyto_location='1p36',
    )


def germline(row_: MorbidmapRow, provenance=None) -> Association:
    provenance = provenance or OmimMorbidmap(
        mapping_key=row_.mapping_key, phenotype_label=row_.phenotype_label
    )
    return Association(
        phenotype_mim=row_.phenotype_mim,
        gene_mim=row_.gene_mim,
        predicate=Predicate.GERMLINE_MUTATION,
        provenance=provenance,
    )


def review(*rows, entries=(), removed=(), associations=None) -> list[ReviewCase]:
    if associations is None:
        associations = [germline(r) for r in rows if r.phenotype_mim]
    return review_cases(associations, rows, entries, removed)


def test_unremarkable_association_is_not_reviewed():
    assert review(row('Wobbly knee syndrome')) == []


def test_digenic_and_somatic_labels():
    info = '(Phenotype: 900001 Wobbly knee, digenic, somatic, ), (Map key: 3), (Gene: 900100)'
    assert review(row('Wobbly knee, digenic, somatic')) == [
        ReviewCase(code=1, value=info),
        ReviewCase(code=3, value=info),
    ]


def test_self_referential_cases_are_numbered():
    """Modelled on OMIM:182280: the phenotype is itself listed as the gene of a phenotype that
    has no MIM number."""
    cases = review(
        row('Wobbly knee', '900001', '900100'),
        row('Wobbly-knee', None, '900001', MappingKey.LINKAGE),
        row('Creaky elbow', '900002', '900101'),
        row('Creaky-elbow', None, '900002', MappingKey.LINKAGE),
    )
    assert [case.value for case in cases] == [
        '1: (Phenotype: 900001 Wobbly knee, ), (Map key: 3), (Gene: 900100)',
        '1: (Phenotype: Wobbly-knee), (Map key: 2), (Gene: 900001)',
        '2: (Phenotype: 900002 Creaky elbow, ), (Map key: 3), (Gene: 900101)',
        '2: (Phenotype: Creaky-elbow), (Map key: 2), (Gene: 900002)',
    ]


def test_phenotype_of_the_wrong_entry_type():
    gene = EntryTitles(
        mim='900001',
        prefix=EntryPrefix.GENE,
        preferred_title='wobbly knee gene',
        preferred_symbols=(),
        alternative=NONE,
        included=NONE,
    )
    [case] = review(row('Wobbly knee'), entries=[gene])
    assert (case.code, case.value.split(',')[0]) == (4, '(Phenotype MIM type GENE)')
    [case] = review(row('Wobbly knee'), removed=[RemovedRow(mim='900001', moved_to=())])
    assert (case.code, case.value.split(',')[0]) == (5, '(Phenotype MIM type OBSOLETE)')


def test_curated_association_is_not_reviewed():
    somatic = row('Wobbly knee, somatic')
    protected = CuratorProtected(
        orcid=None,
        protection_type='manually reviewed',
        mondo_id='MONDO:9900001',
        hgnc_id='HGNC:99901',
        omim=None,
    )
    assert review(somatic, associations=[germline(somatic, protected)]) == []


def test_render():
    assert render_review([ReviewCase(code=3, value='a "quoted" value')]) == (
        'classCode\tclassShortName\tvalue\n3\tD2G: somatic\t"a ""quoted"" value"\n'
    )
