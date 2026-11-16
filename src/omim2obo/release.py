"""Everything one release is rendered from."""

from dataclasses import dataclass
from pathlib import Path

from omim2obo.associations import Association, classify
from omim2obo.genes import approved_symbols, hgnc_mappings
from omim2obo.parse.curation import parse_exclusions, parse_protected
from omim2obo.parse.omim import (
    Mim2GeneRow,
    MorbidmapRow,
    PhenotypicSeries,
    RemovedRow,
    TitleRow,
    parse_genemap2,
    parse_mim2gene,
    parse_mim_titles,
    parse_morbidmap,
    parse_phenotypic_series,
)
from omim2obo.parse.reference import (
    Mapping,
    PubmedReference,
    parse_entry_mappings,
    parse_hgnc,
    parse_known_capitalizations,
    parse_pubmed_refs,
    parse_sssom,
)
from omim2obo.titles import EntryTitles, entry_titles


@dataclass(frozen=True, slots=True, kw_only=True)
class Release:
    version: str
    """The release date, as YYYY-MM-DD."""

    entries: tuple[EntryTitles, ...] = ()
    removed: tuple[RemovedRow, ...] = ()
    mim2gene: tuple[Mim2GeneRow, ...] = ()
    symbols: tuple[tuple[str, str], ...] = ()
    """(gene MIM, HGNC symbol)"""

    hgnc_mappings: tuple[Mapping, ...] = ()
    series: tuple[PhenotypicSeries, ...] = ()
    morbidmap: tuple[MorbidmapRow, ...] = ()
    associations: tuple[Association, ...] = ()
    references: tuple[PubmedReference, ...] = ()
    entry_mappings: tuple[Mapping, ...] = ()
    """OMIM's own mappings of its entries, to UMLS and Orphanet."""

    mondo_mappings: tuple[Mapping, ...] = ()


@dataclass(frozen=True, slots=True, kw_only=True)
class Inputs:
    """Every file a build reads."""

    mim_titles: Path
    morbidmap: Path
    mim2gene: Path
    genemap2: Path
    phenotypic_series: Path
    protected: Path
    exclusions: Path
    known_capitalizations: Path
    hgnc: Path
    pubmed_refs: Path
    entry_mappings: Path
    mondo_mappings: tuple[Path, ...]


def read_release(version: str, inputs: Inputs) -> Release:
    def lines(path: Path) -> list[str]:
        return path.read_text().splitlines()

    titles = parse_mim_titles(lines(inputs.mim_titles))
    capitalizations = parse_known_capitalizations(lines(inputs.known_capitalizations))
    mim2gene = parse_mim2gene(lines(inputs.mim2gene))
    symbols = approved_symbols(mim2gene, parse_genemap2(lines(inputs.genemap2)))
    hgnc = hgnc_mappings(symbols, parse_hgnc(lines(inputs.hgnc)))
    morbidmap = parse_morbidmap(lines(inputs.morbidmap))
    associations = classify(
        morbidmap,
        parse_protected(lines(inputs.protected)),
        parse_exclusions(lines(inputs.exclusions)),
    )
    mondo = []
    for path in inputs.mondo_mappings:
        with path.open() as file:
            mondo += parse_sssom(file)

    return Release(
        version=version,
        entries=tuple(
            entry_titles(row, capitalizations) for row in titles if isinstance(row, TitleRow)
        ),
        removed=tuple(row for row in titles if isinstance(row, RemovedRow)),
        mim2gene=tuple(mim2gene),
        symbols=tuple(symbols.items()),
        hgnc_mappings=tuple(hgnc),
        series=tuple(parse_phenotypic_series(lines(inputs.phenotypic_series))),
        morbidmap=tuple(morbidmap),
        associations=tuple(associations),
        references=tuple(parse_pubmed_refs(lines(inputs.pubmed_refs))),
        entry_mappings=tuple(parse_entry_mappings(lines(inputs.entry_mappings))),
        mondo_mappings=tuple(mondo),
    )
