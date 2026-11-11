"""Everything one release is rendered from."""

from dataclasses import dataclass

from omim2obo.associations import Association
from omim2obo.parse.omim import Mim2GeneRow, MorbidmapRow, PhenotypicSeries, RemovedRow
from omim2obo.parse.reference import Mapping, PubmedReference
from omim2obo.titles import EntryTitles


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
