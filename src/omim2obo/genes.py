"""Which HGNC gene an OMIM gene entry is."""

from collections.abc import Iterable

from omim2obo.parse.omim import Genemap2Row, Mim2GeneRow
from omim2obo.parse.reference import HgncGene, Mapping


def approved_symbols(
    mim2gene: Iterable[Mim2GeneRow], genemap2: Iterable[Genemap2Row]
) -> dict[str, str]:
    """Each entry's HGNC symbol: from mim2gene, else genemap2. Where the two disagree, neither."""
    symbols = {row.mim: row.hgnc_symbol for row in mim2gene if row.hgnc_symbol}
    from_genemap2 = {row.mim: row.approved_symbol for row in genemap2 if row.approved_symbol}
    for mim, symbol in from_genemap2.items():
        if symbols.setdefault(mim, symbol) != symbol:
            del symbols[mim]
    return symbols


def hgnc_mappings(symbols: dict[str, str], hgnc: Iterable[HgncGene]) -> list[Mapping]:
    """OMIM names a gene by symbol only; HGNC's own table turns that into an identifier."""
    id_of = {gene.symbol: gene.hgnc_id for gene in hgnc}
    return [
        Mapping(
            subject_id=f'OMIM:{mim}',
            predicate_id='skos:exactMatch',
            object_id=id_of[symbol],
            mapping_justification='semapv:UnspecifiedMatching',
        )
        for mim, symbol in symbols.items()
        if symbol in id_of
    ]
