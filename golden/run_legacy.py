"""Run the legacy pipeline offline, against inputs that are already on disk.

Executed by capture.sh from the root of an exported legacy tree. Two things stand between
`python -m omim2obo --use-cache` and a frozen-input build, and both are handled here rather than by
editing the legacy code:

- `--use-cache` expects the `.tsv` twins of OMIM's `.txt` files, which are only written on download.
- The PubMed/UMLS/Orphanet cache refresh calls the OMIM API unconditionally. It is replaced with a
  no-op, so the build uses `data/mappings.tsv` and `data/pubmed-refs.tsv` exactly as committed.
"""

import omim2obo.parsers.omim_txt_parser as omim_txt_parser
from omim2obo.main import omim2obo

for name in ('mimTitles', 'mim2gene', 'genemap2', 'morbidmap', 'phenotypicSeries'):
    omim_txt_parser.convert_txt_to_tsv(name)

omim_txt_parser.update_cache__pubmed_refs_and_mappings = lambda *args, **kwargs: None

omim2obo(use_cache=True)
