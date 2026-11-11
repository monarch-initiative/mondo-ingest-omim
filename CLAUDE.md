# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

**IMPORTANT: Always read README.md first for detailed context, curator instructions, and the most up-to-date project information before making any changes or answering questions about this codebase.**

## Overview

This repository transforms OMIM (Online Mendelian Inheritance in Man) data into semantic web formats for integration into the Mondo Disease Ontology. It downloads OMIM data files, parses them, creates disease-gene associations with appropriate relationship ontology predicates, and generates RDF/OWL ontology files.

## Common Commands

The release build runs in the ODK container: prefix `make` with `sh run.sh`. Development outside it
uses uv.

- `make all` - Fetch the latest inputs, then build every release artefact
- `make omim.owl` - Build: `python -m omim2obo build` writes `omim.nt` and the TSVs, `robot convert`
  turns `omim.nt` into `omim.owl`. Rebuilt when an input changes; missing downloads are
  fetched. `python -m omim2obo build --help` lists every input file
- `make fetch` - Fetch everything again: `fetch-omim` (keys in `.env`), `fetch-mondo`, `fetch-hgnc`,
  `fetch-cache` (the PubMed/UMLS/Orphanet cache in `data/`)
- `make clean` - Remove what the build wrote
- `uv run pytest` (or `make test` in the container) - Tests, all on synthetic fixtures
- `uv run ruff format` and `uv run ruff check` - Formatter and linter; both are checked in CI
- `make golden` - Rebuild the pre-rewrite pipeline's outputs into `golden/out/` for comparison

## Architecture

The code is in `src/omim2obo`. It replaced a pipeline that built one large RDF graph and queried
it back out with SPARQL; "legacy" in comments means that pipeline (last at commit `098e956`).

Records in, records out, pure functions between them; OWL is one rendering at the end.

1. **Parse** (`parse/`): each input file becomes a list of frozen dataclasses. `parse/omim.py` for
   OMIM's downloads, `parse/curation.py` for the curated TSVs, `parse/reference.py` for mappings,
   PubMed references and HGNC.
2. **Classify** (`associations.py`): morbidmap rows plus the protected and exclusion files become
   `Association` records, each with a predicate and a typed provenance (`OmimMorbidmap`,
   `CuratorProtected`, `CuratorExcluded`). The mapping-key rules are in `predicates.py`:
   - Key 2 → RO:0003303, key 4 → RO:0003304, when it is the phenotype's only row
   - Key 3, only row, definitive label → RO:0004003 (rendered with its inverse RO:0004013)
   - Otherwise, or if the phenotype is excluded → RO:0003302
   - A protected pair is always RO:0004003
3. **Project** (`project.py`): associations and Mondo's mappings become `MondoGeneRow`s.
4. **Render** (`render/owl.py`, `render/tsv.py`): `Release` (`release.py`) gathers everything a build
   is rendered from; `render_owl` turns it into an rdflib graph, the TSV renderers into text.
   `review.py` lists associations for curators to check.

Also: `titles.py` cleans OMIM's ALL CAPS titles, `genes.py` resolves gene entries to HGNC,
`prefixes.py` is the one CURIE prefix map, `fetch.py` / `omim_client.py` / `api_entry.py` /
`config.py` talk to OMIM.

Code marked `QUIRK` reproduces legacy behaviour that looks wrong but is needed for output parity; do
not fix one in passing. Files carried over from it keep their `git blame` history
(`.git-blame-ignore-revs` lists the reformatting commits).

### Curator Override System

Three key TSV files allow manual curation interventions:
- `data/exclusions-disease-gene.tsv`: Exclude specific disease-gene associations
- `data/protected-disease-gene.tsv`: Protect associations from automated removal
- `data/known_capitalizations.tsv`: Override automatic capitalization

## OMIM Entry Types (`entry_prefix.py`)
- `*` (GENE) and `+` (GENE_AND_PHENOTYPE): gene entries, subclasses of SO:0000704, biolink:Gene
- `#` (PHENOTYPE) and `%` (PHENOTYPE_OR_LOCUS): biolink:Disease
- no prefix (OTHER): excluded with MONDO:excludeTrait
- `^`: removed or moved; a deprecated class pointing to its replacements

## Configuration Files
- `.env` (from `.env.example`): OMIM download and API keys. Never read or print it.
- `pyproject.toml`, `uv.lock`: dependencies, pytest and ruff settings. Dependency floors match what
  the ODK container ships, because the release build uses its system Python.
- `Makefile`: build targets. GNU make 4.3+. Every target lists its real prerequisites: do not reach for
  `-B`. `.PHONY:` goes directly above each phony target
- The command line is `click` (not argparse; typer is not in the ODK container)
- `data/metadata.sssom.yml`: metadata for `omim.sssom.tsv`, still built by `robot` + `sssom parse`

## Notes
- OMIM's downloads are licensed and gitignored; never commit them. Tests use synthetic fixtures (MIM
  numbers 9xxxxx).
- Ontology IRI `http://purl.obolibrary.org/obo/mondo/omim.owl`; version IRI
  `http://purl.obolibrary.org/obo/mondo/releases/{YYYY-MM-DD}/omim.owl`, from `build --version`.
- The graph is written as N-Triples so that no prefixes reach `omim.owl`, keeping its text
  comparable between releases.
