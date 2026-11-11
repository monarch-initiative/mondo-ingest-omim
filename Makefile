export PYTHONPATH := src

# A recipe that fails must not leave a half-written target behind
.DELETE_ON_ERROR:

TSVS = mondo-omim-genes.robot.tsv disease-gene-relationships-qc.tsv \
	mondo-omim-susceptibility-subset.robot.tsv review.tsv

# What a build reads. Each can be pointed elsewhere, e.g. at a frozen snapshot:
#   make omim.owl OMIM_DATA=omim-data VERSION=2026-01-01
# OMIM_DATA: directory of OMIM's licensed .txt downloads. They are never committed
OMIM_DATA ?= data
MONDO_MAPPINGS ?= mondo_exactmatch_omim.sssom.tsv mondo_exactmatch_omimps.sssom.tsv
# VERSION: the release date, YYYY-MM-DD. Today's if not given
VERSION ?=

OMIM_NAMES = mimTitles.txt morbidmap.txt mim2gene.txt genemap2.txt phenotypicSeries.txt
OMIM_FILES = $(addprefix $(OMIM_DATA)/,$(OMIM_NAMES))
HGNC = data/hgnc/hgnc_complete_set.txt
CURATED = data/protected-disease-gene.tsv data/exclusions-disease-gene.tsv data/known_capitalizations.tsv
API_CACHE = data/mappings.tsv data/pubmed-refs.tsv

# MAIN COMMANDS / GOALS ------------------------------------------------------------------------------------------------
# all: Fetch the latest inputs, then build every release artefact from them
.PHONY: all
all:
	$(MAKE) fetch
	$(MAKE) omim.owl omim.sssom.tsv

# build: One run writes omim.nt and the TSVs together (a grouped target, so GNU make 4.3 or later).
# Rebuilt when any input changes
omim.nt $(TSVS) &: $(OMIM_FILES) $(MONDO_MAPPINGS) $(HGNC) $(CURATED) $(API_CACHE)
	python3 -m omim2obo build \
		--mim-titles $(OMIM_DATA)/mimTitles.txt \
		--morbidmap $(OMIM_DATA)/morbidmap.txt \
		--mim2gene $(OMIM_DATA)/mim2gene.txt \
		--genemap2 $(OMIM_DATA)/genemap2.txt \
		--phenotypic-series $(OMIM_DATA)/phenotypicSeries.txt \
		$(addprefix --mondo-mappings ,$(MONDO_MAPPINGS)) \
		$(if $(VERSION),--version $(VERSION))

# Create OWL artefact: the graph in functional syntax
omim.owl: omim.nt
	robot convert -i $< -f ofn -o $@

# Create mapping artefact
omim.sssom.tsv: omim.json
	sssom parse omim.json -I obographs-json -m data/metadata.sssom.yml -o omim.sssom.tsv

# Only a step on the way to omim.sssom.tsv; make removes it afterwards
.INTERMEDIATE: omim.json
omim.json: omim.owl
	robot convert -i $< -o omim.json

# FETCHING -------------------------------------------------------------------------------------------------------------
# make cannot tell when a file on someone else's server has changed, so everything that comes over
# HTTP has a target that always runs. Some of what they write is committed (HGNC, the API cache), so
# expect a diff afterwards.
# fetch: All of them
.PHONY: fetch
fetch:
	$(MAKE) fetch-omim fetch-mondo fetch-hgnc
	$(MAKE) fetch-cache

# fetch-omim: OMIM's five downloadable files, into data/. Needs the download key in .env
.PHONY: fetch-omim
fetch-omim:
	python3 -m omim2obo fetch files

# fetch-cache: PubMed references and UMLS / Orphanet mappings for entries OMIM has changed since the
# last run, into data/ (committed). Needs the API key in .env, and OMIM's files
.PHONY: fetch-cache
fetch-cache: | data/mimTitles.txt data/morbidmap.txt
	python3 -m omim2obo fetch cache

# fetch-mondo: Mondo's mappings to OMIM entries and to phenotypic series
MONDO_DOWNLOADS = mondo_exactmatch_omim.sssom.tsv mondo_exactmatch_omimps.sssom.tsv
.PHONY: fetch-mondo
fetch-mondo:
	for f in $(MONDO_DOWNLOADS); do \
		wget "http://purl.obolibrary.org/obo/mondo/mappings/$$f" -O $$f.tmp && mv $$f.tmp $$f || exit 1; \
	done

# fetch-hgnc: HGNC's complete gene set (committed)
.PHONY: fetch-hgnc
fetch-hgnc:
	@echo "Downloading HGNC complete set..."
	@mkdir -p data/hgnc
	wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 --tries=3 \
		--continue --progress=bar:force \
		"https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt" -O $(HGNC).tmp
	@echo "Verifying download..."
	@if [ ! -s $(HGNC).tmp ]; then \
		echo "ERROR: Downloaded file is empty!"; \
		rm -f $(HGNC).tmp; \
		exit 1; \
	fi
	@line_count=$$(wc -l < $(HGNC).tmp | tr -d ' '); \
	if [ $$line_count -lt 40000 ]; then \
		echo "ERROR: Downloaded file has only $$line_count lines (expected >40000)"; \
		rm -f $(HGNC).tmp; \
		exit 1; \
	fi
	@if ! head -1 $(HGNC).tmp | grep -q "hgnc_id"; then \
		echo "ERROR: Downloaded file doesn't have expected header"; \
		rm -f $(HGNC).tmp; \
		exit 1; \
	fi
	@mv $(HGNC).tmp $(HGNC)
	@echo "✓ Successfully downloaded and verified HGNC file ($$line_count lines)"

# A download that is missing altogether is fetched when something needs it
$(addprefix data/,$(OMIM_NAMES)) &:
	$(MAKE) fetch-omim

$(MONDO_DOWNLOADS) &:
	$(MAKE) fetch-mondo

$(HGNC):
	$(MAKE) fetch-hgnc

# clean: Remove everything the build writes. Downloaded inputs are left alone
.PHONY: clean
clean:
	rm -f omim.nt omim.owl omim.json omim.sssom.tsv $(TSVS)

# SETUP / INSTALLATION -------------------------------------------------------------------------------------------------
# The ODK container already has everything the build imports. This is for working outside it.
.PHONY: install
install:
	uv sync

# QA / TESTING ---------------------------------------------------------------------------------------------------------
.PHONY: test
test:
	 python3 -m pytest tests

.PHONY: lint
lint:
	ruff format --check
	ruff check

# Capture golden outputs from the legacy pipeline as it was before the rewrite. Needs Docker and OMIM's downloads
# GOLDEN_OMIM_DATA: the frozen copy of OMIM's downloads to build it from
GOLDEN_OMIM_DATA ?= omim-data

.PHONY: golden
golden:
	OMIM_DATA=$(GOLDEN_OMIM_DATA) golden/capture.sh

# HELP -----------------------------------------------------------------------------------------------------------------
.PHONY: help
help:
	@echo "----------------------------------------"
	@echo "	Command reference: OMIM"
	@echo "----------------------------------------"
	@echo "all"
	@echo "Creates all release artefacts.\n"
	@echo "fetch"
	@echo "Downloads every input again and refreshes the API cache.\n"
	@echo "fetch-omim, fetch-mondo, fetch-hgnc, fetch-cache"
	@echo "One source at a time.\n"
	@echo "omim.owl"
	@echo "Builds the main release artefact; redone when an input changes.\n"
	@echo "omim.sssom.tsv"
	@echo "Creates an SSSOM TSV of OMIM terms.\n"
	@echo "install"
	@echo "Installs Python requirements with uv, for working outside the ODK container.\n"
	@echo "clean"
	@echo "Removes the built artefacts.\n"
	@echo "test"
	@echo "Runs unit tests.\n"
	@echo "lint"
	@echo "Checks formatting and lint rules.\n"
