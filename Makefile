export PYTHONPATH := src

.PHONY: all fetch help install test lint cleanup golden

TSVS = mondo-omim-genes.robot.tsv disease-gene-relationships-qc.tsv \
	mondo-omim-susceptibility-subset.robot.tsv review.tsv

# MAIN COMMANDS / GOALS ------------------------------------------------------------------------------------------------
# all: Fetch the latest inputs, then build every release artefact from them
all: fetch
	$(MAKE) -B omim.owl omim.sssom.tsv

# fetch: Download OMIM's files (needs the keys in .env), Mondo's mappings and HGNC, and refresh the
# cache of PubMed references and UMLS / Orphanet mappings in data/
fetch:
	$(MAKE) -B mondo_exactmatch_omim.sssom.tsv mondo_exactmatch_omimps.sssom.tsv data/hgnc/hgnc_complete_set.txt
	python3 -m omim2obo fetch

# build: From files already on disk. One run writes omim.nt and the TSVs together (a grouped target,
# so GNU make 4.3 or later)
# BUILD_ARGS: e.g. "--omim omim-data --version 2026-01-01"; see python3 -m omim2obo build --help
omim.nt $(TSVS) &:
	python3 -m omim2obo build $(BUILD_ARGS)

data/hgnc/hgnc_complete_set.txt:
	@echo "Downloading HGNC complete set..."
	@mkdir -p data/hgnc
	wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 --tries=3 \
		--continue --progress=bar:force \
		"https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt" -O $@.tmp
	@echo "Verifying download..."
	@if [ ! -s $@.tmp ]; then \
		echo "ERROR: Downloaded file is empty!"; \
		rm -f $@.tmp; \
		exit 1; \
	fi
	@line_count=$$(wc -l < $@.tmp | tr -d ' '); \
	if [ $$line_count -lt 40000 ]; then \
		echo "ERROR: Downloaded file has only $$line_count lines (expected >40000)"; \
		rm -f $@.tmp; \
		exit 1; \
	fi
	@if ! head -1 $@.tmp | grep -q "hgnc_id"; then \
		echo "ERROR: Downloaded file doesn't have expected header"; \
		rm -f $@.tmp; \
		exit 1; \
	fi
	@mv $@.tmp $@
	@echo "✓ Successfully downloaded and verified HGNC file ($$line_count lines)"

omim.sssom.tsv: omim.json
	sssom parse omim.json -I obographs-json -m data/metadata.sssom.yml -o omim.sssom.tsv
	make cleanup

mondo_exactmatch_omim.sssom.tsv:
	wget "http://purl.obolibrary.org/obo/mondo/mappings/mondo_exactmatch_omim.sssom.tsv" -O $@

mondo_exactmatch_omimps.sssom.tsv:
	wget "http://purl.obolibrary.org/obo/mondo/mappings/mondo_exactmatch_omimps.sssom.tsv" -O $@

# More commands / goals  -----------------------------------------------------------------------------------------------
# Create mapping artefact(s)
omim.json: omim.owl
	robot convert -i $< -o omim.json

# Create OWL artefact: the graph in functional syntax
omim.owl: omim.nt
	robot convert -i $< -f ofn -o $@

cleanup:
	@rm -f omim.json

# SETUP / INSTALLATION -------------------------------------------------------------------------------------------------
# The ODK container already has everything the build imports. This is for working outside it.
install:
	uv sync

# QA / TESTING ---------------------------------------------------------------------------------------------------------
test:
	 python3 -m pytest tests

lint:
	ruff format --check
	ruff check

# Capture golden outputs from the legacy pipeline as it was before the rewrite. Needs Docker and OMIM's downloads
# OMIM_DATA: where OMIM's licensed .txt downloads are kept. They are never committed.
OMIM_DATA ?= omim-data

golden:
	OMIM_DATA=$(OMIM_DATA) golden/capture.sh

# HELP -----------------------------------------------------------------------------------------------------------------
help:
	@echo "----------------------------------------"
	@echo "	Command reference: OMIM"
	@echo "----------------------------------------"
	@echo "all"
	@echo "Creates all release artefacts.\n"
	@echo "fetch"
	@echo "Downloads the inputs and refreshes the API cache.\n"
	@echo "omim.owl"
	@echo "Builds the main release artefact from files on disk.\n"
	@echo "omim.sssom.tsv"
	@echo "Creates an SSSOM TSV of OMIM terms.\n"
	@echo "install"
	@echo "Installs Python requirements with uv, for working outside the ODK container.\n"
	@echo "test"
	@echo "Runs unit tests.\n"
	@echo "lint"
	@echo "Checks formatting and lint rules.\n"
