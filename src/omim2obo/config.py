"""Configuration"""
from pathlib import Path

from dotenv import dotenv_values

ROOT_DIR = Path(__file__).resolve().parents[2]  # the repository root
DATA_DIR = ROOT_DIR / 'data'
ENV_PATH = ROOT_DIR / '.env'
DISEASE_GENE_PROTECTED_PATH = DATA_DIR / 'protected-disease-gene.tsv'
PUBMED_REFS_PATH = DATA_DIR / 'pubmed-refs.tsv'
MAPPINGS_PATH = DATA_DIR / 'mappings.tsv'
CACHE_LAST_UPDATED_PATH = DATA_DIR / 'cache-last-updated.txt'
CACHE_INCOMPLETENESS_INDICATOR_PATH = DATA_DIR / 'initial-cache-incomplete.txt'

CONFIG = dotenv_values(ENV_PATH)
