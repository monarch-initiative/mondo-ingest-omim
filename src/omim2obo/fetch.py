"""Fetching what a build needs from OMIM: its downloadable files, and for each entry the PubMed
references and UMLS and Orphanet mappings that only its API gives."""

import os
from datetime import datetime, timedelta
from pathlib import PosixPath

import pandas as pd
import requests

from omim2obo.api_entry import ORPHANET, UMLS, get_mapped_ids, get_pubs
from omim2obo.config import (
    CACHE_INCOMPLETENESS_INDICATOR_PATH,
    CACHE_LAST_UPDATED_PATH,
    CONFIG,
    DATA_DIR,
    DISEASE_GENE_PROTECTED_PATH,
    MAPPINGS_PATH,
    PUBMED_REFS_PATH,
)
from omim2obo.omim_client import OmimClient
from omim2obo.parse.curation import parse_protected
from omim2obo.parse.omim import parse_mim_titles, parse_morbidmap

MIM_FILES = ('mimTitles', 'genemap2', 'mim2gene', 'phenotypicSeries', 'morbidmap')


def fetch():
    """Bring data/ up to date with OMIM. Needs the keys in .env."""
    for file_name in MIM_FILES:
        download_mim_file(file_name)
    update_cache__pubmed_refs_and_mappings()


def download_mim_file(file_name: str):
    """Retrieve OMIM downloadable text file from the OMIM download server"""
    file_name = file_name if file_name.endswith('.txt') else file_name + '.txt'
    mim_file_path: PosixPath = DATA_DIR / file_name

    print(f'Downloading {file_name} from OMIM...')
    url = f'https://data.omim.org/downloads/{CONFIG["DOWNLOAD_KEY"]}/{file_name}'
    # todo: This doesn't work for genemap2.txt. But does the previous URL work? If so, why not just
    #  use that?
    if file_name == 'mim2gene.txt':
        url = f'https://omim.org/static/omim/data/{file_name}'
    resp = requests.get(url)
    if resp.status_code == 200:
        text = resp.text
        if not text.startswith('<!DOCTYPE html>'):
            # Save
            with open(mim_file_path, 'w') as fout:
                fout.write(text)
        else:
            raise RuntimeError('Unexpected response: ' + text)
    else:
        msg = 'Response from server: ' + resp.text
        # LOG.warning(msg)
        # with open(mim_file, 'r') as fin:
        #     lines = fin.readlines()
        # LOG.warning('Failed to retrieve mimTitles.txt. Using the cached file.')
        raise RuntimeError(msg)


def get_all_phenotype_mims() -> set[str]:
    """Get all phenotype MIM numbers"""
    # Protected phenotypes count too, whether or not OMIM lists them in morbidmap.
    with (
        open(DATA_DIR / 'morbidmap.txt') as morbidmap,
        open(DISEASE_GENE_PROTECTED_PATH) as protected,
    ):
        p_mims = {row.phenotype_mim for row in parse_morbidmap(morbidmap)}
        p_mims |= {row.phenotype_mim for row in parse_protected(protected)}
    p_mims.discard(None)
    return p_mims


def _read_cached_entry_df(path: str) -> pd.DataFrame:
    """Read a dataframe containing MIM entries data cached from the OMIM API"""
    return pd.read_csv(
        path,
        sep='\t',
        true_values=['True'],
        false_values=['False'],
        dtype={
            'is_phenotype': bool,
            'mim': str,
            'umls_ids': str,
            'orphanet_ids': str,
            'date_fetched': str,
        },
    )


def update_cache__pubmed_refs_and_mappings(phenotypes_only_for_cache_init=False, overwrite=False):
    """Update cache for MIM entries (pubmed refs & mappings) if cache is not complete or there is
    possibly new data
    """
    # Load existing data
    mims_phenos: set[str] = get_all_phenotype_mims()
    mappings_df_cached = (
        _read_cached_entry_df(MAPPINGS_PATH)
        if os.path.exists(MAPPINGS_PATH) and not overwrite
        else pd.DataFrame()
    )
    pubmed_df_cached = (
        _read_cached_entry_df(PUBMED_REFS_PATH)
        if os.path.exists(PUBMED_REFS_PATH) and not overwrite
        else pd.DataFrame()
    )

    # Fetch
    client = OmimClient(api_key=CONFIG['API_KEY'])
    # - Fetch everything if no cache or cache incomplete
    if not os.path.exists(CACHE_LAST_UPDATED_PATH) or os.path.exists(
        CACHE_INCOMPLETENESS_INDICATOR_PATH
    ):
        print('Cache for pubmed references and mappings is incomplete.')
        # Determine MIMs to fetch
        # - Get all MIMs
        if phenotypes_only_for_cache_init:
            mims_all = mims_phenos
        else:
            with open(DATA_DIR / 'mimTitles.txt') as titles:
                mims_all = {row.mim for row in parse_mim_titles(titles)}
        mims_all.discard('')
        # - Get cached MIMs
        mims_cached: set[str] = (
            set(mappings_df_cached['mim']) if len(mappings_df_cached) > 0 else set()
        )
        mims_cached |= set(pubmed_df_cached['mim']) if len(pubmed_df_cached) > 0 else set()
        # Fetch
        results: list[dict] = client.fetch(
            ids=list(mims_all - mims_cached), update_cache_metadata=True, seed_run=True
        )
    # - Else fetch new data if available
    else:
        print('Checking for recently updated MIMs.')
        # Applies precaution: -1 day from last fetched: to account for MIMs updated later in the day
        with open(CACHE_LAST_UPDATED_PATH) as f:
            last_updated_str = f.readline().strip()
        last_updated: datetime = datetime.strptime(last_updated_str, '%Y-%m-%d')
        last_updated = last_updated - timedelta(days=1)
        results: list[dict] = client.fetch(since_date=last_updated, update_cache_metadata=True)

    # Save
    if len(results) == 0:
        return
    # - Create dataframes from fetched data
    mappings_rows: list[dict] = []
    pubmed_rows: list[dict] = []
    for entry in results:
        mim = str(entry['mimNumber'])
        mappings = get_mapped_ids(entry)
        common_data = {
            'mim': str(mim),
            'is_phenotype': str(mim) in mims_phenos,
            'date_fetched': datetime.now().strftime('%Y-%m-%d'),
        }
        mappings_rows.append(
            {
                **common_data,
                **{
                    'umls_ids': '|'.join(mappings[UMLS]),
                    'orphanet_ids': '|'.join(mappings[ORPHANET]),
                },
            }
        )
        pubmed_rows.append(
            {
                **common_data,
                **{
                    'pmid_refs': '|'.join(get_pubs(entry)),
                },
            }
        )
    mappings_df_new = pd.DataFrame(mappings_rows)
    pubmed_df_new = pd.DataFrame(pubmed_rows)
    # Update cache & save
    # - remove old data from cache if new data has been fetched
    mappings_df_cached_del_old = mappings_df_cached[
        ~mappings_df_cached['mim'].isin(mappings_df_new['mim'])
    ]
    pubmed_df_cached_del_old = pubmed_df_cached[~pubmed_df_cached['mim'].isin(pubmed_df_new['mim'])]
    # - concat & save
    mappings_df = pd.concat(
        [mappings_df_cached_del_old, mappings_df_new], ignore_index=True
    ).sort_values(['mim'])
    pubmed_df = pd.concat([pubmed_df_cached_del_old, pubmed_df_new], ignore_index=True).sort_values(
        ['mim']
    )
    mappings_df.to_csv(MAPPINGS_PATH, sep='\t', index=False)
    pubmed_df.to_csv(PUBMED_REFS_PATH, sep='\t', index=False)
