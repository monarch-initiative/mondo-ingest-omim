"""Reading an entry as OMIM's API returns it."""

from collections import defaultdict

ORPHANET, UMLS = 'orphanet', 'umls'


def get_pubs(entry) -> list[str]:
    """Get pubmed information from an OMIM entry"""
    result = []
    for rlst in entry.get('referenceList', []):
        if 'pubmedID' in rlst['reference']:
            result.append(str(rlst['reference']['pubmedID']))
    return result


def get_mapped_ids(entry) -> dict[str, list[str]]:
    """Get mapped IDs from an OMIM entry"""
    external_links = entry.get('externalLinks', {})
    result = defaultdict(list)
    if 'orphanetDiseases' in external_links:
        for item in (
            external_links['orphanetDiseases'].strip().split(';;;')
        ):  # Separates multiple Orphanet entries, e.g. '634475;;31692;;Mosaic NF2-related schwannomatosis;;;637;;183;;Full NF2-related schwannomatosis'
            result[ORPHANET].append(
                item.split(';;')[0].strip()
            )  # e.g. '634475;;31692;;Mosaic NF2-related schwannomatosis'
    if 'umlsIDs' in external_links:
        result[UMLS] = external_links['umlsIDs'].split(',')
    return result
