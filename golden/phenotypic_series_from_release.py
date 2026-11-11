"""Reconstruct phenotypicSeries.txt from a released omim.owl, for accounts not entitled to it.

    gh release download -R monarch-initiative/omim -p omim.owl
    python3 golden/phenotypic_series_from_release.py omim.owl > omim-data/phenotypicSeries.txt

The pipeline takes three things from that file: each series' number, its title, and its member MIMs.
All three are in the released OWL, so the file can be rebuilt to give the same axioms. The third
column of member rows, a phenotype label the pipeline never reads, is filled with a placeholder.
"""

import re
import sys
from collections import defaultdict

label_re = re.compile(
    r'^AnnotationAssertion\(rdfs:label <https://omim\.org/phenotypicSeries/(PS\d+)> '
    r'"((?:[^"\\]|\\.)*)"(?:\^\^xsd:string)?\)$'
)
member_re = re.compile(
    r'^SubClassOf\(<https://omim\.org/entry/(\d+)> <https://omim\.org/phenotypicSeries/(PS\d+)>\)$'
)
decl_re = re.compile(r'^Declaration\(Class\(<https://omim\.org/phenotypicSeries/(PS\d+)>\)\)$')
labels, members, declared = {}, defaultdict(list), set()
for line in open(sys.argv[1]):
    line = line.rstrip('\n')
    if m := label_re.match(line):
        labels[m[1]] = re.sub(r'\\(.)', r'\1', m[2])
    elif m := member_re.match(line):
        members[m[2]].append(m[1])
    elif m := decl_re.match(line):
        declared.add(m[1])
assert declared == set(labels), declared ^ set(labels)
assert set(members) <= declared
print('# Reconstructed from the monarch-initiative/omim release omim.owl; not an OMIM download')
print('# Phenotypic Series Number\tMIM Number\tPhenotype')
for ps in sorted(labels):
    print(f'{ps}\t{labels[ps]}')
    for mim in sorted(members[ps]):
        print(f'{ps}\t{mim}\t-')
print(f'{len(labels)} series, {sum(map(len, members.values()))} memberships', file=sys.stderr)
