"""OMIM Types"""

from enum import Enum


class EntryPrefix(Enum):
    """From OMIM docs: What do the symbols preceding a MIM number represent?

    *: An asterisk (*) before an entry number indicates a **gene**.

    #: A number symbol (#) before an entry number indicates that it is a **descriptive**
    entry, usually of a phenotype, and does not represent a unique locus. The
    reason for the use of the number symbol is given in the first paragraph of
    the entry. Discussion of any gene(s) related to the phenotype resides in
    another entry(ies) as described in the first paragraph.

    +: A plus sign (+) before an entry number indicates that the entry contains the
    description of a gene of known **sequence** and a **phenotype**.

    %: A percent sign (%) before an entry number indicates that the entry describes
    a confirmed mendelian **phenotype** or phenotypic locus for which the underlying
    molecular basis is not known.

    '': No symbol before an entry number generally indicates a **description** of a phenotype
    for which the mendelian basis, although suspected, has not been clearly established
    or that the separateness of this phenotype from that in another entry is unclear.

    ^: A caret (^) before an entry number means the entry no longer exists because
    it was **removed** from the database or moved to another entry as indicated.

    ----
    See also:
    - the description of symbols used in the disorder column of the OMIM Gene Map and Morbid Map.
    - data/symbolic_prefixes.tsv for a tabular version of this information
    """

    GENE = 'Asterisk'  # *
    GENE_AND_PHENOTYPE = 'Plus'  # +
    PHENOTYPE = 'Number Sign'  # #
    PHENOTYPE_OR_LOCUS = 'Percent'  # %
    OTHER = 'NULL'  # ''


BIOLINK_CATEGORY = {
    EntryPrefix.GENE: 'Gene',
    EntryPrefix.GENE_AND_PHENOTYPE: 'Gene',
    EntryPrefix.PHENOTYPE: 'Disease',
    EntryPrefix.PHENOTYPE_OR_LOCUS: 'Disease',
}
