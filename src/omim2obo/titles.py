"""Turning OMIM's ALL CAPS title fields into readable titles and symbols."""
import re
from typing import List, Dict, Set, Tuple, Union

from omim2obo.romanplus import fromRoman, romanNumeralPattern, toRoman


def detect_abbreviations(label: str, capitalization_threshold=0.75) -> List[str]:
    """Detect possible abbreviations / acronyms"""
    # Compile regexp
    # todo: handle several warnings: {1} redundant, {1,} simplified to +
    acronyms_without_periods_compiler = re.compile('[A-Z]{1}[A-Z0-9]{1,}')
    # todo: PyCharm flagged next 2 lines as invalid escape sequence, but this code seems to work? Should double check
    acronyms_with_periods_compiler = re.compile(r'[A-Z]{1}\.([A-Z0-9]\.){1,}')
    title_cased_abbrev_compiler = re.compile(r'[A-Z]{1}[a-zA-Z]{1,}\.')

    # Special threshold-based logic, because incoming data has highly inconsistent capitalization
    # is_largely_uppercase = synonym.upper() == synonym  # too many false positives
    fully_capitalized_count = 0
    words = label.split()
    for word in words:
        if word.upper() == word:
            fully_capitalized_count += 1
    is_largely_uppercase = \
        fully_capitalized_count / len(words) >= capitalization_threshold

    # Detect cases
    if is_largely_uppercase:
        acronyms_without_periods = []  # can't infer because everything was uppercase
    else:
        acronyms_without_periods: List[str] = acronyms_without_periods_compiler.findall(label)
    title_cased_abbrevs: List[str] = title_cased_abbrev_compiler.findall(label)
    acronyms_with_periods: List[str] = acronyms_with_periods_compiler.findall(label)

    return acronyms_with_periods + acronyms_without_periods + title_cased_abbrevs


# todo: rename? It's doing more than cleaning; it's mutating
def cleanup_title(
    title: str,
    replacement_case_method: str = 'lower',  # 'upper', 'title', 'lower', 'capitalize' (=sentence case)
    conjunctions: List[str] = ['and', 'but', 'yet', 'for', 'nor', 'so'],
    little_preps: List[str] = ['at', 'by', 'in', 'of', 'on', 'to', 'up', 'as', 'it', 'or'],
    articles: List[str] = ['a', 'an', 'the'],
    word_replacements: Dict[str, str] = None,
) -> str:
    """Reformat the ALL CAPS OMIM labels to something more pleasant to read.

    :param title: A preferred, alternative, or included title.

    1. Converts roman numerals to arabic
    2. Makes the text adhere to the case of `replacement_case_method`, except for supplied
    conjunctions, prepositions, and articles, which will always be lowercased. NOTE: The default for this is 'lower',
    meaning that this operation by default does nothing.

    Assumptions:
    1. All acronyms are capitalized

    todo later's:
    1: Find a pattern for hyphenated types, and maintain acronym capitalization
       e.g. MITF-related melanoma and renal cell carcinoma predisposition syndrome
       e.g. ATP1A3-associated neurological disorder
    2. Make pattern for chromosomes
       agonadism, 46,XY, with intellectual disability, short stature, retarded bone age, and multiple extragenital
         malformations
       Chromosome special formatting capitalization?
       There seems to be special formatting for chromosome refs; they have a comma in the middle, but with no space
       after the comma, though some places I saw on the internet contained a space.
       e.g. "46,XY" in: agonadism, 46,XY, with intellectual disability, short stature, retarded bone age, and multiple
         extragenital malformations
    3. How to find acronym if it is capitalized but only includes char [A-Z], and
        every other char in the string is also capitalized? I don't see a way unless
        checking every word against an explicit dictionary of terms, though there are sure
        to also be (i) acronyms in that dictionary, and (ii) non-acronyms missing from
        that dictionary. And also concern (iii), where to get such an extensive dictionary?
    4. Add "special character" inclusion into acronym regexp. But which special
        chars to include, and which not to include?
    5. Acronym capture extension: case where at least 1 word is not capitalized:
        any word that is fully capitalized might as well be acronym, so long
       as at least 1 other word in the label is not all caps. Maybe not a good rule,
       because there could be some violations, and this probably would not happen
       that often anwyay
        - Not sure what I meant about (5) - joeflack4 2021/09/10
    6. Eponyms: re-capitalize first char?
       e.g.: Balint syndrome, Barre-Lieou syndrome, Wallerian degeneration, etc.
       How to do this? Simply get/create a list of known eponyms? Is this feasible?
    """
    fixedwords = []
    i = 0
    for wrd in title.split():
        i += 1
        # convert the roman numerals to numbers,
        # but assume that the first word is not
        # a roman numeral (this permits things like "X inactivation"
        if i > 1 and re.match(romanNumeralPattern, wrd):
            num = fromRoman(wrd)
            # make the assumption that the number of syndromes are <100
            # this allows me to retain "SYNDROME C"
            # and not convert it to "SYNDROME 100"
            if 0 < num < 100:
                # get the non-roman suffix, if present.
                # for example, IIIB or IVA
                suffix = wrd.replace(toRoman(num), '', 1)
                fixed = ''.join((str(num), suffix))
                wrd = fixed
        wrd = getattr(wrd, replacement_case_method)()
        # replace interior conjunctions, prepositions, and articles with lowercase, always
        if wrd in (conjunctions + little_preps + articles) and i != 1:
            wrd = wrd.lower()
        if word_replacements:
            wrd = word_replacements.get(wrd, wrd)
        fixedwords.append(wrd)
    label_newcase = ' '.join(fixedwords)

    return label_newcase


def recapitalize_acronyms_in_title(title: str, known_abbrevs: Set[str] = None, capitalization_threshold=0.75) -> str:
    """Re-capitalize acronyms / words based on information contained w/in original label

    todo: If title has been used on cleanup_title() using a replacement_case_method other than the non-default 'lower',
     then the .replace() operation will not work. To solve, this (a) capture the replacement_case_method used and
     pass that here, or (b) duplicate the .replace() line and call it on alternative casing variations (.title() and
     capitalize() (=sentence case)), (c) possibly just compare to word.lower() instead of 'word.
    todo: (more important): It's probable that .split(' ') is not enough to cover all cases. Should also run the check
     by splitting on other characters. E.g. consider the following potential cases: "TITLE (ACRONYM)",
     "TITLE: ACRONYM1&ACRONYM2", "TITLE/ACRONYM" or "TITLE ACRONYM/ACRONYM", "TITLE {ACRONYM1,ACRONYM2}",
     "TITLE[ACRONYM]",  "TITLE-ACRONYM", or less likely cases such as "TITLE_ACRONYM", "TITLE.ACRONYM". There are quite
      a few different combos of special char usage that could theoretically arise. It might be possible for thisthat to
      utilize the regular expressions in detect_abbreviations(), and substitute in the acronym in the place of the [A-Z]
      part. It is also possible to improve detect_abbreviations() by considering some of thes eother possible example
      cases above.
    """
    inferred_abbrevs: Set[str] = set(detect_abbreviations(title, capitalization_threshold))
    abbrevs: Set[str] = known_abbrevs.union(inferred_abbrevs)
    if not abbrevs:
        return title
    title2_words: List[str] = []
    for word in title.split():
        abbrev_match = False
        for abbrev in abbrevs:
            if abbrev.lower() == word:
                title2_words.append(abbrev)
                abbrev_match = True
                break
        if not abbrev_match:
            title2_words.append(word)
    title2 = ' '.join(title2_words)
    return title2


def recapitalize_acronyms_in_titles(
    titles: Union[str, List[str]], known_abbrevs: Set[str] = None, capitalization_threshold=0.75
) -> Union[str, List[str]]:
    """Re-capitalize acronyms in a list of titles"""
    if isinstance(titles, str):
        return recapitalize_acronyms_in_title(titles, known_abbrevs, capitalization_threshold)
    return [recapitalize_acronyms_in_title(title, known_abbrevs, capitalization_threshold) for title in titles]


def remove_included_and_formerly_suffixes(title: str) -> str:
    """Remove ', INCLUDED' and ', FORMERLY' suffixes from a title"""
    for suffix in ['FORMERLY', 'INCLUDED']:
        title = re.sub(r',\s*' + suffix, '', title, re.IGNORECASE)
    return title


def separate_former_titles_and_symbols(
    titles: List[str], symbols: List[str]
) -> Tuple[List[str], List[str], List[str], List[str]]:
    """Separate current title/symbols from deprecated (marked 'former') ones"""
    former_titles = [x for x in titles if ', FORMERLY' in x.upper()]
    former_symbols = [x for x in symbols if ', FORMERLY' in x.upper()]
    current_titles = [x for x in titles if ', FORMERLY' not in x.upper()]
    current_symbols = [x for x in symbols if ', FORMERLY' not in x.upper()]
    return current_titles, current_symbols, former_titles, former_symbols


def clean_alt_and_included_titles(
    titles: List[str], symbols: List[str], word_replacements: Dict[str, str]
) -> Tuple[List[str], List[str]]:
    """Remove ', INCLUDED' and ', FORMERLY' suffixes from titles/symbols & misc title reformatting"""
    # remove ', included' and ', formerly', if present
    titles2 = [remove_included_and_formerly_suffixes(x) for x in titles]
    symbols2 = [remove_included_and_formerly_suffixes(x) for x in symbols]
    # additional reformatting for titles
    titles3 = [cleanup_title(x, word_replacements=word_replacements) for x in titles2]
    return titles3, symbols2


def parse_title_symbol_pairs(title_symbol_pairs_str: str) -> Tuple[List[str], List[str]]:
    """Parses a string containing title-symbol pairs.

    :param title_symbol_pairs_str: A string representing title-symbol pairs.
    Format:
    - Pairs are separated by ';;'
    - Within each pair:
        - The first element is always a title
        - Optionally followed by zero or more symbols, separated by ';'

    Examples:
      Positional semantics:
        Title1;Symbol1;Symbol2;;Title2;;Title3;Symbol3
      Alternative Title(s); symbol(s):
        ACROCEPHALOSYNDACTYLY, TYPE V; ACS5;; ACS V;; NOACK SYNDROME
      Included Title(s); symbols:
        CRANIOFACIAL-SKELETAL-DERMATOLOGIC DYSPLASIA, INCLUDED
    """
    titles: List[str] = []
    symbols: List[str] = []
    title_symbol_pairs: List[str] = title_symbol_pairs_str.split(';;')
    for pair_str in title_symbol_pairs:
        pair: List[str] = [x.strip() for x in pair_str.split(';')]
        titles.append(pair[0])
        symbols.extend(pair[1:])
    return titles, symbols


def get_alt_and_included_titles_and_symbols(
    title_symbol_pair_str, word_replacements: Dict[str, str]
) -> Tuple[List[str], List[str], List[str], List[str]]:
    """Separates different types of titles/symbols, and cleans them."""
    titles: List[str] = []
    symbols: List[str] = []
    former_titles: List[str] = []
    former_symbols: List[str] = []
    if title_symbol_pair_str:
        titles, symbols = parse_title_symbol_pairs(title_symbol_pair_str)
        titles, symbols, former_titles, former_symbols = separate_former_titles_and_symbols(titles, symbols)
        titles, symbols = clean_alt_and_included_titles(titles, symbols, word_replacements)
        former_titles, former_symbols = clean_alt_and_included_titles(former_titles, former_symbols, word_replacements)
    return titles, symbols, former_titles, former_symbols
