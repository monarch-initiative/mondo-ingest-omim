"""Characterisation tests: these pin what the title functions do today, on invented titles.

QUIRK marks behaviour that is probably unintended.
"""
import pytest

from omim2obo.parsers.omim_entry_parser import (
    cleanup_title, detect_abbreviations, get_alt_and_included_titles_and_symbols, parse_title_symbol_pairs,
    recapitalize_acronyms_in_title, recapitalize_acronyms_in_titles, remove_included_and_formerly_suffixes,
    separate_former_titles_and_symbols)

EPONYMS = {'behcet': 'Behcet', 'danlos': 'Danlos'}


@pytest.mark.parametrize('title, cleaned', [
    ('WOBBLY KNEE SYNDROME', 'wobbly knee syndrome'),
    ('Mixed Case Title', 'mixed case title'),
    ('WOBBLY KNEE SYNDROME 1A', 'wobbly knee syndrome 1a'),
    ('WOBBLY KNEE SYNDROME, X-LINKED', 'wobbly knee syndrome, x-linked'),
])
def test_cleanup_lowercases(title, cleaned):
    assert cleanup_title(title, word_replacements=EPONYMS) == cleaned


def test_cleanup_collapses_whitespace():
    assert cleanup_title('WOBBLY  KNEE   SYNDROME', word_replacements=EPONYMS) == 'wobbly knee syndrome'


@pytest.mark.parametrize('title, cleaned', [
    ('WOBBLY KNEE SYNDROME, TYPE IV', 'wobbly knee syndrome, type 4'),
    ('WOBBLY KNEE SYNDROME XCIX', 'wobbly knee syndrome 99'),
    ('WOBBLY KNEE SYNDROME IIIB', 'wobbly knee syndrome 3b'),
    ('WOBBLY KNEE SYNDROME IA', 'wobbly knee syndrome 1a'),
])
def test_cleanup_converts_roman_numerals_keeping_a_suffix_letter(title, cleaned):
    assert cleanup_title(title, word_replacements=EPONYMS) == cleaned


@pytest.mark.parametrize('title, cleaned', [
    ('X INACTIVATION, FAMILIAL SKEWED', 'x inactivation, familial skewed'),
    ('IV DRIP SYNDROME', 'iv drip syndrome'),
])
def test_cleanup_never_reads_the_first_word_as_a_roman_numeral(title, cleaned):
    assert cleanup_title(title, word_replacements=EPONYMS) == cleaned


@pytest.mark.parametrize('title, cleaned', [
    ('WOBBLY KNEE SYNDROME C', 'wobbly knee syndrome c'),      # 100
    ('WOBBLY KNEE SYNDROME CD', 'wobbly knee syndrome cd'),    # 400
    ('VITAMIN D DEFICIENCY', 'vitamin d deficiency'),          # 500
    ('WOBBLY KNEE SYNDROME M', 'wobbly knee syndrome m'),      # 1000
    ('DI GEORGE SYNDROME', 'di george syndrome'),              # 501
])
def test_cleanup_leaves_roman_numerals_of_100_or_more(title, cleaned):
    assert cleanup_title(title, word_replacements=EPONYMS) == cleaned


@pytest.mark.parametrize('title, cleaned', [
    ('WOBBLY KNEE SYNDROME X', 'wobbly knee syndrome 10'),
    ('WOBBLY KNEE SYNDROME L', 'wobbly knee syndrome 50'),
    ('WOBBLY KNEE SYNDROME I', 'wobbly knee syndrome 1'),
])
def test_cleanup_QUIRK_single_letters_below_100_are_roman_numerals(title, cleaned):
    assert cleanup_title(title, word_replacements=EPONYMS) == cleaned


@pytest.mark.parametrize('title, cleaned', [
    ('WOBBLY KNEE SYNDROME, TYPE IV, AUTOSOMAL', 'wobbly knee syndrome, type iv, autosomal'),
    ('WOBBLY KNEE SYNDROME (VI)', 'wobbly knee syndrome (vi)'),
    ('WOBBLY KNEE SYNDROME VI; WKS6', 'wobbly knee syndrome vi; wks6'),
])
def test_cleanup_QUIRK_roman_numeral_touching_punctuation_is_not_converted(title, cleaned):
    assert cleanup_title(title, word_replacements=EPONYMS) == cleaned


def test_cleanup_restores_known_capitalisations():
    assert cleanup_title('BEHCET SYNDROME', word_replacements=EPONYMS) == 'Behcet syndrome'


@pytest.mark.parametrize('title, cleaned', [
    ('EHLERS-DANLOS SYNDROME', 'ehlers-danlos syndrome'),
    ('SYNDROME OF DANLOS, TYPE II', 'syndrome of danlos, type 2'),
])
def test_cleanup_QUIRK_known_capitalisations_apply_only_to_whole_words(title, cleaned):
    assert cleanup_title(title, word_replacements=EPONYMS) == cleaned


@pytest.mark.parametrize('label, found', [
    ('wobbly knee syndrome WKS', ['WKS']),
    ('Wobbly knee syndrome, WKS1-related', ['WKS1']),
    ('ABC ABC def ghi jkl', ['ABC', 'ABC']),
    ('Wobbly St. knee', ['St.']),
    ('wobbly knee', []),
])
def test_detect_abbreviations(label, found):
    assert detect_abbreviations(label) == found


@pytest.mark.parametrize('label', ['WOBBLY KNEE SYNDROME WKS', 'A B c d'])
def test_detect_abbreviations_finds_none_when_most_words_are_capitals(label):
    assert detect_abbreviations(label) == []


def test_detect_abbreviations_QUIRK_dotted_acronym_yields_only_its_last_letter():
    assert detect_abbreviations('syndrome of U.S.A. origin') == ['A.']


@pytest.mark.parametrize('title, symbols, recapitalized', [
    ('wobbly knee syndrome wks', {'WKS'}, 'wobbly knee syndrome WKS'),
    ('dna repair wobble', {'DNA'}, 'DNA repair wobble'),
    ('wobbly knee syndrome 2 wks2', {'WKS', 'WKS2'}, 'wobbly knee syndrome 2 WKS2'),
    ('wobbly knee syndrome wks', set(), 'wobbly knee syndrome wks'),
])
def test_recapitalize_restores_the_entrys_symbols(title, symbols, recapitalized):
    assert recapitalize_acronyms_in_title(title, symbols) == recapitalized


def test_recapitalize_keeps_capitals_already_in_the_title():
    assert recapitalize_acronyms_in_title('the DNA wobble of wks', {'WKS'}) == 'the DNA wobble of WKS'


@pytest.mark.parametrize('title', ['wobbly knee syndrome, wks, familial', 'wks1-related wobbly knee'])
def test_recapitalize_QUIRK_symbol_touching_punctuation_is_not_restored(title):
    assert recapitalize_acronyms_in_title(title, {'WKS', 'WKS1'}) == title


def test_recapitalize_accepts_one_title_or_several():
    assert recapitalize_acronyms_in_titles('wobble wks', {'WKS'}) == 'wobble WKS'
    assert recapitalize_acronyms_in_titles(['wobble wks', 'knee'], {'WKS'}) == ['wobble WKS', 'knee']


@pytest.mark.parametrize('field, titles, symbols', [
    ('ONLY TITLE', ['ONLY TITLE'], []),
    ('KNEE WOBBLE;; GENU INSTABILE; GI', ['KNEE WOBBLE', 'GENU INSTABILE'], ['GI']),
    ('A; S1; S2;; B;; C; S3', ['A', 'B', 'C'], ['S1', 'S2', 'S3']),
    ('TITLE ;  SYM ;;  T2', ['TITLE', 'T2'], ['SYM']),
])
def test_title_symbol_pairs(field, titles, symbols):
    assert parse_title_symbol_pairs(field) == (titles, symbols)


@pytest.mark.parametrize('title, stripped', [
    ('WOBBLY KNEE, INCLUDED', 'WOBBLY KNEE'),
    ('WOBBLY KNEE, FORMERLY', 'WOBBLY KNEE'),
    ('WOBBLY KNEE,INCLUDED', 'WOBBLY KNEE'),
    ('A, FORMERLY, INCLUDED', 'A'),
    ('WOBBLY KNEE', 'WOBBLY KNEE'),
])
def test_suffixes_are_removed(title, stripped):
    assert remove_included_and_formerly_suffixes(title) == stripped


def test_suffix_removal_QUIRK_is_case_sensitive():
    # re.IGNORECASE is passed where re.sub takes `count`
    assert remove_included_and_formerly_suffixes('wobbly knee, included') == 'wobbly knee, included'


@pytest.mark.parametrize('title, stripped', [
    ('A, INCLUDED, INCLUDED, INCLUDED', 'A, INCLUDED'),
    ('A, FORMERLY, FORMERLY, FORMERLY', 'A, FORMERLY'),
])
def test_suffix_removal_QUIRK_stops_after_two(title, stripped):
    assert remove_included_and_formerly_suffixes(title) == stripped


@pytest.mark.parametrize('title, stripped', [
    ('WOBBLY, INCLUDED KNEE', 'WOBBLY KNEE'),
    ('KNEE, INCLUDEDNESS', 'KNEENESS'),
])
def test_suffix_removal_QUIRK_is_not_anchored_to_the_end(title, stripped):
    assert remove_included_and_formerly_suffixes(title) == stripped


def test_former_titles_and_symbols_are_separated_from_current_ones():
    assert separate_former_titles_and_symbols(
        ['KNEE WOBBLE', 'OLD WOBBLE, FORMERLY'], ['KW', 'OW, FORMERLY'],
    ) == (['KNEE WOBBLE'], ['KW'], ['OLD WOBBLE, FORMERLY'], ['OW, FORMERLY'])


def test_field_is_split_sorted_and_cleaned():
    field = ('WOBBLY KNEE SYNDROME, TYPE II, INCLUDED; WKS2, INCLUDED'
             ';; OLD WOBBLE, FORMERLY; OW, FORMERLY'
             ';; KNEE WOBBLE')
    assert get_alt_and_included_titles_and_symbols(field) == (
        ['wobbly knee syndrome, type 2', 'knee wobble'], ['WKS2'], ['old wobble'], ['OW'])


def test_empty_field_has_no_titles():
    assert get_alt_and_included_titles_and_symbols('') == ([], [], [], [])


def test_field_QUIRK_mixed_case_formerly_is_filed_as_former_but_keeps_its_suffix():
    assert get_alt_and_included_titles_and_symbols('Wobbly knee, formerly; wk, formerly') == (
        [], [], ['wobbly knee, formerly'], ['wk, formerly'])
