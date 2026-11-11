"""Parsers: one function per input file, each turning lines of text into immutable records.

This is the only layer that sees raw text, and so the only one that checks it. A record that has
been constructed is well-formed; nothing downstream validates again.
"""
