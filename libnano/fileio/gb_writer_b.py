# Copyright (C) 2014-2018. Nick Conway & Ben Pruitt; Wyss Institute
# Copyright (C) 2023 Nick Conway & Ben Pruitt;
# See LICENSE.TXT for full GPLv2 license.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
'''
libnano.fileio.gb_writer_b
~~~~~~~~~~~~~~~~~~~~~~~~~~

http://www.insdc.org/files/feature_table.html

section 4.1 shows a 79 character feature table format, as well as the LOCUS
being 80 characters but BioPython files output an 80 character version

However section 4.3 Data item positions shows

    22-80              location

as the 80th character is in play.

This code sticks with 79 character limit
'''
import io
import sys
from typing import (
    Dict,
    Union,
)

import _io  # type: ignore

import libnano.helpers.textwrap as textwrap

NEWLINE_BYT: bytes = b'\r\n' if sys.platform == 'win32' else b'\n'
STR_T = Union[str, bytes]


def write(
        fd: _io.TextIOWrapper,
        d: Dict,
        order_qualifiers: bool = False,
) -> None:
    '''Write dictionary to Genbank file

    Args:
        fd: Genbank filedescriptor
        d: Dictionary to write
        order_qualifiers: If True, write features in order

    '''
    d_info = d[b'info']
    writeLocus(fd, d_info)
    writeDefinition(fd, d_info)
    writeAccession(fd, d_info)
    writeVersion(fd, d_info)
    writeDBLINK(fd, d_info)
    writeKeywords(fd, d_info)
    writeSource(fd, d_info)
    i = 0
    for reference in d_info[b'references']:
        i += 1
        writeReference(fd, reference, i)
    writeComment(fd, d_info)
    writeFeatures(fd, d, order_qualifiers)
    writeOrigin(fd, d)


def write_file(
        filepath: str,
        d: Dict,
        order_qualifiers: bool = False,
) -> None:
    '''Write dictionary to Genbank file

    Args:
        filepath: Genbank filepath
        d: Dictionary to write
        order_qualifiers: If True, write features in order

    '''
    with io.open(filepath, 'wb') as fd:
        write(fd, d, order_qualifiers=order_qualifiers)


def multiline(
        main_str: bytes,
        indent_str: bytes,
        lim: int = 67,
) -> bytes:
    '''
    Args:
        main_str: String to write
        indent_str: Line indent string
        lim: Line width limit to wrap at

    Returns:
        ``main_str`` in multiline format
    '''
    o_list = [main_str[i:i + lim] for i in range(0, len(main_str), lim)]
    return indent_str.join(o_list)


def multiline_spaces(
        main_str: bytes,
        indent_str: bytes,
        lim: int = 67,
) -> bytes:
    '''
    Args:
        main_str: string to join with ``indent_str``
        indent_str: string to lead lines with
        lim: line width

    Returns:
        ``indent_str`` joined with a wrapped ``main_str``
    '''
    main_list = textwrap.wrap(
        main_str,
        width=lim,
        drop_whitespace=True,
        break_on_hyphens=False,
    )
    return indent_str.join(main_list)


def writeLocus(
        fd: _io.TextIOWrapper,
        d: Dict,
) -> None:
    '''Write a locus to a file descritor

    Locus in sample file is not valid

    Args:
        fd: File descriptor
        d: Locus dictionary to write keyed by bytes
    '''
    locus_str = b'LOCUS       '
    name = d[b'name']
    molecule_type = d[b'molecule_type']
    form = d[b'form']
    if form is None:
        form = b''
    elif form == b'linear':
        form = b'linear  '

    gb_division = d[b'gb_division']
    if gb_division is None:
        gb_division = b''

    mod_date = d[b'mod_date']
    if mod_date is None:
        mod_date = b''

    out_str = (
        b'%s%-16s %11d bp    %-6s  %8s %s %s%s' %
        (
            locus_str, name, d[b'length'],
            molecule_type, form, gb_division,
            mod_date, NEWLINE_BYT,
        )
    )
    fd.write(out_str)


def writeDefinition(
        fd: _io.TextIOWrapper,
        d: Dict,
) -> None:
    '''Write a definition to a file descritor

    Args:
        fd: File descriptor
        d: Definition dictionary to write keyed by bytes
    '''
    definition_str = b'DEFINITION  '
    indent_str = b'%s            ' % (NEWLINE_BYT)
    out_list = [
        definition_str,
        multiline_spaces(d[b'definition'], indent_str),
        NEWLINE_BYT,
    ]
    fd.write(b''.join(out_list))


def writeAccession(
        fd: _io.TextIOWrapper,
        d: Dict,
) -> None:
    '''Write an accession to a file descritor

    Args:
        fd: File descriptor
        d: Accession dictionary to write keyed by bytes
    '''
    accession_str = b'ACCESSION   %s%s'
    fd.write(accession_str % (d[b'accession'], NEWLINE_BYT))


def writeVersion(
        fd: _io.TextIOWrapper,
        d: Dict,
) -> None:
    '''Write a version to a file descritor

    Args:
        fd: File descriptor
        d: Version dictionary to write keyed by bytes
    '''
    version = d[b'version']
    if version is not None:
        gi = d[b'GI']
        fd.write(b'VERSION     %s  GI:%s%s' % (version, gi, NEWLINE_BYT))


def writeDBLINK(
        fd: _io.TextIOWrapper,
        d: Dict,
) -> None:
    '''Write a DBLink to a file descritor

    Args:
        fd: File descriptor
        d: DBLink dictionary to write keyed by bytes
    '''
    if b'dblink' in d:
        if d[b'dblink'] is not None:
            fd.write(b'DBLINK      %s%s' % (d[b'dblink'], NEWLINE_BYT))


def writeKeywords(
        fd: _io.TextIOWrapper,
        d: Dict,
) -> None:
    '''Write a Keywords to a file descritor

    Args:
        fd: File descriptor
        d: Keywords dictionary to write keyed by bytes
    '''
    keywords_str = b'KEYWORDS    %s%s'
    fd.write(keywords_str % (d[b'keywords'], NEWLINE_BYT))


def writeSource(
        fd: _io.TextIOWrapper,
        d: Dict,
) -> None:
    '''Write a Source to a file descritor

    Args:
        fd: File descriptor
        d: Source dictionary to write keyed by bytes
    '''
    source_str = b'SOURCE      '
    organism_str = b'%s  ORGANISM  ' % (NEWLINE_BYT)
    indent_str = b'%s            ' % (NEWLINE_BYT)
    org = d[b'organism']
    out_list = [
        source_str, d[b'source'],
        organism_str,
        org[0],
    ]
    if org[1] is not None:
        out_list += [
            indent_str,
            multiline_spaces(org[1], indent_str),
        ]
    out_list += [NEWLINE_BYT]
    fd.write(b''.join(out_list))


def writeReference(
        fd: _io.TextIOWrapper,
        ref: Dict,
        i: int,
) -> None:
    '''Write a Reference to a file descritor

    Args:
        fd: File descriptor
        ref: Reference dictionary to write keyed by bytes
        i: index for the reference
    '''
    assert (ref[b'r_index'] == i)
    reference_str = b'REFERENCE   '
    authors_str = b'  AUTHORS   '
    title_str = b'  TITLE     '
    indent_str = b'%s            ' % (NEWLINE_BYT)
    journal_str = b'  JOURNAL   '
    pubmed_str = b'  PUBMED    '
    if ref[b'start_idx'] is not None:
        idx_str = b'  (bases %d to %d)' % (ref[b'start_idx'], ref[b'end_idx'])
    else:
        idx_str = b''
    out_list = [
        reference_str,
        b'%d%s%s' % (i, idx_str, NEWLINE_BYT),
        authors_str,
        multiline_spaces(ref[b'authors'], indent_str), NEWLINE_BYT,
        title_str,
        multiline_spaces(ref[b'title'], indent_str), NEWLINE_BYT,
        journal_str,
        multiline_spaces(ref[b'journal_info'], indent_str), NEWLINE_BYT,
    ]
    if ref[b'pubmed'] is not None:
        out_list += [pubmed_str, ref[b'pubmed'], NEWLINE_BYT]

    fd.write(b''.join(out_list))


def writeComment(
        fd: _io.TextIOWrapper,
        d: Dict,
) -> None:
    '''Write a Comment to a file descritor

    Args:
        fd: File descriptor
        d: Comment dictionary to write keyed by bytes
    '''
    comment_str = b'COMMENT     '
    indent_str = b'%s            ' % NEWLINE_BYT
    # indent_str_gb_asm = b'            '
    if b'comment' in d:
        comment = d[b'comment']
        if isinstance(comment, list):
            out_list = [
                comment_str,
                multiline_spaces(comment[0], indent_str),
                indent_str, indent_str, indent_str,
                indent_str.join(comment[1]), NEWLINE_BYT,
            ]
        else:
            out_list = [
                comment_str, multiline_spaces(
                    comment, indent_str,
                ), NEWLINE_BYT,
            ]
        fd.write(b''.join(out_list))


def writeFeatures(
        fd: _io.TextIOWrapper,
        d: Dict,
        order_qualifiers: bool,
) -> None:
    '''Write Features to a file descritor

    Args:
        fd: File descriptor
        d: Features dictionary to write keyed by bytes
        order_qualifiers: If True, order by the `'qualifiers'` key
    '''
    feature_header = b'FEATURES             Location/Qualifiers%s' % (
        NEWLINE_BYT
    )
    feature_type_prefix_str = b'     '
    indent_str = b'%s                     ' % (NEWLINE_BYT)
    feature_type_field_size = 16
    fd.write(feature_header)
    for feature in d[b'features']:
        ftype = feature[b'type']
        location_str = feature[b'location']
        out_list = [
            feature_type_prefix_str,
            ftype, (feature_type_field_size - len(ftype)) * b' ',
            location_str, NEWLINE_BYT,
        ]
        if order_qualifiers:
            quals = feature[b'qualifiers']
            for key in sorted(quals.keys()):
                value_list = quals[key]
                if not isinstance(value_list, list):
                    value_list = [value_list]
                # assumes value_list is all the same type
                for value in sorted(value_list):
                    # codon start is a weird edge case
                    if key in (b'codon_start', b'transl_table'):
                        out_list += [
                            b'                     /%s=%d%s' %
                            (key, value, NEWLINE_BYT),
                        ]
                    elif key == b'transl_except':
                        out_list += [
                            b'                     /%s=%s%s' %
                            (key, value, NEWLINE_BYT),
                        ]
                    elif key == b'translation':
                        qualifier_str = b"/%s=\"%s\"" % (key, value)
                        qual_list = multiline(
                            qualifier_str, indent_str, lim=58,
                        )  # or 58
                        out_list += [
                            b'                     ' +
                            qual_list, NEWLINE_BYT,
                        ]
                    else:
                        if value is None:
                            qualifier_str = b'/%s%s' % (key, NEWLINE_BYT)
                        else:
                            qualifier_str = b"/%s=\"%s\"%s" % (
                                key, value, NEWLINE_BYT,
                            )
                        qual_list = multiline_spaces(
                            qualifier_str, indent_str, lim=58,
                        )  # or 58
                        out_list += [
                            b'                     ' +
                            qual_list, NEWLINE_BYT,
                        ]
        else:
            for key, value_list in feature[b'qualifiers'].items():
                if not isinstance(value_list, list):
                    value_list = [value_list]
                for value in value_list:
                    # codon start is a weird edge case
                    if key in (b'codon_start', b'transl_table'):
                        out_list += [
                            b'                     /%s=%d%s' %
                            (key, value, NEWLINE_BYT),
                        ]
                    elif key == b'transl_except':
                        out_list += [
                            b'                     /%s=%s%s' %
                            (key, value, NEWLINE_BYT),
                        ]
                    elif key == b'translation':
                        qualifier_str = b"/%s=\"%s\"" % (key, value)
                        qual_list = multiline(
                            qualifier_str, indent_str, lim=58,
                        )  # or 58
                        out_list += [
                            b'                     ' +
                            qual_list, NEWLINE_BYT,
                        ]
                    else:
                        if value is None:
                            qualifier_str = b'/%s%s' % (key, NEWLINE_BYT)
                        else:
                            qualifier_str = b"/%s=\"%s\"%s" % (
                                key, value, NEWLINE_BYT,
                            )
                        qual_list = multiline_spaces(
                            qualifier_str, indent_str, lim=58,
                        )  # or 58
                        out_list += [
                            b'                     ' +
                            qual_list, NEWLINE_BYT,
                        ]
        fd.write(b''.join(out_list))


def writeOrigin(
        fd: _io.TextIOWrapper,
        d: Dict,
) -> None:
    '''Write Origin to a file descritor

    Args:
        fd: File descriptor
        d: Origin dictionary to write keyed by bytes
    '''
    i = 0
    format_string = b'%9d %s %s %s %s %s %s%s'
    origin_header = b'ORIGIN      %s' % (NEWLINE_BYT)
    fd.write(origin_header)
    origin = d[b'seq']
    count = len(origin) - 60
    while i < count:
        str_tup = ((i + 1,) +
                   tuple([origin[j:j + 10] for j in range(i, i + 60, 10)]) +
                   (NEWLINE_BYT,))
        fd.write(format_string % (str_tup))
        i += 60
    last_str = origin[i:]
    last_str = last_str + (60 - (len(last_str))) * b' '
    last_list = [last_str[j:j + 10] for j in range(0, 60, 10)]
    last_list = [x.rstrip() for x in last_list if x[0] != b' '[0]]
    str_tup = (i + 1,) + tuple(last_list)
    format_string = b'%9d' + b' %s' * len(last_list) + NEWLINE_BYT
    fd.write(format_string % (str_tup))
    fd.write(b'//%s' % (NEWLINE_BYT))
