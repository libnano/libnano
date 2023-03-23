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
libnano.fileio.gb_writer
~~~~~~~~~~~~~~~~~~~~~~~~

http://www.insdc.org/files/feature_table.html

section 4.1 shows a 79 character feature table format, as well as the LOCUS
being 80 characters but BioPython files output an 80 character version

However section 4.3 Data item positions shows

    22-80              location

as the 80th character is in play.

This code sticks with 79 character limit
'''
import io
import textwrap

import _io  # type: ignore


def write(
        fd: _io.TextIOWrapper,
        d: dict,
        order_qualifiers: bool = False,
) -> None:
    '''Write dictionary to Genbank file

    Args:
        fd: Genbank filedescriptor
        d: Dictionary to write
        order_qualifiers: If True, write features in order

    '''
    d_info = d['info']
    writeLocus(fd, d_info)
    writeDefinition(fd, d_info)
    writeAccession(fd, d_info)
    writeVersion(fd, d_info)
    writeDBLINK(fd, d_info)
    writeKeywords(fd, d_info)
    writeSource(fd, d_info)
    i = 0
    for reference in d_info['references']:
        i += 1
        writeReference(fd, reference, i)
    writeComment(fd, d_info)
    writeFeatures(fd, d, order_qualifiers)
    writeOrigin(fd, d)


def write_file(
        filepath: str,
        d: dict,
        order_qualifiers: bool = False,
) -> None:
    '''Write dictionary to Genbank file

    Args:
        filepath: Genbank filepath
        d: Dictionary to write
        order_qualifiers: If True, write features in order

    '''
    with io.open(filepath, 'w', encoding='utf-8') as fd:
        write(fd, d, order_qualifiers=order_qualifiers)


def multiline(
        main_str: str,
        indent_str: str,
        lim: int = 67,
) -> str:
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
        main_str: str,
        indent_str: str,
        lim: int = 67,
) -> str:
    '''
    Args:
        main_str: string to join with ``indent_str``
        indent_str: string to lead lines with
        lim: line width

    Returns:
        ``indent_str`` joined with a wrapped ``main_str``
    '''
    main_list = textwrap.wrap(
        main_str, lim, drop_whitespace=True,
        break_on_hyphens=False,
    )
    return indent_str.join(main_list)


def writeLocus(
        fd: _io.TextIOWrapper,
        d: dict,
) -> None:
    '''Write a locus to a file descritor

    Locus in sample file is not valid

    Args:
        fd: File descriptor
        d: Locus dictionary to write keyed by bytes
    '''
    locus_str = 'LOCUS       '
    name = d['name']
    molecule_type = d['molecule_type']
    form = d['form']
    if form is None:
        form = ''
    elif form == 'linear':
        form = 'linear  '

    gb_division = d['gb_division']
    if gb_division is None:
        gb_division = ''

    mod_date = d['mod_date']
    if mod_date is None:
        mod_date = ''

    out_str = locus_str + '%-16s %11d bp    %-6s  %8s %s %s\n' %\
        (name, d['length'], molecule_type, form, gb_division, mod_date)
    fd.write(out_str)


def writeDefinition(
        fd: _io.TextIOWrapper,
        d: dict,
) -> None:
    '''Write a definition to a file descritor

    Args:
        fd: File descriptor
        d: Definition dictionary to write keyed by bytes
    '''
    definition_str = 'DEFINITION  '
    indent_str = '\n            '
    out_list = [
        definition_str,
        multiline_spaces(d['definition'], indent_str),
        '\n',
    ]
    fd.write(''.join(out_list))


def writeAccession(
        fd: _io.TextIOWrapper,
        d: dict,
) -> None:
    '''Write an accession to a file descritor

    Args:
        fd: File descriptor
        d: Accession dictionary to write keyed by bytes
    '''
    accession_str = 'ACCESSION   '
    fd.write(''.join([accession_str, d['accession'], '\n']))


def writeVersion(
        fd: _io.TextIOWrapper,
        d: dict,
) -> None:
    '''Write a version to a file descritor

    Args:
        fd: File descriptor
        d: Version dictionary to write keyed by bytes
    '''
    version_str = 'VERSION     '
    version = d['version']
    if version is not None:
        gi = d['GI']
        fd.write(''.join([version_str, '%s  GI:%s' % (version, gi), '\n']))


def writeDBLINK(
        fd: _io.TextIOWrapper,
        d: dict,
) -> None:
    '''Write a DBLink to a file descritor

    Args:
        fd: File descriptor
        d: DBLink dictionary to write keyed by bytes
    '''
    if 'dblink' in d:
        if d['dblink'] is not None:
            fd.write('DBLINK      %s\n' % (d['dblink']))


def writeKeywords(
        fd: _io.TextIOWrapper,
        d: dict,
) -> None:
    '''Write a Keywords to a file descritor

    Args:
        fd: File descriptor
        d: Keywords dictionary to write keyed by bytes
    '''
    keywords_str = 'KEYWORDS    '
    fd.write(''.join([keywords_str, d['keywords'], '\n']))


def writeSource(
        fd: _io.TextIOWrapper,
        d: dict,
) -> None:
    '''Write a Source to a file descritor

    Args:
        fd: File descriptor
        d: Source dictionary to write keyed by bytes
    '''
    source_str = 'SOURCE      '
    organism_str = '\n  ORGANISM  '
    indent_str = '\n            '
    org = d['organism']
    out_list = [
        source_str, d['source'],
        organism_str,
        org[0],
    ]
    if org[1] is not None:
        out_list += [
            indent_str,
            multiline_spaces(org[1], indent_str),
        ]
    out_list += ['\n']
    fd.write(''.join(out_list))


def writeReference(
        fd: _io.TextIOWrapper,
        ref: dict,
        i: int,
) -> None:
    '''Write a Reference to a file descritor

    Args:
        fd: File descriptor
        ref: Reference dictionary to write keyed by bytes
        i: index for the reference
    '''
    assert (ref['r_index'] == i)
    reference_str = 'REFERENCE   '
    authors_str = '  AUTHORS   '
    title_str = '  TITLE     '
    indent_str = '\n            '
    journal_str = '  JOURNAL   '
    pubmed_str = '  PUBMED    '
    if ref['start_idx'] is not None:
        idx_str = '  (bases %d to %d)' % (ref['start_idx'], ref['end_idx'])
    else:
        idx_str = ''
    out_list = [
        reference_str,
        '%d%s\n' % (i, idx_str),
        authors_str,
        multiline_spaces(ref['authors'], indent_str), '\n',
        title_str,
        multiline_spaces(ref['title'], indent_str), '\n',
        journal_str,
        multiline_spaces(ref['journal_info'], indent_str), '\n',
    ]
    if ref['pubmed'] is not None:
        out_list += [pubmed_str, ref['pubmed'], '\n']

    fd.write(''.join(out_list))


def writeComment(
        fd: _io.TextIOWrapper,
        d: dict,
) -> None:
    '''Write a Comment to a file descritor

    Args:
        fd: File descriptor
        d: Comment dictionary to write keyed by bytes
    '''
    comment_str = 'COMMENT     '
    indent_str = '\n            '
    # indent_str_gb_asm = '            '
    if 'comment' in d:
        comment = d['comment']
        if isinstance(comment, list):
            out_list = [
                comment_str,
                multiline_spaces(comment[0], indent_str),
                indent_str, indent_str, indent_str,
                indent_str.join(comment[1]), '\n',
            ]
        else:
            out_list = [
                comment_str, multiline_spaces(
                    comment, indent_str,
                ), '\n',
            ]
        fd.write(''.join(out_list))


def writeFeatures(
        fd: _io.TextIOWrapper,
        d: dict,
        order_qualifiers: bool,
) -> None:
    '''Write Features to a file descritor

    Args:
        fd: File descriptor
        d: Features dictionary to write keyed by bytes
        order_qualifiers: If True, order by the `'qualifiers'` key
    '''
    feature_header = 'FEATURES             Location/Qualifiers\n'
    feature_type_prefix_str = '     '
    indent_str = '\n                     '
    feature_type_field_size = 16
    fd.write(feature_header)
    for feature in d['features']:
        ftype = feature['type']
        location_str = feature['location']
        out_list = [
            feature_type_prefix_str,
            ftype, (feature_type_field_size - len(ftype)) * ' ',
            location_str, '\n',
        ]
        if order_qualifiers:
            quals = feature['qualifiers']
            for key in sorted(quals.keys()):
                value_list = quals[key]
                if not isinstance(value_list, list):
                    value_list = [value_list]
                # assumes value_list is all the same type
                for value in sorted(value_list):
                    # codon start is a weird edge case
                    if key in ('codon_start', 'transl_table'):
                        out_list += [
                            '                     /%s=%d\n' %
                            (key, value),
                        ]
                    elif key == 'transl_except':
                        out_list += [
                            '                     /%s=%s\n' %
                            (key, value),
                        ]
                    elif key == 'translation':
                        qualifier_str = "/%s=\"%s\"" % (key, value)
                        qual_list = multiline(
                            qualifier_str, indent_str, lim=58,
                        )  # or 58
                        out_list += ['                     ' + qual_list, '\n']
                    else:
                        if value is None:
                            qualifier_str = '/%s\n' % (key)
                        else:
                            qualifier_str = "/%s=\"%s\"\n" % (key, value)
                        qual_list = multiline_spaces(
                            qualifier_str, indent_str, lim=58,
                        )  # or 58
                        out_list += ['                     ' + qual_list, '\n']
        else:
            for key, value_list in feature['qualifiers'].items():
                if not isinstance(value_list, list):
                    value_list = [value_list]
                for value in value_list:
                    # codon start is a weird edge case
                    if key in ('codon_start', 'transl_table'):
                        out_list += [
                            '                     /%s=%d\n' %
                            (key, value),
                        ]
                    elif key == 'transl_except':
                        out_list += [
                            '                     /%s=%s\n' %
                            (key, value),
                        ]
                    elif key == 'translation':
                        qualifier_str = "/%s=\"%s\"" % (key, value)
                        qual_list = multiline(
                            qualifier_str, indent_str, lim=58,
                        )  # or 58
                        out_list += ['                     ' + qual_list, '\n']
                    else:
                        if value is None:
                            qualifier_str = '/%s\n' % (key)
                        else:
                            qualifier_str = "/%s=\"%s\"\n" % (key, value)
                        qual_list = multiline_spaces(
                            qualifier_str, indent_str, lim=58,
                        )  # or 58
                        out_list += ['                     ' + qual_list, '\n']
        fd.write(''.join(out_list))


def writeOrigin(
        fd: _io.TextIOWrapper,
        d: dict,
) -> None:
    '''Write Origin to a file descritor

    Args:
        fd: File descriptor
        d: Origin dictionary to write keyed by bytes
    '''
    i = 0
    format_string = '%9d %s %s %s %s %s %s\n'
    origin_header = 'ORIGIN      \n'
    fd.write(origin_header)
    origin = d['seq']
    count = len(origin) - 60
    while i < count:
        str_tup = (i + 1,) + tuple([
            origin[j:j + 10]
            for j in range(i, i + 60, 10)
        ])
        fd.write(format_string % str_tup)
        i += 60
    last_str = origin[i:]
    last_str = last_str + (60 - (len(last_str))) * ' '
    last_list = [last_str[j:j + 10] for j in range(0, 60, 10)]
    last_list = [x.rstrip() for x in last_list if x[0] != ' ']
    str_tup = (i + 1,) + tuple(last_list)
    format_string = '%9d' + ' %s' * len(last_list) + '\n'
    fd.write(format_string % str_tup)
    fd.write('//\n')


# if __name__ == '__main__':
#     import filecmp
#     import os.path as op

#     import gb_reader as gbr
#     path = op.dirname(op.dirname(op.dirname(op.abspath(__file__))))
#     fn = op.join(path, 'tests', 'test_data', 'mds42_full.gb')
#     fn_out = op.join(path, 'tests', 'test_data', 'sample_out.gb')
#     gbd = gbr.parse(fn, is_ordered=True)
#     write_file(fn_out, gbd)
#     print(filecmp.cmp(fn, fn_out))
