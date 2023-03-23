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
libnano.fileio.gb_reader_b
~~~~~~~~~~~~~~~~~~~~~~~~~~

http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html
http://www.insdc.org/documents/feature_table.html

All keys are native strings, as are values, except the origin which is always
a python 3 byte string (not unicode)

'''
import io
import re
import sys
from collections import OrderedDict
from typing import (
    Any,
    Dict,
    List,
    Union,
)


def _bytes(x: Union[str, bytes]) -> bytes:
    return x.encode('utf8') if isinstance(x, str) else x


NEWLINE_STR: str = '\r\n' if sys.platform == 'win32' else '\n'
NEWLINE_BYT: bytes = b'\r\n' if sys.platform == 'win32' else b'\n'


def parse(
        filepath: str,
        is_ordered: bool = False,
) -> Dict:
    '''Parse a Genbank file to a dictionary

    Args:
        filepath: Filepath to read Genbank file from

        is_ordered: If True will retain the order of the qualifiers

    Returns:
        Parsed Genbank file dictionary

    '''
    d: Dict[bytes, Any] = {b'info': {}}
    d_info = d[b'info']
    with io.open(filepath, 'rb') as fd:
        raw = fd.read()
    start, _, origin = raw.partition(b'ORIGIN')
    start, _, features = start.partition(
        b'FEATURES             Location/Qualifiers%s' % (NEWLINE_BYT),
    )

    parseLocus(start, d_info)
    parseDefinition(start, d_info)
    parseAccession(start, d_info)
    parseVersion(start, d_info)
    parseDBLink(start, d_info)
    parseKeywords(start, d_info)
    parseSource(start, d_info)
    parseOrganism(start, d_info)

    d_info[b'references'] = parseReference(start)
    _, _, comment = start.partition(b'COMMENT     ')
    parseComment(d_info, comment)
    d[b'features'] = parseFeatures(features, is_ordered)  # type: ignore
    d[b'seq'] = parseOrigin(origin)
    return d


def parseComment(
        d: Dict,
        comment: bytes,
) -> None:
    '''Parse comment into dictionary

    Args:
        d: Output dictionary
        comment: Raw byte comment string to parse

    '''
    if comment != b'':
        # get rid of ApE empty comment
        if comment.startswith(b'%sCOMMENT     ' % (NEWLINE_BYT)):
            comment = comment[13:]
        idx_genome_asm_data = -1
        genome_asm_data_newline_count = 0
        lines = comment.split(NEWLINE_BYT)
        lines = [line.strip() for line in lines]
        # print(lines)
        # handle ##Genome-Assembly-Data-START## edge case
        for i, line in enumerate(lines):
            if line == b'':
                genome_asm_data_newline_count += 1
            elif genome_asm_data_newline_count == 2:
                idx_genome_asm_data = i
                genome_asm_data_newline_count = 0
            else:
                genome_asm_data_newline_count = 0
        # end for
        if idx_genome_asm_data < 0:
            d[b'comment'] = b' '.join(lines)
        else:
            d[b'comment'] = [
                b' '.join(lines[:idx_genome_asm_data - 2]),
                lines[idx_genome_asm_data:-1],
            ]


_RE_LOCUS: List[str] = [
    r'^LOCUS',                                  # field
    r' +(?P<name>[\w|.]+)',                     # name
    r' +(?P<length>[0-9]+) bp',                 # sequence length
    r'(?: +(?P<stranded>[a-z]{2})-)?',          # opt: ss, ds, ms
    r' *(?P<molecule_type>[a-z|A-Z|-|]{2,6})',  # molecule type
    r' +(?P<form>[\w]{6,8})?',                  # linear or circular
    r' +(?P<gb_division>[a-z|A-Z]{3})?',        # Genbank division
    r' +(?P<mod_date>[0-9]+-[A-Z]+-[0-9]+)',    # modification date
    r'.*%s' % (NEWLINE_STR),                    # match line end
]
RE_LOCUS: bytes = _bytes(''.join(_RE_LOCUS))

_RE_DEFINITION: List[str] = [
    r'^DEFINITION',                         # field
    ' +(?P<definition>(?:.*%s)(?: .*%s)*)' % (
        NEWLINE_STR,
        NEWLINE_STR,
    ),  # look ahead assertion for multiline
]
RE_DEFINITION: bytes = _bytes(''.join(_RE_DEFINITION))

_RE_ACCESSION: List[str] = [
    r'^ACCESSION',                  # field
    r' +(?P<accession>[\w|.]*)'     # look ahead assertion for multiline
    r'.*',                          # match line end
    NEWLINE_STR,
]
RE_ACCESSION: bytes = _bytes(''.join(_RE_ACCESSION))

_RE_VERSION: List[str] = [
    r'^VERSION',                # field
    r' +(?P<version>[\w|.]+)',  # version
    r' +GI:(?P<GI>[\w|.]+)'     # gi field
    r'.*',                      # match line end
    NEWLINE_STR,
]
RE_VERSION: bytes = _bytes(''.join(_RE_VERSION))

RE_DBLINK = (
    r'^DBLINK +(?P<dblink>[\w|:| |.]+)'.encode('utf8') + NEWLINE_BYT
)

_RE_KEYWORDS: List[str] = [
    r'^KEYWORDS',
    r' +(?P<keywords>[\w|.]*)'
    r'.*',
    NEWLINE_STR,
]
RE_KEYWORDS: bytes = _bytes(''.join(_RE_KEYWORDS))

_RE_SOURCE: List[str] = [
    r'^SOURCE',
    r' +(?P<source>.*)',
    NEWLINE_STR,
]
RE_SOURCE: bytes = _bytes(''.join(_RE_SOURCE))

_RE_ORGANISM: List[str] = [
    r'^  ORGANISM',                          # field
    r'(?: +(?P<organism0>(?:.*%s))?' % NEWLINE_STR,
    r'(?: +(?P<organism1>(?:.*%s)(?: .*%s)*))?)' % (
        NEWLINE_STR,
        NEWLINE_STR,
    ),  # multiline
]
RE_ORGANISM: bytes = _bytes(''.join(_RE_ORGANISM))


RE_COMP_LOCUS = re.compile(
    RE_LOCUS,
    flags=re.M,
)


def parseLocus(
        raw: bytes,
        d_out: Dict,
) -> None:
    '''Parse locus into dictionary

    Args:
        raw: Raw byte string to parse
        d: Output dictionary

    '''
    m = re.match(RE_COMP_LOCUS, raw)
    if m is not None:
        d = m.groupdict()
        d['length'] = int(d['length'])
        for k, v in d.items():
            d_out[_bytes(k)] = v


RE_COMP_DEFINITION = re.compile(
    RE_DEFINITION,
    flags=re.M,
)


def parseDefinition(
        raw: bytes,
        d_out: Dict,
) -> None:
    '''Parse definition into dictionary

    Args:
        raw: Raw byte string to parse
        d: Output dictionary

    '''
    m = re.search(RE_COMP_DEFINITION, raw)
    if m is None:
        d_out[b'definition'] = None
    else:
        d = m.groupdict()
        if d['definition'] is not None:
            temp_l = d['definition'].split(NEWLINE_BYT)
            temp_l = [x.strip() for x in temp_l]
            d_out[b'definition'] = b' '.join(temp_l)[:-1]
        else:
            d_out[b'definition'] = None


RE_COMP_ACCESSION = re.compile(RE_ACCESSION, flags=re.M)


def parseAccession(
        raw: bytes,
        d_out: Dict,
) -> None:
    '''Parse accession into dictionary

    Args:
        raw: Raw byte string to parse
        d_out: Output dictionary

    '''
    m = re.search(RE_COMP_ACCESSION, raw)
    if m is None:
        d_out[b'accession'] = None
    else:
        d = m.groupdict()
        d_out[b'accession'] = d['accession']


RE_COMP_VERSION = re.compile(RE_VERSION, flags=re.M)


def parseVersion(
        raw: bytes,
        d_out: Dict,
) -> None:
    '''Parse version into dictionary

    Args:
        raw: Raw byte string to parse
        d_out: Output dictionary

    '''
    m = re.search(RE_COMP_VERSION, raw)
    if m is None:
        d_out[b'version'] = None
    else:
        d = m.groupdict()
        d_out[b'version'] = d['version']
        d_out[b'GI'] = d['GI']


RE_COMP_DBLINK = re.compile(RE_DBLINK, flags=re.M)


def parseDBLink(
        raw: bytes,
        d_out: Dict,
) -> None:
    '''Parse DB link into dictionary

    Args:
        raw: Raw byte string to parse
        d_out: Output dictionary

    '''
    m = re.search(RE_COMP_DBLINK, raw)
    if m is None:
        d_out[b'dblink'] = None
    else:
        d = m.groupdict()
        d_out[b'dblink'] = d['dblink']


RE_COMP_KEYWORDS = re.compile(RE_KEYWORDS, flags=re.M)


def parseKeywords(
        raw: bytes,
        d_out: Dict,
) -> None:
    '''Parse keywords into dictionary

    Args:
        raw: Raw byte string to parse
        d_out: Output dictionary

    '''
    m = re.search(RE_COMP_KEYWORDS, raw)
    if m is None:
        d_out[b'keywords'] = None
    else:
        d = m.groupdict()
        d_out[b'keywords'] = d['keywords']


RE_COMP_SOURCE = re.compile(RE_SOURCE, flags=re.M)


def parseSource(
        raw: bytes,
        d_out: Dict,
) -> None:
    '''Parse source into dictionary

    Args:
        raw: Raw byte string to parse
        d_out: Output dictionary

    '''
    m = re.search(RE_COMP_SOURCE, raw)
    if m is None:
        d_out[b'source'] = None
    else:
        d = m.groupdict()
        d_out[b'source'] = d['source']


RE_COMP_ORGANISM = re.compile(RE_ORGANISM, flags=re.M)


def parseOrganism(
        raw: bytes,
        d_out: Dict,
) -> None:
    '''Parse organism into dictionary

    Args:
        raw: Raw byte string to parse
        d_out: Output dictionary

    '''
    m = re.search(RE_COMP_ORGANISM, raw)
    if m is None:
        d_out[b'organism'] = [None, None]
    else:
        d = m.groupdict()

        temp_l = d['organism0'].split(NEWLINE_BYT)
        temp_l = [x.strip() for x in temp_l]
        org0 = b' '.join(temp_l)[:-1]

        org1 = None
        if d['organism1'] is not None:
            temp_l = d['organism1'].split(NEWLINE_BYT)
            temp_l = [x.strip() for x in temp_l]
            org1 = b' '.join(temp_l)[:-1]

        d_out[b'organism'] = [org0, org1]


'''
REFERENCE   1  (bases 1 to 5028)
  AUTHORS   Torpey,L.E., Gibbs,P.E., Nelson,J. and Lawrence,C.W.
  TITLE     Cloning and sequence of REV7, a gene whose function is required for
            DNA damage-induced mutagenesis in Saccharomyces cerevisiae
  JOURNAL   Yeast 10 (11), 1503-1509 (1994)
  PUBMED    7871890
'''
re_reference: List[str] = [
    r'^REFERENCE',
    (
        r' +(?P<r_index>[0-9]+)(?: +\(bases '
        r'(?P<start_idx>[0-9]+) to (?P<end_idx>[0-9]+)\)){0,1}'
    ),
    r'.*',
    NEWLINE_STR,
    r'^  AUTHORS',
    r' +(?P<authors>.+)',
    NEWLINE_STR,
    r'^  TITLE',                                                        # field
    r' +(?P<title>(?:.*%s)(?: .*%s)*)' % (NEWLINE_STR, NEWLINE_STR),    # mline
    r'^  JOURNAL',
    r' +(?P<journal_info>.+%s(?: {12}.+%s)*)' % (NEWLINE_STR, NEWLINE_STR),
    r'(?:^  PUBMED +(?P<pubmed>[0-9]+)%s){0,1}' % (NEWLINE_STR),
]
RE_REFERENCE = _bytes(''.join(re_reference))
RE_COMP_REF = re.compile(RE_REFERENCE, flags=re.M)


def parseReference(
        raw: bytes,
) -> List[dict]:
    '''Parse reference into a list of dictionaries

    Args:
        raw: Raw byte string to parse

    Returns:
        List of dictionaries for the reference

    '''
    ref_list = []

    for m in re.finditer(RE_COMP_REF, raw):
        d_temp: Dict[bytes, Any] = {}
        d = m.groupdict()
        temp_l = d['title'].split(NEWLINE_BYT)
        temp_l = [x.strip() for x in temp_l]
        d_temp[b'title'] = b' '.join(temp_l)[:-1]

        temp_l = d['journal_info'].split(NEWLINE_BYT)
        temp_l = [x.strip() for x in temp_l]
        d_temp[b'journal_info'] = b' '.join(temp_l)[:-1]

        d_temp[b'r_index'] = int(d['r_index'])
        if d['start_idx'] is not None:
            d_temp[b'start_idx'] = int(d['start_idx'])
        else:
            d_temp[b'start_idx'] = None
        if d['end_idx'] is not None:
            d_temp[b'end_idx'] = int(d['end_idx'])
        else:
            d_temp[b'end_idx'] = None
        d_temp[b'authors'] = d['authors']
        d_temp[b'pubmed'] = d['pubmed']
        ref_list.append(d_temp)
    # end for
    return ref_list


def addMultivalue(
        d: Dict,
        key: bytes,
        val: Any,
) -> None:
    '''Add a multivalue item to the dictionary

    Args:
        d: Output dictionary
        key: key to add a multivalue for
        value: New value to add
    '''
    if key in d:
        old_val = d[key]
        if isinstance(old_val, list):
            old_val.append(val)
        else:
            d[key] = [old_val, val]
    else:
        d[key] = val


'''
see section 3.4 Location
'''
_RE_FEATURE: List[str] = [
    r'^ {5}(?P<feature_key>[\w]+)',
    r' +(?P<location>.+)',
    NEWLINE_STR,
    '(?P<qualifiers>(?:^ {21}.*%s)*)' % (NEWLINE_STR),
]
RE_FEATURE: bytes = _bytes(''.join(_RE_FEATURE))
RE_COMP_FEATURE = re.compile(
    RE_FEATURE,
    flags=re.M,
)


# Qualifers can have tags with /'s in the value so it's tough to escape them
# for now we need to split on "                     /"

QUOTE_BYTE: int = b'\"'[0]


def parseFeatures(
        raw: bytes,
        is_ordered: bool = False,
) -> List[dict]:
    '''
    Parse features into a list of dictionaries

    Args:
        raw: Raw byte string to parse
        is_ordered: If True, parse into an ``OrderedDictionary``

    Returns:
        list of dictionaries for the features
    '''
    features_list = []

    for feature_match in re.finditer(RE_COMP_FEATURE, raw):
        feature = feature_match.groupdict()
        if 'qualifiers' not in feature:
            print(feature)
            raise IOError('bad feature')

        d: Dict[bytes, Any] = {
            b'type': feature['feature_key'],
            b'location': feature['location'],
            # 'partials': (feature['partial5'], feature['partial3']),
            b'qualifiers': OrderedDict() if is_ordered else {},
        }
        qs = d[b'qualifiers']
        # prevent splitting on </tags>
        qs_list = feature['qualifiers'].split(b'                     /')
        for qualifier in qs_list[1:]:
            # heal line breaks
            '''
            Need to address the multi value key problem
            i.e. more than one value in the list of qualifiers
            '''
            q_list = qualifier.split(b'=')
            key = q_list[0]
            yes_val = True
            try:
                q_list = q_list[1].split(NEWLINE_BYT)
                if q_list[-1] == b'':
                    q_list.pop()  # remove ending '' item
            except (IndexError, SyntaxError):
                q_list = [b'']
                yes_val = False
            q_list = [x.lstrip() for x in q_list]
            is_str = True
            if key == b'translation':
                temp = b''.join(q_list)
            elif key in (b'codon_start', b'transl_table'):
                is_str = False
                temp = b' '.join(q_list)
            else:
                temp = b' '.join(q_list)
            if yes_val and temp[0] == QUOTE_BYTE and temp[-1] == QUOTE_BYTE:
                value_to_add = temp[1:-1] if is_str else int(temp[1:-1])
            elif not yes_val:
                value_to_add = None
            else:
                value_to_add = temp if is_str else int(temp)
            addMultivalue(
                qs,
                key,
                value_to_add,
            )
        features_list.append(d)
    # end for
    return features_list


def parseOrigin(raw: bytes) -> bytes:
    '''
    Parse origin into a string

    Args:
        raw: Raw byte string to parse

    Returns:
        Origin string
    '''
    out_list = []

    all_lines = raw.split(NEWLINE_BYT)
    start = 1 if all_lines[0].strip() == b'' else 0
    for line in all_lines[start:-1]:
        temp = line.split()
        out_list += temp[1:]
    seq = b''.join(out_list)

    assert (seq.isalpha())
    return seq


if __name__ == '__main__':
    import os.path as op
    path = op.dirname(op.dirname(op.dirname(op.abspath(__file__))))
    fn = op.join(path, 'tests', 'test_data', 'failed.gb')

    def main():
        d = parse(fn)
        return d

    print(main())
