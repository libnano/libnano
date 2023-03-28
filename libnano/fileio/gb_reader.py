# Copyright (C) 2023. Nick Conway;
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
libnano.fileio.gb_reader
~~~~~~~~~~~~~~~~~~~~~~~~~~

http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html
http://www.insdc.org/documents/feature_table.html

All keys are native strings, as are values, except the origin which is always
a python 3 byte string (not unicode)

'''
import io
import re
from collections import OrderedDict
from typing import (
    Any,
    Dict,
    List,
)


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
    d: Dict[str, Any] = {'info': {}}
    d_info = d['info']
    with io.open(filepath, 'r', encoding='utf-8') as fd:
        raw = fd.read()
    start, _, origin = raw.partition('ORIGIN')
    start, _, features = start.partition(
        'FEATURES             Location/Qualifiers\n',
    )

    d_info.update(parse_locus(start))
    d_info.update(parse_definition(start))
    d_info.update(parse_accession(start))
    d_info.update(parse_version(start))
    d_info.update(parse_db_link(start))
    d_info.update(parse_keywords(start))
    d_info.update(parse_source(start))
    d_info.update(parse_organism(start))

    d_info['references'] = parse_reference(start)
    _, _, comment = start.partition('COMMENT     ')
    parse_comment(d_info, comment)
    d['features'] = parse_features(features, is_ordered)
    d['seq'] = parse_origin(origin)
    return d


def parse_comment(
        d: Dict,
        comment: str,
) -> None:
    '''Parse comment into dictionary

    Args:
        d: Output dictionary
        comment: Raw byte comment string to parse

    '''
    if comment != '':
        # get rid of ApE empty comment
        if comment.startswith('\nCOMMENT     '):
            comment = comment[13:]
        idx_genome_asm_data = -1
        genome_asm_data_newline_count = 0
        lines = comment.split('\n')
        lines = [line.strip() for line in lines]
        # print(lines)
        # handle ##Genome-Assembly-Data-START## edge case
        for i, line in enumerate(lines):
            if line == '':
                genome_asm_data_newline_count += 1
            elif genome_asm_data_newline_count == 2:
                idx_genome_asm_data = i
                genome_asm_data_newline_count = 0
            else:
                genome_asm_data_newline_count = 0

        if idx_genome_asm_data < 0:
            d['comment'] = ' '.join(lines)
        else:
            d['comment'] = [
                ' '.join(lines[:idx_genome_asm_data - 2]),
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
    r'.*\n',                                    # match line end
]
RE_LOCUS: str = ''.join(_RE_LOCUS)

_RE_DEFINITION: List[str] = [
    '^DEFINITION',                          # field
    # look ahead assertion for multiline
    ' +(?P<definition>(?:.*\n)(?: .*\n)*)',
]
RE_DEFINITION: str = ''.join(_RE_DEFINITION)

_RE_ACCESSION: List[str] = [
    r'^ACCESSION',                          # field
    # look ahead assertion for multiline
    r' +(?P<accession>[\w|.]*)'
    r'.*\n',                                # match line end
]
RE_ACCESSION: str = ''.join(_RE_ACCESSION)

_RE_VERSION: List[str] = [
    r'^VERSION',                # field
    r' +(?P<version>[\w|.]+)',  # version
    r' +GI:(?P<GI>[\w|.]+)'     # gi field
    r'.*\n',                    # match line end
]

RE_DBLINK: str = r'^DBLINK +(?P<dblink>[\w|:| |.]+)\n'

RE_VERSION: str = ''.join(_RE_VERSION)

_RE_KEYWORDS: List[str] = [
    r'^KEYWORDS',
    r' +(?P<keywords>[\w|.]*)'
    r'.*\n',
]
RE_KEYWORDS = ''.join(_RE_KEYWORDS)

_RE_SOURCE: List[str] = [
    r'^SOURCE',
    r' +(?P<source>.*)',
    r'\n',
]
RE_SOURCE: str = ''.join(_RE_SOURCE)

_RE_ORGANISM: List[str] = [
    r'^  ORGANISM',                          # field
    r'(?: +(?P<organism0>(?:.*\n))?',
    r'(?: +(?P<organism1>(?:.*\n)(?: .*\n)*))?)',  # multiline
]
RE_ORGANISM: str = ''.join(_RE_ORGANISM)

RE_PREAMBLE: str = (
    RE_LOCUS +
    RE_DEFINITION +
    RE_ACCESSION +
    RE_VERSION +
    RE_DBLINK +
    RE_KEYWORDS +
    RE_SOURCE +
    RE_ORGANISM
)
RE_COMP_PREABLE = re.compile(
    RE_PREAMBLE,
    flags=re.M,
)

RE_COMP_LOCUS = re.compile(
    RE_LOCUS,
    flags=re.M,
)


def parse_locus(
        raw: str,
) -> Dict[str, Any]:
    '''Parse locus into dictionary

    Args:
        raw: Raw string to parse
    '''
    m = re.match(RE_COMP_LOCUS, raw)
    if m is not None:
        d = m.groupdict()
        d['length'] = int(d['length'])
        return d
    else:
        return {}


RE_COMP_DEFINITION = re.compile(RE_DEFINITION, flags=re.M)


def parse_definition(
        raw: str,
) -> Dict:
    '''Parse definition into dictionary

    Args:
        raw: Raw string to parse
    '''
    m = re.search(RE_COMP_DEFINITION, raw)
    if m is None:
        return {'definition': None}
    d = m.groupdict()
    if d['definition'] is not None:
        temp_l = d['definition'].split('\n')
        temp_l = [x.strip() for x in temp_l]
        d['definition'] = ' '.join(temp_l)[:-1]
    return d


RE_COMP_ACCESSION = re.compile(RE_ACCESSION, flags=re.M)


def parse_accession(
        raw: str,
) -> Dict:
    '''Parse accession into dictionary

    Args:
        raw: Raw string to parse
    '''
    m = re.search(RE_COMP_ACCESSION, raw)
    if m is None:
        return {'accession': None}
    d = m.groupdict()
    return d


RE_COMP_VERSION = re.compile(RE_VERSION, flags=re.M)


def parse_version(
        raw: str,
) -> Dict:
    '''Parse version into dictionary

    Args:
        raw: Raw string to parse
    '''
    m = re.search(RE_COMP_VERSION, raw)
    if m is None:
        # print("Version none")
        return {'version': None}
    d = m.groupdict()
    # print("version", d)
    return d


RE_COMP_DBLINK = re.compile(
    RE_DBLINK,
    flags=re.M,
)


def parse_db_link(
        raw: str,
) -> Dict:
    '''Parse DB Link into dictionary

    Args:
        raw: Raw string to parse
    '''
    m = re.search(RE_COMP_DBLINK, raw)
    if m is None:
        return {'dblink': None}
    d = m.groupdict()
    return d


RE_COMP_KEYWORDS = re.compile(
    RE_KEYWORDS,
    flags=re.M,
)


def parse_keywords(
        raw: str,
) -> Dict:
    '''Parse keywords into dictionary

    Args:
        raw: Raw string to parse
    '''
    m = re.search(RE_COMP_KEYWORDS, raw)
    if m is None:
        return {'keywords': None}
    d = m.groupdict()
    return d


RE_COMP_SOURCE = re.compile(
    RE_SOURCE,
    flags=re.M,
)


def parse_source(
        raw: str,
) -> Dict:
    '''Parse source into dictionary

    Args:
        raw: Raw string to parse
    '''
    m = re.search(RE_COMP_SOURCE, raw)
    if m is None:
        return {'source': None}
    d = m.groupdict()
    return d


RE_COMP_ORGANISM = re.compile(
    RE_ORGANISM,
    flags=re.M,
)


def parse_organism(
        raw: str,
) -> Dict:
    '''Parse organism into dictionary

    Args:
        raw: Raw string to parse
    '''
    m = re.search(RE_COMP_ORGANISM, raw)
    if m is None:
        return {'organism': [None, None]}
    d = m.groupdict()

    temp_l = d['organism0'].split('\n')
    temp_l = [x.strip() for x in temp_l]
    org0 = ' '.join(temp_l)[:-1]

    org1 = None
    if d['organism1'] is not None:
        temp_l = d['organism1'].split('\n')
        temp_l = [x.strip() for x in temp_l]
        org1 = ' '.join(temp_l)[:-1]

    del d['organism0']
    del d['organism1']

    d['organism'] = [org0, org1]
    return d


'''
REFERENCE   1  (bases 1 to 5028)
  AUTHORS   Torpey,L.E., Gibbs,P.E., Nelson,J. and Lawrence,C.W.
  TITLE     Cloning and sequence of REV7, a gene whose function is required for
            DNA damage-induced mutagenesis in Saccharomyces cerevisiae
  JOURNAL   Yeast 10 (11), 1503-1509 (1994)
  PUBMED    7871890
'''
_RE_REFERENCE: List[str] = [
    r'^REFERENCE',
    (
        r' +(?P<r_index>[0-9]+)(?: +\(bases '
        r'(?P<start_idx>[0-9]+) to (?P<end_idx>[0-9]+)\)){0,1}'
    ),
    r'.*\n',
    r'^  AUTHORS',
    r' +(?P<authors>.+)\n',
    r'^  TITLE',                            # field
    r' +(?P<title>(?:.*\n)(?: .*\n)*)',     # multiline
    r'^  JOURNAL',
    r' +(?P<journal_info>.+\n(?: {12}.+\n)*)',
    r'(?:^  PUBMED +(?P<pubmed>[0-9]+)\n){0,1}',
]
RE_REFERENCE: str = ''.join(_RE_REFERENCE)
RE_COMP_REF = re.compile(
    RE_REFERENCE,
    flags=re.M,
)


def parse_reference(
        raw: str,
) -> List[dict]:
    '''Parse reference into a list of dictionaries

    Args:
        raw: Raw string to parse

    Returns:
        List of dictionaries

    '''
    ref_list = []

    for m in re.finditer(RE_COMP_REF, raw):
        d = m.groupdict()
        temp_l = d['title'].split('\n')
        temp_l = [x.strip() for x in temp_l]
        d['title'] = ' '.join(temp_l)[:-1]

        temp_l = d['journal_info'].split('\n')
        temp_l = [x.strip() for x in temp_l]
        d['journal_info'] = ' '.join(temp_l)[:-1]

        d['r_index'] = int(d['r_index'])
        if d['start_idx'] is not None:
            d['start_idx'] = int(d['start_idx'])
        if d['end_idx'] is not None:
            d['end_idx'] = int(d['end_idx'])
        ref_list.append(d)
    # end for
    return ref_list


def add_multivalue(
        d: Dict,
        key: str,
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
    r' +(?P<location>.+)\n',
    r'(?P<qualifiers>(?:^ {21}.*\n)*)',
]
RE_FEATURE: str = ''.join(_RE_FEATURE)
RE_COMP_FEATURE = re.compile(RE_FEATURE, flags=re.M)


# Qualifers can have tags with /'s in the value so it's tough to escape them
# for now we need to split on "                     /"

def parse_features(
        raw: str,
        is_ordered: bool = False,
) -> List[dict]:
    '''
    Parse features into a list of dictionaries

    Args:
        raw: Raw string to parse
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

        d: Dict[str, Any] = {
            'type': feature['feature_key'],
            'location': feature['location'],
            # 'partials': (feature['partial5'], feature['partial3']),
            'qualifiers': OrderedDict() if is_ordered else {},
        }
        qs = d['qualifiers']
        # prevent splitting on </tags>
        qs_list = feature['qualifiers'].split('                     /')
        for qualifier in qs_list[1:]:
            # heal line breaks
            '''
            Need to address the multi value key problem
            i.e. more than one value in the list of qualifiers
            '''
            q_list = qualifier.split('=')
            key = q_list[0]
            yes_val = True
            try:
                q_list = q_list[1].split('\n')
                if q_list[-1] == '':
                    q_list.pop()  # remove ending '' item
            except (IndexError, SyntaxError):
                q_list = ['']
                yes_val = False
            q_list = [x.lstrip() for x in q_list]
            is_str = True
            if key == 'translation':
                temp = ''.join(q_list)
            elif key in ('codon_start', 'transl_table'):
                is_str = False
                temp = ' '.join(q_list)
            else:
                temp = ' '.join(q_list)
            if yes_val and temp[0] == '\"' and temp[-1] == '\"':
                value_to_add = temp[1:-1] if is_str else int(temp[1:-1])
            elif not yes_val:
                value_to_add = None
            else:
                value_to_add = temp if is_str else int(temp)
            add_multivalue(qs, key, value_to_add)
        features_list.append(d)

    return features_list


def parse_origin(raw: str) -> str:
    '''
    Parse origin into a string

    Args:
        raw: Raw string to parse

    Returns:
        Origin string
    '''
    out_list = []

    all_lines = raw.split('\n')
    start = 1 if all_lines[0].strip() == '' else 0
    for line in all_lines[start:-1]:
        temp = line.split()
        out_list += temp[1:]
    seq = ''.join(out_list)

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
