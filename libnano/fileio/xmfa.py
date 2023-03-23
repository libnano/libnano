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
libnano.fileio.xmfa
~~~~~~~~~~~~~~~~~~~

Lightweight XMFA parser. No dependencies. Files are parsed into simple
namedtuples as documented below.

'''
from __future__ import print_function

import io
import re
from typing import (
    List,
    NamedTuple,
)

_SUB_ALIGNMENT_GROUP = r'(?P<alignments>^>[^=]+)=(?P<notes>.*)$'
_SUB_ALIGNMENT_GROUP_RE = re.compile(
    _SUB_ALIGNMENT_GROUP,
    flags=re.M,
)

_SUB_ALIGNMENT = (
    r'^>\s*(?P<seq_num>\d+):(?P<start_idx>\d+)-'
    r'(?P<end_idx>\d+)\s*(?P<strand>[+|-])\s*'
    r'(?P<notes>[ \S]+)\s*?$(?P<seq>(?:[^=>]+)+)'
)

_SUB_ALIGNMENT_RE = re.compile(
    _SUB_ALIGNMENT,
    flags=re.M,
)

_WS_SPLIT = re.compile(r'\s+')


class SubAlignment(NamedTuple):
    '''SubAlignment class

    '''
    seq_num: int
    start_idx: int
    end_idx: int
    strand: str
    notes: str
    seq: str


class SubAlignmentGroup(NamedTuple):
    '''SubAlignmentGroup class

    '''
    alignments: List[SubAlignment]
    notes: str


def removeWhitespace(raw_string: str) -> str:
    '''Remove all whitespace (including newlines/carriage returns) from
    a string.

    Args:
        raw_string: Raw string to operate on

    Returns:
        Whitespace free version of string
    '''

    return ''.join(re.split(_WS_SPLIT, raw_string))


def parseXMFA(
        xmfa_fp: str,
) -> List[SubAlignmentGroup]:
    '''
    Args:
        xmfa_fp: Filepath to XMFA file

    Returns:
        List of :class:`SubAlignmentGroup` instance
    '''

    sub_alignment_group_comp = _SUB_ALIGNMENT_GROUP_RE
    sub_alignment_comp = _SUB_ALIGNMENT_RE

    with io.open(xmfa_fp, 'r', encoding='utf-8') as xmfa_fd:
        xmfa_data = xmfa_fd.read()

    sub_alignment_groups = re.finditer(sub_alignment_group_comp, xmfa_data)
    sub_alignment_group_list = []

    for sub_alignment_group in sub_alignment_groups:

        sub_alignments = re.finditer(
            sub_alignment_comp,
            sub_alignment_group.group('alignments'),
        )
        sub_alignment_list = []

        for sub_alignment in sub_alignments:

            sub_alignment_list.append(
                SubAlignment(
                    seq_num=int(sub_alignment.group('seq_num')),
                    start_idx=int(sub_alignment.group('start_idx')),
                    end_idx=int(sub_alignment.group('end_idx')),
                    strand=str(sub_alignment.group('strand')),
                    notes=str(sub_alignment.group('notes')),
                    seq=removeWhitespace(sub_alignment.group('seq')),
                ),
            )

        sub_alignment_group_list.append(
            SubAlignmentGroup(
                alignments=sub_alignment_list,
                notes=str(sub_alignment_group.group('notes')),
            ),
        )

    return sub_alignment_group_list


if __name__ == '__main__':
    import os

    DIR = os.path.dirname(os.path.realpath(__file__))
    print()
    # output = parseXMFA(os.path.join(DIR, 'example.xmfa'))
    output = parseXMFA(os.path.join(DIR, '1swap100_output', '1swap100.xmfa'))
    print()
    print(output)
