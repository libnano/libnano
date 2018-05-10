# -*- coding: utf-8 -*-
'''
xmfa | xmfa.py
~~~~~~~~~~~~~~

Lightweight XMFA parser. No dependencies. Files are parsed into simple
namedtuples as documented below.

'''
from __future__ import print_function
import re
from collections import namedtuple
import io
from typing import List

_sub_alignment_group: str = r'(?P<alignments>^>[^=]+)=(?P<notes>.*)$'
_sub_alignment_group_re: '_sre.SRE_Pattern' = re.compile(_sub_alignment_group, flags=re.M)

_sub_alignment: str= (  r'^>\s*(?P<seq_num>\d+):(?P<start_idx>\d+)-'
                        r'(?P<end_idx>\d+)\s*(?P<strand>[+|-])\s*'
                        r'(?P<notes>[ \S]+)\s*?$(?P<seq>(?:[^=>]+)+)')

_sub_alignment_re: '_sre.SRE_Pattern' = re.compile(_sub_alignment, flags=re.M)

_ws_split: '_sre.SRE_Pattern' = re.compile('\s+')

SubAlignmentGroup: namedtuple = namedtuple( 'SubAlignmentGroup',
                                            ['alignments', 'notes'])

SubAlignment: namedtuple = namedtuple('SubAlignment',
                                        [   'seq_num',
                                            'start_idx',
                                            'end_idx',
                                            'strand',
                                            'notes',
                                            'seq'])


def removeWhitespace(raw_string: str) -> str:
    ''' Remove all whitespace (including newlines/carriage returns) from
    a string.
    '''

    return ''.join(re.split(_ws_split, raw_string))

def parseXMFA(xmfa_fp: str) -> List[SubAlignmentGroup]:
    """
    """

    sub_alignment_group_comp = _sub_alignment_group_re
    sub_alignment_comp = _sub_alignment_re

    with io.open(xmfa_fp, 'r', encoding='utf-8') as xmfa_fd:
        xmfa_data = xmfa_fd.read()

    sub_alignment_groups = re.finditer(sub_alignment_group_comp, xmfa_data)
    sub_alignment_group_list = []

    for sub_alignment_group in sub_alignment_groups:

        sub_alignments = re.finditer(sub_alignment_comp,
                                     sub_alignment_group.group('alignments'))
        sub_alignment_list = []

        for sub_alignment in sub_alignments:

            sub_alignment_list.append(
                SubAlignment(
                    seq_num=        int(sub_alignment.group('seq_num')),
                    start_idx=      int(sub_alignment.group('start_idx')),
                    end_idx=        int(sub_alignment.group('end_idx')),
                    strand=         str(sub_alignment.group('strand')),
                    notes=          str(sub_alignment.group('notes')),
                    seq=            removeWhitespace(sub_alignment.group('seq'))
            ))

        sub_alignment_group_list.append(
            SubAlignmentGroup(
                alignments=         sub_alignment_list,
                notes=              str(sub_alignment_group.group('notes'))
        ))

    return sub_alignment_group_list


if __name__ == '__main__':
    import os

    DIR = os.path.dirname(os.path.realpath(__file__))
    print()
    # output = parseXMFA(os.path.join(DIR, 'example.xmfa'))
    output = parseXMFA(os.path.join(DIR, '1swap100_output', '1swap100.xmfa'))
    print()
    print(output)

