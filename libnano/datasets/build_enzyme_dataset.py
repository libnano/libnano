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
libnano.datasets.build_enzyme_dataset
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Includes methods for parsing and visualizing restriction enzyme cut sites
from the NEB Rebase list found at:

http://rebase.neb.com/rebase/link_allenz

Restriction enzyme information will be serialized and stored in a JSON file
with entries in the following format:

{
  'name': <name>
  'cutsites': [{  # List of cutsites (some enzymes have more than one)
      'core_seq': [fwd, revcomp],  # No ambiguous bases outside of core
      'full_seq': [fwd, revcomp],  # Includes ambiguous bases
      'core_regex': [fwd, revcomp],  # Reg. exps. to match `core_seq`
      'full_regex': [fwd, revcomp],  # Reg. exps. to match `full_seq`
      'cutsite_idxs' [[fwd], [rev]],  # Indexed relative to base 0 of exp. seq.
      'shorthand': <NEB format>  # NEB-format shorthand for cutsite
  }],
  'availability': <NEB format> # NEB-format shorthand for avail.
}

Here are some notes on NEB Rebase cutsite/shorthand formatting:

    If the recognition sequence is symmetric, the cutsite will be represented
    with a caret:

        e.g., EcoRI -- NEB shorthand: G^AATTC
              Cutting:
                GAATTC    ->    AATTC
                CTTAAG    ->   CTTAA

    If the recognition sequence is asymmetric, the cutsite will be represented
    with parenthetical notation (indexing relative to end):

        e.g., BseYI -- NEB shorthand: CCCAGC(-5/-1)
              Cutting:
                CCCAGC    ->     CCAGC
                GGGTCG    ->    GGGTC

        e.g., SapI -- NEB shorthand: GCTCTTC(1/4)
              Cutting:
                GCTCTTCNNNN    ->    GCTCTTCN
                GCTCTTCNNNN    ->    GCTCTTCNNNN

        e.g., Bsp24I -- NEB shorthand: (8/13)GACNNNNNNTGG(12/7)
              Cutting:
                  NNNNNNNNNNNNNGACNNNNNNTGGNNNNNNNNNNNN    ->
                  NNNNNNNNNNNNNCTGNNNNNNACCNNNNNNNNNNNN    ->

                       NNNNNNNNGACNNNNNNTGGNNNNNNNNNNNN
                  NNNNNNNNNNNNNCTGNNNNNNACCNNNNNNN

'''
from __future__ import print_function

import json
import os.path as op
import re
import sys
import urllib.request as urllib2
from typing import (  # Union,
    Any,
    Dict,
    List,
    NamedTuple,
    Optional,
    Tuple,
)

try:
    import libnano
except (ImportError, ModuleNotFoundError):
    LIBNANO_PATH = op.dirname(
        op.dirname(
            op.dirname(op.abspath(__file__)),
        ),
    )
    sys.path = [LIBNANO_PATH] + sys.path
from libnano.helpers.jsonbytes import base_decode_dict
from libnano.seqstr import reverse_complement as _rc  # type: ignore

LOCAL_DIR: str = op.dirname(
    op.realpath(__file__),
)

REGEX_BASE_LUT: Dict[str, str] = {
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'T',
    'R': '[AG]',    # purine
    'Y': '[CG]',    # pyrimidine
    'K': '[GT]',    # keto
    'M': '[AC]',    # amino
    'S': '[GC]',    # strong bonds
    'W': '[AT]',    # weak bonds
    'B': '[GTC]',   # all but A
    'D': '[GAT]',   # all but C
    'H': '[ACT]',   # all but G
    'V': '[GCA]',   # all but T
    'N': 'N',       # leave N's put
}

NEB_ALPHABET = ''.join(REGEX_BASE_LUT.keys())

# Checks for symmetric + internal cutsite annotation
NEB_I_SHORTHAND_RE = re.compile(
    r'^[\^' + NEB_ALPHABET + ']+$',
)

# Pulls out indices and sequence from a parenthetical cutsite annotation
NEB_P_SHORTHAND_RE = re.compile(
    r'^(?:\((?P<startidx>[\d|-]+)/(?P<startrevidx>[\d|-]+)\))?'
    '(?P<seq>[' + NEB_ALPHABET + ']+)'
    r'(?:\((?P<endidx>[\d|-]+)/(?P<endrevidx>[\d|-]+)\))?$',
)


class EnzymeTuple(NamedTuple):
    name: str
    prototype: str
    microorganism: str
    source: str
    recognition_sequence: str
    methylation_site: str
    commercial_availability: str
    references: List[int]


def conditional_int(
        i: Optional[str],
) -> Optional[int]:
    '''
    Args:
        i: string or None

    Returns:
        integer version of the string or None
    '''
    return int(i) if i is not None else None


def seq_to_regex(
        seq: str,
) -> str:
    '''Convert sequence to regex

    Args:
        seq: Sequence to Convert

    Returns:
        Regex expression to match the sequence against
    '''
    if seq is None:
        return None
    regex_str = ''
    prev_b = ''
    num_n = 0
    for b in seq:
        if prev_b == 'N' and b != 'N':
            regex_str += '[ATGC]{%d}' % num_n
            regex_str += REGEX_BASE_LUT[b]
            num_n = 0
        elif b == 'N':
            num_n += 1
        else:
            regex_str += REGEX_BASE_LUT[b]
            pass
        prev_b = b
    if prev_b == 'N':
        regex_str += '[ATGC]{%d}' % num_n
    return regex_str


def get_rebase_list(
) -> Tuple[int, List[EnzymeTuple], List[str]]:
    '''Get the latest list of restriction enzymes from NEB Rebase.

        The enzyme list is a list of tuples with the following fields as
    desribed on NEB Rebase:

    <ENZYME NAME>
    <PROTOTYPE>
    <MICROORGANISM>
    <SOURCE>
    <RECOGNITION SEQUENCE>
    <METHYLATION SITE>
    <COMMERCIAL AVAILABILITY>
    <REFERENCES>

    The reference list is indexed by reference number (i.e., reference 1 can
    be found at index 1)

    Returns:
        Tuple of the form::
            (
                <Database version>,
                <a list of restriction enzymes>,
                <the respective references>,
            )

    '''
    enzyme_list_url = 'http://rebase.neb.com/rebase/link_allenz'
    enzyme_list_re = r''.join(r'<%d>(.*)\n' % n for n in range(1, 9))
    reference_list_re = r'(\d+)[.][ \t]+(.*)\n'
    request_headers = {
        'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.11 ' +
                      '(KHTML, like Gecko) Chrome/23.0.1271.64 Safari/537.11',
        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,' +
                  '*/*;q=0.8',
        'Accept-Charset': 'ISO-8859-1,utf-8;q=0.7,*;q=0.3',
        'Accept-Encoding': 'none',
        'Accept-Language': 'en-US,en;q=0.8',
        'Connection': 'keep-alive',
    }

    req = urllib2.Request(
        enzyme_list_url,
        headers=request_headers,
    )
    res = urllib2.urlopen(req)
    raw_data = res.read().decode('utf8')
    reference_section = re.search(      # type: ignore
        r'References:(.*)',
        raw_data,
        flags=re.S,
    ).group(1)
    database_version = int(      # type: ignore
        re.search(r'allenz\.(\d+)', raw_data).group(1),  # type: ignore
    )
    enzyme_data_raw = re.findall(
        enzyme_list_re,
        raw_data,
    )
    enzyme_data = [
        EnzymeTuple(
            name=ed[0],
            prototype=ed[1],
            microorganism=ed[2],
            source=ed[3],
            recognition_sequence=ed[4],
            methylation_site=ed[5],
            commercial_availability=ed[6],
            references=ed[7].split(','),
        )
        for ed in enzyme_data_raw
    ]

    references = [None] + [
        ref for _, ref in re.findall(
            reference_list_re, reference_section,
        )
    ]
    return database_version, enzyme_data, references


def process_neb_shorthand(
        neb_shorthand: str,
) -> Dict[str, Any]:
    '''Process NEB shorthand for an enzyme cutsite.

    Args:
        neb_shorthand: a single cutsite represented by NEB shorthand

    Returns:
        Dictionary of the form::

            {
                'core_seq': [fwd, revcomp],
                'full_seq': [fwd, revcomp],
                'core_regex': [fwd, revcomp],
                'full_regex': [fwd, revcomp],
                'cutsite_idxs' [[fwd], [rev]],
                'shorthand': <NEB format>
            }

            where `core_seq`, `full_seq`, `core_regex`, and `full_regex` are
            all lists of two strings (representing fwd and rev seqs), and
            `cutsite_idxs` is a list of lists with the following format::

                [[fwd_idx_1, fwd_idx_2, ...], [rev_idx_1, rev_idx_2, ...]]

            Indexing is from the 0th index of the expanded seq. and the cut
            occurs before the provided index

    Raises:
        ValueError (if the neb_shorthand is not in the proper format)
    '''
    neb_shorthand = neb_shorthand.strip()
    core_seq: List[str] = ['', '']
    full_seq: List[str] = ['', '']
    core_regex: List[str] = ['', '']
    full_regex: List[str] = ['', '']
    cutsite_idxs: List[List[int]] = [[], []]
    # Unknown / uncharacterized recognition sequence
    if neb_shorthand == '?':
        pass
    # Symmetric recognition sequence with internal cutsites
    elif re.match(NEB_I_SHORTHAND_RE, neb_shorthand):
        if '^' in neb_shorthand:
            caret_idx = neb_shorthand.index('^')
            cutsite_idxs = [
                [caret_idx],
                [len(neb_shorthand) - caret_idx - 1],
            ]
        else:
            # NOTE skip this? Consider raising an Exception?
            pass
        b_seq = neb_shorthand.replace('^', '')
        core_seq = [b_seq, _rc(b_seq)]
        full_seq = core_seq
    else:
        match = re.match(
            NEB_P_SHORTHAND_RE,
            neb_shorthand,
        )
        # Asymmetric cutsite or cutsite outside of recognition sequence
        if match is not None:
            si, sri, ei, eri = (
                conditional_int(match.group('startidx')),
                conditional_int(match.group('startrevidx')),
                conditional_int(match.group('endidx')),
                conditional_int(match.group('endrevidx')),
            )
            seq = match.group('seq')
            core_seq = [seq, _rc(seq)]
            # Assymetric internal cutsite

            # max(si, sri) is None doesn't work in Python 3
            # not sure I am translating logic correctly
            if ei is not None and eri is not None:
                if (
                    (si is None and sri is None) and
                    max(ei, eri) < 1
                ):
                    full_seq = core_seq
                # Cutsite outside of recognition sequence
                else:
                    # IMPORTANT: Check if si and sri were found to be None by
                    # the regex to create the the sequence strings
                    si_temp: int = 0 if si is None else si
                    sri_temp: int = 0 if sri is None else sri

                    front_offset: int = max(si_temp, sri_temp, 0)
                    ex_seq = 'N' * front_offset + seq + 'N' * max(ei, eri, 0)
                    full_seq = [ex_seq, _rc(ex_seq)]
                    if si is not None and sri is not None:
                        cutsite_idxs[0].append(front_offset - si)
                        cutsite_idxs[1].append(front_offset - sri)
                    if ei is not None and eri is not None:
                        cutsite_idxs[0].append(ei + front_offset + len(seq))
                        cutsite_idxs[1].append(eri + front_offset + len(seq))
        else:
            raise ValueError(
                'Improperly formatted NEB shorthand: %s' %
                neb_shorthand,
            )
    # Build regular expressions
    if core_seq is not None:
        core_regex = [seq_to_regex(core_seq[0]), seq_to_regex(core_seq[1])]
        full_regex = [seq_to_regex(full_seq[0]), seq_to_regex(full_seq[1])]
    output_dict = {
        'core_seq': core_seq,
        'full_seq': full_seq,
        'core_regex': core_regex,
        'full_regex': full_regex,
        'cutsite_idxs': cutsite_idxs,
        'shorthand': neb_shorthand,
    }
    return output_dict


def process_enzyme_record(
        record_tuple: EnzymeTuple,
) -> Dict[str, Any]:
    '''Process a raw record tuple from NEB Rebase.

    Recall the raw record tuple indexing is as follows::

        Index       Contents
        ~~~~~       ~~~~~~~~~~~~~~~~~~~~~~~~~~
          0         <ENZYME NAME>
          1         <PROTOTYPE>
          2         <MICROORGANISM>
          3         <SOURCE>
          4         <RECOGNITION SEQUENCE>
          5         <METHYLATION SITE>
          6         <COMMERCIAL AVAILABILITY>
          7         <REFERENCES>

    Returns:
        a dictionary in the following format (see module docstring for more
            information)::

            {
              'name': <name>
              'cutsites': [{
                  'core_seq': [fwd, revcomp],
                  'full_seq': [fwd, revcomp],
                  'core_regex': [fwd, revcomp],
                  'full_regex': [fwd, revcomp],
                  'cutsite_idxs' [[fwd], [rev]],
                  'shorthand': <NEB format>,
                  'availability', # NEB-format shorthand for avail.
              }]
            }
    '''
    name = record_tuple[0]
    neb_shorthand = record_tuple[4].strip()
    cutsites = []
    # Check for multiple sequences
    if ',' in neb_shorthand:
        neb_shorthands = neb_shorthand.split(',')
        for ns in neb_shorthands:
            cutsites.append(
                process_neb_shorthand(ns),
            )
    else:
        cutsites.append(
            process_neb_shorthand(neb_shorthand),
        )
    enzyme_rec = {
        'name': name,
        'cutsites': cutsites,
        'availability:': record_tuple[6],
    }
    return enzyme_rec


def qc_enzyme_dataset(
        enzyme_dataset_by_name: Dict,
) -> None:
    '''Quick QC of an enzyme dataset, against expected records for common
    enzymes.

    Args:
        enzyme_dataset_by_name: Expected dataset dictionary to check

    Raises:
        ValueError: Mismatch between expected value
    '''
    expected_records = {
        'EcoRI': {
            'cutsites': [
                {
                    'shorthand': 'G^AATTC',
                    'core_seq': ['GAATTC', 'GAATTC'],
                    'full_regex': ['GAATTC', 'GAATTC'],
                    'cutsite_idxs': [[1], [5]],
                    'core_regex': ['GAATTC', 'GAATTC'],
                    'full_seq': ['GAATTC', 'GAATTC'],
                },
            ],
            'availability:': 'BCFIJKMNOQRSUVXY',
            'name': 'EcoRI',
        },
        'BsaI': {
            'cutsites': [
                {
                    'shorthand': 'GGTCTC(1/5)',
                    'core_seq': ['GGTCTC', 'GAGACC'],
                    'full_regex': ['GGTCTC[ATGC]{5}', '[ATGC]{5}GAGACC'],
                    'cutsite_idxs': [[7], [11]],
                    'core_regex': ['GGTCTC', 'GAGACC'],
                    'full_seq': ['GGTCTCNNNNN', 'NNNNNGAGACC'],
                },
            ],
            'availability:': 'N',
            'name': 'BsaI',
        },
        'BaeI': {
            'cutsites': [
                {
                    'shorthand': '(10/15)ACNNNNGTAYC(12/7)',
                    'core_seq': ['ACNNNNGTAYC', 'GRTACNNNNGT'],
                    'full_regex': [
                        '[ATGC]{15}AC[ATGC]{4}GTA[CG]C[ATGC]{12}',
                        '[ATGC]{12}G[AG]TAC[ATGC]{4}GT[ATGC]{15}',
                    ],
                    'cutsite_idxs': [[5, 38], [0, 33]],
                    'core_regex': [
                        'AC[ATGC]{4}GTA[CG]C',
                        'G[AG]TAC[ATGC]{4}GT',
                    ],
                    'full_seq': [
                        'NNNNNNNNNNNNNNNACNNNNGTAYCNNNNNNNNNNNN',
                        'NNNNNNNNNNNNGRTACNNNNGTNNNNNNNNNNNNNNN',
                    ],
                },
            ],
            'availability:': 'N',
            'name': 'BaeI',
        },
    }
    if isinstance(
        list(enzyme_dataset_by_name.keys())[0],
        str,
    ):
        def coerce_b(s):
            return s
    else:
        def coerce_b(s):
            return s.encode('utf-8')

        expected_records = base_decode_dict(
            expected_records,
        )
    for name, rec in expected_records.items():
        new_rec = enzyme_dataset_by_name[name]
        for i, cutsite in enumerate(rec[coerce_b('cutsites')]):
            for k, v in cutsite.items():    # type: ignore
                new_value = new_rec[coerce_b('cutsites')][i][k]
                if not new_value == v:
                    raise ValueError(
                        f'Mismatch between expected value ({v}) '
                        f'and value ({new_value}) for key {k} for enzyme '
                        f'{name}',
                    )


def update_enzyme_dataset():
    '''Update the enzyme dataset if a new version is available.

    The top-level JSON file structure is as follows::

        {
            'rebase_version': <integer version # of Rebase database>,
            'enzyme_data: <dict of enzyme data keyed by enzyme name>,
        }
    '''
    needs_update = True
    current_version = -1
    dataset_fp = op.join(
        LOCAL_DIR,
        'enzyme_dataset.json',
    )
    # Get the latest data from NEB Rebase
    latest_version, enzyme_data, _ = get_rebase_list()
    # Check the version of the current dataset
    try:
        with open(dataset_fp) as fd:
            current_dataset = json.load(fd)
            current_version = current_dataset['rebase_version']
            if latest_version == current_version:
                needs_update = False
                needs_update = True
    except (IOError, OSError):
        pass
    if needs_update:
        enzyme_recs = list(
            map(
                process_enzyme_record,
                enzyme_data,
            ),
        )
        enzyme_data_by_name = {
            rec['name']: rec
            for rec in enzyme_recs
        }
        qc_enzyme_dataset(
            enzyme_data_by_name,
        )
        print(
            'New NEB Rebase version %d supercedes verion %d, dataset '
            'loaded and validated: %s' % (
                latest_version, current_version,
                dataset_fp,
            ),
        )
        with open(dataset_fp, 'w') as fd:
            json.dump(
                {
                    'rebase_version': latest_version,
                    'enzyme_data': enzyme_data_by_name,
                }, fd,
            )
    else:
        print(
            'Latest NEB Rebase version (%d) matches latest dataset, no '
            'update necessary' % current_version,
        )


if __name__ == '__main__':
    update_enzyme_dataset()
    # version, res, _ = get_rebase_list()
    # print(version)
    # processed_recs = list(map(process_enzyme_record, res))
    # print(REGEX_BASE_LUT)


# ~~~~~~~~~~~ Needs to be updated for new cutsite indexing scheme ~~~~~~~~~~~ #

# from blessings import Terminal

# term = Terminal()

# def visualize_cut_site(seq, cutsites, term_bg='black'):
#     if term_bg == 'black':
#         colorA = term.black_on_white
#         colorB = term.black_on_green
#     else:
#         colorA = term.black_on_bright_black
#         colorB = term.black_on_green

#     compseq = seqstr.complement(seq)

#     if isinstance(cutsites, int):   # one cutsite
#         i = cutsites
#         sl = len(seq)
#         if i != 0:
#             print('5\'-' + \
#                     colorA(seq[:i]) + \
#                     colorB(seq[i:]) + \
#                     '-5\'')
#             print('3\'-' + \
#                     colorA(compseq[:sl-i]) + \
#                     colorB(compseq[sl-i:]) + \
#                     '-5\'')
#         else:
#             print('5\'-' + \
#                     colorA(seq) + \
#                     '-5\'')
#             print('3\'-' + \
#                     colorB(compseq) + \
#                     '-5\'')
#     elif isinstance(cutsites, tuple): # tuple
#         max3p = max5p = 0

#         if cutsites[0] is not None:
#             max5p = max(cutsites[0])
#             Ns5p_fwd = cutsites[0][0]
#             Ns3p_rev = cutsites[0][1]
#             if Ns5p_fwd < Ns3p_rev:
#                 idx5p_fwd = Ns3p_rev - Ns5p_fwd
#                 idx3p_rev = 0
#             else:
#                 idx5p_fwd = 0
#                 idx3p_rev = Ns5p_fwd - Ns3p_rev
#         if cutsites[1] is not None:
#             max3p = max(cutsites[1])
#             Ns3p_fwd = cutsites[1][0]
#             Ns5p_rev = cutsites[1][1]

#             if Ns3p_fwd < Ns5p_rev:
#                 idx3p_fwd = max5p + len(seq) + Ns3p_fwd
#                 idx5p_rev = len(seq) + max5p + max3p
#             else:
#                 idx3p_fwd = len(seq) + max5p + max3p
#                 idx5p_rev = max5p + len(seq) + Ns5p_rev

#         fwdseq = 'N'*max5p + seq + 'N'*max3p
#         revseq = 'N'*max5p + compseq + 'N'*max3p

#         if idx5p_fwd !=  0:
#             print('5\'-' + \
#                     colorA(fwdseq[:idx5p_fwd]) + \
#                     colorB(fwdseq[idx5p_fwd:]) + \
#                     '-5\'')
#             print('3\'-' + \
#                     colorB(revseq[idx3p_rev:idx5p_rev]) + \
#                     colorA(revseq[idx5p_rev:]) + \
#                     '-5\'')
#         else:
#             print('5\'-' + \
#                     colorB(fwdseq[idx5p_fwd:idx3p_fwd]) + \
#                     colorA(fwdseq[idx3p_fwd:]) + \
#                     '-5\'')
#             print('3\'-' + \
#                     colorA(revseq[:idx3p_rev]) + \
#                     colorB(revseq[idx3p_rev:]) + \
#                     '-5\'')
#     elif isinstance(cutsites, list): # list
#         max3p = max(cutsites)
#         Ns3p_fwd = cutsites[0]
#         Ns5p_rev = cutsites[1]

#         fwdseq = seq + 'N'*max3p
#         revseq = compseq + 'N'*max3p

#         if Ns3p_fwd < Ns5p_rev:
#             idx3p_fwd = len(seq) + Ns3p_fwd
#             if Ns5p_rev < 0:
#                 idx5p_rev = len(seq) + Ns5p_rev
#             else:
#                 idx5p_rev = len(revseq)

#             print('5\'-' + \
#                     colorA(fwdseq[:idx3p_fwd]) + \
#                     colorB(fwdseq[idx3p_fwd:]) + \
#                     '-5\'')
#             print('3\'-' + \
#                     colorA(revseq[:idx5p_rev]) + \
#                     colorB(revseq[idx5p_rev:]) + \
#                     '-5\'')
#         else:
#             if Ns3p_fwd < 0:
#                 idx3p_fwd = len(seq) + Ns3p_fwd
#             else:
#                 idx3p_fwd = len(fwdseq)
#             idx5p_rev = len(seq) + Ns5p_rev
#             print('5\'-' + \
#                     colorA(fwdseq[:idx3p_fwd]) +
#                     colorB(fwdseq[idx3p_fwd:]) +
#                     '-5\'')
#             print('3\'-' + \
#                     colorA(revseq[:idx5p_rev]) + \
#                     colorB(revseq[idx5p_rev:]) + \
#                     '-5\'')
#     else:  # None
#         print('5\'-' + \
#                 colorA(seq) + \
#                 '-5\'')
#         print('3\'-' + \
#                 colorA(compseq) + \
#                 '-5\'')
