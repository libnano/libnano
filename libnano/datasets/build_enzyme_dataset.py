# -*- coding: utf-8 -*-
'''libnano.datasets.build_enzyme_dataset

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
import re
import os
import sys
import urllib.request as urllib2
from typing import (
    Dict,
    List,
    Tuple
)

from libnano import seqstr
from libnano.helpers.jsonbytes import base_decode_dict

LOCAL_DIR: str = os.path.dirname(os.path.realpath(__file__))

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
    'N': 'N'        # leave N's put
}

NEB_ALPHABET: str = ''.join(REGEX_BASE_LUT.keys())

# Checks for symmetric + internal cutsite annotation
NEB_I_SHORTHAND_RE: '_sre.SRE_Pattern' = re.compile('^[\^' + NEB_ALPHABET + ']+$')

# Pulls out indices and sequence from a parenthetical cutsite annotation
NEB_P_SHORTHAND_RE: '_sre.SRE_Pattern'  = re.compile(
    r'^(?:\((?P<startidx>[\d|-]+)/(?P<startrevidx>[\d|-]+)\))?'
     '(?P<seq>[' + NEB_ALPHABET + ']+)'
    r'(?:\((?P<endidx>[\d|-]+)/(?P<endrevidx>[\d|-]+)\))?$')

rc = seqstr.reverseComplement
condInt = lambda i: int(i) if i is not None else None

def seqToRegex(seq: str) -> str:
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


def getRebaseList() -> Tuple[int, List, List]:
    ''' Get the latest list of restriction enzymes from NEB Rebase.

    Returns:
        Database version, a list of restriction enzymes, and
        the respective references.

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

    '''
    ENZYME_LIST_URL = 'http://rebase.neb.com/rebase/link_allenz'
    ENZYME_LIST_RE = r''.join(r'<%d>(.*)\n' % n for n in range(1,9))
    REFERENCE_LIST_RE = r'(\d+)[.][ \t]+(.*)\n'
    REQUEST_HEADERS = {
        'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.11 ' + \
                      '(KHTML, like Gecko) Chrome/23.0.1271.64 Safari/537.11',
        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,' + \
                  '*/*;q=0.8',
        'Accept-Charset': 'ISO-8859-1,utf-8;q=0.7,*;q=0.3',
        'Accept-Encoding': 'none',
        'Accept-Language': 'en-US,en;q=0.8',
        'Connection': 'keep-alive'
    }

    req = urllib2.Request(ENZYME_LIST_URL, headers=REQUEST_HEADERS)
    res = urllib2.urlopen(req)
    raw_data = res.read().decode('utf8')
    reference_section = re.search(r'References:(.*)', raw_data,
                                  flags=re.S).group(1)
    database_version = int(re.search(r'allenz\.(\d+)', raw_data).group(1))
    enzyme_data = re.findall(ENZYME_LIST_RE, raw_data)
    references = [None] + [ref for ref_num, ref in re.findall(
                           REFERENCE_LIST_RE, reference_section)]
    return database_version, enzyme_data, references


def processNebShorthand(neb_shorthand: str) -> dict:
    ''' Process NEB shorthand for an enzyme cutsite.

    Args:
        neb_shorthand (str): a single cutsite represented by NEB shorthand

    Returns:
        dictionary of the form::

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
    core_seq = [None, None]
    full_seq = [None, None]
    core_regex = [None, None]
    full_regex = [None, None]
    cutsite_idxs = [[], []]
    # Unknown / uncharacterized recognition sequence
    if neb_shorthand == '?':
        pass
    # Symmetric recognition sequence with internal cutsites
    elif re.match(NEB_I_SHORTHAND_RE, neb_shorthand):
        if '^' in neb_shorthand:
            caret_idx = neb_shorthand.index('^')
            cutsite_idxs = [[caret_idx, None],
                            [len(neb_shorthand)-caret_idx-1, None]]
        else:
            pass
        b_seq = neb_shorthand.replace('^', '')
        core_seq = [b_seq, rc(b_seq)]
        full_seq = core_seq
    else:
        match = re.match(NEB_P_SHORTHAND_RE, neb_shorthand)
        # Asymmetric cutsite or cutsite outside of recognition sequence
        if match is not None:
            si, sri, ei, eri = (condInt(match.group('startidx')),
                                condInt(match.group('startrevidx')),
                                condInt(match.group('endidx')),
                                condInt(match.group('endrevidx')))
            seq = match.group('seq')
            core_seq = [seq, rc(seq)]
            # Assymetric internal cutsite

            # max(si, sri) is None doesn't work in Python 3
            # not sure I am translating logic correctly
            if (si is None and sri is None) and max(ei, eri) < 1:
                full_seq = core_seq
            # Cutsite outside of recognition sequence
            else:
                if si is None:
                    si = 0
                if sri is None:
                    sri = 0
                front_offset = max(si, sri, 0)
                ex_seq = 'N' * max(si, sri, 0) + seq + 'N' * max(ei, eri, 0)
                full_seq = [ex_seq, rc(ex_seq)]
                if max(si, sri) is not None:
                    cutsite_idxs[0].append(front_offset - si)
                    cutsite_idxs[1].append(front_offset - sri)
                if max(ei, eri) is not None:
                    cutsite_idxs[0].append(ei + front_offset + len(seq))
                    cutsite_idxs[1].append(eri + front_offset + len(seq))
        else:
            raise ValueError('Improperly formatted NEB shorthand: %s' %
                             neb_shorthand)
    # Build regular expressions
    if core_seq is not None:
        core_regex = [seqToRegex(core_seq[0]), seqToRegex(core_seq[1])]
        full_regex = [seqToRegex(full_seq[0]),
                          seqToRegex(full_seq[1])]
    output_dict = {
        'core_seq': core_seq,
        'full_seq': full_seq,
        'core_regex': core_regex,
        'full_regex': full_regex,
        'cutsite_idxs': cutsite_idxs,
        'shorthand': neb_shorthand
    }
    return output_dict


def processEnzymeRecord(record_tuple: tuple) -> dict:
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
    neb_shorthand = record_tuple[4].strip()
    cutsites = []
    # Check for multiple sequences
    if ',' in neb_shorthand:
        neb_shorthands = neb_shorthand.split(',')
        for ns in neb_shorthands:
            cutsites.append(processNebShorthand(ns))
    else:
        cutsites.append(processNebShorthand(neb_shorthand))
    enzyme_rec = {
        'name': record_tuple[0],
        'cutsites': cutsites,
        'availability:': record_tuple[6]
    }
    return enzyme_rec


def qcEnzymeDataset(enzyme_dataset_by_name: dict):
    '''Quick QC of an enzyme dataset, against expected records for common
    enzymes.
    '''
    expected_records = {
    'EcoRI': {
        'cutsites': [
            {
                'shorthand': 'G^AATTC',
                'core_seq': ['GAATTC', 'GAATTC'],
                'full_regex': ['GAATTC', 'GAATTC'],
                'cutsite_idxs': [[1, None], [5, None]],
                'core_regex': ['GAATTC', 'GAATTC'],
                'full_seq': ['GAATTC', 'GAATTC']
            }
        ],
        'availability:': 'BCFIJKMNOQRSUVXY',
        'name': 'EcoRI'
    },
    'BsaI': {
        'cutsites': [
            {
                'shorthand': 'GGTCTC(1/5)',
                'core_seq': ['GGTCTC', 'GAGACC'],
                'full_regex': ['GGTCTC[ATGC]{5}', '[ATGC]{5}GAGACC'],
                'cutsite_idxs': [[0, 7], [0, 11]],
                'core_regex': ['GGTCTC', 'GAGACC'],
                'full_seq': ['GGTCTCNNNNN', 'NNNNNGAGACC']
            }
        ],
        'availability:': 'N',
        'name': 'BsaI'
    },
    'BaeI': {
        'cutsites': [
            {
                'shorthand': '(10/15)ACNNNNGTAYC(12/7)',
                'core_seq': ['ACNNNNGTAYC', 'GRTACNNNNGT'],
                'full_regex': ['[ATGC]{15}AC[ATGC]{4}GTA[CG]C[ATGC]{12}',
                                   '[ATGC]{12}G[AG]TAC[ATGC]{4}GT[ATGC]{15}'],
                'cutsite_idxs': [[5, 38], [0, 33]],
                'core_regex': ['AC[ATGC]{4}GTA[CG]C', 'G[AG]TAC[ATGC]{4}GT'],
                'full_seq': ['NNNNNNNNNNNNNNNACNNNNGTAYCNNNNNNNNNNNN',
                             'NNNNNNNNNNNNGRTACNNNNGTNNNNNNNNNNNNNNN']
            }
        ],
    'availability:': 'N',
    'name': 'BaeI'
    }
    }
    if isinstance(list(enzyme_dataset_by_name.keys())[0], str):
        coerce_b = lambda s: s
    else:
        coerce_b = lambda s: s.encode('utf-8')
        expected_records = base_decode_dict(expected_records)
    for name, rec in expected_records.items():
        new_rec = enzyme_dataset_by_name[name]
        for i, cutsite in enumerate(rec[coerce_b('cutsites')]):
            for k, v in cutsite.items():
                new_value = new_rec[coerce_b('cutsites')][i][k]
                if not new_value == v:
                    raise ValueError('Mismatch between expected value (%s) '
                                     'and value (%s) for key %s for enzyme '
                                     '%s' % (v, new_value, k, name))
# end def

def updateEnzymeDataset():
    '''Update the enzyme dataset if a new version is available.

    The top-level JSON file structure is as follows::

        {
            'rebase_version': <integer version # of Rebase database>,
            'enzyme_data: <dict of enzyme data keyed by enzyme name>,
        }
    '''
    needs_update = True
    current_version = -1
    dataset_fp = os.path.join(LOCAL_DIR, 'enzyme_dataset.json')
    # Get the latest data from NEB Rebase
    latest_version, enzyme_data, _ = getRebaseList()
    # Check the version of the current dataset
    try:
        with open(dataset_fp) as fd:
            current_dataset = json.load(fd)
            current_version = current_dataset['rebase_version']
            if latest_version == current_version:
                needs_update = False
    except (IOError, OSError):
        pass
    if needs_update:
        enzyme_recs = list(map(processEnzymeRecord, enzyme_data))
        enzyme_data_by_name = {rec['name']: rec for rec in enzyme_recs}
        qcEnzymeDataset(enzyme_data_by_name)
        print('New NEB Rebase version %d supercedes verion %d, dataset '
              'loaded and validated: %s' % (latest_version, current_version,
                                            dataset_fp))
        with open(dataset_fp, 'w') as fd:
            json.dump({'rebase_version': latest_version,
                       'enzyme_data': enzyme_data_by_name}, fd)
    else:
        print('Latest NEB Rebase version (%d) matches latest dataset, no '
              'update necessary' % current_version)


if __name__ == '__main__':
    updateEnzymeDataset()
    # version, res, _ = getRebaseList()
    # print(version)
    # processed_recs = list(map(processEnzymeRecord, res))
    # print(REGEX_BASE_LUT)


# ~~~~~~~~~~~ Needs to be updated for new cutsite indexing scheme ~~~~~~~~~~~ #

# from blessings import Terminal

# term = Terminal()

# def visualizeCutSite(seq, cutsites, term_bg='black'):
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
# # end def
