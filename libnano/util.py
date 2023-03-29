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
libnano.util
~~~~~~~~~~~~

Helpful functions and data structures

'''

import random
from collections import Counter
from typing import (
    Dict,
    List,
    Optional,
)

from libnano.seqstr import reverse_complement  # type: ignore

# ~~~~~~~~~~~~~~~~~~~~~~~~ DNA Sequence Manipulation ~~~~~~~~~~~~~~~~~~~~~~~~ #


def random_dna_seq(length: int) -> str:
    '''
    Args:
        length: length of output DNA sequence

    Returns:
        Random DNA sequence of length ``length``
    '''
    return ''.join([random.choice('ATGC') for x in range(length)])


def random_rna_seq(length: int) -> str:
    '''
    Args:
        length: length of output RNA sequence

    Returns:
        Random RNA sequence of length ``length``
    '''
    return ''.join([random.choice('AUGC') for x in range(length)])


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Filters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ (from nucleic) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
def gc_percent(
        seq: str,
        gc_min: float,
        gc_max: float,
) -> bool:
    '''
    gc_min and gc_max should be decimal percentages (i.e., .4 for 40%)

    Args:
        seq: Sequence to operate on
        gc_min: GC minimum percentage (i.e., .4 for 40%)
        gc_max: GC maximum percentage (i.e., .4 for 40%)

    Returns:
        True if seq has a gc percentage betweeen gc_min and gc_max
        (inclusive)

    '''
    gc = seq.count('G') + seq.count('C')
    return gc_min <= gc / len(seq) <= gc_max


def gc_run(
        seq: str,
        run_length: int,
) -> bool:
    '''
    Args:
        seq: Sequence to operate on
        run_length: maximum GC run limit

    Returns:
        ``True`` of seq has a maximum GC run length <= run_length

    '''
    lrun = 0
    for b in seq:
        if b in 'GC':
            lrun += 1
            if lrun > run_length:
                return False
        else:
            lrun = 0
    return True


def at_run(
        seq: str,
        run_length: int,
) -> bool:
    '''
    Args:
        seq: Sequence to operate on
        run_length: maximum AT run limit

    Returns:
        ``True`` of seq has a maximum AT run length <= run_length

    '''
    lrun = 0
    for b in seq:
        if b in 'AT':
            lrun += 1
            if lrun > run_length:
                return False
        else:
            lrun = 0
    return True


def homopol_run(
        seq: str,
        run_length: int,
) -> bool:
    '''
    Args:
        seq: Sequence to operate on
        run_length: maximum homopolymer run limit

    Returns:
        ``True`` of seq has a maximum homopolymer run length <= run_length

    '''
    prev = ''
    lrun = 1
    for b in seq:
        if b == prev:
            lrun += 1
            if lrun > run_length:
                return False
        else:
            lrun = 1
        prev = b
    return True


def max_run(
        seq: str,
        max_a: Optional[float] = None,
        max_t: Optional[float] = None,
        max_g: Optional[float] = None,
        max_c: Optional[float] = None,
        max_at: Optional[float] = None,
        max_gc: Optional[float] = None,
) -> bool:
    '''
    Args:
        seq: Sequence to operate on
        max_a: Max A run limimt
        max_t: Max T run limimt
        max_g: Max G run limimt
        max_c: Max C run limimt
        max_at: Max AT run limimt
        max_gc: Max GC run limimt

    Returns:
        ``True`` of ``seq`` has maximum A, T, G, C, AT, or GC run <= a
            provided maximum.

    '''
    prev = ''
    lrun = 0
    gcrun = 0
    atrun = 0
    for b in seq:
        if b == prev:
            lrun += 1
            if max_a and b == 'A' and lrun > max_a:
                return False
            if max_t and b == 'T' and lrun > max_t:
                return False
            if max_g and b == 'G' and lrun > max_g:
                return False
            if max_c and b == 'C' and lrun > max_c:
                return False
        else:
            lrun = 0
        if b in 'GC':
            gcrun += 1
            if max_gc and gcrun > max_gc:
                return False
            atrun = 0
        elif b in 'AT':
            atrun += 1
            if max_at and atrun > max_at:
                return False
            gcrun = 0
    return True

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Word dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class SeqSet:
    '''SeqSet class

    Attributes:
        seqs: List of sequences in the set
        max_run: Maximum run guard length when screening
    '''

    def __init__(
            self,
            max_run: int,
    ):
        '''
        Args:
            max_run: Maximum run guard length when screening
        '''
        self.seqs: List[str] = []
        self.max_run = max_run

    def add_seq(self, seq: str) -> None:
        '''
        Args:
            seq: Sequence to add to the set
        '''
        self.seqs.append(seq)

    def reset(self) -> None:
        '''Clear the sequence set
        '''
        self.seqs = []

    def check_against_words(
            self,
            seq: str,
    ) -> bool:
        '''
        Args:
            seq: sequence to check against ``self.seqs``

        Returns:
            True on pass, False otherwise
        '''
        for seq2 in self.seqs:
            run_len = 0
            for idx, b in enumerate(seq):
                if b == seq2[idx]:
                    run_len += 1
                    if run_len > self.max_run:
                        return False
                else:
                    run_len = 0
        return True


def generate_word_list(
        seq: str,
        word_size: int,
        include_rc: bool = True,
) -> List[str]:
    '''
    Args:
        seq: Sequence to window for the word list
        word_size: Length of the word windows
        include_rc: If True, include reverse complements in the word lsit

    Returns:
        Word list of strings created from ``seq``

    '''
    word_list = [
        seq[idx:idx + word_size].upper() for idx in
        range(len(seq) - word_size)
    ]
    if include_rc:
        seq_rc = reverse_complement(seq)
        word_list += [
            seq_rc[idx:idx + word_size].upper() for idx in
            range(len(seq_rc) - word_size)
        ]
    return word_list


def generate_word_dict(
        seq: str,
        word_size: int,
        include_rc: bool = True,
) -> Dict:
    '''
    Args:
        seq: Sequence to window for the word list
        word_size: Length of the word windows
        include_rc: If True, include reverse complements in the word lsit

    Returns:
        Dictionary of words to counts of the word list of strings created
        from ``seq``

    '''
    word_list = generate_word_list(
        seq,
        word_size,
        include_rc,
    )
    return Counter(word_list)


def duplicate_word_count(
        seq: str,
        word_size: int,
        include_rc: bool = True,
) -> int:
    '''
    Args:
        seq: Sequence to window for the word list
        word_size: Length of the word windows
        include_rc: If True, include reverse complements in the word lsit

    Returns:
        Number of words with incidence more than 1 in the word set

    '''
    word_list = generate_word_list(seq, word_size, include_rc)
    return len([k for k, v in Counter(word_list).items() if v > 1])


# ~~~~~~~~~~~~~~~~~~~~~~~~ Random sequence generation ~~~~~~~~~~~~~~~~~~~~~~~ #

def _build_base_list(
        probs: Dict[str, int],
) -> List[str]:
    '''
    Args:
        probs: Dictionary of base to probability

    Returns:
        List of each `base` key in ``probs`` times the frequency value

    '''
    base_list = []
    for base, freq in probs.items():
        base_list += [base] * freq
    return base_list


def weighted_choice(
        weights: List[int],
) -> int:
    '''
    Args:
        weights: List of weights

    Returns:
        Choice index generated from the ``weights`` list.  -1 is returned on
        failure

    '''
    choice = random.random() * sum(weights)
    for i, w in enumerate(weights):
        choice -= w
        if choice < 0:
            return i
    return -1


def random_probability_sequence(
        length: int,
        probs: Optional[Dict[str, int]] = None,
) -> str:
    '''
    Args:
        length: Length of generated sequence
        probs: Probability per base in the format::
            {'A':25, 'T':25, 'G':25, 'C':25} (25% liklihood of each)

    Returns:
        Random sequence of length ``length``

    '''
    if probs is None:
        probs = {'A': 25, 'T': 25, 'G': 25, 'C': 25}
    weights = list(probs.values())
    bases = list(probs.keys())
    return ''.join([
        bases[weighted_choice(weights)]
        for x in range(length)
    ])


def _check_seq_against_word_list(
        seq: str,
        word_list: List[str],
) -> bool:
    '''
    Args:
        seq: Sequece to check against the ``word_list``
        word_list: List of strings to check ``seq`` against

    Returns:
        ``True`` if no intersection with word_list, else ``False``

    '''
    word_size = len(word_list[0])
    sub_seqs = [
        seq[idx:idx + word_size].upper() for idx in
        range(len(seq) - word_size)
    ]
    if len(list(set(sub_seqs) & set(word_list))) > 0:  # Length of intersection
        return False
    else:
        return True


def random_probability_sequence_avoid_word_list(
        length: int,
        flank_left: str = '',
        flank_right: str = '',
        probs: Optional[Dict[str, int]] = None,
        word_list: Optional[List[str]] = None,
        try_count_limit: int = 0,
) -> str:
    '''Generate a random sequence screened against homopolymer runs of length
    3 and screen against an optional ``word_list``

    Args:
        length: Length of generated sequence
        flank_left: Left flanking sequence
        flank_right: Right flanking sequence
        probs: Probability per base in the format::
            {'A':25, 'T':25, 'G':25, 'C':25} (25% liklihood of each)
        word_list: List of strings to check ``seq`` against

    Returns:
        Random sequence screened successfully against ``word_list``

    '''
    try_count = 0
    while True:
        raw_rand = random_probability_sequence(
            length,
            probs=probs,
        )
        if homopol_run(raw_rand, 3):
            if word_list:
                if _check_seq_against_word_list(
                    seq=f'{flank_left}{raw_rand}{flank_right}',
                    word_list=word_list,
                ):
                    return raw_rand
            else:
                return raw_rand
        if try_count_limit:
            try_count += 1
            if try_count >= try_count_limit:
                break
    return ''
