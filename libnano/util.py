"""
    libnano.util
    ~~~~~~~~~~~

    Helpful functions and data structures

"""

import random
import sys
import math

from itertools import chain
from collections import Counter

_PY3 = sys.version_info[0] == 3

if _PY3:
    maketrans = str.maketrans
else:
    from string import maketrans


# ~~~~~~~~~~~~~~~~~~~~~~~~ DNA Sequence Manipulation ~~~~~~~~~~~~~~~~~~~~~~~~ #

_DNAcomp = maketrans('ACGTacgt','TGCATGCA')

def reverseComplement(seq):
    return seq.translate(_DNAcomp)[::-1]
# end def

rc = reverseComplement

def complement(seq):
    return seq.translate(_DNAcomp)
# end def

def randomSeq(length):
    return ''.join([random.choice('ATGC') for x in range(length)])
# end def


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Filters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ (from nucleic) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
def gc_percent(seq, gc_min, gc_max):
    ''' Return True if seq has a gc percentage betweeen gc_min and gc_max
    (inclusive)

    gc_min and gc_max should be decimal percentages (i.e., .4 for 40%)

    '''
    gc = seq.count('G') + seq.count('C')
    return gc_min <= gc/len(seq) <= gc_max


def gc_run(seq, run_length):
    ''' Return True of seq has a maximum GC run length <= run_length
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


def at_run(seq, run_length):
    ''' Return True of seq has a maximum AT run length <= run_length
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


def homopol_run(seq, run_length):
    ''' Return True of seq has a maximum homopolymer run length <= run_length
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


def max_run(seq, max_a=None, max_t=None, max_g=None, max_c=None,
                 max_at=None, max_gc=None):
    ''' Return True of seq has maximum A, T, G, C, AT, or GC run <= a provided
    maximum.
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

class SeqSet(object):
    def __init__(self, max_run):
        self.seqs = []
        self.max_run = max_run

    def addSeq(self, seq):
        self.seqs.append(seq)

    def reset(self):
        self.seqs = []
        self.word_list = []

    def checkAgainstWords(self, seq):
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


def genWordList(seq, word_size, include_rc=True):
    word_list = [seq[idx:idx+word_size].upper() for idx in
                 range(len(seq) - word_size)]
    if include_rc:
        seq_rc = reverseComplement(seq)
        word_list += [seq_rc[idx:idx+word_size].upper() for idx in
                      range(len(seq_rc) - word_size)]
    return word_list
# end def

def genWordDict(seq, word_size, include_rc=True):
    word_list = genWordList(seq, word_size, include_rc)
    return Counter(word_list)
# end def

def duplicateWordCount(seq, word_size, include_rc=True):
    word_list = genWordList(seq, word_size, include_rc)
    return len([k for k,v in Counter(word_list).items() if v>1])
# end def


# ~~~~~~~~~~~~~~~~~~~~~~~~ Random sequence generation ~~~~~~~~~~~~~~~~~~~~~~~ #

def _buildBaseList(probs):
    base_list = []
    # print("the probs", probs)
    for base, freq in probs.items():
        base_list += [base] * freq
    return base_list

def weighted_choice(weights):
    choice = random.random() * sum(weights)
    for i, w in enumerate(weights):
        choice -= w
        if choice < 0:
            return i
# end def

def randSeqInv(add_length, probs, start_length=None):
    """ Probs format {'A':25, 'T':25, 'G':25, 'C':25}
    (25% liklihood of each)
    or
    Probs format {'A':0.25, 'T':0.25, 'G':0.25, 'C':0.25}

    Computes the inverse distribution from the probs

    start_length is the length that the original weights
    was calculated from and if not None will rebalance the
    distribution as the sequence is built up
    """
    seq = ''
    L = start_length
    weights = list(probs.values())
    bases = list(probs.keys())
    # bases = ['A', 'T', 'G', 'C']
    invweights = [0,0,0,0] # initialize
    for x in range(add_length):
        enum_weights = [x for x in enumerate(weights)]
        # sort by weight
        enum_weights_sorted = sorted(enum_weights, key=lambda x: x[1])
        # swap high and low
        invweights[invweight_sorted[0][0]] = invweight_sorted[3][1]
        invweights[invweight_sorted[3][0]] = invweight_sorted[0][1]
        invweights[invweight_sorted[1][0]] = invweight_sorted[2][1]
        invweights[invweight_sorted[2][0]] = invweight_sorted[1][1]

        base_idx = weighted_choice(invweights)
        seq += bases[base_idx]

        if L is not None:
            # recalculate weights
            weights[base_idx] = (weights[base_idx]*L+1)/(L+1)
            for b in bases:
                if b == base_idx:
                    continue
                else:
                    weights[b] = weights[b]*L/(L+1)
# end def

def randSeq(length, probs=None):
    """ Probs format {'A':25, 'T':25, 'G':25, 'C':25} (25% liklihood of each)
    """
    if probs is None:
        probs = {'A':25, 'T':25, 'G':25, 'C':25}
    weights = list(probs.values())
    bases = list(probs.keys())
    return ''.join([bases[weighted_choice(weights)] for x in range(length)])
# end def

def _checkAgainstWords(seq, word_list):
    word_size = len(word_list[0])
    sub_seqs = [seq[idx:idx+word_size].upper() for idx in
                 range(len(seq) - word_size)]
    if len(list(set(sub_seqs) & set(word_list))) > 0: # Length of intersection
        return False
    else:
        return True
# end def

def randSeqAvoidWords(length, flank_left='', flank_right='', probs=None,
                      word_list=[]):
    while True:
        raw_rand = randSeq(length, probs=probs)
        if homopol_run(raw_rand, 3):
            if _checkAgainstWords(flank_left + raw_rand + flank_right,
                                  word_list=word_list):
                return raw_rand
# end def

# TODO: add filtering capabilies a la nucleic (avoid GC runs etc)
