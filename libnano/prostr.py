# -*- coding: utf-8 -*-
'''libnano.prostr
~~~~~~~~~~~~~~~~~~~

Methods to translate from DNA/RNA to amino acid sequence and vice versa

'''

import copy
import bisect
from collections import Counter
import math
import random
import sys
from typing import (
    Dict
)

from libnano.datasets import dataset_container

# Modified from http://eli.thegreenplace.net/2010/01/22/
#                       weighted-random-generation-in-python
class _WeightedRandomGenerator(object):
    '''
    '''
    def __init__(self, weights, values, rng=None):
        self.totals = []
        self.values = values
        self.rng = rng or random
        running_total = 0

        for w in weights:
            running_total += w
            self.totals.append(running_total)
        self.total = running_total

    def next(self):
        rnd = self.rng.random() * self.total
        return self.values[bisect.bisect_right(self.totals, rnd)]

    def __call__(self):
        return self.next()


class ReverseTranslator(object):
    ''' Basic reverse translator that randomly assigns codons for each amino
    acid.
    '''

    def __init__(self, aa_seq: str, na_type: str = 'DNA', codon_lut: dict = None):
        '''
        Args:
            aa_seq (str): primary amino acid sequence in single letter code
            na_type (str): whether to return DNA-space or RNA-space sequences
            codon_lut (Optional[dict]): overrides the standard lut

        '''
        self.aa_seq = aa_seq
        self.na_type = na_type
        rng = random.Random(aa_seq)

        codon_lut = codon_lut or dataset_container.AA_TO_DNA
        self.codon_gens = self._generateRandomCodonGens(codon_lut, na_type,
                                                        rng)
    # end def

    def _generateRandomCodonGens(self, codon_lut: dict, na_type: str, rng):
        # Need to wrap lambda function for list scoping
        def _gen(codon_list):
            return lambda: rng.choice(codon_list)
        codon_gens = {}
        for aa in set(list(self.aa_seq)):
            codon_list = (codon_lut[aa] if na_type == 'DNA' else
                          [c.replace('T', 'U') for c in codon_lut[aa]])
            codon_gens[aa] = _gen(codon_list)
        return codon_gens
    # end def

    def revTranslate(self) -> str:
        '''Generate a random reverse translation

        Returns:
            str: reverse translation of `self.aa_seq`
        '''
        trans = ''.join([self.codon_gens[aa]() for aa in self.aa_seq])
        return trans
    # end def
# end class


class WeightedReverseTranslator(ReverseTranslator):
    ''' Base class for providing reverse translations of amino acid sequences
    '''

    def __init__(self, aa_seq, organism, na_type='DNA', freq_threshold=0.0):
        '''
        Args:
            aa_seq (str): primary amino acid sequence in single letter code
            organism (str): organism name as per codon freq dataset
            na_type (Optional[str]): whether to return DNA-space or RNA-space sequences
            freq_threshold (Optional[float]): threshold
        '''
        self.aa_seq = aa_seq
        self.organism = organism
        self.na_type = na_type
        rng = random.Random(aa_seq)

        if isinstance(aa_seq, str):
            coerce_type = lambda s: s.decode('utf-8')
            freq_lut = \
                dataset_container.CODON_FREQ_DATASET_U[coerce_type(
                                                       organism)]
        else:
            coerce_type = lambda s: s.encode('utf-8')
            freq_lut = dataset_container.CODON_FREQ_DATASET[coerce_type(
                                                            organism)]
        self.h_freq_codon_lut = self._findHighestFreqCodons(freq_lut)
        freq_lut = self._applyFreqThreshold(freq_lut, freq_threshold)
        self.codon_gens = self._generateWeightedCodonGens(freq_lut,
                                                          na_type,
                                                          rng)
    # end def

    def _applyFreqThreshold(self, freq_lut: dict, freq_threshold: float) -> dict:
        '''
        Args:
            freq_lut (dict):
            freq_threshold (float):
        '''
        n_freq_lut = dict()
        for aa, codon_dict in freq_lut.items():
            t_codon_dict = {codon: freq for codon, freq in codon_dict.items()
                            if freq > freq_threshold}
            if len(t_codon_dict) == 0:
                raise ValueError('Frequency threshold is too high for '
                                 'amino acid %s' % aa)
            n_freq_lut[aa] = t_codon_dict
        return n_freq_lut
    # end def

    def _findHighestFreqCodons(self, freq_lut: dict) -> Dict[str, str]:
        '''
        Args:
            freq_lut (dict): Look up table for a an amino acid to get the codons
                used by frequency for an organism

        Returns:
            dict: Look up table keyed by amino acid sequence returning the
                highest frequency codon for that amino acid sequence based on
                the `freq_lut`

        '''
        n_freq_lut = {}
        for aa, codon_dict in freq_lut.items():
            h_codon = sorted(codon_dict.items(), key=lambda r: r[1])[-1][0]
            n_freq_lut[aa] = h_codon
        return n_freq_lut
    # end def

    def _generateWeightedCodonGens(self, freq_lut: dict, na_type: str, rng) -> dict:
        '''
        Args:
            freq_lut (dict):
            na_type (str):
            rng (object): random number generator

        Returns:
            dict: generated codons
        '''
        codon_gens = {}
        for aa in set(list(self.aa_seq)):
            aa_codon_freqs = freq_lut[aa].items()
            codons = [r[0] for r in aa_codon_freqs]
            if na_type == 'RNA':
                codons = [c.replace('T', 'U') for c in codons]
            freqs = [r[1] for r in aa_codon_freqs]
            codon_gens[aa] = _WeightedRandomGenerator(freqs, codons, rng)
        return codon_gens
    # end def

    def revTranslate(self, highest_freq_only: bool = False) -> str:
        '''Generate a reverse translation from a frequency table.
        Use with seqscreen module to determine sequence acceptability
        returns codon sequence in DNA

        Args:
            highest_freq_only (Optional[bool]):

        Returns:
            string
        '''
        if highest_freq_only:
            trans = ''.join(map(
                self.h_freq_codon_lut.__getitem__, self.aa_seq
            ))
        else:
            trans = ''.join([self.codon_gens[aa]() for aa in self.aa_seq])
        return trans
    # end def
# end class


def dnaToAA(seq: str) -> str:
    ''' Translate a DNA sequence into a primary amino acid sequence.

    Uses the first reading frame and assumes that the sequence contains
    complete codons. It is the user's responsibility to provide a string in
    the appropriate reading frame and of the appropriate length.

    Args:
        seq (str): DNA sequence

    Returns:
        str: converted DNA sequence to an amino acid sequence
    '''
    return ''.join([dataset_container.DNA_TO_AA[seq[i:i + 3]] for i in
                    range(0, len(seq), 3)])
# end def

def rnaToAA(seq: str) -> str:
    ''' Translate an RNA sequence into a primary amino acid sequence.

    Uses the first reading frame and assumes that the sequence contains
    complete codons. It is the user's responsibility to provide a string in
    the appropriate reading frame and of the appropriate length.

    Args:
        seq (str): RNA sequence

    Returns:
        str: converted RNA sequence to an amino acid sequence
    '''
    return ''.join([dataset_container.RNA_TO_AA[seq[i:i + 3]] for i in
                    range(0, len(seq), 3)])
# end def
