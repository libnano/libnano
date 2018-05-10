# -*- coding: utf-8 -*-
import unittest
import json
import itertools

import _setup
from libnano.prostr import (
    dnaToAA,
    rnaToAA,
    ReverseTranslator,
    WeightedReverseTranslator
)
from libnano.datasets import dataset_container

class TestProStr(unittest.TestCase):

    def setUp(self):
        pass
    # end def

    def test_ReverseTranslator(self):
        aa_seq = ''.join(dataset_container.AA_TO_DNA.keys())
        # DNA random
        a = ReverseTranslator(aa_seq)
        self.assertIsNotNone(a)
        for _ in range(1000):
            rt = a.revTranslate()
            self.assertEqual(dnaToAA(rt), aa_seq)

        # RNA random
        a = ReverseTranslator(aa_seq, na_type='RNA')
        self.assertIsNotNone(a)
        for _ in range(1000):
            self.assertEqual(rnaToAA(a.revTranslate()), aa_seq)

        # Weighted by organism codon freqs
        for organism in dataset_container.CODON_FREQ_DATASET.keys():
            # DNA, weighted
            a = WeightedReverseTranslator(aa_seq, organism)
            self.assertIsNotNone(a)
            for _ in range(1000):
                self.assertEqual(dnaToAA(a.revTranslate()), aa_seq)
            # RNA, weighted
            a = WeightedReverseTranslator(aa_seq, organism, na_type='RNA')
            self.assertIsNotNone(a)
            for _ in range(1000):
                self.assertEqual(rnaToAA(a.revTranslate()), aa_seq)
    # end def

    def test_dnaToAA(self):
        self.assertEqual(dnaToAA('TTTTTCAGT'), 'FFS')
        with self.assertRaises(KeyError):
            dnaToAA('TTTTTCAGTA')  # seq % 3 != 0
        with self.assertRaises(KeyError):
            dnaToAA('TTTTTCAGX')  # Invalid character
        # Exhaustive (test translation of every codon)
        dna_seq = ''
        aa_seq = ''
        for aa, codons in dataset_container.AA_TO_DNA.items():
            dna_seq += ''.join(codons)
            aa_seq += aa * len(codons)
        self.assertEqual(dnaToAA(dna_seq), aa_seq)
    # end def

    def test_rnaToAA(self):
        self.assertEqual(rnaToAA('UUUUUCAGU'), 'FFS')
        with self.assertRaises(KeyError):
            rnaToAA('UUUUUCAGUU')  # seq % 3 != 0
        with self.assertRaises(KeyError):
            rnaToAA('UUUUUCAGX')  # Invalid character
        # Exhaustive (test translation of every codon)
        rna_seq = ''
        aa_seq = ''
        for aa, codons in dataset_container.AA_TO_RNA.items():
            rna_seq += ''.join(codons)
            aa_seq += aa * len(codons)
        self.assertEqual(rnaToAA(rna_seq), aa_seq)
# end class


if __name__ == '__main__':
    unittest.main(verbosity=2)